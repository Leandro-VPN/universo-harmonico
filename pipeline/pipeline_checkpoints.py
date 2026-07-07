# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PIPELINE COM CHECKPOINTS — processa TODOS os 330M, em blocos             ║
# ║                                                                            ║
# ║  Arquitetura:                                                             ║
# ║    1. Processa N_BLOCO trials por vez (default 10M)                      ║
# ║    2. Acumula estatísticas do espectro de Rayleigh (soma de cos/sin)     ║
# ║       que são ADITIVAS -- somar por blocos = somar tudo de uma vez       ║
# ║    3. Salva checkpoint em disco após cada bloco                          ║
# ║    4. Se a sessão cair, rodar de novo retoma do último checkpoint        ║
# ║    5. RAM de pico ~ tamanho de 1 bloco, não dos 330M inteiros            ║
# ║                                                                            ║
# ║  Para rodar do zero: apague a pasta CHECKPOINT_DIR                       ║
# ║  Para retomar: só rode o script de novo, sem mudar nada                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

import zipfile
import numpy as np
import h5py
import gc
import json
import os
import time
import matplotlib.pyplot as plt
from scipy.stats import chi2
from google.colab import drive

drive.mount('/content/drive')
DRIVE = "/content/drive/MyDrive/Colab Notebooks/Dados Bell"
ZIP_A = f"{DRIVE}/19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat.zip"
ZIP_B = f"{DRIVE}/19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat.zip"
HDF5  = f"{DRIVE}/19_45_CH_pockel_100kHz.run.nolightconeshift.dat.compressed.build.hdf5"

CHECKPOINT_DIR = f"{DRIVE}/checkpoints_bell"
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

BIN_PS = 78.125
CANAL_SYNC   = 6
CANAIS_CLICK = [2, 4]
SENTINEL_TT  = np.uint32(0xFFFFFFFF)
SENTINEL_SET = np.int8(-1)

N_BLOCO     = 10_000_000   # trials por bloco -- ajuste conforme a RAM disponível
CHUNK_BYTES = 24 * 5_000_000

COMBOS = {
    'XX': (1, 1, 8.6, 178),
    'XY': (1, 2, 29.7, 15),
    'YX': (2, 1, 29.7, 15),
    'YY': (2, 2, 50.8, 5),
}
K_MAX_GLOBAL = max(n*5 for _, _, _, n in COMBOS.values())  # maior k_max necessário


def log_mem(label=""):
    import psutil
    rss = psutil.Process(os.getpid()).memory_info().rss / 1e9
    print(f"  [RAM] {label}: {rss:.2f} GB")


# ── Estado persistente em disco ──────────────────────────────────────────────
PROGRESSO_PATH = f"{CHECKPOINT_DIR}/progresso.json"
ACUMULADO_PATH = f"{CHECKPOINT_DIR}/acumulado.npz"


def carregar_progresso():
    if os.path.exists(PROGRESSO_PATH):
        with open(PROGRESSO_PATH) as f:
            return json.load(f)
    return {"bloco_atual": 0, "byte_offset_a": 0, "byte_offset_b": 0,
            "trial_offset_a": 0, "trial_offset_b": 0,
            "ultimo_sync_a": None, "ultimo_sync_b": None,
            "concluido": False}


def salvar_progresso(estado):
    with open(PROGRESSO_PATH, 'w') as f:
        json.dump(estado, f)


def carregar_acumulado():
    """
    Acumuladores: para cada combo, soma_cos[k], soma_sin[k], N total.
    Como cos/sin somam linearmente, blocos diferentes apenas se somam.
    """
    if os.path.exists(ACUMULADO_PATH):
        d = np.load(ACUMULADO_PATH)
        return {nome: {
                    'soma_cos': d[f'{nome}_cos'].copy(),
                    'soma_sin': d[f'{nome}_sin'].copy(),
                    'N': int(d[f'{nome}_N'])
                } for nome in COMBOS}
    return {nome: {
                'soma_cos': np.zeros(K_MAX_GLOBAL + 1),
                'soma_sin': np.zeros(K_MAX_GLOBAL + 1),
                'N': 0
            } for nome in COMBOS}


def salvar_acumulado(acc):
    kwargs = {}
    for nome, d in acc.items():
        kwargs[f'{nome}_cos'] = d['soma_cos']
        kwargs[f'{nome}_sin'] = d['soma_sin']
        kwargs[f'{nome}_N']   = d['N']
    np.savez(ACUMULADO_PATH, **kwargs)


# ── Extração de UM bloco, retomando de onde parou no arquivo ───────────────
def extrair_bloco(zip_path, byte_offset, sync_offset, n_trials_bloco):
    """
    Lê o arquivo zip a partir de byte_offset (posição em bytes dentro do
    stream descomprimido) e extrai até preencher n_trials_bloco trials
    completos (ou até o arquivo acabar).

    Retorna:
      tt_rel, setting, sync_abs  -- arrays do bloco (tamanho <= n_trials_bloco)
      novo_byte_offset           -- onde parar na próxima chamada
      novo_sync_offset           -- timetag do último sync visto (para
                                     continuidade entre blocos)
      fim_arquivo                -- True se chegou ao fim do zip
    """
    tt_rel   = np.full(n_trials_bloco, SENTINEL_TT, dtype=np.uint32)
    setting  = np.full(n_trials_bloco, SENTINEL_SET, dtype=np.int8)
    sync_abs = np.full(n_trials_bloco, -1, dtype=np.int64)

    trial_local = -1
    sync_atual_tt = sync_offset
    fim_arquivo = False
    bytes_lidos_total = 0

    with zipfile.ZipFile(zip_path, 'r') as zf:
        fname = zf.namelist()[0]
        with zf.open(fname) as f:
            f.seek(byte_offset)
            while trial_local < n_trials_bloco - 1:
                raw = f.read(CHUNK_BYTES)
                if not raw:
                    fim_arquivo = True
                    break
                n = len(raw) // 24
                # se sobrou um resto de bytes (não múltiplo de 24), devolve
                resto = len(raw) - n * 24
                if resto > 0:
                    raw = raw[:-resto]
                    f.seek(f.tell() - resto)  # devolve o resto pro próximo read

                bytes_lidos_total += len(raw)

                data = np.frombuffer(raw, dtype=np.uint64, count=n*3).reshape(n, 3)
                ch = data[:, 0].astype(np.int64)
                tt = data[:, 1].astype(np.int64)
                del data

                is_sync  = (ch == CANAL_SYNC)
                is_click = np.isin(ch, CANAIS_CLICK)
                idx_sync  = np.flatnonzero(is_sync)
                idx_click = np.flatnonzero(is_click)

                if idx_sync.size == 0 and idx_click.size == 0:
                    del ch, tt, is_sync, is_click, idx_sync, idx_click
                    continue

                relevantes = np.union1d(idx_sync, idx_click)
                parou_no_meio = False
                for pos, i in enumerate(relevantes):
                    if is_sync[i]:
                        trial_local += 1
                        sync_atual_tt = tt[i]
                        if trial_local < n_trials_bloco:
                            sync_abs[trial_local] = sync_atual_tt
                        if trial_local >= n_trials_bloco - 1:
                            # bloco cheio -- registra onde paramos no CHUNK
                            # para recalcular o byte_offset exato
                            parou_no_meio = True
                            break
                    else:
                        if (trial_local >= 0 and trial_local < n_trials_bloco
                                and tt_rel[trial_local] == SENTINEL_TT):
                            delta = tt[i] - sync_atual_tt
                            if 0 <= delta < SENTINEL_TT:
                                tt_rel[trial_local] = np.uint32(delta)
                                setting[trial_local] = np.int8(0 if ch[i] == 2 else 1)

                del ch, tt, is_sync, is_click, idx_sync, idx_click, relevantes

                if parou_no_meio:
                    break

            byte_offset_final = f.tell()

    trials_obtidos = trial_local + 1
    if trials_obtidos < n_trials_bloco:
        tt_rel   = tt_rel[:trials_obtidos]
        setting  = setting[:trials_obtidos]
        sync_abs = sync_abs[:trials_obtidos]

    return tt_rel, setting, sync_abs, byte_offset_final, sync_atual_tt, fim_arquivo


# ── Acumula o espectro de Rayleigh de um bloco no total ────────────────────
def acumular_rayleigh(difs, acc_nome, acc):
    N = len(difs)
    if N == 0:
        return
    ks = np.arange(1, K_MAX_GLOBAL + 1)
    for k in ks:
        phi = 2 * np.pi * (difs % k) / k
        acc[acc_nome]['soma_cos'][k] += np.sum(np.cos(phi))
        acc[acc_nome]['soma_sin'][k] += np.sum(np.sin(phi))
    acc[acc_nome]['N'] += N


def finalizar_espectro(acc_nome, acc, n_prev):
    N = acc[acc_nome]['N']
    if N == 0:
        return None
    ks = np.arange(1, K_MAX_GLOBAL + 1)
    C = acc[acc_nome]['soma_cos'][1:] / N
    S = acc[acc_nome]['soma_sin'][1:] / N
    Z = 2 * N * (C**2 + S**2)
    V = np.sqrt(Z / (2 * N))
    p = chi2.sf(Z, df=2)
    return {'ks': ks, 'V': V, 'p': p, 'Z': Z, 'N': N}


# ── LOOP PRINCIPAL: processa blocos até cobrir TODO o arquivo ──────────────
print("="*70)
print("INICIANDO / RETOMANDO PROCESSAMENTO EM BLOCOS")
print("="*70)

estado = carregar_progresso()
acc = carregar_acumulado()

print(f"Bloco atual: {estado['bloco_atual']}")
print(f"Trials já processados -- Alice: {estado['trial_offset_a']:,} | "
      f"Bob: {estado['trial_offset_b']:,}")
for nome in COMBOS:
    print(f"  Acumulado {nome}: N={acc[nome]['N']:,}")

if estado['concluido']:
    print("\n✓ Processamento já estava completo. Pulando para análise final.")
else:
    while not estado['concluido']:
        print(f"\n--- BLOCO {estado['bloco_atual']} ---")
        t0 = time.time()

        tt_a, set_a, sync_a, bo_a, last_sync_a, fim_a = extrair_bloco(
            ZIP_A, estado['byte_offset_a'], estado['ultimo_sync_a'], N_BLOCO)
        tt_b, set_b, sync_b, bo_b, last_sync_b, fim_b = extrair_bloco(
            ZIP_B, estado['byte_offset_b'], estado['ultimo_sync_b'], N_BLOCO)

        n_bloco_real = min(len(tt_a), len(tt_b))
        print(f"  Alice: {len(tt_a):,} trials | Bob: {len(tt_b):,} trials | "
              f"usando {n_bloco_real:,} em comum")
        log_mem("após extração do bloco")

        if n_bloco_real == 0:
            print("  Bloco vazio -- fim dos dados.")
            estado['concluido'] = True
            salvar_progresso(estado)
            break

        tt_a, set_a, sync_a = tt_a[:n_bloco_real], set_a[:n_bloco_real], sync_a[:n_bloco_real]
        tt_b, set_b, sync_b = tt_b[:n_bloco_real], set_b[:n_bloco_real], sync_b[:n_bloco_real]

        # acumula coincidências por combo neste bloco
        for nome, (sa, sb, theta_rel, n_prev) in COMBOS.items():
            mask = (
                (tt_a != SENTINEL_TT) & (tt_b != SENTINEL_TT) &
                (set_a == sa - 1) & (set_b == sb - 1)
            )
            idx = np.flatnonzero(mask)
            if idx.size > 0:
                dt = (sync_a[idx].astype(np.int64) - sync_b[idx].astype(np.int64)
                      + tt_a[idx].astype(np.int64) - tt_b[idx].astype(np.int64))
                acumular_rayleigh(dt, nome, acc)
                del dt
            del mask, idx

        del tt_a, set_a, sync_a, tt_b, set_b, sync_b
        gc.collect()
        log_mem("após acumular e liberar bloco")

        # atualiza estado e salva checkpoint
        estado['bloco_atual']     += 1
        estado['byte_offset_a']    = bo_a
        estado['byte_offset_b']    = bo_b
        estado['trial_offset_a']  += n_bloco_real
        estado['trial_offset_b']  += n_bloco_real
        estado['ultimo_sync_a']    = int(last_sync_a) if last_sync_a is not None else None
        estado['ultimo_sync_b']    = int(last_sync_b) if last_sync_b is not None else None
        estado['concluido']        = fim_a or fim_b

        salvar_progresso(estado)
        salvar_acumulado(acc)

        dt_bloco = time.time() - t0
        print(f"  Bloco concluído em {dt_bloco:.1f}s. "
              f"Checkpoint salvo (trials totais: {estado['trial_offset_a']:,}).")

        if fim_a or fim_b:
            print(f"\n  Fim do arquivo alcançado "
                  f"(Alice fim={fim_a}, Bob fim={fim_b}).")


# ── ANÁLISE FINAL — usa TODO o acumulado de todos os blocos ────────────────
print("\n" + "="*70)
print("ANÁLISE FINAL — ESPECTRO ACUMULADO DE TODOS OS BLOCOS")
print("="*70)

resultados = {}
for nome, (sa, sb, theta_rel, n_prev) in COMBOS.items():
    res = finalizar_espectro(nome, acc, n_prev)
    if res is None:
        print(f"\n{nome}: sem dados.")
        continue
    print(f"\n--- {nome} (θ_rel≈{theta_rel}°, n previsto={n_prev}) ---")
    print(f"  N total de coincidências: {res['N']:,}")
    resultados[nome] = {'res': res, 'n_prev': n_prev, 'theta_rel': theta_rel}
    for k in [n_prev, n_prev*2, n_prev*3]:
        if k <= len(res['ks']):
            print(f"    k={k:3d}  V={res['V'][k-1]:.6f}  Z={res['Z'][k-1]:.1f}  p={res['p'][k-1]:.2e}")

# ── plots ────────────────────────────────────────────────────────────────
n_validos = [n for n in resultados]
if n_validos:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Espectro Rayleigh — TODOS os blocos acumulados',
                  fontsize=13, fontweight='bold')
    for idx, nome in enumerate(COMBOS.keys()):
        ax = axes[idx // 2, idx % 2]
        if nome in resultados:
            r = resultados[nome]['res']
            n_prev = resultados[nome]['n_prev']
            ax.semilogy(r['ks'], r['V'], 'b-', linewidth=0.8, alpha=0.8)
            mults = np.arange(n_prev, len(r['ks'])+1, n_prev)
            if len(mults) > 0:
                ax.scatter(mults, r['V'][mults-1], color='red', s=50, zorder=5,
                          label=f'múltiplos de {n_prev}')
            ax.set_title(f"{nome}  N={r['N']:,}  (n previsto={n_prev})")
            ax.legend(fontsize=8)
        else:
            ax.set_title(f"{nome}  (sem dados)")
        ax.set_xlabel('Harmônico k')
        ax.set_ylabel('V(k)')
        ax.grid(alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig(f'{CHECKPOINT_DIR}/espectro_final_completo.png', dpi=150, bbox_inches='tight')
    plt.show()

print("\n" + "="*70)
print("CONCLUÍDO -- todos os dados disponíveis foram processados")
print("="*70)
print(f"Checkpoints salvos em: {CHECKPOINT_DIR}")
log_mem("final")
