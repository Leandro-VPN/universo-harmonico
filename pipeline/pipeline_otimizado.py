# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  PIPELINE OTIMIZADO EM MEMÓRIA — Shalm 2015, dados completos              ║
# ║                                                                            ║
# ║  Mudanças principais vs versão anterior:                                  ║
# ║    1. NaN eliminado -> sentinela inteira (int32/int8 em vez de float64)  ║
# ║    2. Timetags RELATIVOS ao sync, gravados direto como uint32             ║
# ║    3. HDF5 lido e comparado em chunks -- nunca carregado inteiro          ║
# ║    4. del + gc.collect() explícitos após cada etapa pesada                ║
# ║    5. Streaming real do zip -- nunca materializa o arquivo descomprimido  ║
# ║                                                                            ║
# ║  RAM esperada: ~3-4 GB de pico (vs ~14-18 GB da versão anterior)          ║
# ╚══════════════════════════════════════════════════════════════════════════╝

import zipfile
import numpy as np
import h5py
import gc
import time
import matplotlib.pyplot as plt
from scipy.stats import chi2
from google.colab import drive

drive.mount('/content/drive')
DRIVE = "/content/drive/MyDrive/Colab Notebooks/Dados Bell"
ZIP_A = f"{DRIVE}/19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat.zip"
ZIP_B = f"{DRIVE}/19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat.zip"
HDF5  = f"{DRIVE}/19_45_CH_pockel_100kHz.run.nolightconeshift.dat.compressed.build.hdf5"

BIN_PS = 78.125
CANAL_SYNC   = 6
CANAIS_CLICK = [2, 4]

# Janela máxima entre um sync e o próximo (em bins). O laser dispara a
# ~99.14 kHz -> período ~10.087 µs -> ~129,100 bins. Usamos uint32 para o
# timetag relativo, que comporta até ~4.29 bilhões -- enorme folga.
SENTINEL_TT  = np.uint32(0xFFFFFFFF)   # "sem clique neste trial"
SENTINEL_SET = np.int8(-1)             # "sem setting"

CHUNK_BYTES = 24 * 5_000_000  # 5M registros por chunk (~120 MB de leitura bruta)

COMBOS = {
    'XX': (1, 1, 8.6, 178),
    'XY': (1, 2, 29.7, 15),
    'YX': (2, 1, 29.7, 15),
    'YY': (2, 2, 50.8, 5),
}


def log_mem(label=""):
    """Mostra o uso de RAM atual do processo (requer psutil, já vem no Colab)."""
    import psutil, os
    rss = psutil.Process(os.getpid()).memory_info().rss / 1e9
    print(f"  [RAM] {label}: {rss:.2f} GB")


# ── PASSO 1: contar trials (varredura leve, só conta syncs) ─────────────────
def contar_syncs(zip_path):
    total = 0
    with zipfile.ZipFile(zip_path, 'r') as zf:
        fname = zf.namelist()[0]
        with zf.open(fname) as f:
            while True:
                raw = f.read(CHUNK_BYTES)
                if not raw:
                    break
                n = len(raw) // 24
                data = np.frombuffer(raw, dtype=np.uint64, count=n*3).reshape(n, 3)
                total += int(np.sum(data[:, 0] == CANAL_SYNC))
                del data
    return total

print("Contando syncs (passada leve)...")
t0 = time.time()
n_syncs_a = contar_syncs(ZIP_A)
n_syncs_b = contar_syncs(ZIP_B)
N_TRIALS = min(n_syncs_a, n_syncs_b)
print(f"  Alice: {n_syncs_a:,} syncs | Bob: {n_syncs_b:,} syncs")
print(f"  N_TRIALS = {N_TRIALS:,}  ({time.time()-t0:.1f}s)")
log_mem("após contagem")


# ── PASSO 2: extração em UMA passada — sync relativo + clique relativo ─────
def extrair_otimizado(zip_path, max_trials):
    """
    Retorna:
      tt_rel  : uint32[max_trials]  -- timetag do clique relativo ao sync do
                                        próprio trial (SENTINEL se não houve clique)
      setting : int8[max_trials]    -- 0 (canal 2) ou 1 (canal 4), -1 se ausente
      sync_abs: int64[max_trials]   -- timetag absoluto de cada sync (precisa
                                        ficar em int64 pois é usado para alinhar
                                        Alice e Bob entre si)
    Processa o arquivo em UMA única passada por streaming, sem acumular
    listas de chunks na memória (ao contrário da versão anterior que fazia
    `syncs.append(...)` e só concatenava no final).
    """
    tt_rel   = np.full(max_trials, SENTINEL_TT, dtype=np.uint32)
    setting  = np.full(max_trials, SENTINEL_SET, dtype=np.int8)
    sync_abs = np.full(max_trials, -1, dtype=np.int64)

    trial_atual = -1          # índice do último sync processado
    sync_atual_tt = None      # timetag absoluto do sync corrente

    with zipfile.ZipFile(zip_path, 'r') as zf:
        fname = zf.namelist()[0]
        with zf.open(fname) as f:
            while True:
                raw = f.read(CHUNK_BYTES)
                if not raw:
                    break
                n = len(raw) // 24
                data = np.frombuffer(raw, dtype=np.uint64, count=n*3).reshape(n, 3)
                ch = data[:, 0].astype(np.int64)
                tt = data[:, 1].astype(np.int64)
                del data

                is_sync  = (ch == CANAL_SYNC)
                is_click = np.isin(ch, CANAIS_CLICK)

                idx_sync  = np.flatnonzero(is_sync)
                idx_click = np.flatnonzero(is_click)

                if idx_sync.size == 0 and idx_click.size == 0:
                    del ch, tt
                    continue

                relevantes = np.union1d(idx_sync, idx_click)
                for i in relevantes:
                    if is_sync[i]:
                        trial_atual += 1
                        sync_atual_tt = tt[i]
                        if trial_atual < max_trials:
                            sync_abs[trial_atual] = sync_atual_tt
                    else:
                        if (trial_atual >= 0 and trial_atual < max_trials
                                and tt_rel[trial_atual] == SENTINEL_TT):
                            delta = tt[i] - sync_atual_tt
                            if 0 <= delta < SENTINEL_TT:
                                tt_rel[trial_atual] = np.uint32(delta)
                                setting[trial_atual] = np.int8(0 if ch[i] == 2 else 1)

                del ch, tt, is_sync, is_click, idx_sync, idx_click, relevantes

                if trial_atual >= max_trials:
                    break

    return tt_rel, setting, sync_abs


print("\nExtraindo Alice (passada única, streaming)...")
t0 = time.time()
tt_rel_a, set_a, sync_abs_a = extrair_otimizado(ZIP_A, N_TRIALS)
print(f"  cliques válidos: {np.sum(tt_rel_a != SENTINEL_TT):,}  ({time.time()-t0:.1f}s)")
log_mem("após Alice")

print("\nExtraindo Bob (passada única, streaming)...")
t0 = time.time()
tt_rel_b, set_b, sync_abs_b = extrair_otimizado(ZIP_B, N_TRIALS)
print(f"  cliques válidos: {np.sum(tt_rel_b != SENTINEL_TT):,}  ({time.time()-t0:.1f}s)")
log_mem("após Bob")


# ── NOTA SOBRE PERFORMANCE ───────────────────────────────────────────────────
# O loop `for i in relevantes` acima ainda é Python puro e portanto lento em
# 330M de trials (pode levar 30-60 min por arquivo). Se isso for proibitivo,
# use a variante vetorizada abaixo, que evita o loop Python ao custo de mais
# RAM temporária por chunk -- ainda assim muito menor que a versão original.
#
# A ideia da vetorização: para cada chunk, calcular o índice de trial de
# CADA evento de uma vez com searchsorted (como no script original), mas
# usando os tipos compactos (uint32/int8) desde o início, em vez de NaN.
#
# def extrair_vetorizado(zip_path, max_trials):
#     tt_rel   = np.full(max_trials, SENTINEL_TT, dtype=np.uint32)
#     setting  = np.full(max_trials, SENTINEL_SET, dtype=np.int8)
#     sync_abs = np.full(max_trials, -1, dtype=np.int64)
#     trial_offset = 0  # quantos syncs já vimos em chunks anteriores
#     last_sync_tt = None
#
#     with zipfile.ZipFile(zip_path, 'r') as zf:
#         fname = zf.namelist()[0]
#         with zf.open(fname) as f:
#             while True:
#                 raw = f.read(CHUNK_BYTES)
#                 if not raw:
#                     break
#                 n = len(raw) // 24
#                 data = np.frombuffer(raw, dtype=np.uint64, count=n*3).reshape(n, 3)
#                 ch = data[:, 0].astype(np.int64)
#                 tt = data[:, 1].astype(np.int64)
#                 del data
#
#                 syncs_chunk = tt[ch == CANAL_SYNC]
#                 if last_sync_tt is not None:
#                     syncs_busca = np.concatenate(([last_sync_tt], syncs_chunk))
#                 else:
#                     syncs_busca = syncs_chunk
#
#                 # índice de trial de cada evento (relativo ao chunk)
#                 trials_no_chunk = np.searchsorted(syncs_busca, tt, side='right') - 1
#                 trials_globais = trials_no_chunk + trial_offset - (1 if last_sync_tt is not None else 0)
#
#                 valid = (trials_globais >= 0) & (trials_globais < max_trials)
#                 tr = trials_globais[valid]
#                 ev_ch = ch[valid]
#                 ev_tt = tt[valid]
#
#                 # grava syncs
#                 mask_sync = (ev_ch == CANAL_SYNC)
#                 sync_abs[tr[mask_sync]] = ev_tt[mask_sync]
#
#                 # grava cliques (primeiro de cada trial)
#                 mask_click = np.isin(ev_ch, CANAIS_CLICK)
#                 tr_c, tt_c, ch_c = tr[mask_click], ev_tt[mask_click], ev_ch[mask_click]
#                 if len(tr_c) > 0:
#                     _, first_idx = np.unique(tr_c, return_index=True)
#                     tr_final = tr_c[first_idx]
#                     # usa o sync correspondente (já gravado acima ou em chunk anterior)
#                     deltas = tt_c[first_idx] - sync_abs[tr_final]
#                     ok = (deltas >= 0) & (deltas < SENTINEL_TT)
#                     tt_rel[tr_final[ok]] = deltas[ok].astype(np.uint32)
#                     setting[tr_final[ok]] = np.where(ch_c[first_idx][ok] == 2, 0, 1).astype(np.int8)
#
#                 if len(syncs_chunk) > 0:
#                     last_sync_tt = syncs_chunk[-1]
#                     trial_offset += len(syncs_chunk)
#
#                 del ch, tt, syncs_chunk, trials_no_chunk, trials_globais, valid, tr, ev_ch, ev_tt
#
#     return tt_rel, setting, sync_abs


# ── PASSO 3: validação cruzada com HDF5, lida em chunks (sem carregar tudo) ─
print("\n" + "="*60)
print("VALIDAÇÃO CRUZADA COM HDF5 (lido em chunks)")
print("="*60)

CHUNK_HDF5 = 20_000_000
acertos_a, total_comp_a = 0, 0
acertos_b, total_comp_b = 0, 0
contagem_setting3_a = 0

with h5py.File(HDF5, 'r') as f:
    ds_set_a = f['alice/settings']
    ds_set_b = f['bob/settings']
    n_h5 = min(ds_set_a.shape[0], ds_set_b.shape[0], N_TRIALS)

    for start in range(0, n_h5, CHUNK_HDF5):
        end = min(start + CHUNK_HDF5, n_h5)

        s_a_h5 = ds_set_a[start:end]     # carrega só esta fatia
        s_b_h5 = ds_set_b[start:end]

        bruto_a = set_a[start:end]
        bruto_b = set_b[start:end]

        mask_a = (bruto_a >= 0) & np.isin(s_a_h5, [1, 2])
        if mask_a.sum() > 0:
            acertos_a += np.sum(bruto_a[mask_a] == (s_a_h5[mask_a] - 1))
            total_comp_a += mask_a.sum()

        mask_b = (bruto_b >= 0) & np.isin(s_b_h5, [1, 2])
        if mask_b.sum() > 0:
            acertos_b += np.sum(bruto_b[mask_b] == (s_b_h5[mask_b] - 1))
            total_comp_b += mask_b.sum()

        contagem_setting3_a += np.sum(s_a_h5 == 3)

        del s_a_h5, s_b_h5, bruto_a, bruto_b, mask_a, mask_b

if total_comp_a > 0:
    print(f"Alice: concordância = {100*acertos_a/total_comp_a:.2f}%  (N={total_comp_a:,})")
if total_comp_b > 0:
    print(f"Bob:   concordância = {100*acertos_b/total_comp_b:.2f}%  (N={total_comp_b:,})")
print(f"Setting=3 (Alice) em {n_h5:,} trials: {100*contagem_setting3_a/n_h5:.3f}%")

if total_comp_a > 0 and acertos_a/total_comp_a < 0.9:
    print("  ⚠️  ATENÇÃO: baixa concordância Alice — revisar mapeamento de canal.")
if total_comp_b > 0 and acertos_b/total_comp_b < 0.9:
    print("  ⚠️  ATENÇÃO: baixa concordância Bob — revisar mapeamento de canal.")

print("\n⚠️  PARE AQUI e confira a concordância antes de prosseguir.")
print("Se concordância < 90%, NÃO rode o resto — volte e ajuste o mapeamento de canal.\n")

gc.collect()
log_mem("após validação HDF5")


# ── PASSO 4: offset de clock + Δt relativo (tudo em arrays compactos) ──────
print("\n" + "="*60)
print("CALCULANDO OFFSET E COINCIDÊNCIAS")
print("="*60)

offset_sync = sync_abs_a - sync_abs_b
offset_mediano = np.median(offset_sync[offset_sync > -2])  # ignora sentinelas -1,-1
print(f"Offset mediano Alice-Bob: {offset_mediano:,.0f} bins "
      f"({offset_mediano*BIN_PS/1e9:.6f} ms)")

del offset_sync
gc.collect()


def extrair_coincidencias(setting_a, setting_b):
    """Retorna Δt (int64, pequeno array) só para os trials que interessam."""
    mask = (
        (tt_rel_a != SENTINEL_TT) & (tt_rel_b != SENTINEL_TT) &
        (set_a == setting_a) & (set_b == setting_b)
    )
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return np.array([], dtype=np.int64)

    dt = (sync_abs_a[idx].astype(np.int64) - sync_abs_b[idx].astype(np.int64)
          + tt_rel_a[idx].astype(np.int64) - tt_rel_b[idx].astype(np.int64))
    del mask, idx
    return dt


# ── PASSO 5: espectro Rayleigh (igual à versão anterior) ───────────────────
def espectro_rayleigh(difs, k_max=200):
    N = len(difs)
    if N == 0:
        return None
    ks = np.arange(1, k_max + 1)
    V = np.zeros(k_max + 1)
    p = np.ones(k_max + 1)
    Z = np.zeros(k_max + 1)
    for k in ks:
        phi = 2 * np.pi * (difs % k) / k
        C = np.mean(np.cos(phi))
        S = np.mean(np.sin(phi))
        Z[k] = 2 * N * (C**2 + S**2)
        V[k] = np.sqrt(Z[k] / (2 * N))
        p[k] = chi2.sf(Z[k], df=2)
    return {'ks': ks, 'V': V, 'p': p, 'Z': Z, 'N': N}


print("\n" + "="*60)
print("ANÁLISE HARMÔNICA — DADOS COMPLETOS")
print("="*60)

resultados = {}
for nome, (sa, sb, theta_rel, n_prev) in COMBOS.items():
    difs = extrair_coincidencias(sa - 1, sb - 1)
    N = len(difs)
    print(f"\n--- {nome} (θ_rel≈{theta_rel}°, n previsto={n_prev}) ---")
    print(f"  {N:,} coincidências")

    if N < 100:
        print("  Poucas coincidências, pulando.")
        continue

    print(f"  Δt: média={difs.mean():.1f} bins  std={difs.std():.1f} bins")

    res = espectro_rayleigh(difs, k_max=min(200, n_prev*5))
    resultados[nome] = {'res': res, 'n_prev': n_prev, 'theta_rel': theta_rel}

    if res is not None:
        for k in [n_prev, n_prev*2, n_prev*3]:
            if k <= len(res['ks']):
                print(f"    k={k:3d}  V={res['V'][k]:.6f}  Z={res['Z'][k]:.1f}  p={res['p'][k]:.2e}")

    del difs
    gc.collect()

log_mem("após análise espectral")


# ── PASSO 6: plots ───────────────────────────────────────────────────────────
n_validos = [n for n in resultados if resultados[n]['res'] is not None]

if n_validos:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Espectro Rayleigh — dados completos, pipeline otimizado',
                  fontsize=13, fontweight='bold')
    for idx, nome in enumerate(COMBOS.keys()):
        ax = axes[idx // 2, idx % 2]
        if nome in resultados and resultados[nome]['res'] is not None:
            r = resultados[nome]['res']
            n_prev = resultados[nome]['n_prev']
            ax.semilogy(r['ks'], r['V'][1:], 'b-', linewidth=0.8, alpha=0.8)
            mults = np.arange(n_prev, len(r['ks'])+1, n_prev)
            if len(mults) > 0:
                ax.scatter(mults, r['V'][mults], color='red', s=50, zorder=5,
                          label=f'múltiplos de {n_prev}')
            ax.set_title(f"{nome}  N={r['N']:,}  (n previsto={n_prev})")
            ax.legend(fontsize=8)
        else:
            ax.set_title(f"{nome}  (sem dados)")
        ax.set_xlabel('Harmônico k')
        ax.set_ylabel('V(k)')
        ax.grid(alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig('espectro_otimizado.png', dpi=150, bbox_inches='tight')
    plt.show()

print("\n" + "="*60)
print("CONCLUÍDO")
print("="*60)
log_mem("final")
