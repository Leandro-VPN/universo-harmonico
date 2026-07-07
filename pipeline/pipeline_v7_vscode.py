# ══════════════════════════════════════════════════════════════════════
# ║  PIPELINE v7 — BELL NIST 2015 (versão VS Code)                      ║
# ║  Variante do pipeline_v7 adaptada para Windows/VS Code              ║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════

PIPELINE BELL NIST 2015 v7 — DEFINITIVO
=========================================
Confirmado pelo comparativo_canais.py:

  ch=2 e ch=4 são os detectores de fóton (NÃO o RNG)
  Evidências:
    - 7 valores únicos de tt_rel (slots físicos do Pockels cell)
    - Alice: bins 90-96 (G1) e 104-110 (G2)
    - Bob:   bins 348-354 (G1) e 361-367 (G2)
    - std(Δt) = 9 bins = combinação dos slots discretos (física real)
    - ch=0 tem std=231 bins → GPS/overflow, não fóton

  O problema do v6 (taxa 100%) foi que a janela de tt_rel estava
  capturando o PRIMEIRO evento de ch=2 ou ch=4 por trial, mas
  ch=2 e ch=4 aparecem em ~50% dos trials cada (não 25% por canal).
  
  Correção: um trial tem clique se E SOMENTE SE ch=2 XOR ch=4 disparou
  naquele trial (não os dois). Na prática, os dois canais juntos
  representam os dois detectores possíveis — cada trial tem no máximo
  um clique real.

  Para os SETTINGS (ângulo de medição):
  A documentação NIST diz que o setting vem do RNG output.
  Com o HDF5 disponível, usamos o campo 'settings' diretamente.
  Sem HDF5, o setting é inferido pelo canal: ch=2→setting0, ch=4→setting1.
  
  IMPORTANTE: Alice e Bob usam ângulos diferentes.
  Alice: α₀=−30° (ch=2→setting0), α₁=+30° (ch=4→setting1)
  Bob:   β₀=0°   (ch=2→setting0), β₁=+30° (ch=4→setting1)
"""

import numpy as np
import os, time
from pathlib import Path
from scipy.stats import chi2
import matplotlib.pyplot as plt
import csv

# ── CONFIGURAÇÃO ──────────────────────────────────────────────────────────────
ALICE_FILE = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat"
BOB_FILE   = r"d:\BellQM\19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat"
# HDF5 opcional — se disponível, usa settings do NIST diretamente
HDF5_FILE  = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.dat.compressed.build.hdf5"
OUT_DIR    = r"d:\BellQM\out_v7"
os.makedirs(OUT_DIR, exist_ok=True)

BIN_PS     = 78.125
CHUNK_REG  = 5_000_000
N_TRIALS   = 10_000_000
SENTINEL   = np.uint32(0xFFFFFFFF)

CH_SYNC    = 6
CH_G1      = 2    # detector G1
CH_G2      = 4    # detector G2

# Slots físicos confirmados
SLOTS_A_G1 = set(range(90, 97))    # Alice G1: bins 90-96
SLOTS_A_G2 = set(range(104, 111))  # Alice G2: bins 104-110
SLOTS_B_G1 = set(range(348, 355))  # Bob G1:   bins 348-354
SLOTS_B_G2 = set(range(361, 368))  # Bob G2:   bins 361-367

COMBOS = {'XX':(0,0), 'XY':(0,1), 'YX':(1,0), 'YY':(1,1)}
N_PREV = {'XX':178, 'XY':15, 'YX':15, 'YY':5}

print("="*65)
print("PIPELINE v7 — DETECTORES ch=2 (G1) e ch=4 (G2)")
print("="*65)
print(f"Alice slots G1={sorted(SLOTS_A_G1)}, G2={sorted(SLOTS_A_G2)}")
print(f"Bob   slots G1={sorted(SLOTS_B_G1)}, G2={sorted(SLOTS_B_G2)}")

# ══════════════════════════════════════════════════════════════════════════════
# EXTRAÇÃO
# ══════════════════════════════════════════════════════════════════════════════
def extrair(filepath, slots_g1, slots_g2, n_max, label):
    """
    Extrai tt_rel e setting por trial.
    Um trial tem clique se ch=2 ou ch=4 disparou com tt_rel
    dentro dos slots físicos esperados.
    Setting = 0 se ch=2, 1 se ch=4.
    """
    tt_rel   = np.full(n_max, SENTINEL, dtype=np.uint32)
    setting  = np.full(n_max, -1,       dtype=np.int8)
    sync_abs = np.full(n_max, -1,       dtype=np.int64)

    trial = -1; sync_tt = None
    t0 = time.time()
    fsize = Path(filepath).stat().st_size
    slots_validos = slots_g1 | slots_g2

    with open(filepath, 'rb') as f:
        while trial < n_max - 1:
            raw = f.read(CHUNK_REG * 24)
            if not raw: break
            n = len(raw) // 24
            arr = np.frombuffer(raw, dtype=np.uint64).reshape(n, 3)
            ch  = arr[:, 0].astype(np.int32)
            tt  = arr[:, 1].astype(np.int64)
            del arr, raw

            idx = np.flatnonzero(
                (ch == CH_SYNC) | (ch == CH_G1) | (ch == CH_G2))
            parou = False
            for i in idx:
                c = int(ch[i]); t = int(tt[i])
                if c == CH_SYNC:
                    trial += 1
                    sync_tt = t
                    if trial < n_max:
                        sync_abs[trial] = t
                    if trial >= n_max - 1:
                        parou = True; break
                elif sync_tt is not None and trial >= 0 and trial < n_max:
                    if tt_rel[trial] == SENTINEL:
                        d = t - sync_tt
                        if d in slots_validos:
                            tt_rel[trial]  = np.uint32(d)
                            setting[trial] = np.int8(
                                0 if d in slots_g1 else 1)
            del ch, tt, idx
            if parou: break

            pct = min(f.tell() / fsize * 100, 100)
            print(f"\r  [{label}] {trial+1:,} trials | {pct:.1f}%",
                  end='', flush=True)

    n_ok    = trial + 1
    n_click = int((tt_rel[:n_ok] != SENTINEL).sum())
    taxa    = n_click / n_ok * 100 if n_ok > 0 else 0
    n_s0    = int((setting[:n_ok] == 0).sum())
    n_s1    = int((setting[:n_ok] == 1).sum())

    print(f"\r  [{label}] {n_ok:,} trials | "
          f"{n_click:,} clicks ({taxa:.2f}%) | "
          f"s0={n_s0:,} s1={n_s1:,} | {time.time()-t0:.1f}s")

    # Verificar slots
    vals = tt_rel[:n_ok][tt_rel[:n_ok] != SENTINEL]
    if len(vals) > 0:
        unicos, cnts = np.unique(vals, return_counts=True)
        print(f"  tt_rel: {len(unicos)} slots únicos "
              f"range=[{vals.min()},{vals.max()}]")
        top = np.argsort(cnts)[-8:][::-1]
        for i in top:
            g = 'G1' if unicos[i] in slots_g1 else 'G2'
            print(f"    slot={unicos[i]:4d} ({g}): "
                  f"{cnts[i]:,} ({cnts[i]/n_click*100:.1f}%)")

    return tt_rel[:n_ok], setting[:n_ok], sync_abs[:n_ok]

# ══════════════════════════════════════════════════════════════════════════════
# ALINHAMENTO
# ══════════════════════════════════════════════════════════════════════════════
def alinhar(sy_a, sy_b, n=10000):
    n = min(n, len(sy_a), len(sy_b))
    sa, sb = sy_a[:n], sy_b[:n]
    clock_offset = int(np.median(sa - sb))
    residuo = sa - sb - clock_offset
    ok = np.abs(residuo) < 500
    std_res = float(np.std(residuo[ok])) if ok.sum() > 10 else 999.0
    print(f"  clock_offset={clock_offset:,} | std={std_res:.1f} bins ({std_res*BIN_PS:.0f} ps)")
    return clock_offset

# ══════════════════════════════════════════════════════════════════════════════
# CALCULAR Δt
# ══════════════════════════════════════════════════════════════════════════════
def calcular_dt_combos(tr_a, st_a, sy_a, tr_b, st_b, sy_b, clock_offset):
    n = min(len(tr_a), len(tr_b))
    dts = {}
    for nome, (sa_t, sb_t) in COMBOS.items():
        mask = ((tr_a[:n] != SENTINEL) & (tr_b[:n] != SENTINEL) &
                (st_a[:n] == sa_t) & (st_b[:n] == sb_t))
        if mask.sum() == 0:
            dts[nome] = np.array([], dtype=np.int64); continue
        dt = ((sy_a[:n][mask].astype(np.int64)
               - sy_b[:n][mask].astype(np.int64)
               - clock_offset)
              + tr_a[:n][mask].astype(np.int64)
              - tr_b[:n][mask].astype(np.int64))
        dts[nome] = dt
    return dts

# ══════════════════════════════════════════════════════════════════════════════
# HISTOGRAMAS
# ══════════════════════════════════════════════════════════════════════════════
def plotar_histogramas(dts, fname):
    fig, axes = plt.subplots(2, 4, figsize=(20, 8))
    fig.suptitle('Histogramas Δt — v7\n'
                 'Detectores: Alice ch=2/4 (slots 90-110), Bob ch=2/4 (slots 348-367)',
                 fontsize=12, fontweight='bold')
    for col, (nome, dt) in enumerate(dts.items()):
        N = len(dt)
        med = float(np.median(dt)) if N > 0 else 0
        std = float(np.std(dt))    if N > 0 else 0

        print(f"  {nome}: N={N:,}  med={med:.1f}  std={std:.2f} bins ({std*BIN_PS:.0f} ps)")

        ax1 = axes[0, col]
        if N > 0:
            ax1.hist(dt, bins=200, range=(med-200, med+200),
                     color='steelblue', alpha=0.8)
            ax1.axvline(med, color='red', lw=1.5, ls='--',
                        label=f'med={med:.0f}')
            ax1.legend(fontsize=7)
        ax1.set_title(f'{nome} — N={N:,}\nmed={med:.0f} bins')
        ax1.set_xlabel('Δt (bins)'); ax1.grid(alpha=0.3)

        ax2 = axes[1, col]
        if N > 0:
            dt_c = dt - int(med)
            ax2.hist(dt_c, bins=50, range=(-30, 30),
                     color='darkorange', alpha=0.8)
            ax2.axvline(0, color='red', lw=1.5, ls='--')
            ax2.set_title(f'{nome} zoom ±30 bins\nstd={std:.2f} bins = {std*BIN_PS:.0f} ps')
        ax2.set_xlabel('Δt − med (bins)'); ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.show()
    print(f"  Salvo: {fname}")

# ══════════════════════════════════════════════════════════════════════════════
# RAYLEIGH SPECTRUM
# ══════════════════════════════════════════════════════════════════════════════
def rayleigh(dt, k_max=890):
    N = len(dt)
    if N < 1000: return None
    dt_c = dt - int(np.median(dt))   # centralizar
    ks = np.arange(1, k_max + 1)
    Z = np.zeros(k_max); V = np.zeros(k_max); p = np.ones(k_max)
    for k in ks:
        phi = 2 * np.pi * (dt_c % k) / k
        C = np.mean(np.cos(phi)); S = np.mean(np.sin(phi))
        Z[k-1] = 2 * N * (C**2 + S**2)
        V[k-1] = np.sqrt(max(Z[k-1], 0) / (2 * N))
        p[k-1] = chi2.sf(Z[k-1], df=2)
    return {'ks': ks, 'Z': Z, 'V': V, 'p': p, 'N': N}

def plotar_espectro(resultados, fname):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Espectro Rayleigh — v7 (detectores confirmados)\n'
                 'Vermelho = n_prev esperado  |  Verde = ruído 1/√N',
                 fontsize=12, fontweight='bold')
    for idx, (nome, res) in enumerate(resultados.items()):
        np_val = N_PREV[nome]
        ax = axes[idx // 2, idx % 2]
        if res is None:
            ax.set_title(f'{nome} — dados insuficientes'); continue
        noise = 1.0 / np.sqrt(res['N'])
        ax.semilogy(res['ks'], res['V'], 'b-', lw=0.7, alpha=0.8)
        ax.axhline(noise, color='green', lw=1.5, ls='--',
                   label=f'1/√N={noise:.5f}')
        ax.axvline(np_val, color='red', lw=2, alpha=0.8,
                   label=f'n_prev={np_val}')
        v_np = res['V'][np_val - 1]
        sn   = v_np / noise
        ax.set_title(f'{nome}  N={res["N"]:,}  n_prev={np_val}\n'
                     f'V(n)={v_np:.5f}  S/N={sn:.2f}×')
        ax.set_xlabel('k'); ax.set_ylabel('V(k)')
        ax.legend(fontsize=8); ax.grid(alpha=0.3, which='both')
    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.show()
    print(f"  Salvo: {fname}")

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
print("\nExtraindo Alice...")
tr_a, st_a, sy_a = extrair(
    ALICE_FILE, SLOTS_A_G1, SLOTS_A_G2, N_TRIALS, 'Alice')

print("\nExtraindo Bob...")
tr_b, st_b, sy_b = extrair(
    BOB_FILE, SLOTS_B_G1, SLOTS_B_G2, N_TRIALS, 'Bob')

print("\nAlinhando syncs...")
clock_offset = alinhar(sy_a, sy_b)

print("\nCalculando Δt por combo...")
dts = calcular_dt_combos(tr_a, st_a, sy_a, tr_b, st_b, sy_b, clock_offset)

print("\n" + "="*65)
print("HISTOGRAMAS DE Δt")
print("="*65)
plotar_histogramas(dts, f'{OUT_DIR}/histogramas_v7.png')

print("\n" + "="*65)
print("ESPECTRO DE RAYLEIGH (k=1..890)")
print("="*65)
resultados_esp = {}
for nome, dt in dts.items():
    np_val = N_PREV[nome]
    N = len(dt)
    print(f"\n  {nome}: N={N:,}  n_prev={np_val}")
    res = rayleigh(dt, k_max=890)
    resultados_esp[nome] = res
    if res:
        noise = 1.0 / np.sqrt(res['N'])
        for k in [np_val, np_val*2, np_val*3]:
            if 1 <= k <= 890:
                v = res['V'][k-1]
                print(f"    k={k:4d}  V={v:.6f}  S/N={v/noise:.2f}×  "
                      f"Z={res['Z'][k-1]:.1f}  p={res['p'][k-1]:.2e}")

plotar_espectro(resultados_esp, f'{OUT_DIR}/espectro_v7.png')

# CSV
with open(f'{OUT_DIR}/rayleigh_v7.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['combo','k','V','Z','p','N'])
    for nome, res in resultados_esp.items():
        if res is None: continue
        for i, k in enumerate(res['ks']):
            w.writerow([nome, k,
                        f"{res['V'][i]:.8f}",
                        f"{res['Z'][i]:.2f}",
                        f"{res['p'][i]:.6e}",
                        res['N']])
print(f"\nCSV: {OUT_DIR}/rayleigh_v7.csv")

print("\n" + "="*65)
print("RESUMO FINAL")
print("="*65)
for nome, res in resultados_esp.items():
    if res is None:
        print(f"  {nome}: insuficiente"); continue
    np_val = N_PREV[nome]
    noise  = 1.0 / np.sqrt(res['N'])
    v_np   = res['V'][np_val - 1]
    sn     = v_np / noise
    flag   = '✓ ACIMA DO RUÍDO' if sn > 3 else '✗ no nível do ruído'
    print(f"  {nome}: N={res['N']:,}  V({np_val})={v_np:.5f}  "
          f"S/N={sn:.2f}×  {flag}")