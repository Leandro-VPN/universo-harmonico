# ══════════════════════════════════════════════════════════════════════
# ║  COMPARAÇÃO QRNG vs XOR — via HDF5 (laserPulseNumber + phase)       ║
# ║  XOR HDF5 sem phase/LPN: comparação inconclusiva por falta de dados ║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════

COMPARAÇÃO QRNG vs XOR — via HDF5 (CORRIGIDO)
===============================================
Corrigido: laserPulseNumber pode não existir no HDF5 XOR.
Usa clicks diretamente para calcular Δt via phase → bins.
"""

import numpy as np, h5py, os, time
import matplotlib.pyplot as plt, csv
from scipy import stats

HDF5_QRNG = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.dat.compressed.build.hdf5"
HDF5_XOR  = r"d:\BellQM\23_55_CH_pockel_100kHz.run.ClassicalRNGXOR.dat.compressed.build.hdf5"
OUT       = r"d:\BellQM\out_qrng_vs_xor_final"
os.makedirs(OUT, exist_ok=True)

BIN_PS  = 78.125
COMBOS  = {'XX':(1,1),'XY':(1,2),'YX':(2,1),'YY':(2,2)}
N_PREV  = {'XX':178,'XY':15,'YX':15,'YY':5}
N_PERM  = 1000
LASER_P = 129104
G_PER_B = 360.0 / LASER_P

rng = np.random.default_rng(42)

# ── Inspecionar estrutura do HDF5 ────────────────────────────────────────────
def inspecionar_hdf5(path, label):
    print(f"\n{'='*55}")
    print(f"ESTRUTURA: {label}")
    print(f"{'='*55}")
    def visitar(nome, obj):
        if isinstance(obj, h5py.Dataset):
            print(f"  {nome}: shape={obj.shape} dtype={obj.dtype}")
    with h5py.File(path, 'r') as f:
        f.visititems(visitar)

inspecionar_hdf5(HDF5_QRNG, 'QRNG')
inspecionar_hdf5(HDF5_XOR,  'XOR')

# ── Carregar e calcular Δt de um HDF5 ────────────────────────────────────────
def dts_from_hdf5(path, label):
    print(f"\n{'='*55}")
    print(f"Carregando {label}...")
    with h5py.File(path, 'r') as f:
        # campos obrigatórios
        a_set = f['alice/settings'][:]
        b_set = f['bob/settings'][:]
        a_clk = f['alice/clicks'][:]
        b_clk = f['bob/clicks'][:]

        # phase — pode não existir; usa zeros se ausente
        if 'alice/phase' in f and 'bob/phase' in f:
            a_ph_raw = f['alice/phase'][:]
            b_ph_raw = f['bob/phase'][:]
            tem_phase = True
        else:
            tem_phase = False
            print(f"  ⚠️  phase ausente em {label} — usando laserPulseNumber se disponível")

        # laserPulseNumber — opcional
        if 'alice/laserPulseNumber' in f and 'bob/laserPulseNumber' in f:
            a_lpn = f['alice/laserPulseNumber'][:]
            b_lpn = f['bob/laserPulseNumber'][:]
            tem_lpn = True
        else:
            tem_lpn = False

    N = len(a_set)
    print(f"  Trials: {N:,}")
    print(f"  Settings únicos Alice: {np.unique(a_set)}")
    print(f"  Settings únicos Bob:   {np.unique(b_set)}")
    print(f"  tem_phase={tem_phase}  tem_lpn={tem_lpn}")

    # Coincidências
    coinc  = (a_clk > 0) & (b_clk > 0)
    N_c    = coinc.sum()
    print(f"  Coincidências: {N_c:,} ({N_c/N*100:.4f}%)")

    if N_c == 0:
        print(f"  ⚠️  Sem coincidências em {label}")
        return {nome: np.array([], dtype=np.int64) for nome in COMBOS}

    # Índices dos trials com coincidência
    c_idx  = np.flatnonzero(coinc)
    sa_all = a_set[c_idx]
    sb_all = b_set[c_idx]

    # Calcular Δt
    if tem_phase:
        # Converte phase (graus) → bins de laser
        a_trials_idx = np.flatnonzero(a_clk > 0)
        b_trials_idx = np.flatnonzero(b_clk > 0)

        # Para cada coincidência, acha o índice no array de phase
        ap = np.searchsorted(a_trials_idx, c_idx)
        bp = np.searchsorted(b_trials_idx, c_idx)

        # Garante que o mapeamento é válido
        ap = np.clip(ap, 0, len(a_ph_raw) - 1)
        bp = np.clip(bp, 0, len(b_ph_raw) - 1)

        # Verifica correspondência exata
        valid = ((ap < len(a_ph_raw)) & (bp < len(b_ph_raw)))
        ap = ap[valid]; bp = bp[valid]
        sa_all = sa_all[valid]; sb_all = sb_all[valid]

        a_bins = np.round(a_ph_raw[ap] / G_PER_B).astype(np.int64)
        b_bins = np.round(b_ph_raw[bp] / G_PER_B).astype(np.int64)
        dt_all = a_bins - b_bins

    elif tem_lpn:
        # Fallback: usa laserPulseNumber como proxy de tempo
        a_trials_idx = np.flatnonzero(a_clk > 0)
        b_trials_idx = np.flatnonzero(b_clk > 0)
        ap = np.searchsorted(a_trials_idx, c_idx)
        bp = np.searchsorted(b_trials_idx, c_idx)
        ap = np.clip(ap, 0, len(a_lpn) - 1)
        bp = np.clip(bp, 0, len(b_lpn) - 1)
        valid = ((ap < len(a_lpn)) & (bp < len(b_lpn)))
        ap = ap[valid]; bp = bp[valid]
        sa_all = sa_all[valid]; sb_all = sb_all[valid]
        dt_all = a_lpn[ap].astype(np.int64) - b_lpn[bp].astype(np.int64)
        print(f"  (usando laserPulseNumber como Δt)")

    else:
        # Sem phase nem LPN: Δt = 0 (só conta coincidências)
        print(f"  ⚠️  Sem phase nem LPN — Δt = índice de trial")
        dt_all = np.zeros(len(sa_all), dtype=np.int64)

    # Separar por combo
    dts = {}
    for nome, (sa_t, sb_t) in COMBOS.items():
        mask = (sa_all == sa_t) & (sb_all == sb_t)
        dt   = dt_all[mask]
        dts[nome] = dt
        if len(dt) > 0:
            print(f"  {nome}: N={len(dt):,}  "
                  f"med={np.median(dt):.1f}  std={np.std(dt):.2f} bins")
        else:
            print(f"  {nome}: sem dados")

    return dts

dts_qrng = dts_from_hdf5(HDF5_QRNG, 'QRNG')
dts_xor  = dts_from_hdf5(HDF5_XOR,  'XOR')

# ── ACF ──────────────────────────────────────────────────────────────────────
def acf_lag(x, lag):
    x = x.astype(np.float64)
    mu = np.mean(x); var = np.var(x)
    if var == 0 or len(x) <= lag: return 0.0
    N = len(x)
    return float(np.mean((x[:N-lag] - mu) * (x[lag:] - mu)) / var)

def acf_full(x, max_lag=300):
    x = x.astype(np.float64)
    mu = np.mean(x); var = np.var(x)
    if var == 0: return np.zeros(max_lag)
    N = len(x)
    return np.array([
        np.mean((x[:N-l] - mu) * (x[l:] - mu)) / var
        for l in range(1, max_lag + 1)
    ])

# ── Comparação ────────────────────────────────────────────────────────────────
print("\n" + "="*65)
print("COMPARAÇÃO ACF — QRNG vs XOR")
print("="*65)

fig, axes = plt.subplots(3, 4, figsize=(22, 14))
fig.suptitle(
    'QRNG vs XOR Clássico — ACF dos Δt (via HDF5)\n'
    'Azul=QRNG  Vermelho=XOR  '
    'Idênticos → hardware  |  Diferentes → física quântica',
    fontsize=12, fontweight='bold')

resultados = []

for col, nome in enumerate(COMBOS):
    np_val = N_PREV[nome]
    dt_q   = dts_qrng.get(nome, np.array([]))
    dt_x   = dts_xor.get(nome,  np.array([]))
    N_ac   = min(len(dt_q), len(dt_x), 100_000)

    if N_ac < 100:
        print(f"  {nome}: dados insuficientes "
              f"(QRNG={len(dt_q)}, XOR={len(dt_x)})")
        for ax in axes[:, col]:
            ax.set_title(f'{nome} — sem dados suficientes')
        continue

    acf_q = acf_lag(dt_q[:N_ac], np_val)
    acf_x = acf_lag(dt_x[:N_ac], np_val)

    # Permutações
    print(f"\n  {nome} — {N_PERM} permutações...")
    t0 = time.time()
    null_q = np.array([
        acf_lag(rng.permutation(dt_q[:N_ac]), np_val)
        for _ in range(N_PERM)
    ])
    null_x = np.array([
        acf_lag(rng.permutation(dt_x[:N_ac]), np_val)
        for _ in range(N_PERM)
    ])
    print(f"    ({time.time()-t0:.1f}s)")

    mu_q = null_q.mean(); std_q = null_q.std()
    mu_x = null_x.mean(); std_x = null_x.std()
    z_q  = (acf_q - mu_q) / std_q if std_q > 0 else 0
    z_x  = (acf_x - mu_x) / std_x if std_x > 0 else 0
    std_comb = np.sqrt(std_q**2 + std_x**2)
    z_diff   = (acf_q - acf_x) / std_comb if std_comb > 0 else 0
    ic95     = 2 / np.sqrt(N_ac)

    print(f"    QRNG: ACF={acf_q:.5f}  z={z_q:.1f}σ")
    print(f"    XOR:  ACF={acf_x:.5f}  z={z_x:.1f}σ")
    print(f"    Δ={acf_q-acf_x:+.5f}  z_diff={z_diff:.2f}σ")

    resultados.append({
        'combo': nome, 'np_val': np_val,
        'acf_q': acf_q, 'acf_x': acf_x,
        'z_q': z_q, 'z_x': z_x, 'z_diff': z_diff,
        'delta': acf_q - acf_x, 'N_ac': N_ac
    })

    # ACF completa
    acf_q_f = acf_full(dt_q[:N_ac])
    acf_x_f = acf_full(dt_x[:N_ac])
    lags = np.arange(1, 301)

    # Row 0: ACF completa
    ax = axes[0, col]
    ax.plot(lags, acf_q_f, 'b-', lw=1.0, label='QRNG', alpha=0.9)
    ax.plot(lags, acf_x_f, 'r-', lw=1.0, label='XOR',  alpha=0.8)
    ax.axhline(0,     color='k',     lw=0.5)
    ax.axhline( ic95, color='green', lw=1, ls='--', label='IC95')
    ax.axhline(-ic95, color='green', lw=1, ls='--')
    ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7,
               label=f'n={np_val}')
    ax.set_title(
        f'{nome} — ACF(Δt)\n'
        f'QRNG({np_val})={acf_q:.4f}  XOR={acf_x:.4f}  '
        f'Δ={acf_q-acf_x:+.4f}', fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_xlabel('lag')

    # Row 1: distribuições nulas
    ax = axes[1, col]
    bins = np.linspace(
        min(null_q.min(), null_x.min()),
        max(null_q.max(), null_x.max()), 50)
    ax.hist(null_q, bins=bins, color='blue', alpha=0.5, density=True,
            label=f'nula QRNG\n±{std_q:.4f}')
    ax.hist(null_x, bins=bins, color='red',  alpha=0.5, density=True,
            label=f'nula XOR\n±{std_x:.4f}')
    ax.axvline(acf_q, color='blue', lw=3, label=f'QRNG={acf_q:.4f}')
    ax.axvline(acf_x, color='red',  lw=3, label=f'XOR={acf_x:.4f}')
    ax.set_title(
        f'{nome} — Distribuições nulas\n'
        f'z_diff={z_diff:.2f}σ', fontsize=9)
    ax.legend(fontsize=6); ax.grid(alpha=0.3); ax.set_xlabel('ACF')

    # Row 2: diferença QRNG − XOR
    ax = axes[2, col]
    diff = acf_q_f - acf_x_f
    ax.plot(lags, diff, 'purple', lw=1.0, label='QRNG − XOR')
    ax.axhline(0, color='k', lw=0.5)
    ax.axhline( ic95 * np.sqrt(2), color='green', lw=1, ls='--',
                label='IC95 diff')
    ax.axhline(-ic95 * np.sqrt(2), color='green', lw=1, ls='--')
    ax.axvline(np_val, color='red', lw=1.5, alpha=0.7)
    ax.set_title(
        f'{nome} — ACF(QRNG)−ACF(XOR)\n'
        f'diff(lag={np_val})={diff[np_val-1]:+.5f}', fontsize=9)
    ax.legend(fontsize=6); ax.grid(alpha=0.3); ax.set_xlabel('lag')

plt.tight_layout()
saida = f'{OUT}/qrng_vs_xor_final.png'
plt.savefig(saida, dpi=150, bbox_inches='tight')
plt.show()
print(f"\nGráfico salvo: {saida}")

# ── CSV ───────────────────────────────────────────────────────────────────────
with open(f'{OUT}/resultado_final.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['combo', 'np_val', 'N_ac', 'acf_qrng', 'acf_xor',
                'delta', 'z_diff'])
    for r in resultados:
        w.writerow([r['combo'], r['np_val'], r['N_ac'],
                    f"{r['acf_q']:.6f}", f"{r['acf_x']:.6f}",
                    f"{r['delta']:+.6f}", f"{r['z_diff']:.3f}"])
print(f"CSV: {OUT}/resultado_final.csv")

# ── Veredito ──────────────────────────────────────────────────────────────────
print("\n" + "="*65)
print("VEREDICTO FINAL")
print("="*65)

if not resultados:
    print("\n⚠️  Nenhum combo teve dados suficientes para comparação.")
    print("   Verifique se o HDF5 XOR tem o mesmo formato do QRNG.")
else:
    for r in resultados:
        z = r['z_diff']
        if abs(z) < 2:
            v = '✗ IDÊNTICOS — origem no hardware (Pockels cell)'
        elif abs(z) < 5:
            v = '~ MARGINAL — investigar mais'
        else:
            v = '✓ DIFERENTES — sinal ligado ao tipo de RNG'
        print(f"  {r['combo']}: QRNG={r['acf_q']:.5f}  XOR={r['acf_x']:.5f}  "