# ══════════════════════════════════════════════════════════════════════
# ║  ANÁLISE — ΔΔt, FFT SÉRIE TEMPORAL E AUTOCORRELAÇÃO                 ║
# ║  Testes que dependem da ordem: std(ΔΔt), FFT, ACF vs permutado      ║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════


import numpy as np, os, time
from scipy.stats import chi2
from scipy.signal import welch
import matplotlib.pyplot as plt, csv

ALICE_FILE = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat"
BOB_FILE   = r"d:\BellQM\19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat"
OUT        = r"d:\BellQM\out_ddelta"
os.makedirs(OUT, exist_ok=True)

BIN_PS   = 78.125
SENTINEL = np.uint32(0xFFFFFFFF)
COMBOS   = {'XX':(0,0),'XY':(0,1),'YX':(1,0),'YY':(1,1)}
N_PREV   = {'XX':178,'XY':15,'YX':15,'YY':5}

SLOTS_A_G1 = set(range(90, 97));  SLOTS_A_G2 = set(range(104, 111))
SLOTS_B_G1 = set(range(347, 355)); SLOTS_B_G2 = set(range(361, 368))
SLOTS_A    = SLOTS_A_G1 | SLOTS_A_G2
SLOTS_B    = SLOTS_B_G1 | SLOTS_B_G2

rng = np.random.default_rng(42)

# ── Extrair Δt ────────────────────────────────────────────────────────────────
def extrair(filepath, slots_g1, slots_g2, slots_all, n_max, label):
    tt_rel  = np.full(n_max, SENTINEL, dtype=np.uint32)
    setting = np.full(n_max, -1,       dtype=np.int8)
    sync_abs= np.full(n_max, -1,       dtype=np.int64)
    trial=-1; sync_tt=None
    with open(filepath,'rb') as f:
        while trial < n_max-1:
            raw = f.read(5_000_000*24)
            if not raw: break
            n = len(raw)//24
            arr = np.frombuffer(raw,dtype=np.uint64).reshape(n,3)
            ch = arr[:,0].astype(np.int32); tt = arr[:,1].astype(np.int64)
            del arr,raw
            idx = np.flatnonzero((ch==6)|(ch==2)|(ch==4))
            parou=False
            for i in idx:
                c=int(ch[i]); t=int(tt[i])
                if c==6:
                    trial+=1; sync_tt=t
                    if trial<n_max: sync_abs[trial]=t
                    if trial>=n_max-1: parou=True; break
                elif sync_tt is not None and trial>=0 and trial<n_max:
                    if tt_rel[trial]==SENTINEL:
                        d=t-sync_tt
                        if d in slots_all:
                            tt_rel[trial]=np.uint32(d)
                            setting[trial]=np.int8(0 if d in slots_g1 else 1)
            del ch,tt,idx
            if parou: break
    return tt_rel[:trial+1], setting[:trial+1], sync_abs[:trial+1]

print("Extraindo dados...")
t0=time.time()
tr_a,st_a,sy_a = extrair(ALICE_FILE,SLOTS_A_G1,SLOTS_A_G2,SLOTS_A,10_000_000,'Alice')
tr_b,st_b,sy_b = extrair(BOB_FILE,  SLOTS_B_G1,SLOTS_B_G2,SLOTS_B,10_000_000,'Bob')
n = min(len(tr_a),len(tr_b))
clock_offset = int(np.median(sy_a[:5000].astype(np.int64)-sy_b[:5000].astype(np.int64)))
print(f"  clock_offset={clock_offset:,}  ({time.time()-t0:.1f}s)")

dts_orig = {}
for nome,(sa_t,sb_t) in COMBOS.items():
    mask = ((tr_a[:n]!=SENTINEL)&(tr_b[:n]!=SENTINEL)&
            (st_a[:n]==sa_t)&(st_b[:n]==sb_t))
    dt = ((sy_a[:n][mask].astype(np.int64)-sy_b[:n][mask].astype(np.int64)-clock_offset)
          + tr_a[:n][mask].astype(np.int64)-tr_b[:n][mask].astype(np.int64))
    dts_orig[nome] = dt
    print(f"  {nome}: N={len(dt):,}")

# ── Funções de análise ────────────────────────────────────────────────────────
def rayleigh_V(dt, k_max=200):
    N=len(dt); dt_c=dt-int(np.median(dt))
    ks=np.arange(1,k_max+1); V=np.zeros(k_max)
    for k in ks:
        phi=2*np.pi*(dt_c%k)/k
        C=np.mean(np.cos(phi)); S=np.mean(np.sin(phi))
        V[k-1]=np.sqrt(max(2*N*(C**2+S**2),0)/(2*N))
    return ks, V

def fft_serie(dt, max_k=500):
    """FFT da série Δt(i) vs i — detecta periodicidade na sequência."""
    x = dt.astype(np.float64) - np.mean(dt)
    N = len(x)
    # Janela de Hann
    win = np.hanning(N)
    fft = np.fft.rfft(x * win)
    freqs = np.fft.rfftfreq(N)
    pot = np.abs(fft)**2
    # S/N vs fundo mediano
    med_pot = np.median(pot[5:])
    sn = pot / med_pot
    periodos = 1.0 / freqs[1:]
    sn_p = sn[1:]
    # Só períodos de 2 a max_k
    mask = periodos <= max_k
    return periodos[mask], sn_p[mask]

def autocorr(dt, max_lag=300):
    """Autocorrelação normalizada da série Δt(i)."""
    x = dt.astype(np.float64) - np.mean(dt)
    N = len(x)
    result = np.zeros(max_lag)
    var = np.var(x)
    if var == 0: return np.arange(1,max_lag+1), result
    for lag in range(1, max_lag+1):
        result[lag-1] = np.mean(x[:N-lag]*x[lag:]) / var
    return np.arange(1, max_lag+1), result

# ── Análise por combo ─────────────────────────────────────────────────────────
print("\n"+"="*65)
print("ANÁLISE ΔΔt, FFT E AUTOCORRELAÇÃO")
print("="*65)

fig, axes = plt.subplots(len(COMBOS), 4, figsize=(22, 5*len(COMBOS)))
fig.suptitle('ΔΔt, FFT(série) e Autocorrelação — dependem da ORDEM\n'
             'Azul=original  Vermelho=permutado  '
             'Diferença entre eles = estrutura temporal real',
             fontsize=11, fontweight='bold')

for row, (nome, dt) in enumerate(dts_orig.items()):
    np_val = N_PREV[nome]
    N = len(dt)
    noise_r = 1/np.sqrt(N)

    dt_perm = rng.permutation(dt)

    # ── Col 0: ΔΔt distribuição ───────────────────────────────────────────
    ddt      = np.diff(dt.astype(np.int64))
    ddt_perm = np.diff(dt_perm.astype(np.int64))

    ax = axes[row, 0]
    bins = np.linspace(-200, 200, 100)
    ax.hist(ddt,      bins=bins, color='blue',   alpha=0.6, label='original', density=True)
    ax.hist(ddt_perm, bins=bins, color='red',    alpha=0.6, label='permutado', density=True)
    ax.set_title(f'{nome} — ΔΔt distribuição\n'
                 f'std_orig={np.std(ddt):.1f}  std_perm={np.std(ddt_perm):.1f}',
                 fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    ax.set_xlabel('ΔΔt (bins)')

    std_ddt_o = float(np.std(ddt))
    std_ddt_p = float(np.std(ddt_perm))
    print(f"\n  {nome}: N={N:,}")
    print(f"    ΔΔt std: original={std_ddt_o:.2f}  permutado={std_ddt_p:.2f}  "
          f"ratio={std_ddt_o/std_ddt_p:.3f}")

    # ── Col 1: Rayleigh do ΔΔt ───────────────────────────────────────────
    ks_r, V_ddt      = rayleigh_V(ddt)
    ks_r, V_ddt_perm = rayleigh_V(ddt_perm)
    noise_ddt = 1/np.sqrt(len(ddt))

    ax = axes[row, 1]
    ax.semilogy(ks_r, V_ddt,      'b-', lw=0.8, alpha=0.8, label='original')
    ax.semilogy(ks_r, V_ddt_perm, 'r-', lw=0.8, alpha=0.5, label='permutado')
    ax.axhline(noise_ddt, color='black', lw=1, ls='--', label='1/√N')
    ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7, label=f'n={np_val}')
    v_o = V_ddt[np_val-1]; v_p = V_ddt_perm[np_val-1]
    ax.set_title(f'{nome} — Rayleigh(ΔΔt)\n'
                 f'V({np_val}): orig={v_o:.5f} perm={v_p:.5f} '
                 f'ratio={v_o/max(v_p,1e-9):.2f}×',
                 fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3, which='both')
    print(f"    Rayleigh(ΔΔt) V({np_val}): orig={v_o:.5f}  perm={v_p:.5f}  "
          f"ratio={v_o/max(v_p,1e-9):.2f}×")

    # ── Col 2: FFT da série temporal ──────────────────────────────────────
    # Usar subconjunto de 200k para velocidade
    N_fft = min(N, 200_000)
    per_o, sn_o = fft_serie(dt[:N_fft])
    per_p, sn_p = fft_serie(dt_perm[:N_fft])

    ax = axes[row, 2]
    ax.semilogy(per_o, sn_o, 'b-', lw=0.5, alpha=0.7, label='original')
    ax.semilogy(per_p, sn_p, 'r-', lw=0.5, alpha=0.5, label='permutado')
    ax.axhline(3, color='green', lw=1, ls='--', label='S/N=3×')
    ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7, label=f'n={np_val}')
    # pico mais alto no original
    k_top_o = float(per_o[np.argmax(sn_o)])
    sn_top_o = float(sn_o.max())
    k_top_p = float(per_p[np.argmax(sn_p)])
    ax.set_title(f'{nome} — FFT(série) N={N_fft:,}\n'
                 f'orig top: k={k_top_o:.0f} S/N={sn_top_o:.1f}×  '
                 f'perm top: k={k_top_p:.0f}',
                 fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3, which='both')
    ax.set_xlim(1, 500)
    print(f"    FFT(série): orig_top k={k_top_o:.0f} S/N={sn_top_o:.1f}×  "
          f"perm_top k={k_top_p:.0f}")

    # ── Col 3: Autocorrelação ─────────────────────────────────────────────
    lags_o, ac_o = autocorr(dt[:min(N, 50_000)])
    lags_p, ac_p = autocorr(dt_perm[:min(N, 50_000)])
    ic95 = 2/np.sqrt(min(N, 50_000))

    ax = axes[row, 3]
    ax.plot(lags_o, ac_o, 'b-', lw=0.7, alpha=0.8, label='original')
    ax.plot(lags_p, ac_p, 'r-', lw=0.7, alpha=0.5, label='permutado')
    ax.axhline( ic95, color='green', lw=1, ls='--', label=f'IC95 ±{ic95:.4f}')
    ax.axhline(-ic95, color='green', lw=1, ls='--')
    ax.axhline(0, color='black', lw=0.5)
    ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7, label=f'n={np_val}')
    ac_at_nprev = ac_o[np_val-1]
    ax.set_title(f'{nome} — Autocorrelação\n'
                 f'ACF({np_val})={ac_at_nprev:.5f}  IC95=±{ic95:.4f}',
                 fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    print(f"    ACF({np_val})={ac_at_nprev:.5f}  IC95=±{ic95:.4f}  "
          f"{'FORA DO IC' if abs(ac_at_nprev)>ic95 else 'dentro do IC'}")

plt.tight_layout()
saida = f'{OUT}/ddelta_fft_acf.png'
plt.savefig(saida, dpi=130, bbox_inches='tight')
plt.show()
print(f"\nGráfico salvo: {saida}")

# ── CSV com resultados numéricos ──────────────────────────────────────────────
with open(f'{OUT}/ddelta_resultados.csv','w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['combo','N','std_ddt_orig','std_ddt_perm',
                'V_ddt_nprev_orig','V_ddt_nprev_perm',
                'fft_top_k','fft_top_sn','acf_nprev','ic95'])
    for nome, dt in dts_orig.items():
        np_val = N_PREV[nome]; N = len(dt)
        ddt = np.diff(dt.astype(np.int64))
        ddt_p = np.diff(rng.permutation(dt).astype(np.int64))
        _,V_o = rayleigh_V(ddt); _,V_p = rayleigh_V(ddt_p)
        N_fft = min(N, 200_000)
        per_o, sn_o = fft_serie(dt[:N_fft])
        lags, ac = autocorr(dt[:min(N,50_000)])
        w.writerow([nome, N,
                    f"{np.std(ddt):.4f}", f"{np.std(ddt_p):.4f}",
                    f"{V_o[np_val-1]:.6f}", f"{V_p[np_val-1]:.6f}",
                    f"{float(per_o[np.argmax(sn_o)]):.0f}",
                    f"{float(sn_o.max()):.2f}",
                    f"{ac[np_val-1]:.6f}",
                    f"{2/np.sqrt(min(N,50_000)):.6f}"])

print(f"CSV: {OUT}/ddelta_resultados.csv")

print("\n"+"="*65)
print("VEREDICTO FINAL")
print("="*65)
print("""
Se ΔΔt std(orig) ≠ std(perm) → a sequência temporal dos Δt tem estrutura.
Se FFT(série) mostra pico que some na permutação → periodicidade real.
Se ACF(n_prev) fora do IC95 e some na permutação → autocorrelação real.
""")
