# ══════════════════════════════════════════════════════════════════════
# ║  CONTROLES AVANÇADOS — Surrogate Fourier e Circular Shift           ║
# ║  ACF(surr)≈0.97 confirma estrutura no espectro de potência (hardware)║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════


import numpy as np, os, time
from scipy.stats import chi2
import matplotlib.pyplot as plt, csv

ALICE_FILE = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat"
BOB_FILE   = r"d:\BellQM\19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat"
OUT        = r"d:\BellQM\out_controles2"
os.makedirs(OUT, exist_ok=True)

BIN_PS   = 78.125
SENTINEL = np.uint32(0xFFFFFFFF)
COMBOS   = {'XX':(0,0),'XY':(0,1),'YX':(1,0),'YY':(1,1)}
N_PREV   = {'XX':178,'XY':15,'YX':15,'YY':5}

SLOTS_A_G1 = set(range(90,97));  SLOTS_A_G2 = set(range(104,111))
SLOTS_B_G1 = set(range(347,355)); SLOTS_B_G2 = set(range(361,368))
SLOTS_A    = SLOTS_A_G1|SLOTS_A_G2; SLOTS_B = SLOTS_B_G1|SLOTS_B_G2

rng = np.random.default_rng(42)

# ── Extrair Δt ────────────────────────────────────────────────────────────────
def extrair(filepath, sg1, sg2, sa, n_max, label):
    tt_rel=np.full(n_max,SENTINEL,dtype=np.uint32)
    setting=np.full(n_max,-1,dtype=np.int8)
    sync_abs=np.full(n_max,-1,dtype=np.int64)
    trial=-1; sync_tt=None
    with open(filepath,'rb') as f:
        while trial<n_max-1:
            raw=f.read(5_000_000*24)
            if not raw: break
            n=len(raw)//24
            arr=np.frombuffer(raw,dtype=np.uint64).reshape(n,3)
            ch=arr[:,0].astype(np.int32); tt=arr[:,1].astype(np.int64)
            del arr,raw
            idx=np.flatnonzero((ch==6)|(ch==2)|(ch==4))
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
                        if d in sa:
                            tt_rel[trial]=np.uint32(d)
                            setting[trial]=np.int8(0 if d in sg1 else 1)
            del ch,tt,idx
            if parou: break
    return tt_rel[:trial+1],setting[:trial+1],sync_abs[:trial+1]

print("Extraindo dados...")
t0=time.time()
tr_a,st_a,sy_a=extrair(ALICE_FILE,SLOTS_A_G1,SLOTS_A_G2,SLOTS_A,10_000_000,'Alice')
tr_b,st_b,sy_b=extrair(BOB_FILE,  SLOTS_B_G1,SLOTS_B_G2,SLOTS_B,10_000_000,'Bob')
n=min(len(tr_a),len(tr_b))
clock_offset=int(np.median(sy_a[:5000].astype(np.int64)-sy_b[:5000].astype(np.int64)))

dts_orig={}
for nome,(sa_t,sb_t) in COMBOS.items():
    mask=((tr_a[:n]!=SENTINEL)&(tr_b[:n]!=SENTINEL)&
          (st_a[:n]==sa_t)&(st_b[:n]==sb_t))
    dt=((sy_a[:n][mask].astype(np.int64)-sy_b[:n][mask].astype(np.int64)-clock_offset)
        +tr_a[:n][mask].astype(np.int64)-tr_b[:n][mask].astype(np.int64))
    dts_orig[nome]=dt
print(f"  Extraído ({time.time()-t0:.1f}s)")

# ── Funções ───────────────────────────────────────────────────────────────────
def fourier_surrogate(x, rng):
    """Preserva PSD, aleatoriza fase."""
    N = len(x)
    ft = np.fft.rfft(x.astype(np.float64))
    phases = rng.uniform(0, 2*np.pi, len(ft))
    ft_rand = np.abs(ft) * np.exp(1j * phases)
    surr = np.fft.irfft(ft_rand, n=N)
    # Converter de volta para inteiros com a mesma distribuição de ranks
    ranks = np.argsort(np.argsort(surr))
    sorted_x = np.sort(x)
    return sorted_x[ranks]

def acf(x, max_lag=300):
    x = x.astype(np.float64) - np.mean(x)
    var = np.var(x)
    if var == 0: return np.arange(1,max_lag+1), np.zeros(max_lag)
    N = len(x)
    result = np.array([np.mean(x[:N-lag]*x[lag:])/var
                       for lag in range(1, max_lag+1)])
    return np.arange(1, max_lag+1), result

def rayleigh_V(dt, k_max=200):
    N=len(dt); dt_c=dt-int(np.median(dt))
    ks=np.arange(1,k_max+1); V=np.zeros(k_max)
    for k in ks:
        phi=2*np.pi*(dt_c%k)/k
        C=np.mean(np.cos(phi)); S=np.mean(np.sin(phi))
        V[k-1]=np.sqrt(max(2*N*(C**2+S**2),0)/(2*N))
    return ks,V

# ── Análise ───────────────────────────────────────────────────────────────────
print("\n"+"="*65)
print("CONTROLES AVANÇADOS")
print("="*65)

fig, axes = plt.subplots(len(COMBOS), 5, figsize=(26, 5*len(COMBOS)))
fig.suptitle('Controles avançados: Fourier surrogate, Circular shift, Δ³t\n'
             'Azul=original  Vermelho=permutado  '
             'Verde=surrogate  Laranja=circ.shift  Roxo=Δ³t',
             fontsize=10, fontweight='bold')

resultados = []

for row, (nome, dt) in enumerate(dts_orig.items()):
    np_val = N_PREV[nome]
    N = len(dt)
    ic95 = 2/np.sqrt(min(N, 50_000))
    noise = 1/np.sqrt(N)

    print(f"\n  {nome}: N={N:,}  n_prev={np_val}")

    # Gerar controles
    dt_perm  = rng.permutation(dt)
    dt_surr  = fourier_surrogate(dt, rng)
    dt_shift = np.roll(dt, N//4)   # shift de N/4

    # ΔΔt de cada um
    ddt      = np.diff(dt.astype(np.int64))
    ddt_perm = np.diff(dt_perm.astype(np.int64))
    ddt_surr = np.diff(dt_surr.astype(np.int64))
    ddt_shft = np.diff(dt_shift.astype(np.int64))

    # Δ³t
    d3t      = np.diff(ddt)
    d3t_perm = np.diff(ddt_perm)

    print(f"    ΔΔt std: orig={np.std(ddt):.2f}  perm={np.std(ddt_perm):.2f}  "
          f"surr={np.std(ddt_surr):.2f}  shift={np.std(ddt_shft):.2f}")

    # ACF do ΔΔt original
    lags, ac_ddt      = acf(ddt[:50_000])
    lags, ac_ddt_perm = acf(ddt_perm[:50_000])
    lags, ac_ddt_surr = acf(ddt_surr[:50_000])
    lags, ac_ddt_shft = acf(ddt_shft[:50_000])

    # ACF do Δt original
    N_ac = min(N, 50_000)
    lags_dt, ac_dt      = acf(dt[:N_ac])
    lags_dt, ac_dt_perm = acf(dt_perm[:N_ac])
    lags_dt, ac_dt_surr = acf(dt_surr[:N_ac])
    lags_dt, ac_dt_shft = acf(dt_shift[:N_ac])

    ac_np_orig  = ac_dt[np_val-1]
    ac_np_perm  = ac_dt_perm[np_val-1]
    ac_np_surr  = ac_dt_surr[np_val-1]
    ac_np_shft  = ac_dt_shft[np_val-1]

    print(f"    ACF(Δt, n_prev): orig={ac_np_orig:.5f}  perm={ac_np_perm:.5f}  "
          f"surr={ac_np_surr:.5f}  shift={ac_np_shft:.5f}  IC95=±{ic95:.5f}")

    # Rayleigh(ΔΔt)
    ks, V_ddt_o = rayleigh_V(ddt)
    ks, V_ddt_p = rayleigh_V(ddt_perm)
    ks, V_ddt_s = rayleigh_V(ddt_surr)

    v_o=V_ddt_o[np_val-1]; v_p=V_ddt_p[np_val-1]; v_s=V_ddt_s[np_val-1]
    print(f"    Rayleigh(ΔΔt, n_prev): orig={v_o:.5f}  perm={v_p:.5f}  surr={v_s:.5f}")

    # Δ³t ACF
    lags3, ac_d3 = acf(d3t[:50_000])
    ac_d3_np = ac_d3[np_val-1]
    print(f"    ACF(Δ³t, n_prev)={ac_d3_np:.5f}")

    resultados.append({
        'combo':nome, 'N':N, 'np_val':np_val, 'ic95':ic95,
        'std_ddt_orig':np.std(ddt), 'std_ddt_perm':np.std(ddt_perm),
        'std_ddt_surr':np.std(ddt_surr), 'std_ddt_shft':np.std(ddt_shft),
        'acf_np_orig':ac_np_orig, 'acf_np_perm':ac_np_perm,
        'acf_np_surr':ac_np_surr, 'acf_np_shft':ac_np_shft,
    })

    # ── Gráficos ──────────────────────────────────────────────────────────
    # Col 0: ACF do Δt com todos os controles
    ax = axes[row, 0]
    ax.plot(lags_dt, ac_dt,      'b-',  lw=1.0, alpha=0.9, label='original')
    ax.plot(lags_dt, ac_dt_perm, 'r-',  lw=0.7, alpha=0.6, label='permutado')
    ax.plot(lags_dt, ac_dt_surr, 'g-',  lw=0.7, alpha=0.6, label='surrogate')
    ax.plot(lags_dt, ac_dt_shft, 'orange', lw=0.7, alpha=0.6, label='circ.shift')
    ax.axhline( ic95, color='k', lw=1, ls='--')
    ax.axhline(-ic95, color='k', lw=1, ls='--')
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7, label=f'n={np_val}')
    ax.set_title(f'{nome} — ACF(Δt)\norig({np_val})={ac_np_orig:.4f}  '
                 f'surr={ac_np_surr:.4f}  shift={ac_np_shft:.4f}', fontsize=9)
    ax.legend(fontsize=6); ax.grid(alpha=0.3)
    ax.set_xlabel('lag'); ax.set_xlim(0, 300)

    # Col 1: ACF do ΔΔt
    ax = axes[row, 1]
    ax.plot(lags, ac_ddt,      'b-', lw=1.0, alpha=0.9, label='orig ΔΔt')
    ax.plot(lags, ac_ddt_perm, 'r-', lw=0.7, alpha=0.6, label='perm ΔΔt')
    ax.plot(lags, ac_ddt_surr, 'g-', lw=0.7, alpha=0.6, label='surr ΔΔt')
    ax.plot(lags, ac_ddt_shft, 'orange', lw=0.7, alpha=0.6, label='shift ΔΔt')
    ax.axhline( ic95, color='k', lw=1, ls='--')
    ax.axhline(-ic95, color='k', lw=1, ls='--')
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7)
    ax.set_title(f'{nome} — ACF(ΔΔt)', fontsize=9)
    ax.legend(fontsize=6); ax.grid(alpha=0.3)
    ax.set_xlabel('lag'); ax.set_xlim(0, 300)

    # Col 2: Rayleigh(ΔΔt)
    ax = axes[row, 2]
    ax.semilogy(ks, V_ddt_o, 'b-', lw=0.8, alpha=0.9, label='orig')
    ax.semilogy(ks, V_ddt_p, 'r-', lw=0.7, alpha=0.5, label='perm')
    ax.semilogy(ks, V_ddt_s, 'g-', lw=0.7, alpha=0.5, label='surr')
    noise_ddt = 1/np.sqrt(len(ddt))
    ax.axhline(noise_ddt, color='k', lw=1, ls='--', label='1/√N')
    ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7, label=f'n={np_val}')
    ax.set_title(f'{nome} — Rayleigh(ΔΔt)\n'
                 f'orig={v_o:.5f}  perm={v_p:.5f}  surr={v_s:.5f}', fontsize=9)
    ax.legend(fontsize=6); ax.grid(alpha=0.3, which='both')
    ax.set_xlabel('k')

    # Col 3: ACF(Δ³t)
    ax = axes[row, 3]
    lags3_p, ac_d3_p = acf(d3t_perm[:50_000])
    ax.plot(lags3, ac_d3,   'b-', lw=1.0, alpha=0.9, label='orig Δ³t')
    ax.plot(lags3_p, ac_d3_p,'r-', lw=0.7, alpha=0.6, label='perm Δ³t')
    ax.axhline( ic95, color='k', lw=1, ls='--')
    ax.axhline(-ic95, color='k', lw=1, ls='--')
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7, label=f'n={np_val}')
    ax.set_title(f'{nome} — ACF(Δ³t)\n'
                 f'orig({np_val})={ac_d3_np:.5f}  IC95=±{ic95:.5f}', fontsize=9)
    ax.legend(fontsize=6); ax.grid(alpha=0.3)
    ax.set_xlabel('lag'); ax.set_xlim(0, 300)

    # Col 4: comparação numérica clara
    ax = axes[row, 4]
    labels = ['orig', 'perm', 'surr', 'shift']
    acf_vals = [ac_np_orig, ac_np_perm, ac_np_surr, ac_np_shft]
    colors   = ['blue','red','green','orange']
    bars = ax.bar(labels, acf_vals, color=colors, alpha=0.8)
    ax.axhline( ic95, color='k', lw=1.5, ls='--', label=f'IC95 ±{ic95:.4f}')
    ax.axhline(-ic95, color='k', lw=1.5, ls='--')
    ax.axhline(0, color='k', lw=0.5)
    for bar, val in zip(bars, acf_vals):
        ax.text(bar.get_x()+bar.get_width()/2, val+0.005*np.sign(val),
                f'{val:.4f}', ha='center', fontsize=9, fontweight='bold')
    ax.set_title(f'{nome} — ACF(Δt, n={np_val})\npor tipo de controle',
                 fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3, axis='y')
    ax.set_ylabel('ACF')

plt.tight_layout()
saida = f'{OUT}/controles_avancados.png'
plt.savefig(saida, dpi=130, bbox_inches='tight')
plt.show()
print(f"\nGráfico salvo: {saida}")

# ── CSV ───────────────────────────────────────────────────────────────────────
with open(f'{OUT}/controles_avancados.csv','w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['combo','N','np_val','ic95',
                'std_ddt_orig','std_ddt_perm','std_ddt_surr','std_ddt_shft',
                'acf_np_orig','acf_np_perm','acf_np_surr','acf_np_shft'])
    for r in resultados:
        w.writerow([r['combo'],r['N'],r['np_val'],f"{r['ic95']:.6f}",
                    f"{r['std_ddt_orig']:.4f}",f"{r['std_ddt_perm']:.4f}",
                    f"{r['std_ddt_surr']:.4f}",f"{r['std_ddt_shft']:.4f}",
                    f"{r['acf_np_orig']:.6f}",f"{r['acf_np_perm']:.6f}",
                    f"{r['acf_np_surr']:.6f}",f"{r['acf_np_shft']:.6f}"])
print(f"CSV: {OUT}/controles_avancados.csv")

print("\n"+"="*65)
print("VEREDICTO FINAL")
print("="*65)
print(f"""
Critérios:
  ACF(orig) >> IC95 e ACF(perm) ≈ 0  → estrutura temporal real
  ACF(surr) ≈ ACF(orig)              → estrutura vem do espectro de potência
  ACF(shift) ≈ ACF(orig)             → estrutura sobrevive ao shift
  std(ΔΔt_orig) << std(ΔΔt_perm)    → Δt consecutivos são correlacionados
""")
for r in resultados:
    orig = r['acf_np_orig']; perm = r['acf_np_perm']
    surr = r['acf_np_surr']; shft = r['acf_np_shft']
    ic   = r['ic95']
    if abs(orig) > ic and abs(perm) < ic:
        if abs(surr) > ic:
            diag = '✓ estrutura real (sobrevive ao surrogate)'
        else:
            diag = '~ espectro de potência (não sobrevive ao surrogate)'
        if abs(shft) > ic:
            diag += ' + shift'
    else:
        diag = '✗ sem estrutura além do ruído'
    print(f"  {r['combo']}: ACF orig={orig:.4f} perm={perm:.4f} "
          f"surr={surr:.4f} shift={shft:.4f}  → {diag}")
