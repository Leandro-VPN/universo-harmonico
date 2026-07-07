# ══════════════════════════════════════════════════════════════════════
# ║  DISTRIBUIÇÃO NULA (versão VS Code)                                 ║
# ║  Variante Windows do distribuicao_nula.py                           ║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════

======================================
Gera a distribuição nula do ACF(n_prev) via 1000 permutações
e calcula z-score e p-valor para o valor original.
"""

import numpy as np, os, time
import matplotlib.pyplot as plt, csv
from scipy import stats

ALICE_FILE = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat"
BOB_FILE   = r"d:\BellQM\19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat"
OUT        = r"d:\BellQM\out_nula"
os.makedirs(OUT, exist_ok=True)

BIN_PS   = 78.125
SENTINEL = np.uint32(0xFFFFFFFF)
COMBOS   = {'XX':(0,0),'XY':(0,1),'YX':(1,0),'YY':(1,1)}
N_PREV   = {'XX':178,'XY':15,'YX':15,'YY':5}
N_PERM   = 1000

SLOTS_A_G1=set(range(90,97));  SLOTS_A_G2=set(range(104,111))
SLOTS_B_G1=set(range(347,355)); SLOTS_B_G2=set(range(361,368))
SLOTS_A=SLOTS_A_G1|SLOTS_A_G2; SLOTS_B=SLOTS_B_G1|SLOTS_B_G2

rng = np.random.default_rng(42)

# ── Extrair ───────────────────────────────────────────────────────────────────
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
tr_b,st_b,sy_b=extrair(BOB_FILE,  SLOTS_B_G1,SLOTS_B_G2,SLOTS_B,10_000_000,'Bob')
n=min(len(tr_a),len(tr_b))
clock_offset=int(np.median(sy_a[:5000].astype(np.int64)-sy_b[:5000].astype(np.int64)))

dts_orig={}
for nome,(sa_t,sb_t) in COMBOS.items():
    mask=((tr_a[:n]!=SENTINEL)&(tr_b[:n]!=SENTINEL)&
          (st_a[:n]==sa_t)&(st_b[:n]==sb_t))
    dt=((sy_a[:n][mask].astype(np.int64)-sy_b[:n][mask].astype(np.int64)-clock_offset)
        +tr_a[:n][mask].astype(np.int64)-tr_b[:n][mask].astype(np.int64))
    dts_orig[nome]=dt
    print(f"  {nome}: N={len(dt):,}")
print(f"  ({time.time()-t0:.1f}s)")

# ── ACF rápida (só um lag) ────────────────────────────────────────────────────
def acf_lag(x, lag):
    """ACF normalizada no lag específico — vetorizada."""
    x = x.astype(np.float64)
    mu = np.mean(x); var = np.var(x)
    if var == 0: return 0.0
    N = len(x)
    return float(np.mean((x[:N-lag]-mu)*(x[lag:]-mu))/var)

# ── Distribuição nula ─────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"DISTRIBUIÇÃO NULA — {N_PERM} permutações por combo")
print(f"{'='*65}")

fig, axes = plt.subplots(2, 4, figsize=(22, 10))
fig.suptitle(f'Distribuição nula do ACF(n_prev) — {N_PERM} permutações\n'
             'Se original estiver fora da distribuição nula → estrutura real',
             fontsize=12, fontweight='bold')

resultados = []

for col, (nome, dt) in enumerate(dts_orig.items()):
    np_val = N_PREV[nome]
    N = len(dt)
    N_ac = min(N, 100_000)  # usar 100k para velocidade
    dt_sub = dt[:N_ac]

    # ACF original
    acf_orig = acf_lag(dt_sub, np_val)

    # Distribuição nula
    t0 = time.time()
    acf_null = np.zeros(N_PERM)
    for i in range(N_PERM):
        dt_p = rng.permutation(dt_sub)
        acf_null[i] = acf_lag(dt_p, np_val)
        if (i+1) % 100 == 0:
            print(f"\r  {nome}: {i+1}/{N_PERM} permutações | "
                  f"{time.time()-t0:.1f}s", end='', flush=True)

    print(f"\r  {nome}: {N_PERM} permutações | {time.time()-t0:.1f}s")

    # Estatísticas
    mu_null  = float(np.mean(acf_null))
    std_null = float(np.std(acf_null))
    z_score  = (acf_orig - mu_null) / std_null if std_null > 0 else float('inf')
    p_value  = float(2 * stats.norm.sf(abs(z_score)))  # two-tailed
    ic95     = 2 / np.sqrt(N_ac)

    print(f"  {nome}: ACF_orig={acf_orig:.5f}  "
          f"null_mean={mu_null:.5f}  null_std={std_null:.5f}  "
          f"z={z_score:.1f}  p={p_value:.2e}")

    resultados.append({
        'combo':nome, 'N':N, 'N_ac':N_ac, 'np_val':np_val,
        'acf_orig':acf_orig, 'null_mean':mu_null, 'null_std':std_null,
        'z':z_score, 'p':p_value, 'ic95':ic95,
        'acf_null_min':acf_null.min(), 'acf_null_max':acf_null.max(),
        'pct_above': float((acf_null >= acf_orig).mean()),
    })

    # ── Histograma da distribuição nula ───────────────────────────────────
    ax1 = axes[0, col]
    ax1.hist(acf_null, bins=50, color='lightcoral', edgecolor='darkred',
             alpha=0.8, density=True, label=f'nula ({N_PERM} perm.)')
    ax1.axvline(acf_orig, color='blue', lw=3, label=f'original={acf_orig:.4f}')
    ax1.axvline(mu_null,  color='red',  lw=1.5, ls='--',
                label=f'média nula={mu_null:.4f}')
    # Gaussiana teórica
    x_g = np.linspace(acf_null.min()-3*std_null,
                      max(acf_orig*1.1, acf_null.max()+3*std_null), 300)
    ax1.plot(x_g, stats.norm.pdf(x_g, mu_null, std_null),
             'k-', lw=1.5, alpha=0.7, label='Gaussiana')
    ax1.set_title(f'{nome} — ACF(lag={np_val})\n'
                  f'z={z_score:.1f}  p={p_value:.2e}\n'
                  f'orig={acf_orig:.4f}  '
                  f'[nula: {mu_null:.4f}±{std_null:.4f}]',
                  fontsize=9)
    ax1.legend(fontsize=7); ax1.grid(alpha=0.3)
    ax1.set_xlabel('ACF')

    # ── ACF completa original vs faixa nula ───────────────────────────────
    ax2 = axes[1, col]
    # Calcular ACF original para lags 1..300
    dt_sub64 = dt_sub.astype(np.float64)
    mu = np.mean(dt_sub64); var = np.var(dt_sub64)
    N_s = len(dt_sub)
    max_lag = 300
    acf_full = np.array([
        float(np.mean((dt_sub64[:N_s-lag]-mu)*(dt_sub64[lag:]-mu))/var)
        for lag in range(1, max_lag+1)
    ])

    # Faixa da distribuição nula em cada lag (usando 20 permutações rápidas)
    acf_null_band = np.zeros((20, max_lag))
    for i in range(20):
        dp = rng.permutation(dt_sub64)
        acf_null_band[i] = np.array([
            float(np.mean((dp[:N_s-lag]-mu)*(dp[lag:]-mu))/var)
            for lag in range(1, max_lag+1)
        ])
    null_lo = np.percentile(acf_null_band, 2.5, axis=0)
    null_hi = np.percentile(acf_null_band, 97.5, axis=0)

    lags = np.arange(1, max_lag+1)
    ax2.fill_between(lags, null_lo, null_hi, alpha=0.3, color='red',
                     label='IC95 nula (20 perm.)')
    ax2.plot(lags, acf_full, 'b-', lw=0.8, alpha=0.9, label='original')
    ax2.axhline( ic95, color='green', lw=1, ls='--', label=f'±IC95={ic95:.4f}')
    ax2.axhline(-ic95, color='green', lw=1, ls='--')
    ax2.axhline(0, color='k', lw=0.5)
    ax2.axvline(np_val, color='purple', lw=1.5, alpha=0.7,
                label=f'n_prev={np_val}')
    ax2.set_title(f'{nome} — ACF completa (lag 1..300)\n'
                  f'Faixa vermelha = IC95 da distribuição nula',
                  fontsize=9)
    ax2.legend(fontsize=6); ax2.grid(alpha=0.3)
    ax2.set_xlabel('lag')

plt.tight_layout()
saida = f'{OUT}/distribuicao_nula.png'
plt.savefig(saida, dpi=150, bbox_inches='tight')
plt.show()
print(f"\nGráfico salvo: {saida}")

# ── CSV ───────────────────────────────────────────────────────────────────────
with open(f'{OUT}/distribuicao_nula.csv','w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['combo','N','N_ac','np_val','acf_orig',
                'null_mean','null_std','z','p_value','ic95',
                'null_min','null_max','pct_above'])
    for r in resultados:
        w.writerow([r['combo'],r['N'],r['N_ac'],r['np_val'],
                    f"{r['acf_orig']:.6f}",f"{r['null_mean']:.6f}",
                    f"{r['null_std']:.6f}",f"{r['z']:.2f}",
                    f"{r['p']:.4e}",f"{r['ic95']:.6f}",
                    f"{r['acf_null_min']:.6f}",f"{r['acf_null_max']:.6f}",
                    f"{r['pct_above']:.6f}"])
print(f"CSV: {OUT}/distribuicao_nula.csv")

print(f"\n{'='*65}")
print("VEREDICTO FINAL — DISTRIBUIÇÃO NULA")
print(f"{'='*65}")
for r in resultados:
    flag = '✓ ESTRUTURA REAL' if abs(r['z'])>10 else ('~ MARGINAL' if abs(r['z'])>3 else '✗ RUÍDO')
    print(f"  {r['combo']}: ACF_orig={r['acf_orig']:.5f}  "
          f"nula=[{r['null_mean']:.5f}±{r['null_std']:.5f}]  "
          f"z={r['z']:.1f}  p={r['p']:.2e}  {flag}")
print(f"""
Interpretação:
  z >> 10 e p ≈ 0  → o valor original é impossível sob a hipótese nula
  Se surrogate mantém ACF alto → a estrutura vem do espectro de potência
  (os slots têm correlação temporal — física do Pockels cell, não harmônica)
  O teste definitivo continua sendo ângulos com razão harmônica p:q.
""")


__________________________________________________________
import numpy as np
import h5py

def process_file(filename, key, sync_channel):
    print("Processing", key, filename)
    dt = np.dtype('u8')
    data = np.memmap(filename, dtype=dt, mode='r')
    data = data.reshape(-1, 3)
    
    booleanCut = data[:,0] == sync_channel
    sync = data[booleanCut, 1]
    nSyncs = booleanCut.sum()
    print('nSyncs', key, nSyncs)

    syncperiod = np.diff(sync).astype(np.int64)

    checklaserperiod = np.floor(syncperiod / (129102/800.) + 0.5)
    badSyncIdx = np.where(checklaserperiod != 800)[0]
    
    badSyncInfo = None
    modelockOffsetInfo = None

    if badSyncIdx.size > 0:
        print('syncperiod', key, checklaserperiod[checklaserperiod != 800])
        print('badSyncIdx', key, badSyncIdx)
        badSyncInfo = np.vstack((badSyncIdx, sync[badSyncIdx], sync[badSyncIdx+1], checklaserperiod[checklaserperiod != 800])).astype('int64')
        print("badSyncInfo", badSyncInfo)

        badSync = badSyncInfo
        modelockbool = np.abs(badSync[-1,:]) < 790
        if modelockbool.sum() > 0:
            modeLockIdx = np.where(modelockbool)[0]
            print(key, 'modelockbool', modelockbool.sum(), badSync[:,modeLockIdx[0]])
            lastSyncIdx = badSync[:,modeLockIdx[0]][0] 
            lastIdx = np.where(booleanCut)[0][lastSyncIdx]+1                     
            modelockOffsetInfo = np.array([lastSyncIdx, lastIdx])
            print('modelockOffsetInfo', modelockOffsetInfo)

        gapbool = np.abs(badSync[-1,:]) > 810
        gapbool = gapbool & (np.abs(badSync[-1,:]) < (100000*800))
        print(key, 'gapbool', gapbool.sum(), badSync[:,gapbool])

        ttagjumpbool = (np.abs(badSync[-1,:]) > (100000*800))
        if ttagjumpbool.sum() > 0:
            ttagIdx = badSync[0,ttagjumpbool]
            print('patching syncperiod at', ttagIdx)
            syncperiod[ttagIdx] = syncperiod[ttagIdx-1]

    # Save to HDF5
    with h5py.File(f'{key}.sync.hdf5', 'w') as hdf:
        hdf.create_group('cuts')
        hdf.create_group('stuff')
        hdf['cuts'].create_group(key)
        hdf['stuff'].create_group(key)
        
        hdf['cuts'][key]['sync'] = booleanCut
        hdf['stuff'][key]['nSyncs'] = nSyncs
        hdf['stuff'][key]['sync'] = sync
        hdf['stuff'][key]['syncPeriod'] = syncperiod
        if badSyncInfo is not None:
            hdf['stuff'][key]['badSyncInfo'] = badSyncInfo
        if modelockOffsetInfo is not None:
            hdf['stuff'][key]['modelockOffsetInfo'] = modelockOffsetInfo

process_file('19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat', 'alice', 0)
process_file('19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat', 'bob', 0)
