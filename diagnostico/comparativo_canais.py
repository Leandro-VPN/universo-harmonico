# ══════════════════════════════════════════════════════════════════════
# ║  DIAGNÓSTICO — COMPARATIVO DE CANAIS ENTRE RUNS                     ║
# ║  Compara estrutura de canais entre diferentes pares Alice/Bob       ║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════


import numpy as np
import os, time
from pathlib import Path
import matplotlib.pyplot as plt

ALICE_FILE = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat"
BOB_FILE   = r"d:\BellQM\19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat"
OUT_DIR    = r"d:\BellQM\out_comparativo"
os.makedirs(OUT_DIR, exist_ok=True)

BIN_PS    = 78.125
CH_SYNC   = 6
N_TRIALS  = 2_000_000
SENTINEL  = np.uint32(0xFFFFFFFF)

print("="*65)
print("TESTE COMPARATIVO: ch=0 vs ch=2 vs ch=4")
print("="*65)

def extrair_todos_canais(filepath, n_max, label):
    """Extrai tt_rel separadamente para ch=0, ch=2, ch=4"""
    tr = {0: np.full(n_max, SENTINEL, dtype=np.uint32),
          2: np.full(n_max, SENTINEL, dtype=np.uint32),
          4: np.full(n_max, SENTINEL, dtype=np.uint32)}
    sync_abs = np.full(n_max, -1, dtype=np.int64)
    trial=-1; sync_tt=None
    t0=time.time()

    with open(filepath,'rb') as f:
        while trial < n_max-1:
            raw = f.read(5_000_000*24)
            if not raw: break
            n = len(raw)//24
            arr = np.frombuffer(raw,dtype=np.uint64).reshape(n,3)
            ch = arr[:,0].astype(np.int32)
            tt = arr[:,1].astype(np.int64)
            del arr,raw
            idx = np.flatnonzero((ch==6)|(ch==0)|(ch==2)|(ch==4))
            parou=False
            for i in idx:
                c=int(ch[i]); t=int(tt[i])
                if c==6:
                    trial+=1; sync_tt=t
                    if trial<n_max: sync_abs[trial]=t
                    if trial>=n_max-1: parou=True; break
                elif sync_tt is not None and trial>=0 and trial<n_max:
                    d=t-sync_tt
                    if 0<d<130000:  # dentro de 1 ciclo completo
                        if c in tr and tr[c][trial]==SENTINEL:
                            tr[c][trial]=np.uint32(d)
            del ch,tt,idx
            if parou: break

    n_ok=trial+1
    print(f"\n  {label}: {n_ok:,} trials extraídos ({time.time()-t0:.1f}s)")
    for c in [0,2,4]:
        vals = tr[c][:n_ok]
        n_click = int((vals!=SENTINEL).sum())
        taxa = n_click/n_ok*100
        v_ok = vals[vals!=SENTINEL]
        if n_click>0:
            unicos = len(np.unique(v_ok))
            rng = f"[{v_ok.min()},{v_ok.max()}]"
            med = np.median(v_ok)
            print(f"    ch={c}: {n_click:>8,} ({taxa:5.2f}%)  "
                  f"tt_rel: {unicos} únicos  range={rng}  med={med:.0f}")
        else:
            print(f"    ch={c}: 0 cliques")
    return tr, sync_abs, n_ok

print("\nExtraindo Alice...")
tr_a, sy_a, na = extrair_todos_canais(ALICE_FILE, N_TRIALS, 'Alice')
print("\nExtraindo Bob...")
tr_b, sy_b, nb = extrair_todos_canais(BOB_FILE,   N_TRIALS, 'Bob')

n = min(na, nb)

# ── Alinhamento dos syncs ────────────────────────────────────────────────────
n_al = min(n, 5000)
sa = sy_a[:n_al]; sb = sy_b[:n_al]
clock_offset = int(np.median(sa-sb))
residuo = sa - sb - clock_offset
std_res = np.std(residuo[np.abs(residuo)<500])
print(f"\n  Clock offset: {clock_offset:,} bins | std residual: {std_res:.1f} bins ({std_res*BIN_PS:.0f} ps)")

# ── Para cada combinação de canal, calcular Δt e sua std ────────────────────
print("\n" + "="*65)
print("COMPARAÇÃO DO Δt POR CANAL")
print("Critério: fóton real → std(Δt) < 5 bins (<390 ps)")
print("="*65)

combos_canal = [
    ('Alice ch=2, Bob ch=2', tr_a[2], tr_b[2]),
    ('Alice ch=4, Bob ch=4', tr_a[4], tr_b[4]),
    ('Alice ch=2, Bob ch=4', tr_a[2], tr_b[4]),  # G1×G2
    ('Alice ch=4, Bob ch=2', tr_a[4], tr_b[2]),  # G2×G1
    ('Alice ch=0, Bob ch=0', tr_a[0], tr_b[0]),
    ('Alice ch=0, Bob ch=2', tr_a[0], tr_b[2]),
    ('Alice ch=0, Bob ch=4', tr_a[0], tr_b[4]),
    ('Alice ch=2, Bob ch=0', tr_a[2], tr_b[0]),
    ('Alice ch=4, Bob ch=0', tr_a[4], tr_b[0]),
]

resultados = []
for nome, ta, tb in combos_canal:
    mask = (ta[:n]!=SENTINEL)&(tb[:n]!=SENTINEL)
    N = mask.sum()
    if N < 100:
        print(f"  {nome}: N={N} — insuficiente")
        continue
    dt = (sy_a[:n][mask].astype(np.int64)
          - sy_b[:n][mask].astype(np.int64)
          - clock_offset
          + ta[:n][mask].astype(np.int64)
          - tb[:n][mask].astype(np.int64))
    med = float(np.median(dt))
    std = float(np.std(dt))
    # Std dentro de ±500 bins do pico (elimina outliers)
    dt_c = dt - int(med)
    std_pico = float(np.std(dt_c[np.abs(dt_c)<500]))
    flag = '← FÓTON REAL?' if std_pico < 5 else ('← INTERMEDIÁRIO' if std_pico < 50 else '')
    print(f"  {nome}")
    print(f"    N={N:,}  med={med:.0f}  std={std:.1f}  std_pico={std_pico:.2f} bins ({std_pico*BIN_PS:.0f} ps)  {flag}")
    resultados.append((nome, dt, N, med, std_pico))

# ── Plotar histogramas comparativos ─────────────────────────────────────────
n_plots = min(len(resultados), 9)
fig, axes = plt.subplots(3, 3, figsize=(18, 14))
fig.suptitle('Δt por combinação de canal — identificando o fóton real\n'
             'std < 5 bins (390 ps) = coincidência física real',
             fontsize=12, fontweight='bold')

for idx, (nome, dt, N, med, std_pico) in enumerate(resultados[:9]):
    ax = axes[idx//3, idx%3]
    dt_c = dt - int(med)
    ax.hist(dt_c, bins=200, range=(-500, 500), color='steelblue', alpha=0.8)
    ax.axvline(0, color='red', lw=1.5, ls='--')
    ax.set_title(f'{nome}\nN={N:,}  std_pico={std_pico:.1f} bins ({std_pico*BIN_PS:.0f} ps)',
                 fontsize=9)
    ax.set_xlabel('Δt − mediana (bins)'); ax.grid(alpha=0.3)

for idx in range(len(resultados), 9):
    axes[idx//3, idx%3].set_visible(False)

plt.tight_layout()
saida = f'{OUT_DIR}/comparativo_canais.png'
plt.savefig(saida, dpi=150, bbox_inches='tight')
plt.show()
print(f"\nGráfico salvo: {saida}")

# ── Veredicto ────────────────────────────────────────────────────────────────
print("\n" + "="*65)
print("VEREDICTO")
print("="*65)
melhor = min(resultados, key=lambda x: x[4]) if resultados else None
if melhor:
    nome, dt, N, med, std_pico = melhor
    print(f"\n  Menor std_pico: '{nome}'")
    print(f"  std = {std_pico:.2f} bins = {std_pico*BIN_PS:.0f} ps")
    if std_pico < 5:
        print(f"\n  ✓ FÓTON REAL CONFIRMADO")
        print(f"  Use esses canais no pipeline definitivo.")
    elif std_pico < 20:
        print(f"\n  ~ PARCIALMENTE COINCIDENTE — investigar mais")
    else:
        print(f"\n  ✗ Nenhuma combinação mostra coincidência física clara.")
        print(f"  O arquivo .sync.hdf5 do NIST pode ser necessário para")
        print(f"  identificar os syncs válidos (booleanCut).")
