# ══════════════════════════════════════════════════════════════════════
# ║  CONTROLE POR PERMUTAÇÃO — Rayleigh(Δt) é artefato de distribuição  ║
# ║  Prova que R(dt)=R(perm(dt)): Rayleigh não detecta ordem temporal   ║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════


import numpy as np, os, time
from scipy.stats import chi2
import matplotlib.pyplot as plt, csv

# ── Carregar os Δt do v7 (já calculados) ──────────────────────────────────────
# Rodar após o v7 ter gerado o CSV
CSV_V7 = r"d:\BellQM\out_v7\rayleigh_v7.csv"
DAT_V7 = r"d:\BellQM\out_v7"
OUT    = r"d:\BellQM\out_controle"
os.makedirs(OUT, exist_ok=True)

# Se não tiver o CSV do v7, recalcular os Δt aqui
# Caso tenha, carregar direto
ALICE_FILE = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat"
BOB_FILE   = r"d:\BellQM\19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat"

BIN_PS   = 78.125
SENTINEL = np.uint32(0xFFFFFFFF)
COMBOS   = {'XX':(0,0),'XY':(0,1),'YX':(1,0),'YY':(1,1)}
N_PREV   = {'XX':178,'XY':15,'YX':15,'YY':5}
K_MAX    = 500   # até 500 é suficiente para o teste

SLOTS_A_G1 = set(range(90, 97))
SLOTS_A_G2 = set(range(104, 111))
SLOTS_B_G1 = set(range(347, 355))
SLOTS_B_G2 = set(range(361, 368))
SLOTS_A    = SLOTS_A_G1 | SLOTS_A_G2
SLOTS_B    = SLOTS_B_G1 | SLOTS_B_G2

print("="*65)
print("CONTROLE POR PERMUTAÇÃO")
print("="*65)

# ── Extrair Δt (igual ao v7) ──────────────────────────────────────────────────
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

print("\nExtraindo dados (igual ao v7)...")
N = 10_000_000
t0=time.time()
tr_a, st_a, sy_a = extrair(ALICE_FILE, SLOTS_A_G1, SLOTS_A_G2, SLOTS_A, N, 'Alice')
tr_b, st_b, sy_b = extrair(BOB_FILE,   SLOTS_B_G1, SLOTS_B_G2, SLOTS_B, N, 'Bob')
n = min(len(tr_a), len(tr_b))
clock_offset = int(np.median(
    sy_a[:5000].astype(np.int64) - sy_b[:5000].astype(np.int64)))
print(f"  clock_offset={clock_offset:,}  ({time.time()-t0:.1f}s)")

# Calcular Δt por combo
dts_orig = {}
for nome,(sa_t,sb_t) in COMBOS.items():
    mask = ((tr_a[:n]!=SENTINEL)&(tr_b[:n]!=SENTINEL)&
            (st_a[:n]==sa_t)&(st_b[:n]==sb_t))
    dt = ((sy_a[:n][mask].astype(np.int64)-sy_b[:n][mask].astype(np.int64)-clock_offset)
          + tr_a[:n][mask].astype(np.int64) - tr_b[:n][mask].astype(np.int64))
    dts_orig[nome] = dt
    print(f"  {nome}: N={len(dt):,}  med={np.median(dt):.1f}  std={np.std(dt):.2f}")

# ── Função Rayleigh ────────────────────────────────────────────────────────────
def rayleigh_V(dt, k_max=K_MAX):
    N = len(dt)
    if N < 100: return None, None
    dt_c = dt - int(np.median(dt))
    ks = np.arange(1, k_max+1)
    V  = np.zeros(k_max)
    for k in ks:
        phi = 2*np.pi*(dt_c%k)/k
        C = np.mean(np.cos(phi)); S = np.mean(np.sin(phi))
        Z = 2*N*(C**2+S**2)
        V[k-1] = np.sqrt(max(Z,0)/(2*N))
    return ks, V

# ── Para cada combo, rodar os 3 controles ────────────────────────────────────
print("\n"+"="*65)
print("RODANDO CONTROLES")
print("="*65)

rng = np.random.default_rng(42)

fig, axes = plt.subplots(len(COMBOS), 4, figsize=(22, 5*len(COMBOS)))
fig.suptitle('Controle por permutação — v7\n'
             'Se original ≈ permutado → artefato de distribuição\n'
             'Se original ≠ permutado → estrutura temporal real',
             fontsize=12, fontweight='bold')

for row, (nome, dt_orig) in enumerate(dts_orig.items()):
    np_val = N_PREV[nome]
    N      = len(dt_orig)
    noise  = 1/np.sqrt(N) if N>0 else 1

    # Original
    ks, V_orig = rayleigh_V(dt_orig)

    # Controle 1: permutação aleatória
    dt_perm = rng.permutation(dt_orig)
    _, V_perm = rayleigh_V(dt_perm)

    # Controle 2: Normal com mesma média e std
    med_v = float(np.median(dt_orig)); std_v = float(np.std(dt_orig))
    dt_norm = np.round(rng.normal(med_v, std_v, N)).astype(np.int64)
    _, V_norm = rayleigh_V(dt_norm)

    # Controle 3: combinação aleatória de slots reais
    slots_a_arr = sorted(SLOTS_A_G1 | SLOTS_A_G2)
    slots_b_arr = sorted(SLOTS_B_G1 | SLOTS_B_G2)
    sa_rand = rng.choice(slots_a_arr, size=N)
    sb_rand = rng.choice(slots_b_arr, size=N)
    dt_slot = sa_rand.astype(np.int64) - sb_rand.astype(np.int64)
    _, V_slot = rayleigh_V(dt_slot)

    # Estatísticas nos pontos de interesse
    print(f"\n  {nome}: N={N:,}  n_prev={np_val}")
    print(f"  {'':20} {'V(n_prev)':>12} {'S/N':>8} {'top k':>8} {'V(top)':>10}")
    for lbl, V in [('original', V_orig), ('permutado', V_perm),
                   ('normal', V_norm), ('slots aleat.', V_slot)]:
        if V is None: continue
        v_np  = V[np_val-1]
        k_top = int(ks[np.argmax(V)])
        v_top = V.max()
        print(f"  {lbl:20} {v_np:>12.6f} {v_np/noise:>8.2f}× "
              f"{k_top:>8} {v_top:>10.6f}")

    # Plots
    for col, (lbl, V, cor) in enumerate([
        ('original', V_orig, 'blue'),
        ('permutado', V_perm, 'red'),
        ('normal σ=std', V_norm, 'orange'),
        ('slots aleat.', V_slot, 'green'),
    ]):
        ax = axes[row, col]
        if V is None: ax.set_title(f'{nome} {lbl} — sem dados'); continue
        ax.semilogy(ks, V, '-', color=cor, lw=0.8, alpha=0.8)
        ax.axhline(noise, color='black', lw=1, ls='--',
                   label=f'1/√N={noise:.5f}')
        ax.axvline(np_val, color='purple', lw=1.5, alpha=0.7,
                   label=f'n_prev={np_val}')
        v_np = V[np_val-1]
        k_top = int(ks[np.argmax(V)])
        ax.set_title(f'{nome} — {lbl}\nV({np_val})={v_np:.5f} '
                     f'S/N={v_np/noise:.1f}×  top_k={k_top}',
                     fontsize=9)
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3, which='both')
        ax.set_xlabel('k')
        if col==0: ax.set_ylabel('V(k)')

plt.tight_layout()
saida = f'{OUT}/controle_permutacao.png'
plt.savefig(saida, dpi=150, bbox_inches='tight')
plt.show()
print(f"\nGráfico salvo: {saida}")

# ── Veredicto ─────────────────────────────────────────────────────────────────
print("\n"+"="*65)
print("VEREDICTO FINAL")
print("="*65)
print("""
Critério:
  Se V_orig(n_prev) ≈ V_perm(n_prev) → artefato de distribuição
  Se V_orig(n_prev) >> V_perm(n_prev) → estrutura temporal real
  Se top_k(orig) == n_prev e top_k(perm) ≠ n_prev → sinal seletivo
""")

for nome, dt_orig in dts_orig.items():
    np_val = N_PREV[nome]
    N = len(dt_orig)
    if N == 0: continue
    noise = 1/np.sqrt(N)
    _, V_o = rayleigh_V(dt_orig)
    _, V_p = rayleigh_V(rng.permutation(dt_orig))
    if V_o is None or V_p is None: continue
    v_o = V_o[np_val-1]; v_p = V_p[np_val-1]
    ratio = v_o/v_p if v_p > 0 else float('inf')
    k_top_o = int(ks[np.argmax(V_o)])
    seletivo = '✓ SELETIVO' if k_top_o == np_val else f'top_k={k_top_o}'
    artefato = '✗ ARTEFATO' if ratio < 2 else ('✓ ESTRUTURA REAL' if ratio > 5 else '~ AMBÍGUO')
    print(f"  {nome}: V_orig={v_o:.5f}  V_perm={v_p:.5f}  "
          f"ratio={ratio:.1f}×  {artefato}  {seletivo}")
