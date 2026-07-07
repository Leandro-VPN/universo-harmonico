# ══════════════════════════════════════════════════════════════════════
# ║  DIAGNÓSTICO — ESTRUTURA TEMPORAL DENTRO DO TRIAL                   ║
# ║  Verifica mutualidade exclusiva ch=2 e ch=4 (detectores vs RNG)     ║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════


import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt

ALICE_FILE = r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat"
BOB_FILE   = r"d:\BellQM\19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat"
OUT_DIR    = r"d:\BellQM\out_estrutura"
os.makedirs(OUT_DIR, exist_ok=True)

BIN_PS   = 78.125
CH_SYNC  = 6
N_SYNCS  = 200_000

print("="*65)
print("ESTRUTURA TEMPORAL DENTRO DO TRIAL")
print("="*65)

def inspecionar_trial(filepath, n_syncs, label):
    """
    Para cada trial, registra TODOS os eventos de ch=2 e ch=4
    e suas posições relativas ao sync.
    """
    # contadores
    n_s0_only = 0   # trial com só ch=2
    n_s1_only = 0   # trial com só ch=4
    n_both    = 0   # trial com ch=2 E ch=4
    n_none    = 0   # trial sem nenhum

    # posições dos dois eventos quando ambos presentes
    pos_first  = []  # tt_rel do primeiro evento
    pos_second = []  # tt_rel do segundo evento
    ch_first   = []  # canal do primeiro evento
    ch_second  = []  # canal do segundo evento

    # separação temporal entre os dois eventos
    delta_entre = []

    # acumular eventos do trial corrente
    trial_events = []   # lista de (tt_rel, canal)
    trial = -1
    sync_tt = None

    with open(filepath, 'rb') as f:
        while trial < n_syncs - 1:
            raw = f.read(5_000_000 * 24)
            if not raw: break
            n = len(raw) // 24
            arr = np.frombuffer(raw, dtype=np.uint64).reshape(n, 3)
            ch = arr[:, 0].astype(np.int32)
            tt = arr[:, 1].astype(np.int64)
            del arr, raw

            idx = np.flatnonzero((ch==6)|(ch==2)|(ch==4))
            for i in idx:
                c = int(ch[i]); t = int(tt[i])

                if c == 6:
                    # processar trial anterior
                    if trial >= 0:
                        ev2 = [(r,c2) for r,c2 in trial_events if c2==2]
                        ev4 = [(r,c2) for r,c2 in trial_events if c2==4]
                        if ev2 and ev4:
                            n_both += 1
                            # ordenar por tt_rel
                            all_ev = sorted(trial_events, key=lambda x: x[0])
                            pos_first.append(all_ev[0][0])
                            pos_second.append(all_ev[1][0])
                            ch_first.append(all_ev[0][1])
                            ch_second.append(all_ev[1][1])
                            delta_entre.append(all_ev[1][0] - all_ev[0][0])
                        elif ev2:
                            n_s0_only += 1
                        elif ev4:
                            n_s1_only += 1
                        else:
                            n_none += 1

                    trial += 1
                    sync_tt = t
                    trial_events = []
                    if trial >= n_syncs - 1:
                        break

                elif sync_tt is not None and trial >= 0:
                    d = t - sync_tt
                    if 0 < d < 130000:
                        trial_events.append((d, c))

            del ch, tt
            if trial >= n_syncs - 1:
                break

    total = n_s0_only + n_s1_only + n_both + n_none
    print(f"\n  {label} ({total:,} trials):")
    print(f"    só ch=2:       {n_s0_only:>8,} ({n_s0_only/total*100:.2f}%)")
    print(f"    só ch=4:       {n_s1_only:>8,} ({n_s1_only/total*100:.2f}%)")
    print(f"    ch=2 E ch=4:   {n_both:>8,} ({n_both/total*100:.2f}%)")
    print(f"    nenhum:        {n_none:>8,} ({n_none/total*100:.2f}%)")

    if n_both > 100:
        pf = np.array(pos_first)
        ps = np.array(pos_second)
        de = np.array(delta_entre)
        cf = np.array(ch_first)
        cs = np.array(ch_second)

        print(f"\n    Quando ambos presentes (N={n_both:,}):")
        print(f"    1° evento: med={np.median(pf):.1f}  range=[{pf.min()},{pf.max()}]")
        print(f"    2° evento: med={np.median(ps):.1f}  range=[{ps.min()},{ps.max()}]")
        print(f"    Δ entre eles: med={np.median(de):.1f}  std={np.std(de):.1f}  range=[{de.min()},{de.max()}]")

        # canal do 1° evento
        print(f"    Canal do 1° evento: ch=2 em {(cf==2).mean()*100:.1f}%  ch=4 em {(cf==4).mean()*100:.1f}%")
        print(f"    Canal do 2° evento: ch=2 em {(cs==2).mean()*100:.1f}%  ch=4 em {(cs==4).mean()*100:.1f}%")

        # Se Δ entre eventos é ~constante → um é setting (regular), outro é fóton
        print(f"\n    Δ entre eventos (separação temporal):")
        vals_de, cnts_de = np.unique(de, return_counts=True)
        top = np.argsort(cnts_de)[-10:][::-1]
        for i in top:
            print(f"      Δ={vals_de[i]:6d} bins ({vals_de[i]*BIN_PS:.0f} ps): {cnts_de[i]:,}")

    return {
        'n_s0': n_s0_only, 'n_s1': n_s1_only,
        'n_both': n_both, 'n_none': n_none,
        'pos_first': np.array(pos_first) if pos_first else np.array([]),
        'pos_second': np.array(pos_second) if pos_second else np.array([]),
        'delta': np.array(delta_entre) if delta_entre else np.array([]),
        'ch_first': np.array(ch_first) if ch_first else np.array([]),
    }

r_a = inspecionar_trial(ALICE_FILE, N_SYNCS, 'Alice')
r_b = inspecionar_trial(BOB_FILE,   N_SYNCS, 'Bob')

# ── Gráficos ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 4, figsize=(20, 8))
fig.suptitle('Estrutura temporal dentro do trial\n'
             'Quando ch=2 E ch=4 ocorrem no mesmo trial: qual é o fóton?',
             fontsize=12, fontweight='bold')

for row, (lbl, r) in enumerate([('Alice', r_a), ('Bob', r_b)]):
    if len(r['pos_first']) == 0:
        continue

    # col 0: posição do 1° evento
    ax = axes[row, 0]
    ax.hist(r['pos_first'], bins=100, color='steelblue', alpha=0.8)
    ax.set_title(f'{lbl} — 1° evento (tt_rel)')
    ax.set_xlabel('bins'); ax.grid(alpha=0.3)

    # col 1: posição do 2° evento
    ax = axes[row, 1]
    ax.hist(r['pos_second'], bins=200, color='darkorange', alpha=0.8)
    ax.set_title(f'{lbl} — 2° evento (tt_rel)')
    ax.set_xlabel('bins'); ax.grid(alpha=0.3)

    # col 2: Δ entre os dois eventos
    ax = axes[row, 2]
    de = r['delta']
    ax.hist(de, bins=200, range=(0, de.max() if len(de)>0 else 1000),
            color='green', alpha=0.8)
    ax.set_title(f'{lbl} — Δ entre 1° e 2° evento')
    ax.set_xlabel('bins'); ax.grid(alpha=0.3)

    # col 3: canal do 1° evento
    ax = axes[row, 3]
    cf = r['ch_first']
    ax.bar([2, 4], [(cf==2).sum(), (cf==4).sum()],
           color=['steelblue','darkorange'], alpha=0.8)
    ax.set_title(f'{lbl} — canal do 1° evento')
    ax.set_xlabel('canal'); ax.set_ylabel('contagem'); ax.grid(alpha=0.3)

plt.tight_layout()
saida = f'{OUT_DIR}/estrutura_trial.png'
plt.savefig(saida, dpi=150, bbox_inches='tight')
plt.show()
print(f"\nGráfico salvo: {saida}")

print("\n" + "="*65)
print("INTERPRETAÇÃO")
print("="*65)
print("""
Se n_both ≈ 50% dos trials → em metade dos trials temos 2 eventos.
Isso significa que ch=2 e ch=4 NÃO são mutuamente exclusivos:
  - Um canal é o SETTING do RNG (dispara em todo trial)
  - Outro canal é o DETECTOR do fóton (~25% dos trials)

Se Δ entre os dois eventos é constante (ex: sempre ~14 bins):
  → O setting chega sempre 14 bins antes/depois do fóton
  → Podemos filtrar pelo Δ para separar os dois

Se n_both ≈ 25% e n_s0 ≈ 25% e n_s1 ≈ 25%:
  → ch=2 e ch=4 são mutuamente exclusivos por trial
  → Cada trial tem no máximo um evento de ch=2 ou ch=4
  → O evento presente É o fóton (setting inferido pelo canal)
""")
