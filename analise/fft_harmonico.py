# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  FFT DA SÉRIE TEMPORAL DE COINCIDÊNCIAS                                  ║
# ║  Modelo: razões de afinação justa (p:q inteiros) → período MMC(p,q)     ║
# ║                                                                          ║
# ║  Predição: se as partículas seguem séries harmônicas naturais,           ║
# ║  a taxa de coincidência por bloco de trials deve mostrar pico de FFT    ║
# ║  no período MMC(p,q) correspondente ao ângulo relativo.                 ║
# ║                                                                          ║
# ║  Controle: mesmo teste com n=MMC±1 e MMC±2 — sinal real só no MMC.     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

import zipfile, numpy as np, gc, os, time
import matplotlib.pyplot as plt
from scipy.stats import chi2, kstest
from math import gcd
from google.colab import drive

drive.mount('/content/drive')

DRIVE   = "/content/drive/MyDrive/Colab Notebooks/Dados Bell"
ZIP_A   = f"{DRIVE}/19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat.zip"
ZIP_B   = f"{DRIVE}/19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat.zip"
OUT_DIR = f"{DRIVE}/fft_harmonico"
os.makedirs(OUT_DIR, exist_ok=True)

BIN_PS       = 78.125
CANAL_SYNC   = 6
CANAIS_CLICK = [2, 4]
SENTINEL_TT  = np.uint32(0xFFFFFFFF)
SENTINEL_SET = np.int8(-1)
CHUNK_BYTES  = 24 * 5_000_000
N_TRIALS     = 10_000_000

# ── razão de afinação justa mais próxima de tan(θ_A/2)/tan(θ_B/2) ──────────
# Para os combos do NIST:
#   XX: θ_A=-30°, θ_B=0°  → tan(0°)=0  → razão indefinida
#   XY: θ_A=-30°, θ_B=+30° → tan(15°)/tan(15°) = 1:1 (uníssono)
#   YX: θ_A=+30°, θ_B=0°  → razão indefinida
#   YY: θ_A=+30°, θ_B=+30° → 1:1 (uníssono)
#
# Para uníssono (1:1): MMC=1 → encontro em todo trial → sem modulação detectável
# O teste mais informativo é buscar a razão efetiva via ângulo RELATIVO θ_rel
# e mapear para a fração justa mais próxima.

def mmc(a, b):
    return a * b // gcd(a, b)

def fracao_justa_proxima(r, max_denom=20):
    """Fração p/q com denominador ≤ max_denom mais próxima de r."""
    melhor_p, melhor_q, melhor_err = 1, 1, float('inf')
    for q in range(1, max_denom+1):
        p = round(r * q)
        if p < 1: p = 1
        err = abs(p/q - r)
        if err < melhor_err:
            melhor_err = err
            melhor_p, melhor_q = p, q
    return melhor_p, melhor_q, melhor_err

# ângulos reais do NIST (em graus)
ANGULOS = {'X': -30.0, 'Y': +30.0, 'REF': 0.0}

# combos e seus ângulos relativos
COMBOS_INFO = {
    'XX': {'tA': -30.0, 'tB':   0.0, 'label': 'Alice X, Bob ref'},
    'XY': {'tA': -30.0, 'tB': +30.0, 'label': 'Alice X, Bob Y'},
    'YX': {'tA': +30.0, 'tB':   0.0, 'label': 'Alice Y, Bob ref'},
    'YY': {'tA': +30.0, 'tB': +30.0, 'label': 'Alice Y, Bob Y'},
}

print("="*70)
print("MAPEAMENTO ÂNGULO → RAZÃO JUSTA → MMC")
print("="*70)
print(f"\n{'Combo':5} {'θ_rel':>8} {'razão ω':>10} {'p:q justa':>10} "
      f"{'MMC':>6} {'erro %':>8} {'θ_Bell(MMC)':>12}")
print("-"*65)

for nome, info in COMBOS_INFO.items():
    tA, tB = np.deg2rad(info['tA']), np.deg2rad(info['tB'])
    wA = np.tan(tA/2) if abs(info['tA']) > 0.01 else 1e-9
    wB = np.tan(tB/2) if abs(info['tB']) > 0.01 else 1e-9

    theta_rel = abs(info['tA'] - info['tB'])

    if abs(wA) < 1e-6 or abs(wB) < 1e-6:
        r = float('inf') if abs(wA) < 1e-6 else 0.0
        info['mmc'] = None
        info['pq'] = None
        print(f"{nome:5} {theta_rel:>8.1f}° {'∞ ou 0':>10} {'—':>10} "
              f"{'—':>6} {'—':>8} {'—':>12}")
    else:
        r = abs(wB / wA)
        p, q, err = fracao_justa_proxima(r, max_denom=30)
        g = gcd(p, q); p //= g; q //= g
        m = mmc(p, q)
        theta_bell = np.degrees(np.arccos(1 - 2/m)) if m >= 2 else 0.0
        info['mmc'] = m
        info['pq']  = (p, q)
        print(f"{nome:5} {theta_rel:>8.1f}° {r:>10.5f} {p:>4}:{q:<5} "
              f"{m:>6} {100*err/r:>8.3f}% {theta_bell:>11.3f}°")

# ── Intervalo de afinação justa para θ_rel de Bell ──────────────────────────
print("\n" + "="*70)
print("INTERVALOS JUSTOS PARA ÂNGULOS RELATIVOS DOS COMBOS")
print("(usando θ_rel como ângulo de Bell → n = round(2/(1-cos(θ_rel))))")
print("="*70)
print(f"\n{'Combo':5} {'θ_rel':>8} {'n=MMC':>8} {'intervalo justo':>20} {'p:q':>8}")
print("-"*55)

intervalos_justos = {
    2: ('oitava', '2:1'),
    3: ('quinta justa', '3:2'),
    4: ('quarta justa', '4:3'),
    5: ('terça menor (6:5)', '6:5'),
    6: ('terça maior (5:4)', '5:4'),
    7: ('sétima harm.', '7:4'),
    9: ('tom maior', '9:8'),
    10: ('tom menor', '10:9'),
    14: ('semitom diat.', '15:14'),
    21: ('semitom crom.', '22:21'),
    35: ('trítono 7:5', '7:5'),
    70: ('trítono 10:7', '10:7'),
}

for nome, info in COMBOS_INFO.items():
    theta_rel = abs(info['tA'] - info['tB'])
    if theta_rel < 0.01:
        n = 1; interv = 'uníssono'; pq_str = '1:1'
    else:
        t = np.deg2rad(theta_rel)
        n = round(2 / (1 - np.cos(t)))
        interv, pq_str = intervalos_justos.get(n, (f'n={n}', '?'))
    info['n_bell'] = n
    print(f"{nome:5} {theta_rel:>8.1f}° {n:>8} {interv:>20} {pq_str:>8}")


# ── Extração ─────────────────────────────────────────────────────────────────
def extrair(zip_path, n_trials):
    tt_rel  = np.full(n_trials, SENTINEL_TT,  dtype=np.uint32)
    setting = np.full(n_trials, SENTINEL_SET, dtype=np.int8)
    trial = -1; sync_tt = 0

    with zipfile.ZipFile(zip_path, 'r') as zf:
        fname = zf.namelist()[0]
        with zf.open(fname) as f:
            while trial < n_trials - 1:
                raw = f.read(CHUNK_BYTES)
                if not raw: break
                n = len(raw) // 24
                data = np.frombuffer(raw, dtype=np.uint64,
                                     count=n*3).reshape(n, 3)
                ch = data[:, 0].astype(np.int64)
                tt = data[:, 1].astype(np.int64)
                del data
                rel = np.flatnonzero(
                    (ch == CANAL_SYNC) | np.isin(ch, CANAIS_CLICK))
                parou = False
                for i in rel:
                    if ch[i] == CANAL_SYNC:
                        trial += 1
                        sync_tt = tt[i]
                        if trial >= n_trials - 1:
                            parou = True; break
                    elif (trial >= 0 and trial < n_trials
                          and tt_rel[trial] == SENTINEL_TT):
                        d = tt[i] - sync_tt
                        if 0 <= d < SENTINEL_TT:
                            tt_rel[trial]  = np.uint32(d)
                            setting[trial] = np.int8(0 if ch[i] == 2 else 1)
                del ch, tt, rel
                if parou: break

    return tt_rel[:trial+1], setting[:trial+1]

print(f"\n{'='*70}")
print(f"EXTRAINDO {N_TRIALS:,} TRIALS")
print(f"{'='*70}")
t0 = time.time()
tr_a, st_a = extrair(ZIP_A, N_TRIALS)
print(f"Alice: {np.sum(tr_a != SENTINEL_TT):,} cliques  ({time.time()-t0:.1f}s)")
t0 = time.time()
tr_b, st_b = extrair(ZIP_B, N_TRIALS)
print(f"Bob:   {np.sum(tr_b != SENTINEL_TT):,} cliques  ({time.time()-t0:.1f}s)")

n = min(len(tr_a), len(tr_b))
tr_a, st_a = tr_a[:n], st_a[:n]
tr_b, st_b = tr_b[:n], st_b[:n]


# ── FFT da série temporal de coincidências ────────────────────────────────────
print(f"\n{'='*70}")
print("FFT DA SÉRIE TEMPORAL DE COINCIDÊNCIAS")
print("Buscando periodicidade no período MMC da razão justa")
print(f"{'='*70}")

# Tamanho do bloco para calcular taxa de coincidência
# Deve ser >> MMC para detectar oscilação, mas não tão grande que perca resolução
# Usar MMC como unidade natural: bloco = 1 trial (série binária 0/1)
# FFT na série binária diretamente

COMBOS_SETTINGS = {
    'XX': (1, 1),
    'XY': (1, 2),
    'YX': (2, 1),
    'YY': (2, 2),
}

fig, axes = plt.subplots(4, 3, figsize=(18, 22))
fig.suptitle('FFT da série temporal de coincidências — modelo de afinação justa\n'
             'Pico no período MMC(p,q) indicaria modulação harmônica natural',
             fontsize=12, fontweight='bold')

resultados_fft = {}

for row, (nome, (sa, sb)) in enumerate(COMBOS_SETTINGS.items()):
    info   = COMBOS_INFO[nome]
    n_bell = info['n_bell']

    mask = ((tr_a != SENTINEL_TT) & (tr_b != SENTINEL_TT) &
            (st_a == sa-1) & (st_b == sb-1))
    serie_binaria = mask.astype(np.float32)   # 1 = coincidência, 0 = não
    N_total = len(serie_binaria)
    N_coinc = int(mask.sum())

    print(f"\n--- {nome} ({info['label']}) ---")
    print(f"  N_trials={N_total:,}  N_coinc={N_coinc:,}  "
          f"taxa={N_coinc/N_total*100:.2f}%  n_bell={n_bell}")

    # ── FFT direta na série binária ──────────────────────────────────────────
    # Remover média (DC) antes da FFT
    serie_centrada = serie_binaria - serie_binaria.mean()
    fft_vals = np.fft.rfft(serie_centrada)
    freqs    = np.fft.rfftfreq(N_total)   # em ciclos por trial
    potencia = np.abs(fft_vals)**2

    # Períodos de interesse: n_bell e vizinhos ±5
    periodos_interesse = list(range(max(2, n_bell-5), n_bell+6))

    # Normalizar pela potência média (excluindo DC e muito baixa freq)
    idx_dc  = 0
    pot_med = np.median(potencia[10:])   # mediana do fundo

    print(f"  Potência mediana do fundo: {pot_med:.2e}")
    print(f"\n  {'Período (trials)':>18} {'Frequência':>12} {'Potência':>12} "
          f"{'S/N':>8} {'nota musical':>15}")
    print("  " + "-"*70)

    for p_int in periodos_interesse:
        # frequência correspondente ao período p_int
        f_target = 1.0 / p_int
        idx = np.argmin(np.abs(freqs - f_target))
        pot = potencia[idx]
        sn  = pot / pot_med if pot_med > 0 else 0
        nota = intervalos_justos.get(p_int, (f'n={p_int}', '?'))[0]
        flag = ' ← ALVO' if p_int == n_bell else (' ← VIZINHO' if abs(p_int-n_bell)<=2 else '')
        print(f"  {p_int:>18} {f_target:>12.6f} {pot:>12.3e} {sn:>8.2f}x  "
              f"{nota:>15}{flag}")

    resultados_fft[nome] = {
        'freqs': freqs, 'potencia': potencia,
        'n_bell': n_bell, 'N_coinc': N_coinc,
        'pot_med': pot_med
    }

    # ── col 1: espectro de potência FFT, escala log ──────────────────────────
    ax1 = axes[row, 0]
    # Mostrar apenas períodos 2 a 500 (frequências 1/500 a 1/2)
    mask_freq = (freqs >= 1/500) & (freqs <= 1/2)
    periodos_plot = 1.0 / freqs[mask_freq]
    pot_plot = potencia[mask_freq]

    ax1.semilogy(periodos_plot, pot_plot, 'b-', lw=0.5, alpha=0.6)
    ax1.axhline(pot_med, color='green', lw=1.5, linestyle='--',
                label=f'fundo mediano')

    # marcar n_bell
    if n_bell >= 2:
        f_alvo = 1.0 / n_bell
        idx_alvo = np.argmin(np.abs(freqs - f_alvo))
        pot_alvo = potencia[idx_alvo]
        ax1.axvline(n_bell, color='red', lw=2, alpha=0.8,
                    label=f'MMC={n_bell}')
        ax1.scatter([n_bell], [pot_alvo], color='red', s=80, zorder=5)

    ax1.set_title(f'{nome} — espectro FFT completo\nN_coinc={N_coinc:,}', fontsize=9)
    ax1.set_xlabel('período (trials)', fontsize=8)
    ax1.set_ylabel('potência', fontsize=8)
    ax1.set_xlim(2, 500)
    ax1.legend(fontsize=7)
    ax1.grid(alpha=0.3, which='both')

    # ── col 2: zoom ao redor do n_bell ──────────────────────────────────────
    ax2 = axes[row, 1]
    zoom_min = max(2, n_bell - 20)
    zoom_max = n_bell + 20
    mask_zoom = (freqs >= 1/zoom_max) & (freqs <= 1/zoom_min) & (freqs > 0)
    if mask_zoom.sum() > 0:
        p_zoom = 1.0 / freqs[mask_zoom]
        ax2.semilogy(p_zoom, potencia[mask_zoom], 'b-', lw=1, alpha=0.8)
        ax2.axhline(pot_med, color='green', lw=1.5, linestyle='--')
        if n_bell >= 2:
            ax2.axvline(n_bell, color='red', lw=2, alpha=0.9,
                        label=f'n_bell={n_bell}')
        # marcar vizinhos
        for dv in [-2,-1,1,2]:
            nv = n_bell + dv
            if zoom_min <= nv <= zoom_max:
                ax2.axvline(nv, color='orange', lw=1, alpha=0.5,
                            linestyle='--')
        ax2.set_xlim(zoom_min, zoom_max)
        ax2.set_title(f'{nome} — zoom ±20 ao redor de n={n_bell}', fontsize=9)
        ax2.set_xlabel('período (trials)', fontsize=8)
        ax2.set_ylabel('potência', fontsize=8)
        ax2.legend(fontsize=7)
        ax2.grid(alpha=0.3, which='both')

    # ── col 3: taxa de coincidência por bloco de n_bell trials ──────────────
    ax3 = axes[row, 2]
    if n_bell >= 2 and N_total // n_bell > 10:
        n_blocos = N_total // n_bell
        taxa_blocos = np.array([
            serie_binaria[i*n_bell:(i+1)*n_bell].mean()
            for i in range(n_blocos)
        ])
        taxa_media = taxa_blocos.mean()
        taxa_std   = taxa_blocos.std()

        ax3.plot(taxa_blocos, 'b-', lw=0.8, alpha=0.7)
        ax3.axhline(taxa_media, color='red', lw=1.5, linestyle='--',
                    label=f'média={taxa_media:.4f}')
        ax3.fill_between(range(len(taxa_blocos)),
                         taxa_media - taxa_std,
                         taxa_media + taxa_std,
                         alpha=0.2, color='red', label=f'±1σ={taxa_std:.4f}')
        ax3.set_title(f'{nome} — taxa coinc. por bloco de {n_bell} trials\n'
                      f'std/média = {taxa_std/taxa_media*100:.2f}%', fontsize=9)
        ax3.set_xlabel('bloco', fontsize=8)
        ax3.set_ylabel('taxa de coincidência', fontsize=8)
        ax3.legend(fontsize=7)
        ax3.grid(alpha=0.3)

        # FFT da série de taxas por bloco
        fft_taxas = np.abs(np.fft.rfft(taxa_blocos - taxa_media))**2
        freq_taxas = np.fft.rfftfreq(len(taxa_blocos))
        print(f"\n  FFT da série de taxas por bloco (cada bloco = {n_bell} trials):")
        pot_med_taxas = np.median(fft_taxas[5:])
        idx_pico = np.argmax(fft_taxas[1:]) + 1
        f_pico = freq_taxas[idx_pico]
        periodo_pico = round(1/f_pico) if f_pico > 0 else 0
        sn_pico = fft_taxas[idx_pico] / pot_med_taxas
        print(f"    Pico em período={periodo_pico} blocos (={periodo_pico*n_bell} trials)  "
              f"S/N={sn_pico:.1f}x")

    del mask, serie_binaria

plt.tight_layout()
saida = f'{OUT_DIR}/fft_serie_temporal.png'
plt.savefig(saida, dpi=150, bbox_inches='tight')
plt.show()
print(f"\nGráfico salvo em: {saida}")
gc.collect()


# ── Teste de Rayleigh com período MMC da afinação justa ──────────────────────
print(f"\n{'='*70}")
print("TESTE DE RAYLEIGH — t_i mod MMC_justa")
print("Fórmula: fase_i = (i mod MMC) / MMC  onde i = índice da coincidência")
print(f"{'='*70}")
print("""
Interpretação do modelo:
  Se as partículas emitem séries harmônicas com razão justa p:q,
  as coincidências deveriam ocorrer com concentração de fase
  relativa ao período MMC(p,q) — ou seja, não uniformemente
  distribuídas dentro de cada ciclo de MMC trials.

  Sinal real: V(MMC) >> V(MMC±1)  e  V(MMC) >> 1/sqrt(N)
""")

print(f"\n{'Combo':5} {'n_bell=MMC':>10} {'N_coinc':>10} {'V(MMC)':>10} "
      f"{'p(MMC)':>12} {'V(MMC-1)':>10} {'V(MMC+1)':>10} {'destaca?':>10}")
print("-"*80)

for nome, (sa, sb) in COMBOS_SETTINGS.items():
    info   = COMBOS_INFO[nome]
    n_bell = info['n_bell']
    if n_bell < 2:
        print(f"{nome:5} {'uníssono':>10} — sem período definido")
        continue

    mask = ((tr_a != SENTINEL_TT) & (tr_b != SENTINEL_TT) &
            (st_a == sa-1) & (st_b == sb-1))
    idx  = np.flatnonzero(mask)
    N    = idx.size
    t_i  = np.arange(N, dtype=np.int64)

    def rayleigh_V(t, k):
        phi = 2*np.pi*(t % k) / k
        C, S = np.mean(np.cos(phi)), np.mean(np.sin(phi))
        Z = 2*len(t)*(C**2+S**2)
        return np.sqrt(max(Z,0)/(2*len(t))), chi2.sf(Z, df=2)

    V_mmc,  p_mmc  = rayleigh_V(t_i, n_bell)
    V_m1,   _      = rayleigh_V(t_i, max(2, n_bell-1))
    V_p1,   _      = rayleigh_V(t_i, n_bell+1)

    viz_med = (V_m1 + V_p1) / 2
    destaca = V_mmc > 3 * viz_med if viz_med > 0 else False
    flag = '✓ SIM' if destaca else '✗ não'

    print(f"{nome:5} {n_bell:>10} {N:>10,} {V_mmc:>10.6f} "
          f"{p_mmc:>12.2e} {V_m1:>10.6f} {V_p1:>10.6f} {flag:>10}")
    del mask, idx, t_i

print(f"\n{'='*70}")
print("INTERPRETAÇÃO FINAL")
print(f"{'='*70}")
print("""
Com os combos XX/XY/YX/YY do NIST e ângulos −30°/0°/+30°:

  XY e YY: θ_rel = 60° → n_bell = 4 (quarta justa 4:3, MMC=12)
  XX e YX: θ_rel = 30° → n_bell = 14 (semitom diatônico 15:14, MMC=210)

  Se V(n_bell) >> V(vizinhos) → sinal do modelo harmônico justo
  Se V(n_bell) ≈ V(vizinhos) → distribuição uniforme → sem modulação

O experimento ideal para testar o trítono (7:5, MMC=35):
  Precisaria de ângulos com tan(θ_A/2)/tan(θ_B/2) = 7/5
  → θ_A ≈ 37.4°, θ_B ≈ 24.8°  (não disponíveis no NIST 2015)

Resultado deste teste + FFT → incluir na seção de resultados negativos
ou positivos do artigo, dependendo do que aparecer nos gráficos.
""")
