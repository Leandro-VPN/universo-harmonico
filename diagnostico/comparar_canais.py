# ══════════════════════════════════════════════════════════════════════
# ║  DIAGNÓSTICO — COMPARAÇÃO DE CANAIS POR ARQUIVO                     ║
# ║  Mapeia canal→período e tt_rel para identificar sync e detectores   ║
# ║  Leandro Miguel — ORCID: 0009-0003-4655-472X                         ║
# ║  github.com/Leandro-VPN/universo-harmonico                           ║
# ══════════════════════════════════════════════════════════════════════


import numpy as np
import struct, os
from pathlib import Path

# ── Coloque aqui todos os arquivos .dat que você tem ─────────────────────────
# Edite os caminhos conforme necessário
ARQUIVOS = {
    # Par principal (já analisado)
    'alice_19_45': r"d:\BellQM\19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat",
    'bob_19_44':   r"d:\BellQM\19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat",

    # Par 21_15 (200ns delay) — se tiver baixado
    # 'alice_21_15': r"d:\BellQM\21_15_CH_pockel_100kHz.run.200nsadditiondelay_lightconeshift.alice.dat",
    # 'bob_21_15':   r"d:\BellQM\21_15_CH_pockel_100kHz.run.200nsadditiondelay_lightconeshift.bob.dat",

    # Par XOR clássico — se tiver
    # 'alice_xor':   r"d:\BellQM\23_55_CH_pockel_100kHz.run.ClassicalRNGXOR.alice.dat",
    # 'bob_xor':     r"d:\BellQM\23_55_CH_pockel_100kHz.run.ClassicalRNGXOR.bob.dat",
}

BIN_PS   = 78.125
N_REG    = 2_000_000   # 2M registros por arquivo (~48 MB)

print("="*70)
print("MAPEAMENTO DE CANAIS POR ARQUIVO")
print(f"Lendo {N_REG:,} registros por arquivo")
print("="*70)

def analisar_arquivo(nome, filepath):
    if not Path(filepath).exists():
        print(f"\n  {nome}: arquivo não encontrado")
        return None

    fsize = Path(filepath).stat().st_size
    with open(filepath, 'rb') as f:
        raw = f.read(N_REG * 24)
    n = len(raw) // 24
    arr = np.frombuffer(raw, dtype=np.uint64).reshape(n, 3)
    ch  = arr[:, 0].astype(np.int64)
    tt  = arr[:, 1].astype(np.int64)
    del arr, raw

    print(f"\n{'─'*60}")
    print(f"  {nome}  ({fsize/1e9:.2f} GB)")
    print(f"  {n:,} registros analisados")

    # Período por canal
    resultado = {}
    print(f"\n  {'Canal':>6} {'N':>10} {'%':>6}  {'Período (bins)':>15}  {'µs':>10}  {'Tipo':>15}")
    print(f"  {'-'*65}")

    sync_canal = None
    click_canais = []

    for c in sorted(np.unique(ch)):
        ni = int((ch==c).sum())
        pct = ni/n*100

        # Período mediano
        tt_c = tt[ch==c]
        if len(tt_c) >= 10:
            diffs = np.diff(tt_c[:10000])
            med = np.median(diffs)
            # filtrar outliers
            diffs_ok = diffs[np.abs(diffs - med) < med*0.2]
            periodo = float(np.median(diffs_ok)) if len(diffs_ok) > 0 else med
            us = periodo * BIN_PS / 1000
        else:
            periodo = 0; us = 0

        # Tipo
        if 120_000 < periodo < 140_000:
            tipo = '← SYNC 100kHz'
            sync_canal = int(c)
        elif 240_000 < periodo < 280_000:
            tipo = '← 2× sync'
        elif pct < 5 and ni > 100:
            tipo = '← GPS/misc'
        elif 20 < pct < 35:
            tipo = '← clique? (25%)'
            click_canais.append(int(c))
        else:
            tipo = ''

        print(f"  {c:>6} {ni:>10,} {pct:>6.2f}%  {periodo:>15.0f}  {us:>10.3f}  {tipo}")

        resultado[int(c)] = {
            'n': ni, 'pct': pct, 'periodo': periodo, 'us': us
        }

    # tt_rel por canal (relativo ao sync)
    if sync_canal is not None:
        print(f"\n  tt_rel por canal (relativo ao sync ch={sync_canal}):")
        sync_tt = None
        rels = {}
        for i in range(min(n, 500_000)):
            c_i = int(ch[i]); t_i = int(tt[i])
            if c_i == sync_canal:
                sync_tt = t_i
            elif sync_tt is not None:
                rel = t_i - sync_tt
                if 0 < rel < 200_000:
                    if c_i not in rels: rels[c_i] = []
                    if len(rels[c_i]) < 2000:
                        rels[c_i].append(rel)

        for c_i in sorted(rels):
            arr_r = np.array(rels[c_i])
            if len(arr_r) < 5: continue
            med = np.median(arr_r)
            mn, mx = arr_r.min(), arr_r.max()
            # Janelas esperadas: Alice ~90-111 bins, Bob ~347-367 bins
            in_alice = 85 <= mn and mx <= 115
            in_bob   = 343 <= mn and mx <= 373
            janela = '← Alice click!' if in_alice else ('← Bob click!' if in_bob else '')
            print(f"    ch={c_i:3d}: N={len(arr_r):>5,}  "
                  f"[{mn:.0f}..{mx:.0f}]  med={med:.1f}  {janela}")

    return {'sync': sync_canal, 'clicks': click_canais, 'canais': resultado}

# ── Rodar para todos os arquivos ─────────────────────────────────────────────
resultados = {}
for nome, filepath in ARQUIVOS.items():
    r = analisar_arquivo(nome, filepath)
    if r: resultados[nome] = r

# ── Resumo ───────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("RESUMO — Pares compatíveis")
print("="*70)
print("""
Um par Alice/Bob é compatível para o teste se:
  1. Ambos têm sync no mesmo canal (ou canal diferente mas período igual)
  2. Ambos têm canal de clique com tt_rel na janela correta
     Alice: bins 85-115  |  Bob: bins 343-373
  3. Taxa de clique ~25% em ambos

Resultado atual:
""")

for nome, r in resultados.items():
    s = r.get('sync')
    c = r.get('clicks', [])
    print(f"  {nome}: sync=ch{s}, clicks={c}")

print("""
Se Alice 19_45 tem sync=ch2 mas Bob 19_44 tem sync=ch6,
os arquivos têm estruturas diferentes — o que é incomum.

Próximo passo recomendado:
  1. Rodar este script com os outros pares disponíveis
  2. Procurar um par onde Alice e Bob têm o mesmo canal de sync
  3. Se nenhum par funcionar, usar o HDF5 para tt_rel de Alice
     e .dat de Bob para tt_rel de Bob — são dados do mesmo experimento
""")
