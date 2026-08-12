#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verificacao_cxb.py -- reproduz TODOS os numeros de
protocols/protocolo_discriminante_CxB_v1.md (Secoes 2, 4 e 5).

NAO gera numeros novos: recalcula os publicados e CONFERE contra o .md.
Divergencias sao reportadas, nunca corrigidas em silencio.

Uso:  python verificacao_cxb.py
Exit code 0 se tudo confere, 1 se ha divergencia.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np

DOC = Path(__file__).resolve().parents[3] / "protocols" / "protocolo_discriminante_CxB_v1.md"

M_ALVO = 0.03
SEED = 2026
N_REP = 200_000   # mesma escala declarada na Secao 5 do protocolo
ALPHA = 0.05


# --------------------------------------------------------------------------
# Funcoes centrais
# --------------------------------------------------------------------------
def n_geom_deg(theta_rel_deg: float) -> float:
    """D1: n_geom(theta) = 2/(1-cos theta), theta em graus."""
    return 2.0 / (1.0 - math.cos(math.radians(theta_rel_deg)))


def razao_tangentes(theta_a_deg: float, theta_b_deg: float) -> float:
    """AII: r = tan(theta_A/2)/tan(theta_B/2)."""
    return math.tan(math.radians(theta_a_deg) / 2.0) / math.tan(math.radians(theta_b_deg) / 2.0)


def n_min_trials(n_alvo: float, m: float = M_ALVO) -> int:
    """Poder estatistico: N_min ~ 100 * n / m^2."""
    return round(100.0 * n_alvo / (m * m))


# --------------------------------------------------------------------------
# Parsing do .md
# --------------------------------------------------------------------------
def _f(s: str) -> float:
    return float(s.strip().replace("°", "").replace(",", ".").replace("*", ""))


def parse_secao2(doc: Path):
    txt = doc.read_text(encoding="utf-8")
    m = re.search(r"##\s*2\.\s*Configura(.*?)(?=\n##\s*3\.)", txt, re.S)
    if not m:
        raise SystemExit("ERRO: Secao 2 nao encontrada")
    cfgs = []
    for ln in m.group(1).splitlines():
        ln = ln.strip()
        if not ln.startswith("|"):
            continue
        c = [x.strip() for x in ln.strip("|").split("|")]
        if len(c) < 8 or not c[0].isdigit():
            continue
        mr = re.match(r"(\d+)\s*:\s*(\d+)", c[2])
        if not mr:      # linha de controle (razao irracional)
            continue
        ang_bloch = re.findall(r"([\d,]+)\s*°", c[4])
        ang_lab = re.findall(r"([\d,]+)\s*°", c[5])
        cfgs.append({
            "idx": int(c[0]), "nome": c[1],
            "p": int(mr.group(1)), "q": int(mr.group(2)),
            "n_mod": _f(c[3]),
            "bloch": (_f(ang_bloch[0]), _f(ang_bloch[1])),
            "lab": (_f(ang_lab[0]), _f(ang_lab[1])),
            "n_geom_pub": _f(c[6]),
            "prediz": "C" if "C" in c[7] else ("B" if "B" in c[7] else "?"),
        })
    return cfgs


def parse_tabela_N(doc: Path):
    txt = doc.read_text(encoding="utf-8")
    m = re.search(r"##\s*4\.(.*?)(?=\n##\s*5\.)", txt, re.S)
    linhas = []
    for ln in m.group(1).splitlines():
        ln = ln.strip()
        if not ln.startswith("|") or "m_alvo" in ln or set(ln) <= set("|-: "):
            continue
        c = [x.strip() for x in ln.strip("|").split("|")]
        if len(c) < 3:
            continue
        try:
            m_alvo = _f(c[0].split("(")[0])
            n_tot = int(c[2].replace(".", "").replace("*", "").strip())
        except ValueError:
            continue
        linhas.append((m_alvo, c[1], n_tot))
    return linhas


def parse_fp_secao5(doc: Path):
    """Extrai P(3^4^5), P(3^4 apenas) e P(qualquer par) da Secao 5."""
    txt = doc.read_text(encoding="utf-8")
    out = {}
    m = re.search(r"P\(3∧4∧5\)\s*=\s*([\d,]+)%", txt)
    if m:
        out["345"] = _f(m.group(1))
    m = re.search(r"P\(3∧4 apenas\)\s*=\s*([\d,]+)%", txt)
    if m:
        out["34"] = _f(m.group(1))
    m = re.search(r'P\("qualquer par\s*\n?de 3"[^)]*\)\s*=\s*([\d,]+)%', txt)
    if not m:
        m = re.search(r'qualquer par[^=]*=\s*([\d,]+)%', txt.replace("\n", " "))
    if m:
        out["par"] = _f(m.group(1))
    return out


# --------------------------------------------------------------------------
def main() -> int:
    print("=" * 76)
    print("VERIFICACAO C x B -- protocolo_discriminante_CxB_v1.md")
    print("=" * 76)
    print(f"documento: {DOC}")
    if not DOC.exists():
        print("ERRO: documento nao encontrado.")
        return 1
    div: list[str] = []

    # ---- Secao 2: os 5 pares angulares ---------------------------------------
    cfgs = parse_secao2(DOC)
    print(f"\n[1] Secao 2 -- {len(cfgs)} configuracoes com razao racional (esperado 5)")
    if len(cfgs) != 5:
        div.append(f"Secao 2: {len(cfgs)} configs parseadas, esperado 5")

    for cf in cfgs:
        tA, tB = cf["bloch"]
        p, q = cf["p"], cf["q"]
        r_calc = razao_tangentes(tA, tB)
        r_alvo = p / q
        err_r = abs(r_calc - r_alvo) / r_alvo * 100.0
        theta_rel = tA - tB
        ngeom_calc = n_geom_deg(theta_rel)
        nmod_calc = p * q
        soma = tA + tB
        lab_ok = (abs(cf["lab"][0] - tA / 2) < 6e-3) and (abs(cf["lab"][1] - tB / 2) < 6e-3)

        print(f"\n    [{cf['idx']}] {cf['nome']}  ({p}:{q})")
        print(f"        Bloch=({tA}°, {tB}°)  lab=({cf['lab'][0]}°, {cf['lab'][1]}°)  "
              f"lab==Bloch/2 (D6): {'OK' if lab_ok else 'DIVERGE'}")
        print(f"        razao tan: calc={r_calc:.6f}  alvo={r_alvo:.6f}  erro={err_r:.5f}%  "
              f"{'OK' if err_r < 0.002 else 'DIVERGE'}   (doc: erro < 0,002%)")
        print(f"        n_mod=p*q={nmod_calc} (doc {cf['n_mod']:.0f})  "
              f"{'OK' if nmod_calc == cf['n_mod'] else 'DIVERGE'}")
        print(f"        theta_rel={theta_rel:.3f}°  n_geom calc={ngeom_calc:.4f}  "
              f"doc={cf['n_geom_pub']:.2f}  "
              f"{'OK' if abs(ngeom_calc - cf['n_geom_pub']) / cf['n_geom_pub'] < 2e-3 else 'DIVERGE'}")
        if "simétrico" in cf["nome"] or "simetrico" in cf["nome"]:
            print(f"        simetria theta_A+theta_B = {soma:.3f}°  "
                  f"{'OK' if abs(soma - 180.0) < 1e-2 else 'DIVERGE'}  (D3': 180°)")
            if abs(soma - 180.0) >= 1e-2:
                div.append(f"config {cf['idx']}: soma={soma:.3f}, esperado 180")
        if "D3-exato" in cf["nome"]:
            print(f"        D3 exato: n_geom({theta_rel:.3f}°)={ngeom_calc:.4f} vs n_mod={nmod_calc}  "
                  f"{'OK' if abs(ngeom_calc - nmod_calc) / nmod_calc < 2e-3 else 'DIVERGE'}")
            if abs(ngeom_calc - nmod_calc) / nmod_calc >= 2e-3:
                div.append(f"config {cf['idx']}: D3 n_geom={ngeom_calc:.4f} vs n_mod={nmod_calc}")
        if err_r >= 0.002:
            div.append(f"config {cf['idx']}: erro razao {err_r:.5f}% >= 0,002%")
        if not lab_ok:
            div.append(f"config {cf['idx']}: angulos de laboratorio != Bloch/2")
        if nmod_calc != cf["n_mod"]:
            div.append(f"config {cf['idx']}: n_mod calc={nmod_calc} pub={cf['n_mod']}")
        if abs(ngeom_calc - cf["n_geom_pub"]) / cf["n_geom_pub"] >= 2e-3:
            div.append(f"config {cf['idx']}: n_geom calc={ngeom_calc:.4f} pub={cf['n_geom_pub']}")

    # ---- Secao 4: tabela de N -------------------------------------------------
    print("\n[2] Secao 4 -- tabela de N (m_alvo=0,03 e escala 1/m^2)")
    # N por config usa o n_geom COMO IMPRESSO no documento (2 casas), nao o exato
    Ns = [n_min_trials(cf["n_geom_pub"]) for cf in cfgs]
    total_003 = sum(Ns)
    for cf, n in zip(cfgs, Ns):
        print(f"    config {cf['idx']}: n_geom={cf['n_geom_pub']:>6.2f} -> N={n:>10,}")
    print(f"    TOTAL (m=0,03) = {total_003:,}")

    tabN = parse_tabela_N(DOC)
    base_pub = None
    for (m_alvo, fator, n_tot) in tabN:
        if abs(m_alvo - 0.03) < 1e-9:
            base_pub = n_tot
    if base_pub is None:
        div.append("Secao 4: linha m=0,03 nao encontrada")
    else:
        ok = (total_003 == base_pub)
        print(f"    documento (m=0,03): {base_pub:,}   {'OK' if ok else 'DIVERGE'}")
        if not ok:
            div.append(f"Secao 4: N total calc={total_003} pub={base_pub}")

    print("    escala 1/m^2:")
    for (m_alvo, fator, n_tot) in tabN:
        esperado = total_003 * (M_ALVO / m_alvo) ** 2
        # documento arredonda para 4 algarismos significativos
        rel = abs(esperado - n_tot) / n_tot
        ok = rel < 2e-4
        print(f"      m={m_alvo:<6} fator doc={fator:<8} N doc={n_tot:>12,}  "
              f"calc={esperado:>14,.0f}  dif rel={rel:.2e}  {'OK' if ok else 'DIVERGE'}")
        if not ok:
            div.append(f"Secao 4: m={m_alvo} N pub={n_tot} calc={esperado:.0f}")

    # ---- Secao 5: taxa de falso-positivo conjunto -----------------------------
    print(f"\n[3] Secao 5 -- falso-positivo conjunto sob nula pura "
          f"({N_REP:,} repeticoes, alpha={ALPHA})")
    rng = np.random.default_rng(SEED)
    sig = rng.random((N_REP, 3)) < ALPHA        # colunas = configs 3, 4, 5
    c3, c4, c5 = sig[:, 0], sig[:, 1], sig[:, 2]
    p345 = float(np.mean(c3 & c4 & c5)) * 100
    p34 = float(np.mean(c3 & c4 & ~c5)) * 100
    ppar = float(np.mean((c3 & c4) | (c3 & c5) | (c4 & c5))) * 100
    # esperado analitico
    a = ALPHA
    e345, e34, epar = a**3 * 100, a**2 * (1 - a) * 100, (3 * a**2 * (1 - a) + a**3) * 100

    pub = parse_fp_secao5(DOC)
    print(f"    {'estatistica':<26} {'doc':>8} {'simulado':>10} {'analitico':>10}  status")
    for chave, nome, sim, ana in [
        ("345", "P(3^4^5)", p345, e345),
        ("34", "P(3^4 apenas)", p34, e34),
        ("par", "P(qualquer par)", ppar, epar),
    ]:
        pv = pub.get(chave)
        # tolerancia = 3 sigma de Monte Carlo sobre N_REP repeticoes
        sigma = math.sqrt(max(ana, 1e-9) / 100 * (1 - ana / 100) / N_REP) * 100
        ok = (pv is not None) and abs(pv - ana) <= 4 * sigma
        print(f"    {nome:<26} {('-' if pv is None else f'{pv:.3f}%'):>8} "
              f"{sim:>9.3f}% {ana:>9.4f}%  "
              f"{'OK (dentro de 4 sigma MC)' if ok else 'DIVERGE'}  [4s={4*sigma:.4f}%]")
        if pv is None:
            div.append(f"Secao 5: valor '{nome}' nao encontrado no documento")
        elif not ok:
            div.append(f"Secao 5: {nome} pub={pv}% vs analitico={ana:.4f}% (fora de 4 sigma MC)")

    print("\n    Nota: os numeros da Secao 5 sao saidas de SIMULACAO (200k repeticoes),")
    print("    portanto tem ruido de Monte Carlo; a conferencia e contra o valor")
    print("    analitico dentro da barra de erro esperada, nao igualdade exata.")

    print("\n" + "=" * 76)
    if div:
        print(f"RESULTADO: {len(div)} DIVERGENCIA(S) ENCONTRADA(S)")
        for d in div:
            print("  - " + d)
        print("\nNENHUM lado foi corrigido automaticamente (regra do pre-registro).")
        return 1
    print("RESULTADO: todos os numeros conferem com o documento.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
