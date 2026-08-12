#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verificacao_a0.py -- reproduz TODOS os numeros de
protocols/axioma_0_identificabilidade_v1.md (Secoes 3, 4.1 e 5).

NAO gera numeros novos: recalcula os publicados e CONFERE contra o .md.
Divergencias sao reportadas, nunca corrigidas em silencio.

Uso:  python verificacao_a0.py
Exit code 0 se tudo confere, 1 se ha divergencia.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np

DOC = Path(__file__).resolve().parents[3] / "protocols" / "axioma_0_identificabilidade_v1.md"
DOC_CXB = Path(__file__).resolve().parents[3] / "protocols" / "protocolo_discriminante_CxB_v1.md"

M_ALVO = 0.03
SEED = 2026
N_REP = 200_000
ALPHA = 0.05
DELTA_THETA = 21.787   # graus, declarado na Secao 3


def n_geom_deg(theta_deg: float) -> float:
    return 2.0 / (1.0 - math.cos(math.radians(theta_deg)))


def razao_tangentes(theta_a_deg: float, theta_b_deg: float) -> float:
    return math.tan(math.radians(theta_a_deg) / 2.0) / math.tan(math.radians(theta_b_deg) / 2.0)


def n_min_trials(n_alvo: float, m: float = M_ALVO) -> int:
    return round(100.0 * n_alvo / (m * m))


def _f(s: str) -> float:
    return float(s.strip().replace("°", "").replace(",", ".").replace("*", ""))


def parse_orientacoes(doc: Path):
    txt = doc.read_text(encoding="utf-8")
    orients = []
    for ln in txt.splitlines():
        ln = ln.strip()
        if not ln.startswith("|"):
            continue
        c = [x.strip() for x in ln.strip("|").split("|")]
        if len(c) < 5 or not c[0].isdigit():
            continue
        mr = re.match(r"(\d+)\s*:\s*(\d+)", c[3])
        if not mr:
            continue
        orients.append({
            "idx": int(c[0]),
            "theta_b": _f(c[1]), "theta_a": _f(c[2]),
            "p": int(mr.group(1)), "q": int(mr.group(2)),
            "n_mod_pub": _f(c[4]),
        })
    return orients


def parse_N(doc: Path):
    txt = doc.read_text(encoding="utf-8")
    pares = re.findall(r"orienta[çc][ãa]o\s*(\d)\s*\(n=(\d+)\):\s*([\d.]+)", txt)
    tot = re.search(r"total:\s*([\d.]+)", txt)
    out = {int(i): (int(n), int(v.replace(".", ""))) for (i, n, v) in pares}
    total = int(tot.group(1).replace(".", "")) if tot else None
    return out, total


def main() -> int:
    print("=" * 76)
    print("VERIFICACAO A0 -- axioma_0_identificabilidade_v1.md")
    print("=" * 76)
    print(f"documento: {DOC}")
    if not DOC.exists():
        print("ERRO: documento nao encontrado.")
        return 1
    div: list[str] = []

    # ---- Secao 3: as 3 orientacoes -------------------------------------------
    ors = parse_orientacoes(DOC)
    print(f"\n[1] Secao 3 -- {len(ors)} orientacoes (esperado 3), Delta_theta fixo = {DELTA_THETA}°")
    if len(ors) != 3:
        div.append(f"Secao 3: {len(ors)} orientacoes parseadas, esperado 3")

    for o in ors:
        tA, tB, p, q = o["theta_a"], o["theta_b"], o["p"], o["q"]
        dtheta = tA - tB
        r_calc = razao_tangentes(tA, tB)
        r_alvo = p / q
        err = abs(r_calc - r_alvo) / r_alvo * 100
        nmod_calc = p * q
        ok_d = abs(dtheta - DELTA_THETA) < 1e-3
        ok_r = err < 0.01
        ok_n = (nmod_calc == o["n_mod_pub"])
        print(f"\n    orientacao {o['idx']}: theta_B={tB}°  theta_A={tA}°   razao {p}:{q}")
        print(f"        Delta_theta = {dtheta:.4f}°  (fixo {DELTA_THETA}°)  {'OK' if ok_d else 'DIVERGE'}")
        print(f"        razao tan calc={r_calc:.6f} alvo={r_alvo:.6f} erro={err:.5f}%  "
              f"{'OK' if ok_r else 'DIVERGE'}")
        print(f"        n_mod = p*q = {nmod_calc}  (doc {o['n_mod_pub']:.0f})  "
              f"{'OK' if ok_n else 'DIVERGE'}")
        if not ok_d:
            div.append(f"orientacao {o['idx']}: Delta_theta={dtheta:.4f} vs {DELTA_THETA}")
        if not ok_r:
            div.append(f"orientacao {o['idx']}: erro de razao {err:.5f}%")
        if not ok_n:
            div.append(f"orientacao {o['idx']}: n_mod calc={nmod_calc} pub={o['n_mod_pub']}")

    # n_geom(Delta_theta) = 28 exato (a orientacao 2 satisfaz D3)
    ng = n_geom_deg(DELTA_THETA)
    m = re.search(r"n_geom\(21,787°\)\s*=\s*([\d,]+)", DOC.read_text(encoding="utf-8"))
    ng_pub = _f(m.group(1)) if m else None
    ok_ng = (ng_pub is not None) and abs(ng - ng_pub) < 5e-3
    print(f"\n    n_geom({DELTA_THETA}°) = {ng:.4f}   documento: {ng_pub}   "
          f"{'OK' if ok_ng else 'DIVERGE'}")
    print(f"    (a orientacao 2, razao 7:4, tem n_mod=28 -> satisfaz D3 exatamente)")
    if not ok_ng:
        div.append(f"n_geom({DELTA_THETA}) calc={ng:.4f} pub={ng_pub}")

    # ---- consistencia com o protocolo CxB (exigida pelo PROTOCOL_LOCK) -------
    print("\n[2] Consistencia cruzada com CxB (LOCK: orientacao 2 == config 3 do CxB)")
    if DOC_CXB.exists():
        txt = DOC_CXB.read_text(encoding="utf-8")
        mm = re.search(r"\|\s*3\s*\|\s*7:4 D3-exato \(ramo 1\)\s*\|[^|]*\|[^|]*\|\s*([\d,]+)°,\s*([\d,]+)°", txt)
        if mm:
            cxb_tA, cxb_tB = _f(mm.group(1)), _f(mm.group(2))
            o2 = [o for o in ors if o["idx"] == 2]
            if o2:
                same = abs(cxb_tA - o2[0]["theta_a"]) < 1e-3 and abs(cxb_tB - o2[0]["theta_b"]) < 1e-3
                print(f"    CxB config 3: ({cxb_tA}°, {cxb_tB}°)   A0 orientacao 2: "
                      f"({o2[0]['theta_a']}°, {o2[0]['theta_b']}°)   {'OK -- identicas' if same else 'DIVERGE'}")
                if not same:
                    div.append("orientacao 2 do A0 != config 3 do CxB (LOCK afirma que sao a mesma)")
        else:
            print("    AVISO: nao consegui parsear a config 3 do CxB para conferir")
    else:
        print("    AVISO: protocolo CxB nao encontrado; conferencia cruzada pulada")

    # ---- Secao 5: tabela de N ------------------------------------------------
    print("\n[3] Secao 5 -- N por orientacao (m_alvo=0,03; N ~ 100*n_mod/m^2)")
    Npub, total_pub = parse_N(DOC)
    total_calc = 0
    for o in ors:
        n_alvo = o["p"] * o["q"]
        N_calc = n_min_trials(n_alvo)
        total_calc += N_calc
        n_doc, N_doc = Npub.get(o["idx"], (None, None))
        ok = (N_doc == N_calc) and (n_doc == n_alvo)
        print(f"    orientacao {o['idx']} (n={n_alvo:>2}): calc={N_calc:>9,}   doc={N_doc:>9,}   "
              f"{'OK' if ok else 'DIVERGE'}")
        if not ok:
            div.append(f"Secao 5: orientacao {o['idx']} N calc={N_calc} pub={N_doc}")
    ok_t = (total_pub == total_calc)
    print(f"    TOTAL: calc={total_calc:,}   doc={total_pub:,}   {'OK' if ok_t else 'DIVERGE'}")
    if not ok_t:
        div.append(f"Secao 5: total calc={total_calc} pub={total_pub}")

    # ---- Secao 4.1: falso-positivo conjunto da matriz 3x3 --------------------
    print(f"\n[4] Secao 4.1 -- falso-positivo CONJUNTO da matriz 3x3 "
          f"({N_REP:,} repeticoes, alpha={ALPHA})")
    print("    Regra correta: UMA calibracao conjunta por permutacao sobre as 9 celulas;")
    print("    criterio primario = as 3 celulas-alvo (diagonal) excedem simultaneamente.")
    rng = np.random.default_rng(SEED)
    celulas = rng.random((N_REP, 3, 3)) < ALPHA        # 3 orientacoes x 3 alvos
    diag = celulas[:, 0, 0] & celulas[:, 1, 1] & celulas[:, 2, 2]
    p_diag = float(np.mean(diag)) * 100
    linha28 = celulas[:, 0, 1] & celulas[:, 1, 1] & celulas[:, 2, 1]
    p_linha = float(np.mean(linha28)) * 100
    esperado = ALPHA**3 * 100
    ingenuo = ALPHA**9 * 100

    txt = DOC.read_text(encoding="utf-8")
    mm = re.search(r"α³\s*≈\s*([\d,]+)%", txt)
    pub = _f(mm.group(1)) if mm else None
    sigma = math.sqrt(esperado / 100 * (1 - esperado / 100) / N_REP) * 100
    ok = (pub is not None) and abs(pub - esperado) < 1e-3
    print(f"    P(diagonal 3 celulas) simulado = {p_diag:.4f}%   "
          f"analitico alpha^3 = {esperado:.4f}%   (+-{4*sigma:.4f}% a 4 sigma MC)")
    print(f"    P(linha n=28, padrao A0-nula) simulado = {p_linha:.4f}%  (mesmo alpha^3)")
    print(f"    documento (Secao 4.1): {pub}%   {'OK' if ok else 'DIVERGE'}")
    print(f"    contraste -- formula INGENUA alpha^9 seria {ingenuo:.3e}%  "
          f"(o documento corretamente NAO usa esta)")
    if pub is None:
        div.append("Secao 4.1: valor alpha^3 nao encontrado no documento")
    elif not ok:
        div.append(f"Secao 4.1: pub={pub}% vs analitico alpha^3={esperado:.4f}%")

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
