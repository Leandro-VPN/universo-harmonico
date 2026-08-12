#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verificacao_p6.py -- reproduz TODOS os numeros de protocols/admissibilidade_AII_D3_v1.md
(Proposicao P6, admissibilidade AII + D3 no regime C).

NAO gera numeros novos: recalcula os publicados e CONFERE contra o .md.
Qualquer divergencia e reportada como DIVERGENCIA, nunca corrigida em silencio.

Uso:  python verificacao_p6.py
Saida: relatorio no stdout; exit code 0 se tudo confere, 1 se ha divergencia.
"""
from __future__ import annotations

import math
import re
import sys
from fractions import Fraction
from pathlib import Path

DOC = Path(__file__).resolve().parents[3] / "protocols" / "admissibilidade_AII_D3_v1.md"
TOL_REL = 1e-3  # tolerancia relativa p/ comparar com valores impressos com 2-4 casas


# --------------------------------------------------------------------------
# Funcoes centrais (as mesmas que os testes pytest cobrem)
# --------------------------------------------------------------------------
def delta_theta_max(r: float) -> float:
    """Passo 1 da prova: maximo de Delta_theta(theta_B) para razao r=p/q.
    Forma fechada: 4*arctan(sqrt(r)) - pi.  Retorna radianos."""
    return 4.0 * math.atan(math.sqrt(r)) - math.pi


def delta_theta_max_numerico(r: float, n: int = 2_000_001) -> float:
    """Verificacao cruzada independente: varredura fina de
    Delta_theta(theta_B) = 2*arctan(r*tan(theta_B/2)) - theta_B."""
    best = -math.inf
    for i in range(1, n):
        tb = math.pi * i / n
        ta = 2.0 * math.atan(r * math.tan(tb / 2.0))
        d = ta - tb
        if d > best:
            best = d
    return best


def n_geom(theta_rel: float) -> float:
    """D1: n_geom(theta) = 2/(1-cos theta)."""
    return 2.0 / (1.0 - math.cos(theta_rel))


def theta_rel_exigido(p: int, q: int) -> float:
    """D2+D3: arccos(1 - 2/(p*q))."""
    return math.acos(1.0 - 2.0 / (p * q))


def n_min(p: int, q: int) -> float:
    """Fronteira de admissibilidade: ((p+q)/(p-q))^2 (forma citada na Prop. P6)."""
    return ((p + q) / (p - q)) ** 2


def admissivel_criterio_fechado(p: int, q: int) -> bool:
    """Criterio fechado da Prop. P6: (p-q)*sqrt(p*q) >= (p+q)."""
    return (p - q) * math.sqrt(p * q) >= (p + q)


def admissivel_comparacao_direta(p: int, q: int) -> bool:
    """Mesma decisao, por comparacao geometrica direta: Delta_theta_max >= theta_rel exigido."""
    return delta_theta_max(p / q) >= theta_rel_exigido(p, q)


def M(q: int, d: int) -> float:
    """Reducao a familia d=p-q (Secao 5): M(q,d) = (d^2-4)*q*(q+d) - d^2.
    Admissivel  <=>  M >= 0."""
    return (d * d - 4) * q * (q + d) - d * d


# --------------------------------------------------------------------------
# Parsing do .md (a fonte de verdade a conferir)
# --------------------------------------------------------------------------
def _num(txt: str) -> float:
    """'26,68' -> 26.68 ; '8,163...' -> 8.163 ; '13,44' -> 13.44"""
    t = txt.strip().replace("*", "").replace("…", "").replace(".", "")
    t = t.replace("...", "").strip()
    t = t.replace(",", ".")
    return float(t)


def parse_tabela_secao7(doc: Path):
    """Extrai as linhas da tabela da Secao 7: (p, q, d_pub, nmod_pub, nmin_pub, admis_pub)."""
    txt = doc.read_text(encoding="utf-8")
    m = re.search(r"##\s*7\.(.*?)(?=\n##\s*8\.)", txt, re.S)
    if not m:
        raise SystemExit("ERRO: nao encontrei a Secao 7 em " + str(doc))
    linhas = []
    for ln in m.group(1).splitlines():
        ln = ln.strip()
        if not ln.startswith("|") or "razão" in ln.lower() or set(ln) <= set("|-: "):
            continue
        cels = [c.strip() for c in ln.strip("|").split("|")]
        if len(cels) < 5:
            continue
        mr = re.match(r"(\d+)\s*:\s*(\d+)", cels[0])
        if not mr:
            continue
        p, q = int(mr.group(1)), int(mr.group(2))
        d_pub = int(_num(cels[1]))
        nmod_pub = _num(cels[2])
        nmin_pub = _num(cels[3])
        admis_pub = "sim" in cels[4].lower()
        linhas.append((p, q, d_pub, nmod_pub, nmin_pub, admis_pub, cels[4]))
    return linhas


# --------------------------------------------------------------------------
def main() -> int:
    print("=" * 74)
    print("VERIFICACAO P6 -- admissibilidade_AII_D3_v1.md")
    print("=" * 74)
    print(f"documento: {DOC}")
    if not DOC.exists():
        print("ERRO: documento nao encontrado.")
        return 1

    divergencias: list[str] = []

    # ---- (1) Delta_theta_max: forma fechada vs varredura numerica -------------
    print("\n[1] Delta_theta_max(r) = 4*arctan(sqrt(r)) - pi")
    print("    verificacao cruzada contra varredura numerica fina:")
    for (p, q) in [(3, 2), (7, 5), (7, 4), (5, 2), (4, 1)]:
        r = p / q
        fechado = delta_theta_max(r)
        numerico = delta_theta_max_numerico(r)
        dif = abs(fechado - numerico)
        ok = dif < 1e-5
        print(f"    {p}:{q}  fechado={math.degrees(fechado):9.5f}deg  "
              f"numerico={math.degrees(numerico):9.5f}deg  |dif|={dif:.2e}  "
              f"{'OK' if ok else 'DIVERGE'}")
        if not ok:
            divergencias.append(f"Delta_theta_max {p}:{q}: fechado vs numerico difere {dif:.2e}")

    # ---- (2) Tabela da Secao 7, linha a linha --------------------------------
    print("\n[2] Tabela da Secao 7 -- conferencia linha a linha contra o .md")
    linhas = parse_tabela_secao7(DOC)
    print(f"    {len(linhas)} linhas encontradas no documento")
    print(f"    {'razao':>6} | {'d':>2} | {'n_mod':>7} | {'n_min pub':>10} | "
          f"{'n_min calc':>10} | {'admis pub':>9} | {'admis calc':>10} | status")
    print("    " + "-" * 88)
    for (p, q, d_pub, nmod_pub, nmin_pub, admis_pub, admis_txt) in linhas:
        d_calc = p - q
        nmod_calc = p * q
        nmin_calc = n_min(p, q)
        adm_fech = admissivel_criterio_fechado(p, q)
        adm_dir = admissivel_comparacao_direta(p, q)

        probs = []
        if d_calc != d_pub:
            probs.append(f"d pub={d_pub} calc={d_calc}")
        if abs(nmod_calc - nmod_pub) > 1e-9:
            probs.append(f"n_mod pub={nmod_pub} calc={nmod_calc}")
        if abs(nmin_calc - nmin_pub) / max(nmin_calc, 1e-12) > TOL_REL:
            probs.append(f"n_min pub={nmin_pub} calc={nmin_calc:.4f}")
        if adm_fech != admis_pub:
            probs.append(f"admissibilidade pub={admis_pub} calc={adm_fech}")
        if adm_fech != adm_dir:
            probs.append("criterio fechado != comparacao geometrica direta")

        status = "OK" if not probs else "DIVERGE"
        print(f"    {p:>3}:{q:<2} | {d_calc:>2} | {nmod_calc:>7} | {nmin_pub:>10.3f} | "
              f"{nmin_calc:>10.4f} | {str(admis_pub):>9} | {str(adm_fech):>10} | {status}")
        for pr in probs:
            divergencias.append(f"Secao 7, linha {p}:{q}: {pr}")

    # margem / deficit citados no texto
    print("\n    valores citados no texto da Secao 7:")
    margem_74 = (7 * 4) / n_min(7, 4)
    deficit_75 = (n_min(7, 5) - 7 * 5) / (7 * 5)
    print(f"      7:4 margem = n_mod/n_min = {margem_74:.4f}x   (texto: 2,08x)  "
          f"{'OK' if abs(margem_74 - 2.08) < 5e-3 else 'DIVERGE'}")
    print(f"      7:5 deficit = (n_min-n_mod)/n_mod = {deficit_75*100:.4f}%  (texto: 2,86%)  "
          f"{'OK' if abs(deficit_75 * 100 - 2.86) < 5e-3 else 'DIVERGE'}")
    if abs(margem_74 - 2.08) >= 5e-3:
        divergencias.append(f"margem 7:4 calc={margem_74:.4f} vs texto 2,08")
    if abs(deficit_75 * 100 - 2.86) >= 5e-3:
        divergencias.append(f"deficit 7:5 calc={deficit_75*100:.4f}% vs texto 2,86%")

    # ---- (3) Familia d = p-q -------------------------------------------------
    print("\n[3] Reducao a familia d=p-q (Secao 5), M(q,d)=(d^2-4)*q*(q+d)-d^2")
    QMAX = 5000
    d1_neg = all(M(q, 1) < 0 for q in range(1, QMAX + 1))
    d2_val = {M(q, 2) for q in range(1, QMAX + 1)}
    d2_const = (d2_val == {-4})
    print(f"    d=1: M(q,1) < 0 para todo q<= {QMAX} ......... {'OK' if d1_neg else 'DIVERGE'}")
    print(f"    d=2: M(q,2) == -4 constante p/ q<= {QMAX} .... {'OK' if d2_const else 'DIVERGE'}"
          f"   (valores distintos observados: {sorted(d2_val)})")
    d3_pos = all(M(1, d) > 0 for d in range(3, 60))
    print(f"    d>=3: M(1,d) > 0 (menor q possivel) .......... {'OK' if d3_pos else 'DIVERGE'}")
    if not d1_neg:
        divergencias.append("d=1: existe q com M(q,1) >= 0")
    if not d2_const:
        divergencias.append(f"d=2: M(q,2) nao e constante -4; valores={sorted(d2_val)}")
    if not d3_pos:
        divergencias.append("d>=3: existe d com M(1,d) <= 0")

    # ---- (4) Varredura exaustiva p,q <= 60 -----------------------------------
    print("\n[4] Varredura exaustiva de razoes irredutiveis p>q>=1, p,q <= 60")
    total = 0
    discord = 0
    for q in range(1, 61):
        for p in range(q + 1, 61):
            if Fraction(p, q).denominator != q:
                continue  # nao irredutivel
            total += 1
            if admissivel_criterio_fechado(p, q) != admissivel_comparacao_direta(p, q):
                discord += 1
    print(f"    razoes testadas: {total}   (documento: 1101)")
    print(f"    discordancias criterio fechado vs comparacao direta: {discord}")
    if total != 1101:
        divergencias.append(f"varredura: {total} razoes vs 1101 no documento")
    if discord != 0:
        divergencias.append(f"varredura: {discord} discordancias (documento afirma 100% concordancia)")
    print(f"    concordancia: {100.0*(total-discord)/total:.2f}%  (documento: 100%)")

    # ---- veredito ------------------------------------------------------------
    print("\n" + "=" * 74)
    if divergencias:
        print(f"RESULTADO: {len(divergencias)} DIVERGENCIA(S) ENCONTRADA(S)")
        for d in divergencias:
            print("  - " + d)
        print("\nNENHUM lado foi corrigido automaticamente (regra do pre-registro).")
        return 1
    print("RESULTADO: todos os numeros conferem com o documento.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
