#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testes das funcoes centrais usadas na verificacao dos protocolos pre-registrados.

Valida contra valores FECHADOS conhecidos (analiticos), independentes dos
documentos .md -- os scripts verificacao_*.py e que conferem contra os .md.

Uso:  python -m pytest analysis/protocols/ -v
"""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent


def _load(nome: str, caminho: Path):
    spec = importlib.util.spec_from_file_location(nome, caminho)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[nome] = mod
    spec.loader.exec_module(mod)
    return mod


p6 = _load("p6", BASE / "admissibilidade" / "verificacao_p6.py")
cxb = _load("cxb", BASE / "cxb" / "verificacao_cxb.py")
a0 = _load("a0", BASE / "a0" / "verificacao_a0.py")


# --------------------------------------------------------------------------
# Delta_theta_max
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,q", [(3, 2), (4, 3), (7, 5), (7, 4), (8, 5), (5, 2), (4, 1)])
def test_delta_theta_max_fechado_vs_numerico(p, q):
    """Forma fechada 4*arctan(sqrt(r))-pi deve bater com varredura numerica."""
    r = p / q
    assert p6.delta_theta_max(r) == pytest.approx(
        p6.delta_theta_max_numerico(r, n=200_001), abs=1e-4
    )


def test_delta_theta_max_identidade_cosseno():
    """cos(Delta_theta_max) = 1 - 2(p-q)^2/(p+q)^2  (Passo 2 da prova de P6)."""
    for (p, q) in [(3, 2), (7, 5), (7, 4), (5, 2), (4, 1), (8, 5)]:
        lhs = math.cos(p6.delta_theta_max(p / q))
        rhs = 1.0 - 2.0 * (p - q) ** 2 / (p + q) ** 2
        assert lhs == pytest.approx(rhs, abs=1e-12)


def test_delta_theta_max_unissono():
    """r=1 (unissono) -> Delta_theta_max = 0."""
    assert p6.delta_theta_max(1.0) == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------------------
# Criterio de admissibilidade (P6)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,q,esperado", [
    (3, 2, False),   # d=1
    (4, 3, False),   # d=1
    (5, 4, False),   # d=1
    (6, 5, False),   # d=1
    (9, 8, False),   # d=1
    (7, 5, False),   # d=2 -- tritono, o caso de fronteira
    (7, 4, True),    # d=3
    (8, 5, True),    # d=3
    (5, 2, True),    # d=3
    (4, 1, True),    # d=3
])
def test_admissibilidade_tabela_secao7(p, q, esperado):
    assert p6.admissivel_criterio_fechado(p, q) is esperado
    assert p6.admissivel_comparacao_direta(p, q) is esperado


def test_criterio_fechado_equivale_comparacao_direta():
    """Varredura exaustiva p,q<=60: os dois criterios devem concordar 100%."""
    from fractions import Fraction
    total = discord = 0
    for q in range(1, 61):
        for p in range(q + 1, 61):
            if Fraction(p, q).denominator != q:
                continue
            total += 1
            if p6.admissivel_criterio_fechado(p, q) != p6.admissivel_comparacao_direta(p, q):
                discord += 1
    assert total == 1101
    assert discord == 0


@pytest.mark.parametrize("q", [1, 2, 5, 50, 500, 5000])
def test_familia_d1_sempre_negativa(q):
    """d=1: M(q,1) = -3q(q+1)-1 < 0 para todo q."""
    assert p6.M(q, 1) < 0
    assert p6.M(q, 1) == -3 * q * (q + 1) - 1


@pytest.mark.parametrize("q", [1, 2, 5, 50, 500, 5000])
def test_familia_d2_constante_menos4(q):
    """d=2: M(q,2) = -4, constante, independente de q."""
    assert p6.M(q, 2) == -4


@pytest.mark.parametrize("d", [3, 4, 5, 10, 20])
def test_familia_d3_ou_mais_admissivel(d):
    """d>=3: M(1,d) > 0 ja no menor q possivel."""
    assert p6.M(1, d) > 0


# --------------------------------------------------------------------------
# n_min / n_geom
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,q,esperado", [
    (3, 2, 25.0),
    (4, 3, 49.0),
    (5, 4, 81.0),
    (6, 5, 121.0),
    (9, 8, 289.0),
    (7, 5, 36.0),
])
def test_n_min_valores_fechados(p, q, esperado):
    """n_min = ((p+q)/(p-q))^2 -- valores exatos para d=1 e d=2."""
    assert p6.n_min(p, q) == pytest.approx(esperado, rel=1e-12)


def test_n_min_identidade_com_n_geom_no_maximo():
    """n_min deve ser exatamente n_geom avaliado em Delta_theta_max."""
    for (p, q) in [(3, 2), (7, 5), (7, 4), (5, 2), (8, 5)]:
        assert p6.n_min(p, q) == pytest.approx(
            p6.n_geom(p6.delta_theta_max(p / q)), rel=1e-10
        )


@pytest.mark.parametrize("theta_deg,esperado", [
    (36.8698976, 10.0),    # 5:2 D3-exato
    (90.0, 2.0),
    (180.0, 1.0),
])
def test_n_geom_valores_conhecidos(theta_deg, esperado):
    assert cxb.n_geom_deg(theta_deg) == pytest.approx(esperado, rel=1e-6)


# --------------------------------------------------------------------------
# N_min (orcamento de trials)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("n_alvo,esperado", [
    (4, 444_444),
    (28, 3_111_111),
    (40, 4_444_444),
    (10, 1_111_111),
    (25, 2_777_778),
])
def test_n_min_trials_m003(n_alvo, esperado):
    """N ~ 100*n/m^2 com m=0,03 -- valores tabelados nos protocolos."""
    assert cxb.n_min_trials(n_alvo, 0.03) == esperado
    assert a0.n_min_trials(n_alvo, 0.03) == esperado


def test_n_min_trials_escala_inverso_m2():
    """N deve escalar como 1/m^2: m/2 -> 4x, m/3 -> 9x, m/6 -> 36x."""
    base = cxb.n_min_trials(28, 0.03)
    assert cxb.n_min_trials(28, 0.015) == pytest.approx(4 * base, rel=1e-6)
    assert cxb.n_min_trials(28, 0.01) == pytest.approx(9 * base, rel=1e-6)
    assert cxb.n_min_trials(28, 0.005) == pytest.approx(36 * base, rel=1e-6)


# --------------------------------------------------------------------------
# Razao de tangentes (AII)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("tA,tB,r_esperado", [
    (101.537, 78.463, 1.5),     # CxB config 1 (3:2 simetrico)
    (105.827, 74.173, 1.75),    # CxB config 2 (7:4 simetrico)
    (54.756, 32.969, 1.75),     # CxB config 3 / A0 orientacao 2
    (147.031, 125.244, 1.75),   # CxB config 4 (ramo 2)
    (65.985, 29.115, 2.5),      # CxB config 5 (5:2)
    (29.253, 7.466, 4.0),       # A0 orientacao 1
    (65.871, 44.084, 1.6),      # A0 orientacao 3
])
def test_razao_tangentes_configs(tA, tB, r_esperado):
    """Todos os pares angulares publicados devem realizar sua razao com erro < 0,01%."""
    r = cxb.razao_tangentes(tA, tB)
    assert r == pytest.approx(r_esperado, rel=1e-4)


def test_simetria_configs_B():
    """Configs 1 e 2 (regra B / D3'): theta_A + theta_B = 180 graus."""
    assert 101.537 + 78.463 == pytest.approx(180.0, abs=1e-9)
    assert 105.827 + 74.173 == pytest.approx(180.0, abs=1e-9)


def test_a0_tres_orientacoes_mesmo_delta_theta():
    """As 3 orientacoes do A0 compartilham exatamente o mesmo Delta_theta."""
    pares = [(29.253, 7.466), (54.756, 32.969), (65.871, 44.084)]
    deltas = [tA - tB for (tA, tB) in pares]
    for d in deltas:
        assert d == pytest.approx(21.787, abs=1e-9)
