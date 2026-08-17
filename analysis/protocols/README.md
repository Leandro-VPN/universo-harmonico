# analysis/protocols/ — verificação numérica do pré-registro

Reproduz **todos** os números publicados nos três documentos de
`protocols/` e os confere contra o texto. Estes scripts não geram números
novos e não fazem análise de dados experimentais — nenhum dado foi coletado
ainda.

## Como rodar

```bash
python analysis/protocols/admissibilidade/verificacao_p6.py
python analysis/protocols/cxb/verificacao_cxb.py
python analysis/protocols/a0/verificacao_a0.py
python -m pytest analysis/protocols/ -v
```

Dependências: `numpy`, `pytest` (só para a suíte de testes).

**Exit code 0** = todos os números do `.md` conferem.
**Exit code 1** = há divergência; cada uma é listada no fim do relatório.
Nenhum dos dois lados (script ou documento) é corrigido automaticamente —
essa é a regra central do pré-registro.

## O que cada script cobre

**`admissibilidade/verificacao_p6.py`**
- Δθ_max(r) = 4·arctan(√r) − π, com verificação cruzada por varredura
  numérica fina (5 razões).
- Tabela da Seção 7 linha a linha: d, n_mod, n_min, veredito de
  admissibilidade — pelos dois critérios (fechado e comparação geométrica).
- Margem 7:4 (2,08×) e déficit 7:5 (2,86%) citados no texto.
- Família d = p−q: M(q,1) < 0 e M(q,2) ≡ −4 para q até 5000; M(1,d) > 0
  para d ≥ 3.
- Varredura exaustiva das 1101 razões irredutíveis com p,q ≤ 60,
  confirmando concordância entre os dois critérios.

**`cxb/verificacao_cxb.py`**
- Os 5 pares angulares da Seção 2, em Bloch e em laboratório fotônico
  (convenção D6: θ_lab = θ_Bloch/2): razão de tangentes, n_mod, n_geom,
  simetria θ_A+θ_B=180° (configs B) e D3 exato (configs C).
- Tabela de N da Seção 4 (m=0,03) e a escala 1/m² para 0,02/0,01/0,005.
- Simulação da taxa de falso-positivo conjunto da Seção 5:
  P(3∧4∧5), P(3∧4 apenas) e a regra antiga "qualquer par".

**`a0/verificacao_a0.py`**
- As 3 orientações da Seção 3 para Δθ = 21,787°, com n_geom(Δθ) = 28.
- Consistência cruzada exigida pelo `PROTOCOL_LOCK`: a orientação 2 do A0
  é literalmente a configuração 3 do C×B (mesmos ângulos).
- N por **poder uniforme** (Seção 5, A0 v1.2): as três orientações recebem
  o mesmo N, fixado pelo maior alvo da matriz com teto inteiro
  (⌈N_min(n=40)⌉ = 4.444.445; total 13.333.335). Confere também os fatores
  de poder relativo (0,14× e 0,10×) que o documento cita como justificativa,
  e imprime o dimensionamento antigo (por alvo primário, v1.1) como contraste.
- Simulação do falso-positivo **conjunto** da matriz 3×3 (Seção 4.1),
  contrastando com a fórmula ingênua α⁹ que o protocolo corretamente evita.

**`test_protocolos.py`** — suíte pytest das funções centrais
(Δθ_max, critério de admissibilidade, n_min, n_geom, N_min, razão de
tangentes) contra valores fechados analíticos, independentes dos `.md`.

## Nota sobre as simulações

Os números de falso-positivo nos documentos são saídas de simulação
(200.000 repetições) e portanto carregam ruído de Monte Carlo. A
conferência é contra o valor **analítico** dentro da barra de erro
esperada (4σ), não igualdade exata — reproduzir o dígito exato exigiria a
semente original, que não faz parte do pré-registro.
