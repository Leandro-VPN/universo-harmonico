# Protocolo Experimental — Discriminação C vs B — v1

*Pré-registro. Objetivo: construir um experimento em que as regras de
ressonância C (D3 atual) e B (D3′, achado ago/2026) façam previsões
opostas e pré-especificadas, resolvendo o fork E4 por dado, não por
escolha teórica. Pré-requisito: admissibilidade_AII_D3_v1.md (P6).
Todos os ângulos e N verificados numericamente nesta sessão.*

## 1. As duas hipóteses em jogo

- **C (D3 atual):** ressonância = comensurabilidade temporal,
  n_mod = n_geom(θ_rel). Admissível apenas para d=p−q≥3 (P6).
- **B (D3′):** ressonância = simetria geométrica, θ_A+θ_B=180° (ponto
  de espalhamento angular máximo). Admissível para QUALQUER razão p:q,
  inclusive as excluídas de C (d≤2, como 3:2).

Nenhuma das duas é assumida verdadeira. O experimento é desenhado para
produzir padrões de resultados mutuamente discriminantes entre C e B,
segundo regras de decisão pré-especificadas (Seção 5) — não para
garantir, por construção estatística, que apenas uma hipótese
"sobreviva". Nada impede logicamente que ambos os padrões apareçam
simultaneamente (configurações totalmente independentes), que nenhum
apareça (nulo, o desfecho mais provável a priori), ou que o padrão seja
ambíguo — são desfechos empíricos genuínos, não excluídos pelo desenho.

## 2. Configurações

| # | Configuração | Razão | n_mod | Bloch (θ_A, θ_B) | Lab fotônico (θ_A, θ_B) | n_geom | Prediz linha sob |
|---|---|---|---|---|---|---|---|
| 1 | 3:2 simétrico | 3:2 | 6 | 101,537°, 78,463° | 50,769°, 39,231° | 25,00 | **B** |
| 2 | 7:4 simétrico | 7:4 | 28 | 105,827°, 74,173° | 52,913°, 37,087° | 13,44 | **B** |
| 3 | 7:4 D3-exato (ramo 1) | 7:4 | 28 | 54,756°, 32,969° | 27,378°, 16,485° | 28,00 | **C** |
| 4 | 7:4 D3-exato (ramo 2) | 7:4 | 28 | 147,031°, 125,244° | 73,516°, 62,622° | 28,00 | **C** |
| 5 | 5:2 D3-exato | 5:2 | 10 | 65,985°, 29,115° | 32,992°, 14,557° | 10,00 | **C** |
| 6 | Controle | irracional | — | soma ≠ 180°; razão ∉ ℚ simples | — | — | nenhuma |

Todas as razões de tangentes verificadas (erro < 0,002% em todos os
casos). Configurações 3 e 4 têm a MESMA razão (7:4) e o MESMO n_mod
alvo (28), mas ângulos diferentes — controle interno cuja lógica exata
está na Seção 3: a concordância de sinal entre 3 e 4 torna improvável
um artefato específico de um único par angular; contudo, como ambas
compartilham a mesma razão 7:4 e o mesmo período n=28, essa
concordância NÃO exclui um artefato dependente da razão ou do período
— motivando a configuração 5.

## 3. Por que estas cinco (não outras) — e o que cada uma controla

Duas classes de artefato ameaçam qualquer linha espectral encontrada:

- **(i) Artefato ligado ao ÂNGULO específico** (ex. birrefringência ou
  deriva de Pockels cell num par de ângulos particular).
- **(ii) Artefato ligado ao PERÍODO n em si**, independente do ângulo
  (ex. divisor de clock, laço de software). Precedente real neste
  projeto: o pico de S/N~22× em n=15 no NIST (Figura 7 do v10, nunca
  reproduzido) era candidato exatamente a este tipo — período fixo, não
  ligado a ângulo.

Cada configuração controla uma classe:

- **3 ∧ 4** (mesma razão 7:4, mesmo n_mod=28, ângulos muito distintos):
  concordância entre elas descarta (i) — seria coincidência extrema um
  artefato de ângulo específico afetar dois pares tão diferentes
  igualmente. **NÃO descarta (ii)**, pois ambas visam o mesmo n=28.
- **5** (razão diferente, 5:2, n_mod=10): concordância de 5 com 3∧4
  descarta (ii) especificamente para n=28 — um artefato de período fixo
  em 28 não explicaria sinal também em n=10, via mecanismo e ângulos
  totalmente diferentes.
- **1, 2 testam B em duas razões distintas** (3:2 e 7:4) — controla (ii)
  para B, mas **não tem análogo ao controle de ângulo que 3∧4 dá**: a
  condição simétrica (θ_A+θ_B=180°) tem solução angular ÚNICA por
  razão, não dois ramos como D3. **Assimetria declarada:** o desenho
  atual dá a C um controle que B não tem nesta versão; um refinamento
  futuro (v2) poderia buscar duas razões com o mesmo n_geom simétrico
  para replicar o controle de ângulo também em B — não implementado
  aqui por escopo.
- **6 é nulo esperado sob as duas hipóteses** — controle negativo.

## 3.1 Reserva metodológica sobre a Seção 5

3∧4 e (3∧5 ou 4∧5) **não são substitutos** — descartam ameaças
diferentes. Uma regra que trata "qualquer par de {3,4,5}" como
equivalente (como uma versão anterior deste protocolo fazia) esconde
essa diferença e permite que a força da evidência dependa de qual par
os dados fizerem significativo — exatamente o tipo de grau de liberdade
que um pré-registro existe para eliminar. A Seção 5 abaixo resolve isso
com uma regra hierárquica, fixada antes da coleta.

## 4. Alvo de sensibilidade (decisão pré-registrada, não herdada)

N calculado para **m_alvo = 0,03**, via N_min ≈ 100·n_geom/m² (nota:
usa n_geom, que diverge de n_mod sob B — ver admissibilidade_AII_D3
§4). Escala como 1/m²:

| m_alvo | fator de escala | N total (5 configs) |
|---|---|---|
| 0,03 (padrão deste protocolo) | 1× | 11.604.444 |
| 0,02 | 2,25× | 26.110.000 |
| 0,01 | 9× | 104.440.000 |
| 0,005 | 36× | 417.760.000 |

A 100 kHz: 0,03→**116s**; 0,01→**~17,4 min**; 0,005→**~69,6 min**. O
experimento é barato em qualquer alvo razoável; **m_alvo deve ser
declarado antes da coleta**, não escolhido depois de ver o resultado.

## 5. Regra de decisão — hierárquica, fixada antes da coleta

Versão inicial testada em simulação (dados sintéticos, teste único por
hipótese, NPERM=15): produziu falsos positivos marginais (S/N 1,55 vs
limiar 1,47; S/N 2,40 vs limiar 1,48). Diagnóstico com NPERM=60 (não a
causa — o diagnóstico) reduziu a taxa para perto do nominal, apontando
viés de estimador de p95 em pequena amostra (catálogo A2) como origem.
Corrigido para exigir significância conjunta — mas "qualquer par"
mistura duas classes de ameaça distintas (Seção 3.1), então a regra
final é hierárquica:

**Para C:**
- **Confirmação plena:** 3 ∧ 4 ∧ 5 todas excedem o limiar. Descarta
  simultaneamente artefato de ângulo (i) e de período (ii). Padrão-ouro.
- **Confirmação primária, não replicada:** 3 ∧ 4 excedem, 5 não.
  Descarta (i), não descarta (ii). Reportar como "sinal consistente com
  D3 em n=28, robusto a ângulo, sem replicação independente de período"
  — NÃO tratar como equivalente à confirmação plena.
- **Padrão inconsistente com C** (reportar, não promover): apenas um de
  {3,4} excede; ou 5 excede sem 3∧4. Um artefato genuíno de D3 deveria
  aparecer nos DOIS ramos de 7:4 (mesma razão) — encontrar em só um
  sugere artefato de ângulo específico, não ressonância.

**Para B:**
- **Confirmação:** 1 ∧ 2 excedem. Descarta (ii) por razão cruzada; NÃO
  há controle de ângulo disponível neste desenho (assimetria declarada,
  Seção 3) — ler com essa ressalva mesmo em caso positivo.

**Verificado por simulação (200.000 repetições sob nula pura, α=5% por
teste):** P(3∧4∧5) = 0,008% | P(3∧4 apenas) = 0,253% | P("qualquer par
de 3", a regra anterior) = 0,712% — quase 90× mais permissiva que a
confirmação plena. A hierarquia explícita evita que a força da
conclusão dependa de qual par os dados fizerem significativo.

Nenhum padrão limpo em nenhum nível → resultado nulo, limite superior
de m combinado sob C e sob B nos ângulos ressonantes (primeira vez,
distinto de toda calibração anterior em ângulos não-ressonantes).

## 6. Protocolo de análise (reaproveita pipeline já validado)

1. Extrair série de concordância (ou minoritária, D4) por configuração.
2. FFT centrada; S/N no bin do n_mod alvo, fundo por mediana local.

### 6.1 Regra do bin único (fecha grau de liberdade de p-hacking)

O teste confirmatório (Seção 5) é realizado EXCLUSIVAMENTE no bin de
frequência correspondente ao n_mod pré-especificado de cada
configuração. Nenhum bin vizinho pode substituir o alvo após inspeção
dos dados; nenhuma "melhor janela ao redor do pico" é permitida na
regra confirmatória.

**Reconciliação com o Axioma IV-d:** o modelo já admite que decoerência
por difusão contínua (IV-d, vs. reset abrupto IV-c) alargaria a linha
(Lorentziana, FWHM≈σ_φ²N/2π — ver estrutura_v3.md, P5), o que tornaria
uma busca de bin único subótima se a física real for do tipo IV-d. Isso
NÃO abre exceção à regra acima. Em vez disso: **o teste confirmatório
deste protocolo assume IV-c (linha estreita)** — é a hipótese mais
simples e a que fixa o orçamento de N da Seção 4. Se a busca de bin
único for nula em todas as cinco configurações, uma análise
EXPLORATÓRIA por filtro casado (grade de larguras, já validada na
calibração do Axioma IV-d) pode ser rodada — mas seu resultado tem
status estritamente não-confirmatório, mesma assimetria já estabelecida
para a linha CURBy ("ausência de excesso limita m; presença é ambígua,
não confirma"). Um resultado positivo apenas na busca exploratória
exigiria um NOVO pré-registro dedicado a IV-d antes de ser tratado como
evidência.

3. Permutação ≥60 para o limiar; nunca limiar analítico (regra
   permanente do projeto).
4. Se múltiplas sessões/runs, Welch com segmento = sessão real, nunca
   concatenação (regra permanente).
5. Aplicar a regra de decisão da Seção 5.
6. Reportar as cinco configurações juntas, independente do resultado —
   inclusive se nulo em tudo (limite superior de m sob C e sob B, nos
   ângulos ressonantes pela primeira vez, não apenas não-ressonantes
   como toda a calibração anterior).

## 7. O que este experimento NÃO decide

Não decide o fork E3 (M1 vs M2) — é ortogonal (ver estrutura_v3.md,
E3/E4). Não decide o Axioma 0 — requer protocolo próprio (ver
axioma_0_identificabilidade.md, a escrever). Um resultado nulo aqui não
falsifica o modelo harmônico inteiro, apenas C e B nos alvos testados
— outras razões admissíveis sob C (8:5, 4:1, ...) permanecem abertas.
