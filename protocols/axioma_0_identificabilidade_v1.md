# Axioma 0 — Identificabilidade e Protocolo de Teste — v1

*Por que dados de Bell publicados não podem testar A0 (com demonstração
empírica na ETH), e por que essa limitação NÃO se estende a um
protocolo desenhado com o próprio Axioma II — que tem assinatura
espectral, imune ao problema estático. Pré-requisitos:
admissibilidade_AII_D3_v1.md, estrutura_v3.md (Predição 2, VI). Todos
os números verificados numericamente antes da redação.*

## 1. O axioma e o que ele custa

**A0:** ângulos medidos numa base absoluta da fonte; a física depende
de (θ_A, θ_B) individualmente, não só de Δθ = θ_A − θ_B. Custo
declarado: quebra a invariância rotacional que a QM tem (E(a,b) =
−cos Δθ depende só da diferença).

## 2. Por que minerar dados publicados não funciona — o caso ETH

Testado nesta investigação (ago/2026): comparar E medido vs. E_QM(ρ)
usando a matriz densidade publicada da ETH, nos 4 correlatores do CHSH.
Resultado (estrutura_v3.md, VI): dispersão de 19,7σ entre os |E| que a
QM ideal previa idênticos — mas totalmente absorvível por desvios
angulares de ~0,8–1,0° (ajuste exato com 4 parâmetros para 4 números).
Precisão estatística equivaleria a conhecer os ângulos a 0,17°; a
publicada é ~1–2°. Gap de 6–10×.

**Nota (resposta de S. Storz, ago/2026):** o gap de visibilidade entre
a tomografia (S=2,1876) e o run de CHSH (S=2,0747) tem explicação
metodológica esperada. A cadeia: (i) tomografia reconstruída com
correção de erro de leitura → métrica de visibilidade mais alta; (ii)
CHSH loophole-free não pode aplicar essa correção aos próprios dados
(reabriria brecha de detecção) → diferença de visibilidade esperada,
não misteriosa; (iii) diferenças residuais entre os 4 correlatores
individuais ainda podem ser explicadas por calibração angular, como
já havíamos escrito. **Conclusão: o resultado estático não identifica
A0** — e é exatamente por isso que este documento propõe, em vez de
tentar fechar essa degenerescência com mais casas decimais, uma
assinatura espectral dependente de orientação (§3–5). A própria frase
de Storz ("pequenas mudanças de calibração podem explicar as
diferenças") é a classe de explicação alternativa que o argumento de
identificabilidade precisava manter em aberto — apoio externo à
prudência da interpretação, não à hipótese A0.

**Por que isso não é falta de decimais — é degenerescência estrutural.**
Erro de calibração angular e efeito de A0 produzem exatamente a mesma
assinatura sobre um conjunto pequeno de correlatores estáticos: um
deslocamento de médias, absorvível por qualquer ajuste com graus de
liberdade suficientes. Prova adicional (verificada, estrutura_v3.md
Predição 2): erros de ângulo/visibilidade são CONSTANTES por combo —
contaminam apenas a componente DC do espectro, nunca vazam para bins de
frequência ≠ 0. Isso não é peculiaridade da ETH: **qualquer teste
baseado em poucos correlatores estáticos, de qualquer experimento
publicado, tem esse mesmo teto.**

## 3. A saída: A0 tem assinatura espectral, não apenas estática

O Axioma II já opera sobre ângulos absolutos: r = tan(θ_A/2)/tan(θ_B/2).
Fixando Δθ = θ_A − θ_B e variando a orientação absoluta (equivalente a
deslizar θ_B ao longo de sua curva, θ_A = θ_B + Δθ), a razão r **não é
constante** — ela percorre um intervalo contínuo de valores conforme a
orientação muda (é exatamente a curva Δθ(θ_B) de
admissibilidade_AII_D3_v1.md, §3). Em pontos especiais dessa curva, r
cruza razões simples p:q, cada uma com seu próprio n_mod = pq.

**Consequência testável:** para o MESMO Δθ, orientações absolutas
diferentes tornam ressonantes razões p:q diferentes — logo preveem
linhas espectrais em n_mod diferentes. Isso não é degenerado com o
erro estático de calibração (§2), sob as premissas do protocolo: um
erro constante por configuração contamina apenas DC, enquanto a
hipótese prevê estrutura em frequência não-nula. **Isso não é
imunidade absoluta** — não protege contra artefato periódico
instrumental (mesma classe já tratada em outros pontos do projeto:
LFSR, divisores de clock), controlado pelos mesmos meios já em uso
(permutação, linhas de exclusão de ambiente). Vale notar uma
assimetria: um artefato de período FIXO e universal (independente de
ângulo) imitaria o padrão de **A0-nula** (mesmo n=28 nas três
orientações) muito mais facilmente do que o padrão de **A0**
(diagonal 4/28/40, exigindo que o artefato produza três picos
diferentes, cada um só na orientação certa) — a diagonal é a
assinatura mais difícil de forjar por coincidência instrumental.

### Exemplo concreto (verificado)

Fixando Δθ = 21,787° (o mesmo da configuração "7:4 D3-exato, ramo 1" já
usada no protocolo C×B), três orientações absolutas distintas:

| Orientação | θ_B (Bloch) | θ_A (Bloch) | Razão ressonante | n_mod sob A0 |
|---|---|---|---|---|
| 1 | 7,466° | 29,253° | 4:1 | **4** |
| 2 | 32,969° | 54,756° | 7:4 | **28** |
| 3 | 44,084° | 65,871° | 8:5 | **40** |

Confirmado: n_geom(21,787°) = 27,9995 ≈ 28 — a orientação 2 satisfaz
D3 exatamente, por construção (é a mesma config já verificada em
admissibilidade_AII_D3_v1.md).

## 4. Três hipóteses, três padrões distintos, pré-especificados

| Hipótese | Previsão em (1) | (2) | (3) |
|---|---|---|---|
| **A0** (AII usa ângulos absolutos, modelo atual) | linha em n=4 | linha em n=28 | linha em n=40 |
| **A0-nula** (ressonância só via Δθ, ignorando AII) | linha em n=28 | linha em n=28 | linha em n=28 |
| **QM pura** (sem modulação harmônica, m=0) | nenhuma linha | nenhuma linha | nenhuma linha |

A0-nula é o competidor honesto que faltava nos documentos anteriores:
um modelo que mantém D3 (ressonância determinada por Δθ, via D1
diretamente) mas **descarta** o Axioma II como mecanismo de seleção de
n_mod. Sob A0-nula, a mesma linha (n=28) apareceria nas três
orientações — rotacionalmente invariante, como a QM padrão exige.

### 4.1 Regra estatística conjunta (9 células, declarada antes da coleta)

A matriz 3×3 (3 orientações × 3 alvos) não é 9 testes independentes
avaliados a α=5% cada — é UM padrão a testar, com significância
calibrada por permutação conjunta, mesma disciplina do protocolo
C×B: embaralhar rótulos de orientação/trial e recomputar a estatística
completa, nunca calibrar célula a célula com limiar analítico.

**Critério confirmatório primário:** as TRÊS células-alvo (diagonal,
para A0; ou a linha n=28 nas três orientações, para A0-nula) excedem
o limiar de permutação conjunta. Sob independência aproximada,
P(3 falsos positivos conjuntos por acaso) ≈ α³ ≈ 0,0125% — já bem
abaixo do 0,25% usado como referência no protocolo C×B; não é
necessário exigir adicionalmente que as 6 células "negativas"
estejam limpas para ter confiança — uma flutuação isolada de 5% entre
6 testes secundários não deve derrubar uma confirmação primária já
sólida por si.
**As 6 células negativas entram como reforço diagnóstico, não como
exigência dura:** reportadas sempre, mas um único excesso marginal
entre elas não invalida o critério primário; um padrão sistemático
entre várias delas seria motivo para reabrir a interpretação.

**O teste é auto-contido:** não depende de conhecer os ângulos
absolutos com precisão EXTERNA (o problema do §2, ~0,1° para
interpretar um correlator estático a posteriori) — depende apenas de
comparar o padrão de linhas ENTRE as três orientações. **Isso não
elimina a necessidade de controle angular tout court**: o experimento
ainda precisa PRODUZIR/CONFIGURAR as três orientações com precisão
suficiente para realizar cada Δθ pretendido (é controle relativo de
bancada, rotineiro em qualquer aparato de Bell, que já define Δθ hoje)
— o que se elimina é a necessidade de SABER o valor absoluto no
referencial do laboratório com precisão externa para interpretar o
resultado. A pergunta não é "qual o ângulo exato", é "o espectro muda
quando eu giro o aparato mantendo Δθ fixo?".

## 5. Protocolo mínimo

**N (m_alvo=0,03, mesma convenção do protocolo C×B):**
orientação 1 (n=4): 444.444 | orientação 2 (n=28): 3.111.111 |
orientação 3 (n=40): 4.444.444 | **total: 7.999.999 ≈ 8M**
(80s a 100 kHz — barato; m_alvo é decisão pré-registrada, escala 1/m²
como nos demais protocolos).

**Regra de decisão (matriz 3×3 completa — nove células, diagonal vs.
linha vs. nada):**

A matriz confirmatória é avaliada nas NOVE células orientação × alvo
(n = 4, 28, 40), não apenas nas três primárias. O N de cada orientação
é determinado pelo seu alvo primário, mas as três frequências são
avaliadas no mesmo espectro de cada orientação — sem coleta adicional,
pois todos os bins já existem no periodograma calculado.

|  | n=4 | n=28 | n=40 |
|---|---|---|---|
| **Orient. 1** (4:1) | ✓ | ✗ | ✗ |
| **Orient. 2** (7:4) | ✗ | ✓ | ✗ |
| **Orient. 3** (8:5) | ✗ | ✗ | ✓ |

- **Evidência confirmatória para A0**, sob o catálogo de sistemáticas
  pré-especificado: padrão diagonal — presença nas três células-alvo E
  ausência pré-especificada nos seis off-diagonais.
- **Evidência confirmatória para A0-nula:** presença consistente em
  n=28 nas três orientações E ausência nos seis bins off-target.
- **QM pura:** nenhuma célula significativa.
- **Padrão misto/inconsistente:** reportar sem promover a nenhuma das
  três — mesma disciplina do protocolo C×B (Seção 5 daquele
  documento).

**Linguagem de reporte (obrigatória):** nenhum desfecho constitui
prova. Um padrão diagonal é fortemente consistente com A0, mas
sistemáticas instrumentais que correlacionem orientação e período
permanecem, em princípio, explicação alternativa (§3, ressalva sobre
artefato periódico universal). Reportar sempre como "evidência
confirmatória sob o catálogo pré-especificado", nunca como "A0
confirmado".

Mesma regra de bin único do protocolo C×B (Seção 6.1 daquele
documento) aplica-se aqui sem alteração: teste confirmatório apenas no
bin pré-especificado por orientação/alvo, sem substituição pós-hoc.

## 6. O que este protocolo decide, e o que não decide

- **Decide, simultaneamente:** (a) se existe modulação harmônica
  alguma (m>0, o próprio C) nestas três configurações específicas, e
  (b) SE existir, se seu n_mod segue o Axioma II (A0) ou é invariante
  rotacional (A0-nula). É um teste conjunto, não dois testes
  separados — um resultado "A0-nula confirmada" também confirma C **no
  sentido operacional definido por este protocolo** (existência de
  modulação nesta configuração específica) — NÃO confirma o modelo
  harmônico como um todo: o protocolo não testa E3 nem E4, nem
  qualquer configuração fora das três aqui especificadas.
- **Não decide:** o fork E3 (M1 vs M2) nem o fork E4/B — as três
  orientações aqui usam exclusivamente configurações D3-exatas (regime
  C). Um teste equivalente sob B (buscando se θ_A+θ_B=180° também
  gera dependência de orientação análoga) fica para v2.
- **Não requer** a metrologia de ~0,1° identificada como faltante no
  §2 — o teste é relativo (compara padrões entre orientações na mesma
  plataforma), não absoluto.

## 7. Nota de integração com o Programa II

Este protocolo e o protocolo_discriminante_CxB_v1.md compartilham
infraestrutura total (mesmo pipeline espectral, mesma regra de bin
único, mesmo m_alvo pré-registrado) mas testam perguntas ortogonais:
C×B pergunta *qual regra de ressonância* seleciona ângulos ressonantes,
dado que a física é local em cada configuração; A0 pergunta se a
própria física muda com orientação absoluta, a Δθ fixo. Podem ser
executados na mesma sessão experimental — as orientações 2 e a
configuração "7:4 D3-exato ramo 1" do protocolo C×B **são a mesma
configuração física**, reaproveitável sem custo adicional.
