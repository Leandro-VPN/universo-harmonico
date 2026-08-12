# Admissibilidade AII + D3 — v1

*Documento puramente formal: caracteriza o domínio em que o Axioma II
(mapeamento angular-espectral) e a Definição 3 (condição de ressonância,
regime C) podem ser satisfeitos simultaneamente. Não defende a hipótese
harmônica nem decide o fork C vs B (E4) — apenas prova uma consequência
matemática condicional. Pré-requisito de Predição 1A
(estrutura_v3.md). Todas as afirmações verificadas numericamente antes
da redação.*

## 1. Objetivo e escopo

Dado o Axioma II e a Definição 3 do modelo harmônico, para quais razões
p:q existe um par de ângulos (θ_A, θ_B) satisfazendo as duas
simultaneamente? A resposta é uma condição fechada, dependente de uma
única variável inteira, e resolve por que a proposta experimental do
v10 (trítono 7:5) precisou ser abandonada — não por erro de aritmética,
mas por uma propriedade estrutural do modelo.

## 2. Definições usadas

- **Axioma II:** tan(θ_A/2)/tan(θ_B/2) = p/q, razão irredutível.
- **D2 (período dinâmico):** n_mod = p·q.
- **D3 (condição de ressonância, regime C):** n_mod = n_geom(θ_rel),
  onde θ_rel = |θ_A − θ_B| e n_geom(θ) = 2/(1−cosθ).

Combinando D2 e D3: o ângulo relativo *exigido* para ressonância é

  θ_rel^exigido = arccos(1 − 2/(pq)).

## 3. O problema matemático

Fixada a razão p/q, o Axioma II parametriza uma família de pares
(θ_A(θ_B), θ_B) com θ_A = 2·arctan[(p/q)·tan(θ_B/2)]. A diferença
Δθ(θ_B) = θ_A − θ_B não é livre — tem um valor máximo. A pergunta de
admissibilidade é: **esse máximo alcança o θ_rel exigido por D3?**

## 4. Proposição P6

> Para razão irredutível p/q (p>q≥1), existe par angular satisfazendo
> simultaneamente AII e D3 (regime C) se e somente se
>
> **n_mod ≥ ((p+q)/(p−q))²**
>
> equivalentemente, sob D2 (n_mod=pq): **(p−q)·√(pq) ≥ (p+q)**.

### Demonstração

**Passo 1 — o máximo de Δθ.** Derivando Δθ(θ_B) e igualando a zero, o
ponto crítico ocorre em tan(θ_B/2) = 1/√r (r=p/q). Usando a identidade
arctan(x)+arctan(1/x)=π/2:

  Δθ_max(r) = 4·arctan(√r) − π.

*(Verificado por duas vias independentes: varredura numérica fina e
cálculo direto — concordância em todos os casos testados.)*

**Passo 2 — forma fechada sem trigonometria.** Seja x = arctan(√r).
Por cos(2x) = (1−r)/(1+r) = (q−p)/(q+p), e cos(Δθ_max) = cos(4x) =
2cos²(2x)−1:

  cos(Δθ_max) = 1 − 2(p−q)²/(p+q)².

**Passo 3 — a condição exigida por D3.** De sin²(θ_rel^exigido/2) =
1/(pq): cos(θ_rel^exigido) = 1 − 2/(pq).

**Passo 4 — comparação.** Δθ_max ≥ θ_rel^exigido ⟺ cos(Δθ_max) ≤
cos(θ_rel^exigido) (cosseno decrescente em [0,π]) ⟺

  2(p−q)²/(p+q)² ≥ 2/(pq) ⟺ (p−q)²·pq ≥ (p+q)² ⟺
  (p−q)·√(pq) ≥ (p+q). ∎

Equivalentemente, isolando n_mod=pq: **pq ≥ ((p+q)/(p−q))²**, a forma
citada na Proposição.

## 5. Redução à família d = p − q

Escrevendo p = q+d:

  M(q,d) = (d²−4)·q(q+d) − d².

**d=1:** M = −3q(q+1) − 1 < 0 para todo q≥1. *Impossível identicamente*
(não assintótico — o termo (d²−4)=−3 é negativo, então M diverge para
−∞ conforme q cresce).

**d=2:** M = 0·q(q+2) − 4 = **−4**, constante, independente de q.
*Impossível para todo q* — o déficit não é "quase zero", é exatamente
−4 sempre.

**d≥3:** (d²−4)≥5>0, então M cresce sem limite com q; no menor q
possível (q=1) já é positivo para d≥3 (checar: d=3,q=1: M=5·1·4−9=11>0).
*Sempre admissível*, com margem crescente.

**Prova de que d=1,2 são impossíveis:** decorre diretamente do sinal de
M acima — não é enumeração de casos, é o sinal de uma expressão fechada.
**Prova de que d≥3 é admissível:** M(1,d) = (d²−4)(1+d) − d² > 0 para
todo d≥3, verificável por substituição direta.

## 6. Verificação numérica

- Varredura exaustiva de todas as razões irredutíveis p,q ≤ 60 (1101
  razões): concordância 100% entre o critério fechado e a comparação
  direta Δθ_max vs θ_rel^exigido.
- Família d=1: M<0 confirmado até q=5000 (identidade algébrica,
  verificação apenas ilustrativa).
- Família d=2: M≡−4 confirmado até q=5000 (idem).

## 7. Exemplos

| razão | d=p−q | n_mod=pq | n_min=((p+q)/(p−q))² | admissível? |
|---|---|---|---|---|
| 3:2 (quinta) | 1 | 6 | 25,00 | **não** |
| 4:3 (quarta) | 1 | 12 | 49,00 | **não** |
| 5:4 (terça M) | 1 | 20 | 81,00 | **não** |
| 6:5 (terça m) | 1 | 30 | 121,00 | **não** |
| 9:8 (tom) | 1 | 72 | 289,00 | **não** |
| 7:5 (trítono) | 2 | 35 | 36,00 | **não** (déficit 2,86%) |
| 7:4 (7ª harm.) | 3 | 28 | 13,44 | **sim** (margem 2,08×) |
| 8:5 (6ª menor) | 3 | 40 | 18,78 | **sim** |
| 5:2 | 3 | 10 | 5,44 | **sim** |
| 4:1 | 3 | 4 | 2,78 | **sim** |

O trítono (d=2) é o caso mais próximo da fronteira de toda a família
excluída — falha por apenas 2,86%, não é exclusão folgada.

## 8. O que P6 NÃO afirma

- **Não afirma que consonâncias musicais sejam fisicamente proibidas.**
  A exclusão vive inteiramente dentro do formalismo do modelo (Axioma
  II + Definição 3); não é um enunciado sobre acústica ou percepção.
- **Não afirma que D3 seja verdadeiro na natureza.** D3 é uma definição
  escolhida (regime C); P6 caracteriza suas consequências, não sua
  veracidade física.
- **Não seleciona entre C e B (E4).** P6 caracteriza inteiramente o
  domínio de C. B (D3′: θ_A+θ_B=180°) é uma regra de ressonância
  diferente, sob a qual toda razão p:q é admissível — ver E4 em
  estrutura_v3.md e o protocolo experimental (Predição 1A) desenhado
  para que C e B produzam previsões opostas e testáveis.

## 9. Consequência operacional

A proposta experimental do v10 (trítono 7:5) é substituída por
**Predição 1A**: um programa de quatro configurações (3:2-simétrico,
7:4-simétrico, 7:4-D3-exato em dois ramos angulares independentes,
5:2-D3-exato) desenhado para que C e B discordem, com regra de decisão
pré-registrada. Ver estrutura_v3.md, Seção V.
