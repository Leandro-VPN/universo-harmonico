# protocols/ — pré-registro dos protocolos experimentais

**Nenhum dado experimental foi coletado ainda.** Esta pasta contém
pré-registro, não resultado: os protocolos, os ângulos, os N e as regras de
decisão foram fixados *antes* de qualquer coleta, exatamente para que a
interpretação dos dados futuros não possa ser ajustada depois de vê-los.

## Os três documentos

- **`admissibilidade_AII_D3_v1.md`** — pré-requisito formal dos outros dois.
  Prova a Proposição P6: dadas as razões p:q do Axioma II, quais admitem
  solução angular que satisfaça também a condição de ressonância D3 (regime C).
  Resultado: só famílias com d = p−q ≥ 3. É por isso que o trítono 7:5
  (d=2), alvo do v10, foi abandonado.

- **`protocolo_discriminante_CxB_v1.md`** — protocolo C×B. Cinco configurações
  angulares desenhadas para que as regras de ressonância **C** (D3,
  comensurabilidade temporal) e **B** (D3′, simetria θ_A+θ_B=180°) façam
  previsões *opostas* e pré-especificadas, resolvendo o fork E4 por dado.
  Inclui controle negativo e regra de decisão hierárquica.

- **`axioma_0_identificabilidade_v1.md`** — protocolo A0. Matriz 3×3
  (três orientações absolutas × três alvos n=4/28/40) mantendo Δθ fixo em
  21,787°. Testa se n_mod depende da orientação **absoluta** (Axioma 0) ou
  só de Δθ (A0-nula, rotacionalmente invariante como a QM exige).

## `PROTOCOL_LOCK_v1.txt`

Congela, num único arquivo legível, os parâmetros globais (m_alvo, NPERM,
regra do bin único) e as regras de decisão de ambos os protocolos. A âncora
permanente deste estado é a **tag** `v1.0-protocolos-pre-registro`, não um
hash escrito dentro do arquivo.

Qualquer alteração posterior em configurações, alvos, N, m_alvo ou regras de
decisão exige *version bump* explícito (v1.1 para mudanças menores, v2 para
mudanças de lógica de decisão) e **nunca** pode ser aplicada retroativamente
a dados já coletados sob este lock.

## Verificação

O código que reproduz todos os números destes documentos está em
`analysis/protocols/`:

```bash
python analysis/protocols/admissibilidade/verificacao_p6.py
python analysis/protocols/cxb/verificacao_cxb.py
python analysis/protocols/a0/verificacao_a0.py
python -m pytest analysis/protocols/ -v
```

Cada script recalcula os números publicados no `.md` correspondente e os
confere linha a linha. Sai com código 0 se tudo confere e 1 se encontra
divergência — divergências são **listadas, nunca corrigidas em silêncio**,
nem no script nem no documento. Ver `analysis/protocols/README.md` para
detalhes do que cada script cobre.
