# Universo Harmônico — Bell NIST 2015
**Leandro Miguel** | ORCID: [0009-0003-4655-472X](https://orcid.org/0009-0003-4655-472X)  
Preprint Zenodo | github.com/Leandro-VPN/universo-harmonico

---

## Sobre este repositório

Scripts Python para análise dos dados públicos do experimento Bell loophole-free do NIST (Shalm et al., PRL 115, 250402, 2015).

**Dataset:** https://s3.amazonaws.com/nist-belltestdata/  
**Arquivos principais:**
- `19_45_CH_pockel_100kHz.run.nolightconeshift.alice.dat.compressed.zip` (Alice)
- `19_44_CH_pockel_100kHz.run.nolightconeshift.bob.dat.compressed.zip` (Bob)
- `19_45_CH_pockel_100kHz.run.nolightconeshift.dat.compressed.build.hdf5` (HDF5 processado)

---

## Mapeamento de canais confirmado

| Canal | Conteúdo | Alice tt_rel | Bob tt_rel |
|-------|----------|-------------|------------|
| ch=6 | SYNC (99.14 kHz) | — | — |
| ch=2 | Detector G1 | bins 90–96 | bins 347–354 |
| ch=4 | Detector G2 | bins 104–110 | bins 361–367 |
| ch=0 | GPS/overflow | descartado | descartado |

ch=2 e ch=4 são **mutuamente exclusivos** — 100% dos trials têm exatamente um ou outro.

---

## Estrutura

```
/pipeline/          Scripts de extração e cálculo de Δt
/analise/           Controles estatísticos e espectros
/diagnostico/       Identificação de canais e estrutura dos dados
README.md
```

---

## Ordem de execução recomendada

### 1. Diagnosticar os arquivos
```bash
python diagnostico/comparar_canais.py        # identifica canais
python diagnostico/estrutura_trial.py        # confirma mutualidade exclusiva
```

### 2. Extrair e calcular Δt
```bash
python pipeline/pipeline_v7.py               # pipeline principal (Windows .dat)
# ou
python pipeline/pipeline_otimizado.py        # versão Google Colab (.dat.zip)
# ou (para todos os 330M trials com retomada):
python pipeline/pipeline_checkpoints.py
```

### 3. Testes estatísticos
```bash
python analise/controle_permutacao.py        # prova Rayleigh(Δt) = artefato
python analise/ddelta_fft_acf.py             # ΔΔt, FFT, ACF vs permutado
python analise/distribuicao_nula.py          # 1000 permutações, z-score
python analise/controles_avancados.py        # surrogate Fourier, circular shift
python analise/fft_harmonico.py              # FFT da série temporal
```

---

## Resultados principais

| # | Resultado | Status |
|---|-----------|--------|
| 1 | Identidade P_Bell(θ)=1/n — erro 0,000000% | ✓ Verificado |
| 2 | Mapa de intervalos musicais via série de quintas | ✓ Verificado |
| 3 | Cenário harmônico inverso Dó+Sol, batimento=0 | ✓ Verificado |
| 4 | Pearson QRNG+ / XOR−, exceção YY | ⚠ Inconclusivo |
| 5 | LFSR 15 bits no QRNG (pico n=15, S/N=22,5×) | ✓ Verificado |
| 6 | ACF(Δt) z=187-202σ — origem hardware (Pockels cell) | ✓ Hardware |
| 7 | Rayleigh(Δt) = artefato de distribuição (ratio=1,0×) | ✓ Negativo documentado |

---

## Dependências

```
numpy scipy matplotlib h5py
```

---

## Licença

Dados: domínio público (NIST).  
Código: MIT.  
Conflito de interesses: nenhum. Financiamento: nenhum.
