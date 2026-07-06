\# Tests



Testes que validam o framework. Não fazem parte do código de produção, mas garantem 

que a refatoração não introduziu bugs em relação às versões originais publicadas.



\## Filosofia



Cada teste responde a uma pergunta clara: \*\*"a versão refatorada produz resultados

equivalentes à versão original?"\*\*



Dois tipos de teste são usados:



\- \*\*Sanity check bit-a-bit\*\*: verifica que, com mesma seed, as duas versões produzem

&#x20; outputs idênticos (até precisão de máquina). Útil quando o algoritmo é puramente

&#x20; determinístico ou quando consumo de `rand` é controlado.



\- \*\*Equivalência estatística\*\*: verifica que, ao longo de N realizações independentes,

&#x20; as duas versões produzem médias e variâncias estatisticamente indistinguíveis

&#x20; (t-test, p > 0.05). Útil quando refatoração altera ordem de consumo de `rand`

&#x20; mas mantém o algoritmo correto.



\## Estrutura



| Arquivo | Tipo | Função testada | Status |

|---|---|---|---|

| `test\_CE\_LMMSE\_unified.m` | bit-a-bit | `lib/functions/CE\_LMMSE.m` | parcial (SIM passa, BD falha por ordem de rand) |

| `test\_CE\_LMMSE\_statistical.m` | estatístico | `lib/functions/CE\_LMMSE.m` | PASS para todos os configs |



\## Como rodar



Abra o projeto MATLAB (`hmimo\_channel\_estimation.prj`) e rode:



```matlab

test\_CE\_LMMSE\_statistical

```



O output mostra média ± desvio-padrão da NMSE para cada (arch, R, mode), 

comparando versão original (`lib/legacy/`) com refatorada (`lib/functions/`).



Cada teste leva \~30-60 segundos com `parfor` (pool de 8 workers).



\## Quando rodar



\- \*\*Antes de deletar uma versão de `lib/legacy/`\*\*: confirma que a refatoração 

&#x20; capturou o algoritmo corretamente.

\- \*\*Após mudança significativa em `lib/functions/`\*\*: verifica que a alteração

&#x20; não quebrou equivalência.

\- \*\*Antes de submeter o paper\*\*: documenta que as versões usadas foram validadas.



\## Convenção de nomenclatura



\- `test\_<função>\_<tipo>.m`

\- `<tipo>` ∈ `{unified, statistical, integration}`



Exemplos futuros:

\- `test\_IDD\_MMSE\_statistical.m` (quando unificar IDD)

\- `test\_pipeline\_statistical.m` (validação end-to-end)



\## Resultados arquivados



Quando um teste passa de forma estável, registre o output no fim deste README

para histórico. Isso documenta que a validação foi feita, mesmo que o teste

seja modificado no futuro.



\### test\_CE\_LMMSE\_statistical — última execução



Data: 2026-01-15  

Configuração: N\_trials = 100, Pw\_dBm = 5

