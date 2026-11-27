# ¿Llegó el rebote? Nowcasting del EMAE con recaudación tributaria

Sistema de *nowcasting* de la actividad económica argentina (EMAE, INDEC) a partir de la recaudación tributaria, que se publica los primeros días hábiles del mes siguiente. El EMAE se publica unos 50 días después del cierre del mes.

Autor: Lic. Santiago Martínez Soler

## Idea

La recaudación agregada (Coparticipación Federal, CFI) es una señal ruidosa de la actividad: mezcla impuestos ligados al ciclo con otros que siguen el calendario fiscal, como los vencimientos anuales de Ganancias y Bienes Personales. El modelo usa solo los tributos ligados a transacciones y resume su información en un factor común (primer componente principal). Los tributos son:

| Tributo | Canal |
|---|---|
| IVA DGI | Consumo y producción interna |
| Débitos y Créditos bancarios | Volumen de transacciones |
| Seguridad Social (aportes + contribuciones) | Empleo y salarios registrados |
| Derechos de Importación | Importaciones de bienes |

Todas las series se deflactan por IPC y se desestacionalizan, con remoción de outliers aditivos. Entre enero 2010 y marzo 2016 se usa el IPC de San Luis (`197.1_NIVEL_GENERAL_2014_0_13`), porque el IPC nacional de ese período estuvo intervenido; desde abril 2016, el IPC del INDEC. El empalme encadena variaciones mensuales.

## Especificaciones

**Modelo mensual:** variación mensual del EMAE desestacionalizado.

`y_t = a + b·F_t + c·π_t + d·y_{t-2} + dummies COVID (mar–jun 2020) + e_t`

**Modelo de 3 meses:** promedio móvil de 3 meses de la variación mensual.

`ȳ3_t = a + b·F̄3_t + c·π̄3_t + d·y_{t-2} + dummies COVID + e_t`

Se usa `y_{t-2}` porque, cuando sale la recaudación del mes *t*, el último EMAE publicado es el de *t−2*.

## Evaluación (backtest en tiempo real, ventana expandida desde 2016)

En cada mes se estima solo con la información disponible en ese momento: los pesos del factor se reestiman en cada paso y el EMAE se usa hasta *t−2*.

**Mensual (RMSE relativo al AR con `y_{t-2}`)**

| Modelo | Sin 2020 | Desde 2023 |
|---|---|---|
| Factor de recaudación | **0,959** | **0,942** |
| CFI agregada | 0,988 | 1,045 |
| Crecimiento cero | 0,986 | 0,963 |

- Clark-West, factor contra AR: p = 0,007 (sin 2020) y p = 0,030 (desde 2023).
- Error absoluto medio: 0,88 pp. Acierta el signo en el 66 % de los meses.

**Promedio móvil de 3 meses (RMSE relativo al modelo que solo usa el mes ya publicado)**

| Modelo | Sin 2020 | Desde 2023 |
|---|---|---|
| Factor de recaudación | **0,935** | **0,908** |
| Crecimiento cero | 1,186 | 1,139 |

- Correlación nowcast–dato: 0,64. Error absoluto medio: 0,37 pp.
- Clark-West contra el modelo que solo usa el mes publicado: p = 0,005 (sin 2020) y p = 0,037 (desde 2023).

**Lectura:** la recaudación de transacciones anticipa la dirección de la actividad y aporta información significativa frente a los benchmarks. Pero la variación mensual del EMAE desestacionalizado tiene mucho ruido (desestacionalización y revisiones). El nowcast mensual tiene bandas amplias, de ±1,1 pp al 68 %, y el de 3 meses es el indicador de ciclo más preciso, con ±0,5 pp.

## Registro de nowcasts

`resultados/registro_nowcast.csv` guarda cada nowcast con la fecha de la corrida y, cuando el INDEC publica, el dato observado y el error. Es la evaluación en tiempo real del modelo.

| Corrida | Mes | Nowcast | Banda 68 % | Publicado | Error |
|---|---|---|---|---|---|
| 10/11/2026 | sep-2026 | −0,13 % | −1,24 a +0,98 | +1,1 % (24/11/2026) | +1,23 pp (dentro de la banda de 90 %) |
| 27/11/2026 | oct-2026 | +0,07 % | −1,02 a +1,16 | pendiente | — |

## Estructura

```
data/
  emae_sa.csv               EMAE desestacionalizado (índice) y fuente de cada dato
  cfi_ipc.csv               CFI nominal e inflación mensual (IPC San Luis hasta mar-2016, INDEC desde abr-2016)
  recaudacion.csv           Recaudación mensual por tributo, millones de $ (datos.gob.ar)
  recaudacion_manual.csv    Último mes cargado a mano desde el informe de ARCA
python/nowcast_emae.py      Pipeline completo (desestacionalización tipo X11)
R/nowcast_emae.R            Pipeline completo (X13-ARIMA-SEATS) + descarga de datos
resultados/                 Backtests, métricas, tests, coeficientes, nowcast y registro de nowcasts publicados
figuras/                    Gráficos
```

## Uso

Python (numpy, pandas, scipy, matplotlib):

```
python python/nowcast_emae.py
```

R (seasonal, sandwich, lmtest, tseries; readxl para actualizar datos), con la raíz del repo como directorio de trabajo:

```r
source("R/nowcast_emae.R")
```

Las dos versiones usan los mismos datos y la misma especificación, y solo difieren en la desestacionalización. La versión R usa X13 con efectos de días hábiles y Pascua, así que sus resultados van a diferir levemente de los de este README, que salen de la versión Python.

## Actualización mensual

1. En `R/nowcast_emae.R`, poner `actualizar_datos <- TRUE`. El script baja la recaudación de datos.gob.ar y el EMAE del INDEC, y guarda cada vintage del EMAE en `data/vintages/`.
2. datos.gob.ar suele tener un mes de rezago respecto del informe de ARCA. Para usar el último mes, cargar sus valores en `data/recaudacion_manual.csv`.
3. Agregar la inflación y la CFI del mes en `data/cfi_ipc.csv`.

## Series de datos.gob.ar

| Serie | ID |
|---|---|
| IVA DGI | `142.3_IVA_D_2001_M_7` |
| Débitos y Créditos | `142.3_CREDI_2001_M_24` |
| Aportes personales | `142.3_APORT_2001_M_18` |
| Contribuciones patronales | `142.3_CONTR_2001_M_25` |
| Derechos de Importación | `142.3_DEREC_2001_M_26` |
| Recaudación DGI total | `172.3_SOTAL_DDGI_M_0_0_12` |

## Limitaciones

- **El EMAE de `data/emae_sa.csv` combina tres vintages:** la serie previa hasta julio de 2025, las variaciones mensuales del informe del INDEC del 21/10/2026 para agosto y septiembre de 2025, y las del informe del 24/11/2026 desde octubre de 2025. Al actualizar desde el INDEC se reemplaza por una vintage única.
- **La evaluación es pseudo tiempo real:** usa las últimas vintages, no los datos tal como se publicaron en cada momento.
- **El último mes en `recaudacion_manual.csv` (octubre 2026) es parcial:** IVA DGI y Seguridad Social están estimados a partir de variaciones informadas en la prensa (ver la columna `fuente`).

## Referencias

- Baffigi, A., Golinelli, R. y Parigi, G. (2004). Bridge models to forecast the euro area GDP. *International Journal of Forecasting*, 20(3), 447–460.
- Clark, T. E. y West, K. D. (2007). Approximately normal tests for equal predictive accuracy in nested models. *Journal of Econometrics*, 138(1), 291–311.
- Diebold, F. X. y Mariano, R. S. (1995). Comparing predictive accuracy. *Journal of Business & Economic Statistics*, 13(3), 253–263.
- Giannone, D., Reichlin, L. y Small, D. (2008). Nowcasting: The real-time informational content of macroeconomic data. *Journal of Monetary Economics*, 55(4), 665–676.
- Harvey, D., Leybourne, S. y Newbold, P. (1997). Testing the equality of prediction mean squared errors. *International Journal of Forecasting*, 13(2), 281–291.
- Newey, W. K. y West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. *Econometrica*, 55(3), 703–708.
- Olivera, J. H. G. (1967). Money, prices and fiscal lags. *Banca Nazionale del Lavoro Quarterly Review*, 20(80), 258–268.
- Tanzi, V. (1977). Inflation, lags in collection, and the real value of tax revenue. *IMF Staff Papers*, 24(1), 154–167.
