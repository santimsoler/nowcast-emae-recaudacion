# ==============================================================================
# NOWCASTING DEL EMAE CON RECAUDACIÓN TRIBUTARIA
# Autor: Santiago Martínez Soler
#
# Versión R: igual que la versión Python, pero desestacionaliza con X13-ARIMA-SEATS
# (días hábiles, Pascua y outliers AO/LS/TC; los AO se restan de la serie s.e.).
# Uso: con la raíz del repo como directorio de trabajo, source("R/nowcast_emae.R")
# ==============================================================================

library(seasonal)
library(sandwich)
library(lmtest)

actualizar_datos <- FALSE   # TRUE: baja recaudación de datos.gob.ar y EMAE del INDEC

dir.create("resultados", showWarnings = FALSE)
dir.create("figuras", showWarnings = FALSE)

inicio_oos        <- as.Date("2016-01-01")
meses_covid       <- as.Date(c("2020-03-01", "2020-04-01", "2020-05-01", "2020-06-01"))
ventanas_covid_3m <- seq(as.Date("2020-03-01"), as.Date("2020-08-01"), by = "month")
excl_m            <- c(meses_covid, as.Date("2020-07-01"))
excl_3            <- seq(as.Date("2020-03-01"), as.Date("2020-09-01"), by = "month")

# ==============================================================================
# 0. ACTUALIZACIÓN OPCIONAL DE DATOS
# ==============================================================================

if (actualizar_datos) {

  # Recaudación por impuesto (datos.gob.ar)
  ids <- c(iva_dgi          = "142.3_IVA_D_2001_M_7",
           debitos_creditos = "142.3_CREDI_2001_M_24",
           aportes          = "142.3_APORT_2001_M_18",
           contribuciones   = "142.3_CONTR_2001_M_25",
           dgi_total        = "172.3_SOTAL_DDGI_M_0_0_12",
           der_import       = "142.3_DEREC_2001_M_26")
  url <- paste0("https://apis.datos.gob.ar/series/api/series/?ids=",
                paste(ids, collapse = ","),
                "&format=csv&limit=5000&start_date=2004-01-01")
  nueva <- read.csv(url, stringsAsFactors = FALSE)
  names(nueva) <- c("fecha", names(ids))
  nueva$fecha <- substr(nueva$fecha, 1, 7)
  write.csv(nueva, "data/recaudacion.csv", row.names = FALSE)
  message("Recaudación actualizada hasta ", tail(nueva$fecha, 1))

  # EMAE desestacionalizado (INDEC, base 2004). Se guarda la vintage con fecha.
  library(readxl)
  dir.create("data/vintages", showWarnings = FALSE)
  arch <- file.path("data/vintages", paste0("emae_", Sys.Date(), ".xls"))
  download.file("https://www.indec.gob.ar/ftp/cuadros/economia/sh_emae_mensual_base2004.xls",
                arch, mode = "wb", quiet = TRUE)
  raw <- as.data.frame(read_excel(arch, sheet = 1, col_names = FALSE, .name_repair = "minimal"))
  meses_es <- c("enero", "febrero", "marzo", "abril", "mayo", "junio", "julio",
                "agosto", "septiembre", "octubre", "noviembre", "diciembre")
  limpiar <- function(z) {
    z <- gsub("[^a-z]", "", tolower(iconv(as.character(z), to = "ASCII//TRANSLIT")))
    sub("^setiembre$", "septiembre", z)
  }
  coinc   <- sapply(1:3, function(j) sum(limpiar(raw[[j]]) %in% meses_es))
  col_mes <- which.max(coinc)
  filas   <- which(limpiar(raw[[col_mes]]) %in% meses_es)
  enc     <- apply(raw[seq_len(min(filas) - 1), , drop = FALSE], 2,
                   function(z) paste(tolower(z), collapse = " "))
  col_sa  <- NA
  for (j in grep("desestacionaliz", enc)) {
    v <- suppressWarnings(as.numeric(raw[filas, j]))
    if (median(v, na.rm = TRUE) > 50) { col_sa <- j; break }
  }
  if (is.na(col_sa)) stop("No encontré la columna desestacionalizada del EMAE")
  sa <- suppressWarnings(as.numeric(raw[filas, col_sa]))
  emae_nuevo <- data.frame(
    fecha   = format(seq(as.Date("2004-01-01"), by = "month", length.out = length(filas)), "%Y-%m"),
    emae_sa = sa,
    fuente  = paste("INDEC vintage", Sys.Date()))
  emae_nuevo <- emae_nuevo[!is.na(emae_nuevo$emae_sa), ]
  write.csv(emae_nuevo, "data/emae_sa.csv", row.names = FALSE)
  message("EMAE actualizado hasta ", tail(emae_nuevo$fecha, 1))
}

# ==============================================================================
# 1. DATOS
# ==============================================================================

leer <- function(archivo) {
  df <- read.csv(file.path("data", archivo), stringsAsFactors = FALSE)
  df$fecha <- as.Date(paste0(df$fecha, "-01"))
  df[order(df$fecha), ]
}

emae <- leer("emae_sa.csv")
cfi  <- leer("cfi_ipc.csv")
rec  <- leer("recaudacion.csv")
man  <- leer("recaudacion_manual.csv")

# Carga manual del último mes (informe ARCA) cuando datos.gob.ar todavía no lo tiene
cols_rec <- setdiff(names(rec), "fecha")
for (k in seq_len(nrow(man))) {
  f <- man$fecha[k]
  if (!(f %in% rec$fecha)) {
    fila <- rec[1, ]
    fila[1, cols_rec] <- NA
    fila$fecha <- f
    rec <- rbind(rec, fila)
  }
  for (cc in intersect(cols_rec, names(man))) {
    v <- suppressWarnings(as.numeric(man[[cc]][k]))
    if (!is.na(v)) rec[rec$fecha == f, cc] <- v
  }
}
rec <- rec[order(rec$fecha), ]
rec$seg_social <- rec$aportes + rec$contribuciones

# IPC encadenado desde las variaciones mensuales
cfi$ipc      <- 100 * cumprod(c(1, 1 + cfi$delta_IPC_mensual[-1] / 100))
cfi$cfi_real <- cfi$CFI / cfi$ipc
rec <- merge(rec, cfi[, c("fecha", "ipc")], by = "fecha")

# ==============================================================================
# 2. DESESTACIONALIZACIÓN X13 (td + Pascua + outliers; se restan los AO)
# ==============================================================================

desestacionalizar <- function(fechas, valores) {
  ok <- !is.na(valores)
  fechas <- fechas[ok]
  valores <- valores[ok]
  s <- ts(valores,
          start = c(as.numeric(format(fechas[1], "%Y")), as.numeric(format(fechas[1], "%m"))),
          frequency = 12)
  aj <- tryCatch(
    seas(s, transform.function = "log", regression.aictest = c("td", "easter"),
         outlier.types = c("ao", "ls", "tc")),
    error = function(e) seas(s, transform.function = "log"))
  sa <- as.numeric(final(aj))
  co <- coef(aj)
  ao <- co[grepl("^AO", names(co))]
  efecto <- rep(0, length(sa))
  mes_en <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  for (k in seq_along(ao)) {
    nm <- sub("^AO", "", names(ao)[k])
    f  <- as.Date(sprintf("%s-%02d-01", substr(nm, 1, 4), match(substr(nm, 6, 8), mes_en)))
    efecto[fechas == f] <- ao[k]
  }
  list(datos = data.frame(fecha = fechas, sa = sa / exp(efecto)), ao = names(ao))
}

dlog <- function(d, nombre) {
  out <- data.frame(fecha = d$fecha, x = c(NA, diff(log(d$sa))))
  names(out)[2] <- nombre
  out
}

# ==============================================================================
# 3. BASE DE ESTIMACIÓN
# ==============================================================================

base <- data.frame(fecha = sort(unique(c(cfi$fecha, emae$fecha))))
base <- merge(base, emae[, c("fecha", "emae_sa")], all.x = TRUE)
base <- base[order(base$fecha), ]
base$y <- c(NA, diff(log(base$emae_sa)))
base <- merge(base, cfi[, c("fecha", "delta_IPC_mensual")], all.x = TRUE)
base$infl <- base$delta_IPC_mensual / 100

outliers <- list()
r_cfi <- desestacionalizar(cfi$fecha, cfi$cfi_real)
outliers$cfi <- r_cfi$ao
base <- merge(base, dlog(r_cfi$datos, "x_cfi"), all.x = TRUE)

impuestos <- c(x_iva = "iva_dgi", x_dyc = "debitos_creditos", x_ss = "seg_social", x_imp = "der_import")
for (x in names(impuestos)) {
  col <- impuestos[[x]]
  r <- desestacionalizar(rec$fecha, rec[[col]] / rec$ipc)
  outliers[[col]] <- r$ao
  base <- merge(base, dlog(r$datos, x), all.x = TRUE)
}
X <- names(impuestos)
base <- base[order(base$fecha), ]

lagv <- function(v, k) c(rep(NA, k), head(v, -k))
base$y_l2 <- lagv(base$y, 2)   # último EMAE publicado cuando sale la recaudación del mes t

dum <- paste0("d_", format(meses_covid, "%Y_%m"))
for (k in seq_along(meses_covid)) base[[dum[k]]] <- as.numeric(base$fecha == meses_covid[k])

base$y3    <- as.numeric(stats::filter(base$y, rep(1 / 3, 3), sides = 1))
base$infl3 <- as.numeric(stats::filter(base$infl, rep(1 / 3, 3), sides = 1))
d3 <- paste0("w_", format(ventanas_covid_3m, "%Y_%m"))
for (k in seq_along(ventanas_covid_3m)) base[[d3[k]]] <- as.numeric(base$fecha == ventanas_covid_3m[k])

cat("\n--- OUTLIERS ADITIVOS DETECTADOS ---\n")
print(outliers)

# ==============================================================================
# 4. HERRAMIENTAS
# ==============================================================================

ajustar <- function(datos, vars, y = "y") {
  usar <- vars[sapply(vars, function(v) var(datos[[v]], na.rm = TRUE) > 0)]
  lm(reformulate(usar, response = y), data = datos)
}

pred_sd <- function(m, nueva) {
  p <- predict(m, newdata = nueva, se.fit = TRUE)
  c(fit = unname(p$fit), sd = unname(sqrt(p$se.fit^2 + summary(m)$sigma^2)))
}

# Primer componente principal de las variaciones estandarizadas (pesos estimados solo con train)
factor_pc <- function(train, cols, excluir) {
  mu  <- colMeans(train[, cols])
  sdv <- apply(train[, cols], 2, sd)
  Z   <- scale(train[, cols], center = mu, scale = sdv)
  Z   <- Z[!(train$fecha %in% excluir), , drop = FALSE]
  w   <- eigen(cov(Z), symmetric = TRUE)$vectors[, 1]
  if (sum(w) < 0) w <- -w
  proy <- function(df) as.numeric(scale(df[, cols, drop = FALSE], center = mu, scale = sdv) %*% w)
  attr(proy, "pesos") <- setNames(w, cols)
  proy
}

f3_serie <- function(proy) as.numeric(stats::filter(proy(base), rep(1 / 3, 3), sides = 1))

diebold_mariano <- function(real, p1, p2, h = 1) {
  d <- (real - p1)^2 - (real - p2)^2
  T <- length(d)
  s <- mean(d) / sqrt(mean((d - mean(d))^2) / T) * sqrt((T + 1 - 2 * h + h * (h - 1) / T) / T)
  c(estadistico = s, p_valor = 2 * (1 - pt(abs(s), T - 1)))
}

clark_west <- function(real, p_chico, p_grande, L = 3) {
  f <- (real - p_chico)^2 - ((real - p_grande)^2 - (p_chico - p_grande)^2)
  T <- length(f)
  u <- f - mean(f)
  S <- sum(u^2)
  for (l in 1:L) S <- S + 2 * (1 - l / (L + 1)) * sum(u[(l + 1):T] * u[1:(T - l)])
  t <- mean(f) / sqrt(S / T / T)
  c(estadistico = t, p_valor = 1 - pnorm(t))
}

metricas <- function(d, cols, base_col) {
  rmse <- sapply(cols, function(k) sqrt(mean((d$real - d[[k]])^2)))
  data.frame(
    modelo        = cols,
    RMSE_pp       = round(100 * rmse, 4),
    MAE_pp        = round(100 * sapply(cols, function(k) mean(abs(d$real - d[[k]]))), 4),
    RMSE_relativo = round(rmse / rmse[base_col], 4),
    correlacion   = round(sapply(cols, function(k) if (sd(d[[k]]) > 0) cor(d$real, d[[k]]) else NA), 4),
    acierto_signo = round(sapply(cols, function(k) mean(sign(d[[k]]) == sign(d$real))), 4),
    sesgo_pp      = round(100 * sapply(cols, function(k) mean(d[[k]] - d$real)), 4),
    row.names = NULL)
}

# ==============================================================================
# 5. BACKTEST MENSUAL (recaudación hasta t, EMAE publicado hasta t-2)
# ==============================================================================

cat("\nBacktest mensual...\n")
hist <- base[complete.cases(base[, c("y", "y_l2", "x_cfi", "infl", X)]), ]
idx  <- which(hist$fecha >= inicio_oos)

oos_m <- data.frame(fecha = hist$fecha[idx], real = hist$y[idx],
                    Factor = NA, CFI = NA, AR_t2 = NA, Media = NA, Cero = 0)
for (j in seq_along(idx)) {
  i  <- idx[j]
  tr <- hist[1:(i - 2), ]
  te <- hist[i, , drop = FALSE]
  oos_m$AR_t2[j] <- predict(ajustar(tr, c("y_l2", dum)), te)
  oos_m$CFI[j]   <- predict(ajustar(tr, c("x_cfi", "infl", "y_l2", dum)), te)
  oos_m$Media[j] <- mean(tr$y)
  proy  <- factor_pc(tr, X, excl_m)
  tr$f1 <- proy(tr)
  te$f1 <- proy(te)
  oos_m$Factor[j] <- predict(ajustar(tr, c("f1", "infl", "y_l2", dum)), te)
}
cols_m <- c("Factor", "CFI", "AR_t2", "Media", "Cero")

# ==============================================================================
# 6. BACKTEST PROMEDIO MÓVIL 3 MESES
# ==============================================================================

cat("Backtest 3 meses...\n")
hist3 <- base[complete.cases(base[, c("y3", "y_l2", "infl3", X)]) & base$fecha >= as.Date("2010-06-01"), ]
idx3  <- which(hist3$fecha >= inicio_oos)

oos_3 <- data.frame(fecha = hist3$fecha[idx3], real = hist3$y3[idx3],
                    Factor_3m = NA, Solo_conocido = NA, Cero = 0)
for (j in seq_along(idx3)) {
  i  <- idx3[j]
  tr <- hist3[1:(i - 2), ]
  te <- hist3[i, , drop = FALSE]
  f3 <- f3_serie(factor_pc(tr, X, excl_3))
  tr$f3 <- f3[match(tr$fecha, base$fecha)]
  te$f3 <- f3[match(te$fecha, base$fecha)]
  oos_3$Factor_3m[j]     <- predict(ajustar(tr[!is.na(tr$f3), ], c("f3", "infl3", "y_l2", d3), "y3"), te)
  oos_3$Solo_conocido[j] <- predict(ajustar(tr, c("y_l2", d3), "y3"), te)
}
cols_3 <- c("Factor_3m", "Solo_conocido", "Cero")

# ==============================================================================
# 7. EVALUACIÓN
# ==============================================================================

anio <- function(f) as.numeric(format(f, "%Y"))
sin2020_m <- oos_m[anio(oos_m$fecha) != 2020, ]
desde23_m <- oos_m[oos_m$fecha >= as.Date("2023-01-01"), ]
sin2020_3 <- oos_3[!(oos_3$fecha %in% seq(as.Date("2020-01-01"), as.Date("2021-02-01"), by = "month")), ]
desde23_3 <- oos_3[oos_3$fecha >= as.Date("2023-01-01"), ]

tablas <- list(
  mensual_sin2020   = metricas(sin2020_m, cols_m, "AR_t2"),
  mensual_desde2023 = metricas(desde23_m, cols_m, "AR_t2"),
  `3m_sin2020`      = metricas(sin2020_3, cols_3, "Solo_conocido"),
  `3m_desde2023`    = metricas(desde23_3, cols_3, "Solo_conocido"))

tests <- list()
for (nm in c("sin 2020", "desde 2023")) {
  d  <- if (nm == "sin 2020") sin2020_m else desde23_m
  d3m <- if (nm == "sin 2020") sin2020_3 else desde23_3
  tests[[length(tests) + 1]] <- c("mensual", nm, "DM Factor vs AR_t2", diebold_mariano(d$real, d$Factor, d$AR_t2))
  tests[[length(tests) + 1]] <- c("mensual", nm, "DM Factor vs Cero",  diebold_mariano(d$real, d$Factor, d$Cero))
  tests[[length(tests) + 1]] <- c("mensual", nm, "DM Factor vs CFI",   diebold_mariano(d$real, d$Factor, d$CFI))
  tests[[length(tests) + 1]] <- c("mensual", nm, "CW AR_t2 vs Factor", clark_west(d$real, d$AR_t2, d$Factor))
  tests[[length(tests) + 1]] <- c("mensual", nm, "CW Cero vs Factor",  clark_west(d$real, d$Cero, d$Factor))
  tests[[length(tests) + 1]] <- c("3 meses", nm, "CW Solo_conocido vs Factor_3m",
                                  clark_west(d3m$real, d3m$Solo_conocido, d3m$Factor_3m))
  tests[[length(tests) + 1]] <- c("3 meses", nm, "CW Cero vs Factor_3m",
                                  clark_west(d3m$real, d3m$Cero, d3m$Factor_3m))
}
tests <- as.data.frame(do.call(rbind, tests), stringsAsFactors = FALSE)
names(tests) <- c("objetivo", "muestra", "test", "estadistico", "p_valor")
tests$estadistico <- round(as.numeric(tests$estadistico), 4)
tests$p_valor     <- round(as.numeric(tests$p_valor), 4)

# ==============================================================================
# 8. ESTIMACIÓN FINAL Y NOWCAST
# ==============================================================================

proy_m  <- factor_pc(hist, X, excl_m)
hist$f1 <- proy_m(hist)
mod_m   <- ajustar(hist, c("f1", "infl", "y_l2", dum))

proy_3  <- factor_pc(hist3, X, excl_3)
f3_full <- f3_serie(proy_3)
hist3$f3 <- f3_full[match(hist3$fecha, base$fecha)]
mod_3   <- ajustar(hist3[!is.na(hist3$f3), ], c("f3", "infl3", "y_l2", d3), "y3")

ultimo <- max(hist$fecha)
pred   <- base[base$fecha > ultimo & complete.cases(base[, c(X, "infl")]), ]
pred$f1 <- proy_m(pred)
y_pub  <- setNames(base$y, as.character(base$fecha))
y_pub  <- y_pub[!is.na(y_pub)]
estimados <- c()
filas <- list()
for (k in seq_len(nrow(pred))) {
  t  <- pred$fecha[k]
  t1 <- seq(t, by = "-1 month", length.out = 2)[2]
  t2 <- seq(t, by = "-1 month", length.out = 3)[3]
  fila <- pred[k, , drop = FALSE]
  fila$y_l2 <- if (as.character(t2) %in% names(y_pub)) y_pub[[as.character(t2)]] else estimados[[as.character(t2)]]
  p <- pred_sd(mod_m, fila)
  estimados[as.character(t)] <- p["fit"]
  filas[[length(filas) + 1]] <- data.frame(mes = format(t, "%Y-%m"), objetivo = "var. mensual",
    estimacion_pct = 100 * p["fit"], inf68_pct = 100 * (p["fit"] - p["sd"]), sup68_pct = 100 * (p["fit"] + p["sd"]),
    inf90_pct = 100 * (p["fit"] - 1.645 * p["sd"]), sup90_pct = 100 * (p["fit"] + 1.645 * p["sd"]))
  if (as.character(t2) %in% names(y_pub) && !(as.character(t1) %in% names(y_pub))) {
    f3row <- base[base$fecha == t, , drop = FALSE]
    f3row$f3 <- f3_full[match(t, base$fecha)]
    f3row$y_l2 <- y_pub[[as.character(t2)]]
    p3 <- pred_sd(mod_3, f3row)
    filas[[length(filas) + 1]] <- data.frame(mes = format(t, "%Y-%m"), objetivo = "prom. movil 3m",
      estimacion_pct = 100 * p3["fit"], inf68_pct = 100 * (p3["fit"] - p3["sd"]), sup68_pct = 100 * (p3["fit"] + p3["sd"]),
      inf90_pct = 100 * (p3["fit"] - 1.645 * p3["sd"]), sup90_pct = 100 * (p3["fit"] + 1.645 * p3["sd"]))
  }
}
nowcast <- do.call(rbind, filas)
nowcast[, 3:7] <- round(nowcast[, 3:7], 2)
rownames(nowcast) <- NULL

# ==============================================================================
# 9. RESULTADOS
# ==============================================================================

hac <- function(m) {
  ct <- coeftest(m, vcov = NeweyWest(m, lag = 3, prewhite = FALSE))
  data.frame(variable = rownames(ct), coef = ct[, 1], se_HAC = ct[, 2], t = ct[, 3], p_valor = ct[, 4],
             row.names = NULL)
}

write.csv(oos_m, "resultados/oos_mensual.csv", row.names = FALSE)
write.csv(oos_3, "resultados/oos_3meses.csv", row.names = FALSE)
for (nm in names(tablas)) write.csv(tablas[[nm]], paste0("resultados/metricas_", nm, ".csv"), row.names = FALSE)
write.csv(tests, "resultados/tests_precision.csv", row.names = FALSE)
write.csv(hac(mod_m), "resultados/coeficientes_mensual.csv", row.names = FALSE)
write.csv(hac(mod_3), "resultados/coeficientes_3meses.csv", row.names = FALSE)
write.csv(data.frame(variable = X, mensual = attr(proy_m, "pesos"), tres_meses = attr(proy_3, "pesos")),
          "resultados/pesos_factor.csv", row.names = FALSE)
write.csv(nowcast, "resultados/nowcast.csv", row.names = FALSE)

cat("\n--- PESOS DEL FACTOR ---\n")
print(round(attr(proy_m, "pesos"), 3))
for (nm in names(tablas)) { cat("\n---", nm, "---\n"); print(tablas[[nm]]) }
cat("\n--- TESTS ---\n")
print(tests)
cat("\n--- MODELO MENSUAL (HAC Newey-West, lag 3) ---\n")
print(coeftest(mod_m, vcov = NeweyWest(mod_m, lag = 3, prewhite = FALSE))[1:4, ])
cat("\n--- MODELO 3 MESES (HAC Newey-West, lag 3) ---\n")
print(coeftest(mod_3, vcov = NeweyWest(mod_3, lag = 3, prewhite = FALSE))[1:4, ])
cat("\n--- DIAGNÓSTICOS MODELO MENSUAL ---\n")
print(tseries::jarque.bera.test(residuals(mod_m)))
print(bptest(mod_m))
print(bgtest(mod_m, order = 3))
cat("\n--- NOWCAST ---\n")
print(nowcast)

# Gráficos
png("figuras/oos_mensual.png", width = 1000, height = 450)
d <- oos_m[oos_m$fecha >= as.Date("2021-01-01"), ]
xb <- barplot(100 * d$real, names.arg = format(d$fecha, "%Y-%m"), col = "grey75", border = NA, las = 2,
              cex.names = 0.6, ylab = "Var. mensual (%)", main = "Nowcast mensual vs EMAE publicado",
              ylim = range(c(100 * d$real, 100 * d$Factor)))
lines(xb, 100 * d$Factor, col = "red3", lwd = 2, type = "o", pch = 20)
abline(h = 0, col = "grey50")
dev.off()

png("figuras/oos_3meses.png", width = 1000, height = 450)
d <- oos_3[oos_3$fecha >= as.Date("2021-03-01"), ]
plot(d$fecha, 100 * d$real, type = "l", lwd = 2, bty = "l", xlab = "", ylab = "Var. mensual promedio (%)",
     main = "Nowcast del promedio móvil de 3 meses vs EMAE publicado")
lines(d$fecha, 100 * d$Factor_3m, col = "red3", lwd = 2, lty = 2)
abline(h = 0, col = "grey50")
legend("bottomleft", c("EMAE, prom. móvil 3m", "Nowcast fuera de muestra"),
       col = c("black", "red3"), lty = c(1, 2), lwd = 2, bty = "n")
dev.off()

png("figuras/nowcast.png", width = 1000, height = 450)
nm_ <- nowcast[nowcast$objetivo == "var. mensual", ]
fx  <- as.Date(paste0(nm_$mes, "-01"))
ult <- tail(hist, 15)
plot(ult$fecha, 100 * ult$y, type = "b", pch = 20, bty = "l", xlab = "", ylab = "Var. mensual (%)",
     xlim = c(min(ult$fecha), max(fx) + 20), ylim = range(c(100 * ult$y, nm_$inf90_pct, nm_$sup90_pct)),
     main = "Nowcast EMAE")
arrows(fx, nm_$inf90_pct, fx, nm_$sup90_pct, code = 3, angle = 90, length = 0.05, col = "pink3")
arrows(fx, nm_$inf68_pct, fx, nm_$sup68_pct, code = 3, angle = 90, length = 0.07, col = "red3", lwd = 2)
points(fx, nm_$estimacion_pct, col = "red3", pch = 18, cex = 1.8)
abline(h = 0, col = "grey50")
dev.off()

cat("\nListo: tablas en /resultados y gráficos en /figuras\n")
