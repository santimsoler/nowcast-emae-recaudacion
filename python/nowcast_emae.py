"""
Nowcasting del EMAE con recaudación tributaria
Autor: Santiago Martínez Soler

Pipeline completo:
  1. Lee los datos de /data (EMAE s.e., CFI e IPC, recaudación por impuesto, carga manual del último mes)
  2. Deflacta y desestacionaliza cada serie (método tipo X11 con detección de outliers aditivos)
  3. Construye el factor común de recaudación (IVA DGI, Débitos y Créditos, Seguridad Social, Derechos de Importación)
  4. Backtest pseudo out-of-sample en tiempo real (ventana expandida desde 2016)
       - Modelo mensual: var. mensual del EMAE
       - Modelo 3 meses: promedio móvil de 3 meses de la var. mensual
  5. Tests de precisión predictiva (Diebold-Mariano con corrección HLN, Clark-West)
  6. Estimación final, diagnósticos y nowcast de los meses sin EMAE publicado
  7. Exporta tablas a /resultados y gráficos a /figuras

Uso:  python python/nowcast_emae.py   (desde la raíz del repo)
Requiere: numpy, pandas, scipy, matplotlib
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA, RES, FIG = (os.path.join(RAIZ, d) for d in ("data", "resultados", "figuras"))
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

INICIO_OOS = "2016-01-01"
FECHA_POST = "2020-02-01"
MESES_COVID = ["2020-03-01", "2020-04-01", "2020-05-01", "2020-06-01"]
VENTANAS_COVID_3M = pd.date_range("2020-03-01", "2020-08-01", freq="MS")


# ==============================================================================
# 1. DATOS
# ==============================================================================
def leer(nombre):
    df = pd.read_csv(os.path.join(DATA, nombre))
    df["fecha"] = pd.to_datetime(df["fecha"] + "-01")
    return df.set_index("fecha").sort_index()


emae = leer("emae_sa.csv")
cfi = leer("cfi_ipc.csv")
rec = leer("recaudacion.csv")
manual = leer("recaudacion_manual.csv")

# carga manual del último mes (informe ARCA) cuando datos.gob.ar todavía no lo tiene
for f in manual.index:
    if f not in rec.index:
        rec.loc[f] = np.nan
    for col in rec.columns:
        v = manual.loc[f, col] if col in manual.columns else np.nan
        if pd.notna(v) and v != "":
            rec.loc[f, col] = float(v)
rec = rec.sort_index()
rec["seg_social"] = rec["aportes"] + rec["contribuciones"]

# IPC encadenado desde las variaciones mensuales
cfi["ipc"] = 100 * np.cumprod(np.r_[1, 1 + cfi["delta_IPC_mensual"].values[1:] / 100])
cfi["cfi_real"] = cfi["CFI"] / cfi["ipc"]
rec = rec.join(cfi[["ipc"]], how="inner")


# ==============================================================================
# 2. DESESTACIONALIZACIÓN (aproximación tipo X11; en R se usa X13-ARIMA-SEATS)
# ==============================================================================
def _ma(x, w):
    w = np.asarray(w, float) / np.sum(w)
    h = len(w) // 2
    out = np.full(len(x), np.nan)
    for i in range(h, len(x) - h):
        out[i] = np.dot(w, x[i - h:i + h + 1])
    return out


HENDERSON13 = np.array([-0.019, -0.028, 0.0, 0.066, 0.147, 0.214, 0.240,
                        0.214, 0.147, 0.066, 0.0, -0.028, -0.019])
W35 = np.convolve(np.ones(3) / 3, np.ones(5) / 5)


def _factores_estacionales(si, n):
    S = np.full(n, np.nan)
    for m in range(12):
        idx = np.arange(m, n, 12)
        v = si[idx]
        for j, i in enumerate(idx):
            lo, hi = max(0, j - 3), min(len(v), j + 4)
            seg, ww = v[lo:hi], W35[(lo - j + 3):(hi - j + 3)]
            ok = ~np.isnan(seg)
            S[i] = np.sum(seg[ok] * ww[ok]) / np.sum(ww[ok])
    cen = pd.Series(S).rolling(13, center=True, min_periods=1).apply(
        lambda z: np.average(z, weights=np.r_[0.5, np.ones(len(z) - 2), 0.5]) if len(z) == 13 else z.mean(),
        raw=True).values
    return S - cen


def desestacionalizar(serie, z_ao=3.5):
    """Devuelve la serie s.e. sin outliers aditivos (AO) y la lista de AO detectados."""
    L = np.log(serie.values)
    n, E = len(L), 24
    g_fin = np.mean(L[-12:] - L[-24:-12])
    g_ini = np.mean(L[12:24] - L[:12])
    adelante = [L[n - 12 + (k % 12)] + g_fin * (k // 12 + 1) for k in range(E)]
    atras = np.array([L[k % 12] - g_ini * (2 - k // 12) for k in range(E)])
    X = np.r_[atras, L, adelante]
    N = len(X)
    ao = np.zeros(N, bool)
    nucleo = slice(E, E + n)
    for _ in range(3):
        T1 = _ma(X, [0.5] + [1] * 11 + [0.5])
        si = X - T1
        si[ao] = np.nan
        S1 = _factores_estacionales(si, N)
        T2 = _ma(X - S1, HENDERSON13)
        si2 = X - T2
        si2[ao] = np.nan
        S2 = _factores_estacionales(si2, N)
        SA = X - S2
        I = SA - _ma(SA, HENDERSON13)
        Ic = I[nucleo]
        ok = ~np.isnan(Ic)
        s = 1.4826 * np.median(np.abs(Ic[ok] - np.median(Ic[ok])))
        nuevo = np.zeros(N, bool)
        nuevo[nucleo] = np.abs(np.nan_to_num(Ic)) > z_ao * s
        if (nuevo == ao).all():
            break
        ao = nuevo
    sa = SA[nucleo] - np.where(ao[nucleo], I[nucleo], 0.0)
    return pd.Series(np.exp(sa), index=serie.index), list(serie.index[ao[nucleo]].strftime("%Y-%m"))


# ==============================================================================
# 3. BASE DE ESTIMACIÓN
# ==============================================================================
b = emae[["emae_sa"]].copy()
b["y"] = np.log(b["emae_sa"]).diff()
b = b.reindex(cfi.index.union(b.index))
b["infl"] = cfi["delta_IPC_mensual"].reindex(b.index) / 100

outliers = {}
sa_cfi, outliers["cfi"] = desestacionalizar(cfi["cfi_real"])
b["x_cfi"] = np.log(sa_cfi).diff().reindex(b.index)

IMPUESTOS = {"x_iva": "iva_dgi", "x_dyc": "debitos_creditos", "x_ss": "seg_social", "x_imp": "der_import"}
for x, col in IMPUESTOS.items():
    s = rec[col].dropna()
    sa, outliers[col] = desestacionalizar(s / rec.loc[s.index, "ipc"])
    b[x] = np.log(sa).diff().reindex(b.index)
X = list(IMPUESTOS)

b["y_l1"] = b["y"].shift(1)
b["y_l2"] = b["y"].shift(2)   # último EMAE publicado cuando sale la recaudación del mes t
DUM = []
for d in MESES_COVID:
    nm = "d_" + d[:7].replace("-", "_")
    DUM.append(nm)
    b[nm] = (b.index == d).astype(float)

b["y3"] = b["y"].rolling(3).mean()
b["infl3"] = b["infl"].rolling(3).mean()
D3 = []
for d in VENTANAS_COVID_3M:
    nm = "w_" + d.strftime("%Y_%m")
    D3.append(nm)
    b[nm] = (b.index == d).astype(float)


# ==============================================================================
# 4. HERRAMIENTAS
# ==============================================================================
class OLS:
    def __init__(self, df, vars_, y="y"):
        self.vars = [v for v in vars_ if df[v].var() > 0]
        self.names = ["const"] + self.vars
        self.X = np.c_[np.ones(len(df)), df[self.vars].values]
        self.yv = df[y].values
        self.XtXi = np.linalg.pinv(self.X.T @ self.X)
        self.b = self.XtXi @ self.X.T @ self.yv
        self.e = self.yv - self.X @ self.b
        self.n, self.k = self.X.shape
        self.s2 = self.e @ self.e / (self.n - self.k)

    def predict(self, row):
        x = np.r_[1, row[self.vars].values.astype(float)]
        return float(x @ self.b), float(np.sqrt(x @ (self.s2 * self.XtXi) @ x))

    def hac(self, L=3):
        u = self.X * self.e[:, None]
        S = u.T @ u
        for l in range(1, L + 1):
            G = u[l:].T @ u[:-l]
            S += (1 - l / (L + 1)) * (G + G.T)
        V = self.XtXi @ S @ self.XtXi
        se = np.sqrt(np.diag(V))
        t = self.b / se
        p = 2 * (1 - stats.t.cdf(np.abs(t), self.n - self.k))
        return pd.DataFrame({"coef": self.b, "se_HAC": se, "t": t, "p_valor": p}, index=self.names)


def factor(train, cols, excluir):
    """Primer componente principal de las variaciones estandarizadas (pesos estimados solo con train)."""
    mu, sd = train[cols].mean(), train[cols].std()
    Z = ((train[cols] - mu) / sd)[~train.index.isin(excluir)].values
    w = np.linalg.eigh(np.cov(Z.T))[1][:, -1]
    w = w if w.sum() > 0 else -w
    return (lambda df: ((df[cols] - mu) / sd).values @ w), w


def diebold_mariano(real, p1, p2, h=1):
    d = ((real - p1) ** 2 - (real - p2) ** 2).values
    T = len(d)
    stat = d.mean() / np.sqrt(np.var(d) / T) * np.sqrt((T + 1 - 2 * h + h * (h - 1) / T) / T)
    return stat, 2 * (1 - stats.t.cdf(abs(stat), T - 1))


def clark_west(real, p_chico, p_grande, L=3):
    f = ((real - p_chico) ** 2 - ((real - p_grande) ** 2 - (p_chico - p_grande) ** 2)).values
    T = len(f)
    u = f - f.mean()
    S = u @ u
    for l in range(1, L + 1):
        S += 2 * (1 - l / (L + 1)) * (u[l:] @ u[:-l])
    t = f.mean() / np.sqrt(S / T / T)
    return t, 1 - stats.norm.cdf(t)


def diagnosticos(m):
    e, n = m.e, len(m.e)
    jb = n / 6 * (stats.skew(e) ** 2 + (stats.kurtosis(e, fisher=False) - 3) ** 2 / 4)
    e2 = e ** 2
    g = np.linalg.lstsq(m.X, e2, rcond=None)[0]
    bp = n * (1 - np.sum((e2 - m.X @ g) ** 2) / np.sum((e2 - e2.mean()) ** 2))
    E = np.column_stack([np.r_[np.zeros(l), e[:-l]] for l in range(1, 4)])
    Z = np.c_[m.X, E]
    g2 = np.linalg.lstsq(Z, e, rcond=None)[0]
    bg = n * (1 - np.sum((e - Z @ g2) ** 2) / np.sum((e - e.mean()) ** 2))
    return pd.DataFrame({
        "estadistico": [jb, bp, bg],
        "p_valor": [1 - stats.chi2.cdf(jb, 2), 1 - stats.chi2.cdf(bp, m.k - 1), 1 - stats.chi2.cdf(bg, 3)]},
        index=["Jarque-Bera", "Breusch-Pagan", "Breusch-Godfrey (3)"])


def metricas(d, cols, base):
    rm = {k: np.sqrt(np.mean((d.real - d[k]) ** 2)) for k in cols}
    return pd.DataFrame({
        "RMSE_pp": {k: 100 * rm[k] for k in cols},
        "MAE_pp": {k: 100 * np.mean(np.abs(d.real - d[k])) for k in cols},
        "RMSE_rel_" + base: {k: rm[k] / rm[base] for k in cols},
        "correlacion": {k: (np.corrcoef(d.real, d[k])[0, 1] if d[k].std() > 0 else np.nan) for k in cols},
        "acierto_signo": {k: np.mean(np.sign(d[k]) == np.sign(d.real)) for k in cols},
        "sesgo_pp": {k: 100 * np.mean(d[k] - d.real) for k in cols}}).round(4)


# ==============================================================================
# 5. BACKTEST MENSUAL (información disponible: recaudación hasta t, EMAE hasta t-2)
# ==============================================================================
EXCL_M = pd.to_datetime(MESES_COVID + ["2020-07-01"])
hist = b.dropna(subset=["y", "y_l2", "x_cfi", "infl"] + X).copy()

MODELOS_M = {
    "AR_t2": ["y_l2"] + DUM,
    "CFI": ["x_cfi", "infl", "y_l2"] + DUM,
}
filas = []
for t in hist.index[hist.index >= INICIO_OOS]:
    i = hist.index.get_loc(t)
    tr = hist.iloc[:i - 1].copy()            # EMAE publicado hasta t-2
    te = hist.loc[t].copy()
    r = {"fecha": t, "real": te.y, "Cero": 0.0, "Media": tr.y.mean()}
    for m, v in MODELOS_M.items():
        r[m] = OLS(tr, v).predict(te)[0]
    proy, _ = factor(tr, X, EXCL_M)
    tr["f1"] = proy(tr)
    te["f1"] = proy(te.to_frame().T.astype(float))[0]
    r["Factor"] = OLS(tr, ["f1", "infl", "y_l2"] + DUM).predict(te)[0]
    filas.append(r)
oos_m = pd.DataFrame(filas).set_index("fecha")
COLS_M = ["Factor", "CFI", "AR_t2", "Media", "Cero"]


# ==============================================================================
# 6. BACKTEST 3 MESES
# ==============================================================================
EXCL_3 = pd.date_range("2020-03-01", "2020-09-01", freq="MS")
hist3 = b.dropna(subset=["y3", "y_l2", "infl3"] + X).copy()
hist3 = hist3[hist3.index >= "2010-06-01"]


def factor3(tr):
    proy, w = factor(tr, X, EXCL_3)
    return pd.Series(proy(b), index=b.index).rolling(3).mean(), w


filas = []
for t in hist3.index[hist3.index >= INICIO_OOS]:
    i = hist3.index.get_loc(t)
    tr = hist3.iloc[:i - 1].copy()
    te = hist3.loc[t].copy()
    f3, _ = factor3(tr)
    trf = tr.assign(f3=f3.reindex(tr.index)).dropna(subset=["f3"])
    te["f3"] = f3.loc[t]
    filas.append({"fecha": t, "real": te.y3, "Cero": 0.0,
                  "Factor_3m": OLS(trf, ["f3", "infl3", "y_l2"] + D3, y="y3").predict(te)[0],
                  "Solo_conocido": OLS(tr, ["y_l2"] + D3, y="y3").predict(te)[0]})
oos_3 = pd.DataFrame(filas).set_index("fecha")
COLS_3 = ["Factor_3m", "Solo_conocido", "Cero"]


# ==============================================================================
# 7. TABLAS DE EVALUACIÓN
# ==============================================================================
sin2020_m = oos_m[oos_m.index.year != 2020]
desde23_m = oos_m[oos_m.index >= "2023-01-01"]
sin2020_3 = oos_3[~oos_3.index.isin(pd.date_range("2020-01-01", "2021-02-01", freq="MS"))]
desde23_3 = oos_3[oos_3.index >= "2023-01-01"]

tablas = {}
for nm, d in [("mensual_sin2020", sin2020_m), ("mensual_desde2023", desde23_m)]:
    tablas[nm] = metricas(d, COLS_M, "AR_t2")
for nm, d in [("3m_sin2020", sin2020_3), ("3m_desde2023", desde23_3)]:
    tablas[nm] = metricas(d, COLS_3, "Solo_conocido")

tests = []
for nm, d in [("sin 2020", sin2020_m), ("desde 2023", desde23_m)]:
    for a, c in [("Factor", "AR_t2"), ("Factor", "Cero"), ("Factor", "CFI")]:
        s, p = diebold_mariano(d.real, d[a], d[c])
        tests.append(["mensual", nm, f"DM {a} vs {c}", s, p])
    for chico in ["AR_t2", "Cero"]:
        s, p = clark_west(d.real, d[chico], d["Factor"])
        tests.append(["mensual", nm, f"CW {chico} vs Factor", s, p])
for nm, d in [("sin 2020", sin2020_3), ("desde 2023", desde23_3)]:
    for chico in ["Solo_conocido", "Cero"]:
        s, p = clark_west(d.real, d[chico], d["Factor_3m"])
        tests.append(["3 meses", nm, f"CW {chico} vs Factor_3m", s, p])
tests = pd.DataFrame(tests, columns=["objetivo", "muestra", "test", "estadistico", "p_valor"]).round(4)


# ==============================================================================
# 8. ESTIMACIÓN FINAL Y NOWCAST
# ==============================================================================
proy_m, w_m = factor(hist, X, EXCL_M)
hist["f1"] = proy_m(hist)
mod_m = OLS(hist, ["f1", "infl", "y_l2"] + DUM)

f3_full, w_3 = factor3(hist3)
hist3f = hist3.assign(f3=f3_full.reindex(hist3.index))
mod_3 = OLS(hist3f, ["f3", "infl3", "y_l2"] + D3, y="y3")

ultimo = hist.index.max()
pred = b[(b.index > ultimo)].dropna(subset=X + ["infl"]).copy()
pred["f1"] = proy_m(pred)
y_pub = b["y"].dropna()
filas = []
estimados = {}
for t in pred.index:
    row = pred.loc[t].copy()
    t2 = t - pd.DateOffset(months=2)
    row["y_l2"] = y_pub.get(t2, estimados.get(t2, np.nan))
    mu, se = mod_m.predict(row)
    sd = np.sqrt(se ** 2 + mod_m.s2)
    estimados[t] = mu
    fila = {"mes": t.strftime("%Y-%m"), "objetivo": "var. mensual",
            "estimacion_pct": 100 * mu, "inf68_pct": 100 * (mu - sd), "sup68_pct": 100 * (mu + sd),
            "inf90_pct": 100 * (mu - 1.645 * sd), "sup90_pct": 100 * (mu + 1.645 * sd)}
    filas.append(fila)
    # promedio móvil 3 meses: aplica cuando t-2 está publicado y t-1 todavía no
    t1 = t - pd.DateOffset(months=1)
    if t2 in y_pub.index and t1 not in y_pub.index:
        r3 = b.loc[t].copy()
        r3["f3"] = f3_full.loc[t]
        r3["y_l2"] = y_pub[t2]
        mu3, se3 = mod_3.predict(r3)
        sd3 = np.sqrt(se3 ** 2 + mod_3.s2)
        filas.append({"mes": t.strftime("%Y-%m"), "objetivo": "prom. móvil 3m",
                      "estimacion_pct": 100 * mu3, "inf68_pct": 100 * (mu3 - sd3), "sup68_pct": 100 * (mu3 + sd3),
                      "inf90_pct": 100 * (mu3 - 1.645 * sd3), "sup90_pct": 100 * (mu3 + 1.645 * sd3)})
nowcast = pd.DataFrame(filas).round(2)


# ==============================================================================
# 9. EXPORTAR
# ==============================================================================
oos_m.to_csv(os.path.join(RES, "oos_mensual.csv"))
oos_3.to_csv(os.path.join(RES, "oos_3meses.csv"))
for nm, t in tablas.items():
    t.to_csv(os.path.join(RES, f"metricas_{nm}.csv"))
tests.to_csv(os.path.join(RES, "tests_precision.csv"), index=False)
mod_m.hac().round(5).to_csv(os.path.join(RES, "coeficientes_mensual.csv"))
mod_3.hac().round(5).to_csv(os.path.join(RES, "coeficientes_3meses.csv"))
diagnosticos(mod_m).round(4).to_csv(os.path.join(RES, "diagnosticos_mensual.csv"))
pd.DataFrame({"mensual": w_m, "3_meses": w_3}, index=X).round(4).to_csv(os.path.join(RES, "pesos_factor.csv"))
nowcast.to_csv(os.path.join(RES, "nowcast.csv"), index=False)

# Gráficos
d = oos_m[oos_m.index >= "2021-01-01"]
fig, ax = plt.subplots(figsize=(10, 4.3))
ax.bar(d.index, 100 * d.real, width=20, color="#a0aec0", label="EMAE publicado")
ax.plot(d.index, 100 * d.Factor, color="#c53030", lw=1.8, marker="o", ms=3, label="Nowcast fuera de muestra")
ax.axhline(0, color="grey", lw=.6)
ax.set_ylabel("Var. mensual (%)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Nowcast mensual vs EMAE publicado")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "oos_mensual.png"), dpi=130)

d = oos_3[oos_3.index >= "2021-03-01"]
fig, ax = plt.subplots(figsize=(10, 4.3))
ax.plot(d.index, 100 * d.real, color="black", lw=1.6, label="EMAE, promedio móvil 3 meses")
ax.plot(d.index, 100 * d.Factor_3m, color="#c53030", lw=1.6, ls="--", label="Nowcast fuera de muestra")
ax.axhline(0, color="grey", lw=.6)
ax.set_ylabel("Var. mensual promedio (%)")
ax.legend(frameon=False, fontsize=8)
ax.set_title("Nowcast del promedio móvil de 3 meses vs EMAE publicado")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "oos_3meses.png"), dpi=130)

nm_ = nowcast[nowcast.objetivo == "var. mensual"]
ult = hist.tail(15)
fx = pd.to_datetime(nm_.mes + "-01")
fig, ax = plt.subplots(figsize=(10, 4.3))
ax.plot(ult.index, 100 * ult.y, "o-", color="black", ms=4, label="EMAE publicado")
ax.errorbar(fx, nm_.estimacion_pct, yerr=[nm_.estimacion_pct - nm_.inf90_pct, nm_.sup90_pct - nm_.estimacion_pct],
            fmt="none", ecolor="#feb2b2", capsize=4)
ax.errorbar(fx, nm_.estimacion_pct, yerr=[nm_.estimacion_pct - nm_.inf68_pct, nm_.sup68_pct - nm_.estimacion_pct],
            fmt="D", color="#c53030", capsize=5, lw=2, label="Nowcast (68% y 90%)")
ax.axhline(0, color="grey", lw=.6)
ax.set_ylabel("Var. mensual (%)")
ax.legend(frameon=False, loc="lower left")
ax.set_title("Nowcast EMAE")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "nowcast.png"), dpi=130)

# Resumen en consola
pd.set_option("display.width", 200)
print("Outliers aditivos detectados:", outliers)
print("\nPesos del factor:\n", pd.DataFrame({"mensual": w_m, "3_meses": w_3}, index=X).round(3))
for nm, t in tablas.items():
    print(f"\n--- {nm} ---\n", t)
print("\n--- Tests ---\n", tests)
print("\n--- Coeficientes modelo mensual (HAC) ---\n", mod_m.hac().round(4).iloc[:4])
print("\n--- Coeficientes modelo 3 meses (HAC) ---\n", mod_3.hac().round(4).iloc[:4])
print("\n--- Diagnósticos modelo mensual ---\n", diagnosticos(mod_m).round(3))
print("\n--- NOWCAST ---\n", nowcast)
