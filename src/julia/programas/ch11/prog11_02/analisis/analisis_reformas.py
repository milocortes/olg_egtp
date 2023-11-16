import pandas as pd 
import os 
import glob
import matplotlib.pyplot as plt 

plt.rcParams["figure.figsize"] = (10,6)

files = glob.glob("MEX/*.csv")
datos_mex = pd.concat([pd.read_csv(f) for f in files], ignore_index = True)
datos_mex = datos_mex[datos_mex["time"].isin([0,140])].reset_index(drop=True)

## Carga diseño experimental
exp_design = pd.read_csv("DSOLG_MEX_experimental_design.csv")

## Quita algunas columnas
quita_cols = ["Omega", "delta", "tauc", "kappa"]
exp_design = exp_design.drop(columns=quita_cols)

## Carga datos calibrados
calibrados = pd.read_csv("https://raw.githubusercontent.com/milocortes/olg_egtp/main/src/julia/programas/ch11/prog11_02/experimentos/inputs/parametros_calibrados/DSOLG_probs_OLG_params_MEX.csv")
## Cambia nombre
vars_calibrados = ['alpha', 'by', 'gamma', 'gy', 'kappa', 'np', 'nu', 'psi_scalar', 'tauc']
calibrados = calibrados[vars_calibrados]
calibrados = calibrados.rename(columns={i : f"{i}_calib" for i in vars_calibrados})


datos_complete = datos_mex.merge(right=exp_design, on = "run_id", how="inner")

datos_complete = datos_complete.merge(right=calibrados,how="cross")

"""
### The expected fiscal imbalances of social security systems derived from the demographic and economic transformations.

"""
TT = 140
JJ = 16
JR = 10

# calculates year at which age ij agent is ij_p
def year(it, ij, ijp):

    year = it + ijp - ij

    if(it == 0 or year <= 0):
        year = 0
    

    if(it == TT or year >= TT):
        year = TT
    

    return year

def calcula_balance_gob(run_id : str, ):

    explora_datos = datos_complete.query(f"run_id == {run_id}")
    explora_datos = explora_datos.set_index("time")

    it = 140 
    
    # last year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)
    
    GG = explora_datos.loc[it, "gy"] * explora_datos.loc[it, "YY"]
    BB = explora_datos.loc[it, "by"] * explora_datos.loc[it, "YY"]

    # Pensiones parámetros
    m =  [1.0, 0.9057308098299159, 0.8203482998751553, 0.7430147299885191, 0.6729713331080577, 0.6095308705282791, 0.5464008614767704,
          0.4859117630099724, 0.4272118122492289, 0.36964503023419026, 0.3127621415467156, 0.2556267019520767, 0.19368194691603097, 0.1246353543574093,
          0.059456368192476856, 0.016790780547617935]

    m = {i+1 : j for i,j in enumerate(m)}
    pen = {i : explora_datos.loc[it, "kappa"]*explora_datos.loc[itm, "INC"] for i in range(JR,JJ+1)}

    PP = 0.0

    for ij in range(JR,JJ+1):
        PP += pen[ij]*m[ij]

    # calculate government expenditure
    expend = GG + (1.0+explora_datos.loc[it, "r"])*explora_datos.loc[it, "BB"]  + (1+explora_datos.loc[it, "np"])*explora_datos.loc[itp, "BB"]
    
    # Calcula brecha entre el gasto público total del gobierno considerando pago de deuda, y el ingreso tributario obtenido con la estructura tributaria del primer
    # estado estacionario  
    brecha = expend / (explora_datos.loc[0,"tauw"]*(explora_datos.loc[it,"w"]*explora_datos.loc[it,"LL"]) + explora_datos.loc[0,"tauc"]*explora_datos.loc[it,"CC"]+explora_datos.loc[0,"taur"]*explora_datos.loc[it,"r"]*explora_datos.loc[it,"AA"] )
    brecha = ((brecha)-1)*100


    # get budget balancing social security contribution
    #brecha_pensiones = explora_datos.loc[it,"PP"]/(explora_datos.loc[0,"taup"]*(explora_datos.loc[it,"w"]*explora_datos.loc[it,"LL"]))
    brecha_pensiones = PP/(explora_datos.loc[0,"taup"]*(explora_datos.loc[it,"w"]*explora_datos.loc[it,"LL"]))
    brecha_pensiones =  (brecha_pensiones-1)*100

    return (run_id, brecha, brecha_pensiones)

desbalances_fiscales = pd.DataFrame([calcula_balance_gob(i) for i in datos_complete.run_id.unique()], columns = ["run_id", "brecha_gasto", "brecha_pensiones"])

datos_complete = datos_complete.merge(right = desbalances_fiscales, on = "run_id", how = "inner")
datos_complete_plot = datos_complete.query("time==140")

import seaborn as sns
import scipy as sp



datos_complete_plot["gy"] = datos_complete_plot["gy"]*100
datos_complete_plot["gy_calib"] = datos_complete_plot["gy_calib"]*100


datos_complete_plot["by"] = datos_complete_plot["by"]*100
datos_complete_plot["by_calib"] = datos_complete_plot["by_calib"]*100

## Correlaciones entre Brecha de Gasto Público y variables : psi_scalar, gy, by

# GP vs gy
x_var_name = 'gy'

g = sns.lmplot(x=x_var_name, y='brecha_gasto', data=datos_complete_plot, line_kws={"color": "red"})
g.fig.set_figwidth(10)
g.fig.set_figheight(6)
g.set(title='Correlación entre Desequilibrio Fiscal (Sin pensiones) y Gasto Público',
      ylabel = "Desequilibrio Fiscal (sin pensiones) \n[Desviación en puntos % del gasto inicial]",
      xlabel = "Gasto Público\n[% del PIB]"
      )

def annotate(data, **kws):
    r, p = sp.stats.pearsonr(datos_complete_plot[x_var_name], data['brecha_gasto'])
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
            transform=ax.transAxes)
g.map_dataframe(annotate)
plt.legend()
plt.savefig("corr_brecha_gy.png")
plt.close()

# GP vs by 
x_var_name = 'by'
g = sns.lmplot(x=x_var_name, y='brecha_gasto', data=datos_complete_plot, line_kws={"color": "red"})
g.fig.set_figwidth(10)
g.fig.set_figheight(6)
g.set(title="Correlación entre Desequilibrio Fiscal (Sin pensiones) y Deuda Pública",
      ylabel = "Desequilibrio Fiscal (sin pensiones)\n[Desviación en puntos % del gasto inicial]",
      xlabel = "Deuda Pública\n[% del PIB]"
      )

def annotate(data, **kws):
    r, p = sp.stats.pearsonr(datos_complete_plot[x_var_name], data['brecha_gasto'])
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
            transform=ax.transAxes)
g.map_dataframe(annotate)
plt.savefig("corr_brecha_by.png")
plt.close()

# GP vs psi_scalar 
x_var_name = 'psi_scalar'
g = sns.lmplot(x=x_var_name, y='brecha_gasto', data=datos_complete_plot, line_kws={"color": "red"})

g.fig.set_figwidth(10)
g.fig.set_figheight(6)
g.set(title="Correlación entre Desequilibrio Fiscal (Sin pensiones) e incrementos en las\nProbabilidades de Supervivencia entre Cohortes",
      ylabel = "Desequilibrio Fiscal (Sin pensiones)\n[Desviación en puntos % del gasto inicial]",
      xlabel = "Aumento en la Probabilidad de Supervivencia entre Cohortes\n[valor real]"
      )

def annotate(data, **kws):
    r, p = sp.stats.pearsonr(datos_complete_plot[x_var_name], data['brecha_gasto'])
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
            transform=ax.transAxes)
g.map_dataframe(annotate)
plt.savefig("corr_brecha_psi.png")
plt.close()


"""
### Changes to social security and social protection programs to minimize these fiscal imbalances.

"""

### Preparamos funciones a minimizar

def brecha_pensiones_min(X) -> float:
    kappa = X[0]

    explora_datos = datos_complete.query(f"run_id == {run_id}")
    explora_datos = explora_datos.set_index("time")


    it = 140 
    
    # last year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    m =  [1.0, 0.9057308098299159, 0.8203482998751553, 0.7430147299885191, 0.6729713331080577, 0.6095308705282791, 0.5464008614767704,
          0.4859117630099724, 0.4272118122492289, 0.36964503023419026, 0.3127621415467156, 0.2556267019520767, 0.19368194691603097, 0.1246353543574093,
          0.059456368192476856, 0.016790780547617935]

    m = {i+1 : j for i,j in enumerate(m)}
    pen = {i : kappa*explora_datos.loc[itm, "INC"] for i in range(JR,JJ+1)}

    PP = 0.0

    for ij in range(JR,JJ+1):
        PP += pen[ij]*m[ij]
    # get budget balancing social security contribution
    brecha_pensiones = PP/(explora_datos.loc[0,"taup"]*(explora_datos.loc[it,"w"]*explora_datos.loc[it,"LL"]))
    brecha_pensiones =  (brecha_pensiones-1)*100

    return abs(brecha_pensiones)

def brecha_gasto_min(X):
    explora_datos = datos_complete.query(f"run_id == {run_id}")
    explora_datos = explora_datos.set_index("time")

    gy, by = X

    it = 140 
    
    # last year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    GG = gy * explora_datos.loc[it, "YY"]
    BB = by * explora_datos.loc[it, "YY"]

    # calculate government expenditure
    expend = GG + (1.0+explora_datos.loc[it, "r"])*explora_datos.loc[it, "BB"]  + (1+explora_datos.loc[it, "np"])*explora_datos.loc[itp, "BB"]
    
    # Calcula brecha entre el gasto público total del gobierno considerando pago de deuda, y el ingreso tributario obtenido con la estructura tributaria del primer
    # estado estacionario  
    brecha = expend / (explora_datos.loc[0,"tauw"]*(explora_datos.loc[it,"w"]*explora_datos.loc[it,"LL"]) + explora_datos.loc[0,"tauc"]*explora_datos.loc[it,"CC"]+explora_datos.loc[0,"taur"]*explora_datos.loc[it,"r"]*explora_datos.loc[it,"AA"] )
    brecha = ((brecha)-1)*100

    return abs(brecha)

##### Mimizamos las funciones para cada run_id

#### Brecha Gasto

### PSO paralelo
from  optimization_algorithms import PSO
import numpy as np 

# 1 - Define the upper and lower bounds of the search space
n_dim =  2              # Number of dimensions of the problem
lb = [0.04562, 0.2536]  # lower bound for the search space
ub =  [0.31842999999999999, 0.6304] # upper bound for the search space
lb = np.array(lb)
ub = np.array(ub)

# 2 - Define the parameters for the optimization
max_iters = 60  # maximum number of iterations

# 3 - Parameters for the algorithm
pop_size = 100   # number of individuals in the population
# Cognitive scaling parameter
α = 8
# Social scaling parameter
β = 8

# velocity inertia
w = 0.5
# minimum value for the velocity inertia
w_min = 0.4
# maximum value for the velocity inertia
w_max = 0.9

# 4 - Define the cost function

# 5 - Run the PSO algorithm
# call DE

acumula_resultados_brecha_gasto = []

for run_id in datos_complete_plot.run_id:

    print(run_id)

    try:
        fitness_values, best_vector ,mejor_valor = PSO(brecha_gasto_min,pop_size, max_iters, lb, ub,α,β,w,w_max,w_min)

        acumula_resultados_brecha_gasto.append(
            [run_id] + best_vector + [mejor_valor]
        )
    except:
        print("fallo")

df_resultados_brecha_gasto = pd.DataFrame(acumula_resultados_brecha_gasto, columns = ["run_id", "gy_minimal", "by_minimal", "brecha_gasto_min_val"])

### Brecha pensiones
# 1 - Define the upper and lower bounds of the search space
n_dim =  1              # Number of dimensions of the problem
lb = [0.3215]  # lower bound for the search space
ub =  [0.9645] # upper bound for the search space
lb = np.array(lb)
ub = np.array(ub)

# 2 - Define the parameters for the optimization
max_iters = 40  # maximum number of iterations

# 3 - Parameters for the algorithm
pop_size = 100   # number of individuals in the population
# Cognitive scaling parameter
α = 8
# Social scaling parameter
β = 8

# velocity inertia
w = 0.5
# minimum value for the velocity inertia
w_min = 0.4
# maximum value for the velocity inertia
w_max = 0.9

# 4 - Define the cost function

# 5 - Run the PSO algorithm
# call DE

acumula_resultados_brecha_pensiones = []

for run_id in datos_complete_plot.run_id:

    print(run_id)

    try:
        fitness_values, best_vector ,mejor_valor = PSO(brecha_pensiones_min,pop_size, max_iters, lb, ub,α,β,w,w_max,w_min)

        acumula_resultados_brecha_pensiones.append(
            [run_id] + best_vector + [mejor_valor]
        )
    except:
        print("fallo")

df_resultados_brecha_pensiones = pd.DataFrame(acumula_resultados_brecha_pensiones, columns = ["run_id", "kappa_minimal", "brecha_pensiones_min_val"])

## Recolecta valores
df_resultados_brecha_gasto = pd.read_csv("df_resultados_brecha_gasto.csv")
df_resultados_brecha_pensiones = pd.read_csv("df_resultados_brecha_pensiones.csv")

datos_complete_plot = pd.concat([datos_complete_plot.set_index("run_id"),df_resultados_brecha_gasto.set_index("run_id"), df_resultados_brecha_pensiones.set_index("run_id")], axis = 1)
datos_complete_plot.reset_index(inplace= True)

datos_complete_plot["gy_minimal"] = datos_complete_plot["gy_minimal"]*100
datos_complete_plot["by_minimal"] = datos_complete_plot["by_minimal"]*100


datos_complete_plot[["brecha_gasto"]].rename(columns = {"brecha_gasto":"Desequilibrio Fiscal (Sin pensiones)"}).plot.kde()
plt.xlabel("Desequilibrio Fiscal (Sin pensiones)\n[Desviación en puntos % del gasto inicial]")
plt.title("Densidad Desequilibrio Fiscal (Sin pensiones)\n[Desviación en puntos % del gasto inicial]")
plt.savefig("densidad_brecha_gasto.png")
plt.close()

datos_complete_plot[["brecha_pensiones"]].rename(columns = {"brecha_presupuesto":"Desequilibrio Fiscal del Presupuesto de Pensiones"}).plot.kde()
plt.xlabel("Desequilibrio Fiscal del Presupuesto de Pensiones\n[Desviación en puntos % del gasto inicial en pensiones]")
plt.title("Densidad Desequilibrio Fiscal del Presupuesto de Pensiones\n[Desviación en puntos % del gasto inicial en pensiones]")
plt.savefig("densidad_brecha_pensiones.png")
plt.close()

datos_complete_plot[["gy_minimal","gy"]].rename(columns={"gy_minimal" : "Gasto Público (minimiza desequilibrio fiscal)", "gy" : "Gasto Público (estado de equilibrio)"}).plot.kde()
plt.axvline(x = datos_complete_plot["gy_calib"].mean(), color = "r", label = "Baseline")
plt.title("Densidades del gasto público de equilibrio vs gasto público que minimiza el desequilibrio fiscal (sin pensiones))\n[n=825]")
plt.xlabel("Gasto Público\n[% del PIB]")
plt.legend()
plt.savefig("densidad_gy_minimal.png")
plt.close()

datos_complete_plot[["by_minimal","by"]].rename(columns={"by_minimal" : "Deuda Pública (minimiza desequilibrio fiscal)", "by" : "Deuda Pública (estado de equilibrio)"}).plot.kde(xlabel = "Deuda Pública\n[% del PIB]")
plt.axvline(x = datos_complete_plot["by_calib"].mean(), color = "r", label = "Baseline")
plt.title("Densidades de la deuda pública de equilibrio vs deuda pública que minimiza el desequilibrio fiscal(sin pensiones)\n[n=825]")
plt.xlabel("Deuda Pública\n[% del PIB]")
plt.legend()
plt.savefig("densidad_by_minimal.png")
plt.close()


datos_complete_plot[["kappa"]].rename(columns = {"kappa":"Tasa de Reemplazo"}).plot.kde()
plt.axvline(x = datos_complete_plot.kappa_minimal.describe().loc["mean"], color = 'r', label = 'Tasa de Reemplazo que corrige desequilibrio')
plt.axvline(x = datos_complete_plot.kappa_calib.describe().loc["mean"], color = 'orange', label = 'Baseline')
plt.title("Densidades de la Tasa de Reemplazo de las ejecuciones\n[n=825]")
plt.xlabel("Tasa de reemplazo")
plt.legend()
plt.savefig("densidad_kappa_minimal.png")
plt.close()


"""
### The impacts in the labor markets derived from the transitions and policy changes.

"""

### Horas trabajadas
mean = datos_complete_plot[[f"l_coh_cohort_{i}" for i in range(1,10)]].describe().T["mean"].to_numpy()
std = datos_complete_plot[[f"l_coh_cohort_{i}" for i in range(1,10)]].describe().T["std"].to_numpy()
x = [20 + 5*i  for i in range(1,10)]
plt.plot(x, mean, 'k', color='#2d4a5a', label = "Media")
horas_baseline = [0.2558668454681908,0.28652452367020537, 0.29055257297256587, 0.243667565216225, 0.19572211872823342, 0.12286360110873185, 0.06577582885098014, 0.024718039237347906, 0.014518216809766511]
plt.plot(x, horas_baseline, color = "red",label = "Baseline")
plt.fill_between(x, mean-std, mean+std,
    alpha=0.5, edgecolor='#2d4a5a', facecolor='#6fcbe2',
    linewidth=2, antialiased=True, label = "STD")
plt.xlabel("Ciclo de vida [años]")
plt.ylabel("Horas trabajadas")
plt.title("Horas Trabajadas para cada Cohorte [n = 825]")
plt.legend()
plt.savefig("horas_trabajadas.png")
plt.close()


### Ingresos
mean = datos_complete_plot[[f"y_coh_cohort_{i}" for i in range(1,10)]].describe().T["mean"].to_numpy()
std = datos_complete_plot[[f"y_coh_cohort_{i}" for i in range(1,10)]].describe().T["std"].to_numpy()
x = [20 + 5*i  for i in range(1,10)]
plt.plot(x, mean, 'k', color='#2d4a5a', label = "Media")
ingresos_baseline =  [0.3203104321369967, 0.49721692380595756, 0.6558581312729471, 0.6368784267541635, 0.6005318411981384, 0.4534704719586792, 0.3109948934385704, 0.16781636886462678, 0.11789336413885833]
plt.plot(x, ingresos_baseline, color = "red", label = "Baseline")
plt.fill_between(x, mean-std, mean+std,
    alpha=0.5, edgecolor='#2d4a5a', facecolor='#6fcbe2',
    linewidth=2, antialiased=True, label = "STD")
plt.xlabel("Ciclo de vida [años]")
plt.ylabel("Ingreso Laboral")
plt.title("Ingreso Laboral para cada Cohorte [n = 825]")
plt.legend()
plt.savefig("ingresos.png")
plt.close()

### Riqueza
mean = datos_complete_plot[[f"a_coh_cohort_{i}" for i in range(1,17)]].describe().T["mean"].to_numpy()
std = datos_complete_plot[[f"a_coh_cohort_{i}" for i in range(1,17)]].describe().T["std"].to_numpy()
x = [20 + 5*i  for i in range(1,17)]
plt.plot(x, mean, 'k', color='#2d4a5a', label = "Media")
riqueza_baseline = [ 0.0, 0.05177910963420018, 0.14828132401030525, 0.30599520416004083, 0.5012241128562037, 0.7443471949845972, 1.025235433340271, 1.2162882888010573,1.3748669949360148,1.4724400540225042,1.6989136853813587,1.8718257675561403,1.9277807935399525,1.78609811686002,1.385376221860643, 0.742005102717571]
plt.plot(x, riqueza_baseline, color = "r", label = "Baseline")
plt.fill_between(x, mean-std, mean+std,
    alpha=0.5, edgecolor='#2d4a5a', facecolor='#6fcbe2',
    linewidth=2, antialiased=True, label = "STD")
plt.xlabel("Ciclo de vida [años]")
plt.ylabel("Riqueza")
plt.title("Riqueza para cada Cohorte [n = 825]")
plt.legend()
plt.savefig("riqueza.png")
plt.close()

"""
### The fiscal effect of policies that promote productive aging and allow exploiting the employment generation for older persons, particularly the care economy.
"""

"""
### The labor market changes generated by migration incentives.
"""

"""
### The key political economy elements increase the chances of successful social security and social protection reforms.

"""

# Declare feature vector and Meta variable
datos_complete_plot["Meta"] = 0
datos_complete_plot.loc[(datos_complete_plot["kappa"] > datos_complete_plot["kappa_calib"].mean()) & ((datos_complete_plot["PP"]/datos_complete_plot["YY"]) < 0.0496),"Meta"] = 1


mi_scatter = plt.scatter(datos_complete_plot["kappa"],(datos_complete_plot["PP"]/datos_complete_plot["YY"])*100, c = datos_complete_plot["Meta"], cmap='Paired')
plt.xlabel("Tasa de reemplazo")
plt.ylabel("Pago a Pensiones \n[% PIB]")
plt.title("Pago a pensiones vs Tasa de reemplazo")
plt.axhline(y = 4.96, color = 'red', linestyle = '--')
plt.axvline(x = datos_complete_plot["kappa_calib"].mean(), color = "r",  linestyle = '--')
plt.legend(handles=mi_scatter.legend_elements()[0], labels=["Desinterés","Interés"])
plt.savefig("puntos_interes.png")
plt.close()

X = datos_complete_plot[["r", "rn", "w", "wn", "p", "KK", "AA", "BB", "LL", "HH", "YY", "CC", "II", "GG", "INC", "BQ", "tauc", "tauw", "taur", "taup", "PP", "alpha", "np", "nu", "gamma", "psi_scalar", "gy", "by"]]
y = datos_complete_plot["Meta"]



# import XGBoost
import xgboost as xgb

# definimos data_dmatrix
data_dmatrix = xgb.DMatrix(data=X,label=y)

# dividimos y y X en los conjuntos de entrenamiento y prueba
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 0)

# Entrenamos el clasificador
from xgboost import XGBClassifier
# definimos hiperparámetros
params = {
            'max_depth': 4,
            'alpha': 10,
            'learning_rate': 1.0,
            'n_estimators':100
        }
# instanciamos el clasificador
xgb_clf = XGBClassifier(**params)

# ajustamos el clasificador en los datos de entrenamiento
xgb_clf.fit(X_train, y_train)

# realizamos el pronóstico en el conjunto de entrenamiento
y_pred = xgb_clf.predict(X_test)

# verificamos accuracy
from sklearn.metrics import  accuracy_score
print('XGBClassifier model accuracy score: {0:0.4f}'.format(accuracy_score(y_test, y_pred)))

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
cm = confusion_matrix(y_test, y_pred, labels=xgb_clf.classes_,normalize='true')
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=xgb_clf.classes_)
disp.plot(cmap=plt.cm.Blues)
plt.title("Matriz de confusión")
plt.savefig("matriz_conf.png")
plt.close()

# Validamos con validación cruzada 
from xgboost import cv

params = {'colsample_bytree': 0.3,'learning_rate': 0.1,
                'max_depth': 5, 'alpha': 10}

xgb_cv = cv(dtrain=data_dmatrix, params=params, nfold=3,
                    num_boost_round=50, early_stopping_rounds=10, metrics="auc", as_pandas=True, seed=123)

# calculamos la importancia de las caracterísiticas
nombres = ["Elasticidad del Capital", "Tasa de contribución\nsistema de pensiones", "Consumo", "Deuda Pública", "Delta probabilidad \nde supervivencia", "Producto"]
nombres.reverse()
xgb.plot_importance(xgb_clf, max_num_features = 10).set_yticklabels(nombres)
plt.savefig("importancia.png")
plt.close()


features_importance = ["alpha", "CC", "taup", "BB", "Meta"]
labels_name = ["Elasticidad del Capital", "Consumo", "Tasa de contribución\nsistema de pensiones", "Deuda Pública"]
g = sns.pairplot(datos_complete_plot[features_importance], hue="Meta", markers=["o", "s"])
g.axes[0][0].set_ylabel("Elasticidad del Capital")
g.axes[-1][0].set_xlabel("Elasticidad del Capital")

g.axes[1][0].set_ylabel("Consumo")
g.axes[-1][1].set_xlabel("Consumo")

g.axes[2][0].set_ylabel("Tasa de contribución\nsistema de pensiones")
g.axes[-1][2].set_xlabel("Tasa de contribución\nsistema de pensiones")

g.axes[3][0].set_ylabel("Deuda Pública")
g.axes[-1][3].set_xlabel("Deuda Pública")

g.fig.suptitle("Pairplot de 4 características con mayor importancia", verticalalignment = "bottom")