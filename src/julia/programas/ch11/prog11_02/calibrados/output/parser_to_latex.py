import pandas as pd
from jinja2 import Environment, FileSystemLoader

import os
import glob 

## Configuración de jinja 
file_loader = FileSystemLoader('.')
env = Environment(loader=file_loader)

template = env.get_template('template')

## Carga resultados 
paises = ["MEX", "CHL", "CRI"]
temas = ["labour_market", "pension_system", "params", "capital_market", "good_market", "gov_accounts"] 

datos_paises = {p : {m :{} for m in ["con_probs", "sin_probs"] } for p in paises }
files = glob.glob("csvs/*/*.csv")


for pais in paises:
    for tema in temas:
        for dato in files:
            if pais in dato and tema in dato:
                if "probs" in dato:
                    datos_paises[pais]["con_probs"][tema] = pd.read_csv(dato) 
                    datos_paises[pais]["con_probs"]["params"] = pd.read_csv(os.path.join("csvs",pais,f"DSOLG_probs_OLG_params_{pais}.csv")) 
                else:
                    datos_paises[pais]["sin_probs"][tema] = pd.read_csv(dato) 

## Carga datos observados
datos_observados = pd.read_csv("../datos/parametros_olg.csv")

#### Construye diccionarios con datos
def get_datos(modelo, pais):
    resultados_DSOLG = {}

    ## Consumo Privado % PIB
    # Modelo
    resultados_DSOLG["consumo_privado_modelo"] = datos_paises[pais][modelo]["good_market"]["C"].iloc[1].round(2)
    # Observado
    resultados_DSOLG["consumo_privado_observado"] = (datos_observados.query(f"countrycode=='{pais}'")["csh_c"]*100).values[0].round(2)

    # Gasto Público % PIB 
    # Modelo
    resultados_DSOLG["gasto_publico_modelo"] = datos_paises[pais][modelo]["good_market"]["G"].iloc[1].round(2)
    # Observado
    resultados_DSOLG["gasto_publico_observado"] = (datos_observados.query(f"countrycode=='{pais}'")["csh_g"]*100).values[0].round(2)

    # Inversión % PIB 
    # Modelo
    resultados_DSOLG["inversion_modelo"] = datos_paises[pais][modelo]["good_market"]["I"].iloc[1].round(2)
    # Observado
    resultados_DSOLG["inversion_observado"] = (datos_observados.query(f"countrycode=='{pais}'")["csh_i"]*100).values[0].round(2)

    # Impuesto al Consumo
    # Modelo
    resultados_DSOLG["tax_consumo_modelo"] = datos_paises[pais][modelo]["gov_accounts"]["TAUC"].iloc[2].round(2)
    # Observado
    resultados_DSOLG["tax_consumo_observado"] = (datos_observados.query(f"countrycode=='{pais}'")["vat_tax"]*100).values[0].round(2)

    # Impuesto al Ingreso
    # Modelo
    resultados_DSOLG["tax_ingreso_modelo"] = datos_paises[pais][modelo]["gov_accounts"]["TAUW"].iloc[2].round(2)

    # Impuesto al Ingreso Medio (Observado)
    # Observado
    resultados_DSOLG["tax_ingreso_observado_min"] = (datos_observados.query(f"countrycode=='{pais}'")["tax_w_ocde_min"]*100).values[0].round(2)
    # Impuesto al Ingreso Máximo (Observado)
    # Observado
    resultados_DSOLG["tax_ingreso_observado_max"] = (datos_observados.query(f"countrycode=='{pais}'")["tax_w_ocde_max"]*100).values[0].round(2)
    # Impuesto al Ingreso Mínimo (Observado)
    # Observado
    resultados_DSOLG["tax_ingreso_observado_medio"] = (datos_observados.query(f"countrycode=='{pais}'")["tax_w_ocde_medio"]*100).values[0].round(2)
    
    # Ingresos por impuesto al consumo
    # Modelo
    resultados_DSOLG["recaudacion_tax_consumo_modelo"] = datos_paises[pais][modelo]["gov_accounts"]["TAUC"].iloc[1].round(2)
    # Observado
    resultados_DSOLG["recaudacion_tax_consumo_observado"] = (datos_observados.query(f"countrycode=='{pais}'")["consump_tax_rev_per_gdp"]*100).values[0].round(2)

    # Ingresos por impuesto al ingreso
    # Modelo
    resultados_DSOLG["recaudacion_tax_ingreso_modelo"] = datos_paises[pais][modelo]["gov_accounts"]["TAUW"].iloc[1].round(2)
    # Observado
    resultados_DSOLG["recaudacion_tax_ingreso_observado"] = (datos_observados.query(f"countrycode=='{pais}'")["income_tax_rev_per_gdp"]*100).values[0].round(2)

    # Endeudamiento Público
    # Modelo
    resultados_DSOLG["endeudamiento_modelo"] = datos_paises[pais][modelo]["gov_accounts"]["B"].iloc[1].round(2)
    # Observado
    resultados_DSOLG["endeudamiento_observado"] = datos_paises[pais][modelo]["gov_accounts"]["B"].iloc[1].round(2)

    # Costo necesitado para mantener el nivel de deuda constante
    # Modelo
    resultados_DSOLG["deuda_constante"] = datos_paises[pais][modelo]["params"]["deuda_necesaria"].values[0].round(2)

    # Tasa de Reemplazo
    # Modelo
    resultados_DSOLG["kappa_modelo"] =  datos_paises[pais][modelo]["params"]["kappa"].values[0].round(2)
    # Observado
    resultados_DSOLG["kappa_observado"] =  datos_paises[pais][modelo]["params"]["kappa"].values[0].round(2)

    # Pagos a pensiones (%PIB)
    # Modelo
    resultados_DSOLG["pago_pensiones_modelo"] = datos_paises[pais][modelo]["pension_system"]["PP"].iloc[1].round(2)
    # Observado
    resultados_DSOLG["pago_pensiones_observado"] = (datos_observados.query(f"countrycode=='{pais}'")["public_expenditure_pensions_perc_gdp"]*100).values[0].round(2)

    # Tasa de crecimiento poblacional
    # Modelo
    resultados_DSOLG["natalidad_modelo"] = (datos_paises[pais][modelo]["params"]["np"]*100).values[0].round(2)
    # Observado
    resultados_DSOLG["natalidad_observada"] = (datos_paises[pais][modelo]["params"]["np"]*100).values[0].round(2)

    return resultados_DSOLG

resultados_to_jinja = {pais : {modelo : get_datos(modelo, pais) for modelo in ["con_probs", "sin_probs"]} for pais in paises}


output = template.render(resultados = resultados_to_jinja)

text_file = open("latex_template/slides.tex", "w")
text_file.write(output)
text_file.close()