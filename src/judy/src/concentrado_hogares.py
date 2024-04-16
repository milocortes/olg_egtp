## Este programa obtiene información de la tabla Gasto-Hogares y concentrado hogares

## Cargamos paquetería
import pandas as pd 
import os 

## Definimos directorios
FILE_PATH = os.getcwd()
DATA_PATH = os.path.abspath(os.path.join(FILE_PATH, "..","datos"))
OUTPUT_PATH = os.path.abspath(os.path.join(FILE_PATH, "..","output"))

CONCENTRADO_PATH = os.path.join(DATA_PATH, "enigh2022_ns_concentradohogar_csv.zip")
GASTO_HOGAR_PATH = os.path.join(DATA_PATH, "enigh2022_ns_gastoshogar_csv.zip")
POB_HOGAR_PATH = os.path.join(DATA_PATH, "enigh2022_ns_poblacion_csv.zip")

## Cargamos datos
gasto_hogar = pd.read_csv(GASTO_HOGAR_PATH)
concentrado_hogar = pd.read_csv(CONCENTRADO_PATH)
poblacion = pd.read_csv(POB_HOGAR_PATH)

## Nos quedamos con las columnas folioviv, foliohog, clave y gasto_tri
gasto_hogar_min = gasto_hogar[["folioviv", "foliohog", "clave", "gasto_tri"]]

## Nos quedamos sólo con las claves que tengan J
gasto_hogar_min = gasto_hogar_min[gasto_hogar_min["clave"].apply(lambda x : "J" in x)]

## Generamos una llave folioviv + foliohog
gasto_hogar_min["llave"] = gasto_hogar_min["folioviv"].astype(str) + "-" +gasto_hogar_min["foliohog"].astype(str)

## Sustituimos los NaN por 0.0
gasto_hogar_min["gasto_tri"] = gasto_hogar_min["gasto_tri"].replace(' ', 0.0)

## Convertimos a flotante la columna gasto_tri
gasto_hogar_min.gasto_tri = gasto_hogar_min.gasto_tri.astype(float)

## Hacemos merge con el concentrado hogar
### Generamos la llave única en concentrado_hogar
concentrado_hogar["llave"] = concentrado_hogar["folioviv"].astype(str) + "-" + concentrado_hogar["foliohog"].astype(str)

### Juntamos concentrado_hogar y las columnas ["clave", "gasto_tri", "llave"] del gasto_hogar_min
concentrado_hogar = concentrado_hogar.merge(right = gasto_hogar_min[["clave", "gasto_tri", "llave"]], how = "left", on = "llave")

## Generamos una llave folioviv + foliohog
poblacion["llave"] = poblacion["folioviv"].astype(str) + "-" + poblacion["foliohog"].astype(str)

poblacion["pop_insabi_tiene"] = 1
poblacion.loc[poblacion["pop_insabi"]==2, "pop_insabi_tiene"] = 0

poblacion_seguridad_social = poblacion[["llave","inst_1", "inst_2", "inst_3", "inst_4"]].set_index("llave").replace(" ",0.0).astype(float).sum(axis = 1).reset_index()

poblacion_seguridad_social = poblacion_seguridad_social.rename(columns = {0:"seg_sog"})

poblacion_seguridad_social.loc[poblacion_seguridad_social["seg_sog"]>0.0, "seg_sog" ] = 1.0

### Identificamos a los hogares donde al menos un integrante tenga seguridad social
hogares_al_menos_un_seg_sog = poblacion_seguridad_social.groupby("llave").sum().reset_index().query("seg_sog >0").llave.to_list()


### Creamos una variable al concentrado con los hogares que tengan al menos un integrante con seguridad social
concentrado_hogar["seguridad_social"] = 0
concentrado_hogar.loc[concentrado_hogar["llave"].isin(hogares_al_menos_un_seg_sog), "seguridad_social"] = 1

concentrado_hogar = concentrado_hogar.drop(columns = "llave")

## Guardamos datos del concentrado de hogares
concentrado_hogar.to_csv(os.path.join(OUTPUT_PATH, "enigh2022_ns_concentradohogar_extendido.csv"), index = False)
