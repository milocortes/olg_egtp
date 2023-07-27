import pandas as pd 
import glob 

experimental_design = pd.read_csv("experimental_design.csv")
add_time = lambda x : x.reset_index().rename(columns = {"index":"time"})

resultados = pd.concat(
                    [add_time(pd.read_csv(f)) for f in glob.glob("csvs/*.csv")], 
                    ignore_index = True
                    )

expe_resultados = resultados.merge(right=experimental_design, how = "left", on = "id_exp")                

expe_resultados.to_csv("expe_resultados_capital_humano_OLG.csv", index = False)