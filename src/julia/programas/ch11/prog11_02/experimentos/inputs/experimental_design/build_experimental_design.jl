using YAML
using OrderedCollections
using DataFrames
using CSV 

using LatinHypercubeSampling

## Carga parámetros calibrados
CALIB_PATH = abspath(joinpath(pwd(), "..","parametros_calibrados"))

## Obtenemos los archivos en el directorio
calib_files = readdir(CALIB_PATH)

## Cargamos la configuración de los limites de muestro de LHS
parameters = YAML.load_file("parametros_sampling_config.yaml"; dicttype=OrderedDict{String,Any})

## Definimos una función para construir los limites de muestreo
function LHS_limits(df_params, param_limits)
    return [(df_params[!,param][1]*bounds[1],df_params[!,param][1]*bounds[2]) for (param, bounds) in param_limits]
end

## Cargamos los parámetros calibrados para cada pais
paises = ["MEX", "CHL", "CRI"]

calib_params_paises = Dict(p => LHS_limits(DataFrame(CSV.File(joinpath(CALIB_PATH, f))),parameters["parametros"]) 
                            for p in paises 
                            for f in calib_files 
                            if occursin(p,f))

## Definimos el tamaño de muestra y la cantidad de parámetros a muestrear
plan = randomLHC(parameters["simulacion"]["LHS_sampling"],length(parameters["parametros"]));

## Obtenemos el muestreo para cada país
parameters_LHC_paises = Dict( pais => DataFrame(scaleLHC(plan, bounds), parameters["parametros"].keys) for (pais, bounds) in calib_params_paises)

## Definimos un dataframe con 3 políticas
#policies = DataFrame((policy=[1, 2, 3])) 

## Hacemos un producto cruz entre el muestreo y las políticas
#for (k,v) in parameters_LHC_paises
#    parameters_LHC_paises[k] = crossjoin(v, policies)
#end


## Guardamos el diseño experimental para cada país
for (k,v) in parameters_LHC_paises
    CSV.write(joinpath(pwd(),"datos" ,"DSOLG_"*k*"_experimental_design.csv"), v)
end

#experimental_design_tup = Tuple.(eachrow(experimental_design))
