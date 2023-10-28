#=##############################################################################
! PROGRAM SOLG_TR_SRV
!
! ## OLG model with survival probabilities
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
=##############################################################################
include("utils_dsolg.jl")

pais = "MEX"
## Carga diseño experimental
EXP_DESIGN_PATH = joinpath(abspath(pwd(), "..", "inputs","experimental_design", "datos"))
experimental_design = DataFrame(CSV.File(joinpath(EXP_DESIGN_PATH,"DSOLG_MEX_experimental_design.csv")));
experimental_design[!,"it"] = 1:size(experimental_design, 1);

experimental_design_tup = Tuple.(eachrow(experimental_design));

# Definimos un diccionario de parámetros calibrados
CALIB_DATA_PATH = EXP_DESIGN_PATH = joinpath(abspath(pwd(), "..", "..", "calibrados", "output", "csvs", pais))
df_calibrados_parametros = DataFrame(CSV.File(joinpath(CALIB_DATA_PATH, "DSOLG_probs_OLG_params_MEX.csv")))
calib_params_names = ["Omega", "alpha", "delta", "np", "nu", "gamma", "gy", "by", "tauc", "kappa"]
calib_params = Dict(param => df_calibrados_parametros[:,param][1] for param in calib_params_names)


#global exec_run_id
#exec_run_id = get_run_id()
exec_run_id = parse(Int,ARGS[1])


# Definimos un diccionario de parámetros muestreados por LHC
stress_params_names = ["Omega", "alpha", "delta", "np", "nu", "gamma", "psi_scalar", "gy", "by", "tauc", "kappa", "it"]

stress_params = Dict{Any,Any}(zip(stress_params_names, experimental_design_tup[exec_run_id]))
stress_params["pais"] = pais 

try
    DSOLG(calib_params, stress_params)
    #global exec_run_id = get_run_id()
catch e
    println("Error en la ejecución ")
    #global exec_run_id = get_run_id()
end
    

