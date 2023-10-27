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

## Carga diseño experimental
experimental_design = DataFrame(CSV.File("/home/milo/Documents/egap/BID/OLG/olg_respaldo/olg_egtp/src/julia/programas/ch11/prog11_02/experimentos/inputs/experimental_design/datos/DSOLG_MEX_experimental_design.csv"));
experimental_design[!,"it"] = 1:size(experimental_design, 1);

experimental_design_tup = Tuple.(eachrow(experimental_design));

# Definimos un diccionario de parámetros calibrados
df_calibrados_parametros = DataFrame(CSV.File("/home/milo/Documents/egap/BID/OLG/olg_respaldo/olg_egtp/src/julia/programas/ch11/prog11_02/calibrados/output/csvs/MEX/DSOLG_probs_OLG_params_MEX.csv"))
calib_params_names = ["Omega", "alpha", "delta", "np", "nu", "gamma", "gy", "by", "tauc", "kappa"]
calib_params = Dict(param => df_calibrados_parametros[:,param][1] for param in calib_params_names)


global exec_run_id
exec_run_id = get_run_id()


while exec_run_id <= 10000
    println(exec_run_id)

    # Definimos un diccionario de parámetros muestreados por LHC
    stress_params_names = ["Omega", "alpha", "delta", "np", "nu", "gamma", "gy", "by", "tauc", "kappa", "policy", "iteracion"]

    stress_params = Dict{Any,Any}(zip(stress_params_names, experimental_design_tup[exec_run_id]))
    stress_params["pais"] = "MEX"

    try
        DSOLG(calib_params, stress_params)
        global exec_run_id = get_run_id()
    catch e
        println("Error en la ejecución ")
        global exec_run_id = get_run_id()
    end
    
end 

