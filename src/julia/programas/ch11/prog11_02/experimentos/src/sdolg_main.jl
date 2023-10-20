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
pais = "MEX"
experimental_design = DataFrame(CSV.File("/home/milo/Documents/egap/BID/OLG/olg_respaldo/olg_egtp/src/julia/programas/ch11/prog11_02/experimentos/inputs/experimental_design/datos/DSOLG_MEX_experimental_design.csv"));
experimental_design[!,"it"] = 1:size(experimental_design, 1);

experimental_design_tup = Tuple.(eachrow(experimental_design));
DSOLG(experimental_design_tup[1]...)


