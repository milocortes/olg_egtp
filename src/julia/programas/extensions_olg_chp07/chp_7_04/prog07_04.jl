#=##############################################################################
# PROGRAM TRHCEG_OLG
#
# ## The OLG model with uncertain survival
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr and Fabian Kindermann
#          contact@ce-fortran.com
#
=##############################################################################

using OffsetArrays
using Printf
#using DataFrames
#using StatsPlots

# model parameters
global TT      = 24
global JJ      = 3
global JR      = 3
global gamma   = 0.5
global egam    = 1.0 - 1.0/gamma
global beta    = 0.9
global alpha   = 0.3
global delta   = 0.0
global tol     = 1e-6
global damp    = 0.30
global itermax = 200

global gy = 0.0

global lsra_on = false
global smopec = false

# model variables

for param = [:w, :r, :wn, :Rn, :p, :tauw, :taur, :tauc, :taup, :by, :kappa, :n_p, :tauk, :KK, :LL, :YY, :AA, :CC, :II, :BB, :GG, :BA, :BF, :TB, :BQ, :Tpen, :TXR, :tax, :eps]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT) ;
end

for param = [:h, :omega_b]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT) ;
end

for param = [:m, :a, :c, :pen, :util, :psi, :beq, :GAM]
    @eval global $param = OffsetArray(zeros(TT+1, 3), 0:TT, 1:JJ);
end

global v = OffsetArray(zeros(TT + 2), -JJ+2:TT)

# calculates year at which age ij agent is ij_p
function year(it, ij, ijp)

    year = it + ijp - ij

    if(it == 0 || year <= 0)
        year = 0
    end

    if(it == TT || year >= TT)
        year = TT
    end

    return year
end 

# initializes variables and government parameters
function initialize()

    # set model parameters
    n_p .= 0.2
    gy = 0.20
    by .= 0.0
    kappa .= 0.0
    eps .= 0.0
    tauk .= 0.0
    tax .= 1
    lsra_on = false
    smopec = false

    # set life probs
   # psi(1, :) = 1.0
    psi[:, 1] .=1.0
   # psi(2, :) = 0.85
    psi[:, 2] .=0.85
   # psi(3, :) = 0.80
    psi[:, 3] .=0.80

    # setup for Table 7.11
    #kappa(1:TT) = 0.5
    kappa[1:TT] .= 0.5 

    # setup for Table 7.12
    #gy = 0d0
    #gy = 0.0

    #psi(2, 2) = 0.90d0
    #psi[2,2] .=0.90

    #psi(2, 3:TT) = 0.95d0
    #psi[3:TT, 2] .=0.95

    #psi(3, 3:4) = 0.85d0
    #psi[3:4, 3] .=0.85

    #psi(3, 4:TT) = 0.9d0
    #psi[4:TT, 3] .=0.90

    # initialize tax rates and pensions
    tauc .= 0.0
    tauw .= 0.0
    taur .= 0.0
    taup .= 0.0
    pen .= 0.0

    # initialize assets, LSRA payments and debt holdings
    a .= 0.0
    v .= 0.0
    beq .= 0.0
    BQ .= 0.0
    GAM .= 0.0
    YY .= 0.0
    BA .= 0.0
    BF .= 0.0
    TB .= 0.0
    TXR .= 0.0

    # human capital profile
    #h(1:JR-1) = 1
    h[JR-1:1] .= 1.0

    #h(JR:JJ) = 0
    h[JJ:JR] .= 0.0

    # set bequest distribution
    #omega_b(1) = 0.5d0
    omega_b[1] = 0.5
    
    #omega_b(2) = 0.5d0
    omega_b[2] = 0.5
    
    #omega_b(3) = 0
    omega_b[3] = 0

    # size of cohorts and bequest distribution in specific year
    for it in 0:TT

        #m(1, it) = 1
        m[it,1] = 1.0

        itm = year(it, 2, 1)

        for ij in 2:JJ
            #m(ij, it) = m(ij-1, itm)*psi(ij, it)/(1 + n_p(it))
            m[it, ij] = m[itm,ij-1]*psi[it, ij]/(1.0 + n_p[it])
        end

        #GAM_total = omega_b(1)
        GAM_total = omega_b[1]

        for ij in 2:JJ
            #GAM_total = GAM_total + omega_b(ij)*m(ij, it)
            GAM_total = GAM_total + omega_b[ij]*m[it, ij]
        end

        for ij in 1:JJ
            #GAM(ij, it) = omega_b(ij)/GAM_total
            GAM(it, ij) = omega_b[ij]/GAM_total
        end
    end

end


