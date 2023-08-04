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

global gy = 0.20

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
    @eval global $param = OffsetArray(zeros(3, TT+1), 1:JJ, 0:TT);
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
    tax .= 1.0
    lsra_on = false
    smopec = false

    # set life probs
    psi[1, :] .= 1.0
    psi[2, :] .= 0.85
    psi[3, :] .= 0.80


    # setup for Table 7.11
    kappa[1:TT] .= 0.5


   # setup for Table 7.12
    #gy = 0d0
    #gy = 0.0

    #psi(2, 2) = 0.90d0
    #psi[2,2] .=0.90

    #psi(2, 3:TT) = 0.95d0
    #psi[2, 3:TT] .=0.95

    #psi(3, 3:4) = 0.85d0
    #psi[3, 3:4] .=0.85

    #psi(3, 4:TT) = 0.9d0
    #psi[3, 4:TT] .=0.90


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
    h[1:JR-1] .= 1.0
    h[JR:JJ] .= 0.0

    # set bequest distribution
    omega_b[1] = 0.5
    omega_b[2] = 0.5
    omega_b[3] = 0.0

    # size of cohorts and bequest distribution in specific year
    for it in 0:TT
        m[1, it] = 1.0
        itm = year(it, 2, 1)
        for ij in 2:JJ
            m[ij, it]= m[ij-1, itm]*psi[ij, it]/(1.0 + n_p[it])
        end

        GAM_total = omega_b[1]
        for ij in 2:JJ
            GAM_total= GAM_total + omega_b[ij]*m[ij, it]
        end

        for ij in 1:JJ
            GAM[ij, it]= omega_b[ij]/GAM_total
        end
    end

end

# subroutine for calculating factor prices in a certain year
function factor_prices(it)

    # factor prices and pension payments in year t
    if(smopec && it > 0)
        r[it] = r[0]
    else
        r[it] = (1.0-tauk[it])/(1.0-eps[it]*tauk[it])*(alpha*(KK[it]/LL[it])^(alpha-1.0)-delta)
    end

    w[it] = (1.0-alpha)*(KK[it]/LL[it])^alpha

    # derive after tax prices
    wn[it] = w[it]*(1.0-tauw[it]-taup[it])
    Rn[it] = 1.0 + r[it]*(1.0-taur[it])
    p[it] = 1.0 + tauc[it]

end 

# present value of resources for household aged ij in year it
function get_W(ij, it)

    # get current value of resources
    get_W = wn[it]*h[ij] + beq[ij, it] + pen[ij, it]

    if(it == 1 && ij > 1)
        get_W = get_W + Rn[it]*a[ij, it] + v[-ij+2]
    end

    if(it >= 1 && ij == 1)
        get_W = get_W + v[it]
    end    

    # iterate over remainder of life span
    PRn = 1.0
    for ijp in ij+1:JJ
        itp = year(it, ij, ijp)
        PRn = PRn*Rn[itp]
        get_W = get_W + (wn[itp]*h[ijp] + beq[ijp, itp] + pen[ijp, itp])/PRn
    end

    return get_W
end 


# marginal consumption for household aged ij in year it
function get_Psi(ij, it)

    get_Psi = 1.0
    PRn = 1.0
    PPs = 1.0

    for ijp in ij+1:JJ
        itp = year(it, ij, ijp)
        PRn = PRn*Rn[itp]
        PPs = PPs*psi[ijp, itp]
        get_Psi = get_Psi + (beta^(ijp-ij)*PPs)^gamma*(p[itp]/PRn/p[it])^(1.0-gamma)
    end
    get_Psi = 1.0/p[it]/get_Psi

    return get_Psi
end 


# subroutine for calculating the optimal consumption path
function get_path(ij, it)

    # determine bequests at age ij
    PRn = 1.0
    beq[ij, it] = GAM[ij, it]*BQ[it]

    # determine consumption path for remainder of lifetime
    for ijp in ij+1:JJ

        # get future and previous year as well as interest factor
        itp = year(it, ij, ijp)
        itm = year(it, ij, ijp-1)
        PRn = PRn*psi[ijp, itp]*Rn[itp]

        # get consumption and bequests
        c[ijp, itp] = (beta^(ijp-ij)*PRn*p[it]/p[itp])^gamma*c[ij, it]
        beq[ijp, itp] = GAM[ijp, itp]*BQ[itp]

        # calculate assets
        a[ijp, itp] = wn[itm]*h[ijp-1] + beq[ijp-1, itm] + pen[ijp-1, itm] + Rn[itm]*a[ijp-1, itm] - p[itm]*c[ijp-1, itm]
        if(itp == 2)
            a[ijp, itp] = a[ijp, itp] + v[-ijp+3]
        end

        if(itp > 2 && ijp == 2)
            a[ijp, itp] = a[ijp, itp] + v[itm]
        end 

    end

end


# subroutine for calculating individual decisions in a certain year
function decisions(it)

    # consumption path for cohort that just entered the economy
    c[1, it] = get_Psi(1, it)*get_W(1, it)
    get_path(1, it)

    # consumption for all other cohorts in year 1 of transition
    if(it == 1)
        for ij in 2:JJ
            c[ij, it] = get_Psi(ij, it)*get_W(ij, it)
            get_path(ij, it)
        end
    end
end

# subroutine for calculating quantities in a certain year
function quantities(it)

    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    if(it == 0)
        GG[it] = gy*YY[it]
    else
        GG[it] = GG[0]
    end

    # aggregate individual decisions
    CC[it] = 0.0
    AA[it] = 0.0
    LL[it] = 0.0
    BQ[it] = 0.0
    
    for ij in 1:JJ
        CC[it] = CC[it] + c[ij, it]*m[ij, it]
        AA[it] = AA[it] + a[ij, it]*m[ij, it]/psi[ij, it]
        LL[it] = LL[it] + h[ij]*m[ij, it]
        BQ[it] = BQ[it] + (1.0-psi[ij, it])*Rn[it]*a[ij, it]*m[ij, it]/psi[ij, it]
    end

    YY[it] = KK[it]^alpha*LL[it]^(1.0-alpha)
    BB[it] = by[itm]*YY[it]

    # derive capital in small open or closed economy
    if(smopec && it > 0)

        KK[it] = LL[it]*((r[it]*(1.0-eps[it]*tauk[it])/(1.0-tauk[it])+delta)/alpha)^(1.0/(alpha-1.0))
        BF[it] = AA[it] - KK[it] - BA[it] - BB[it]
        TB[it] = (1.0+n_p(itp))*BF(itp) - (1.0+r[it])*BF[it]
    else
        KK[it] = damp*(AA[it]-BB[it]-BA[it]) + (1.0-damp)*KK[it]

    end

    II[it] = (1.0+n_p[itp])*KK[itp] - (1.0-delta)*KK[it]

end 

# subroutine for calculating government parameters
function government(it)

    taxrev = zeros(4)

    itp = year(it, 1, 2)

    taxrev[1] = tauc[it]*CC[it]
    taxrev[2] = tauw[it]*w[it]*LL[it]
    taxrev[3] = taur[it]*r[it]*AA[it]
    taxrev[4] = tauk[it]*(YY[it] - w[it]*LL[it] - (delta+eps[it]*r[it])*KK[it])

    # get budget balancing tax rate
    if(tax[it] == 1)
        tauc[it] = ((1.0+r[it])*BB[it] + GG[it] - (taxrev[2] + taxrev[3] + taxrev[4] + (1.0+n_p[itp])*BB[itp]))/CC[it]
    elseif(tax[it] == 2)
        tauw[it] = ((1.0+r[it])*BB[it] + GG[it] - (taxrev[1] + taxrev[4]+(1.0+n_p[itp])*BB[itp]))/(w[it]*LL[it] + r[it]*AA[it])
        taur[it] = tauw[it]
    elseif(tax[it] == 3)
        tauw[it] = ((1.0+r[it])*BB[it] + GG[it] - (taxrev[1] + taxrev[3] +taxrev[4] + (1.0+n_p[itp])*BB[itp]))/(w[it]*LL[it])
    else
        taur[it] = ((1.0+r[it])*BB[it] + GG[it] - (taxrev[1] + taxrev[2] +taxrev[4] + (1.0+n_p[itp])*BB[itp]))/(r[it]*AA[it])
    end

    TXR[it] = sum(taxrev)

    # get budget balancing social security contribution
    pen[JR:JJ, it] .= kappa[it]*w[it]
    Tpen[it] = 0.0

    for ij in JR:JJ
        Tpen[it] = Tpen[it] + pen[ij, it]*m[ij, it]
    end
    taup[it] = Tpen[it]/w[it]/LL[it]

end 

# solves initial steady state using Gauss-Seidel
function get_SteadyState()

    # initial guess for capital
    KK[0] = 1.0
    LL[0] = 1.0

    for iter in 1:itermax

        # get prices, decisions, quantities and taxes
        factor_prices(0)
        decisions(0)
        quantities(0)
        government(0)

        # check for the number of markets in equilibrium
        if(abs(YY[0] - CC[0] - II[0] - GG[0])/YY[0] < tol)
            break
        end
    end

    #=
    if(iter < itermax)then
        write(*,'(a,i4,a,f16.10)')'Iteration: ', iter, ' Diff: ', &
                                  abs(YY(0)-CC(0)-II(0)-GG(0))/YY(0)
        write(*,*)
    else
        write(*, '(/a/)')'!!! No equilibrium found !!!'
    endif
    =#

end 

# solves for transition path using Gauss-Seidel
function get_Transition()

    # initialize values from initial equilibrium
    a[:, 1] .= a[:, 0]
    KK[:] .= KK[0]
    LL[:] .= LL[0]

    nmarket = 0

    for iter in 1:itermax

        # get prices, decisions and quantities
        for it in 1:TT
            factor_prices(it)
        end

        for it in 1:TT
            decisions(it)
        end

        if(lsra_on) 
            lsra()
        end

        for it in 1:TT
            quantities(it)
        end

        for it in 1:TT
            government(it)
        end

        # check for the number of markets in equilibrium
        nmarket = 0
        for it in 1:TT
            if(abs(YY[it] - CC[it] - II[it] - GG[it] - TB[it])/YY[it] < tol*10.0) 
                nmarket = nmarket + 1
        end
    end

        if(nmarket == TT)
        break

    end
end

#=
    if(iter > itermax)then
        write(*, '(/a/)') '!!! No equilibrium found !!!'
    else
        write(*,'(a,i4,a,i4,a,f16.10)')'Iteration: ',iter, &
            ' Markets: ', nmarket,' Diff: ', &
            maxval(abs(YY(:) - CC(:) - II(:) - GG(:) - TB(:))/YY(:))
        write(*,*)
    endif

    =#
end 



# subroutine to compute household utility
function utility(it)

    # for first generation
    util[1, it] = 0.0
    PPs = 1.0
    for ij in 1:JJ
        itp = year(it, 1, ij)
        PPs = PPs*psi[ij, itp]
        util[1, it] = util[1, it] + beta^(ij-1)*PPs*c[ij, itp]^egam/egam
    end

    # for current total popilation if year = 0 or 1
    if(it < 2)
        for ij in 2:JJ
            util[ij, it] = c[ij, it]^egam/egam
            PPs = 1.0
            for ijp in ij+1:JJ
                itp = year(it, ij, ijp)
                PPs = PPs*psi[ijp, itp]
                util[ij, it] = util[ij, it] + beta^(ijp-ij)*PPs*c[ijp, itp]^egam/egam
            end
        end
    end

end

# subroutine for calculating lsra transfers
function lsra()

    # calculate utility for each generation
    for it in 1:TT
        utility(it)
    end

    # transfers to old generations
    BA .= 0.0
    for ij in 2:JJ
        v[-ij+2] = v[-ij+2] + get_W(ij, 1)*((util[ij, 0]/util[ij, 1])^(1.0/egam)-1.0)
        BA[2] = BA[2] + v[-ij+2]*m[ij, 1]
    end

    # long run equilibrium
    PVV  = v[TT]*(1.0+r[TT])/(r[TT]-n_p[TT])
    sum1 = get_W(1, TT)*(1.0+r[TT])/(r[TT]-n_p[TT])
    sum2 = get_W(1, TT)*(util[1, TT]*egam)^(-1.0/egam)*(1.0+r[TT])/(r[TT]-n_p[TT])

    # transition path
    for it in TT-1:-1:1
        itp = year(it, 1, 2)
        PVV = PVV*(1.0+n_p[itp])/(1.0+r[itp]) + v[it]
        sum1 = sum1*(1.0+n_p[itp])/(1.0+r[itp]) + get_W(1, it)
        sum2 = sum2*(1.0+n_p[itp])/(1.0+r[itp]) + get_W(1, it)*(util[1,it]*egam)^(-1.0/egam)
    end

    # calculate ustar for future generations
    ustar = ((sum1 - BA[2] - PVV)/sum2)^egam/egam

    # calculate transfers to future generations and debt of LSRA
    for it in 1:TT
        v[it] = v[it] + get_W(1, it)*((ustar/util[1, it])^(1.0/egam)-1.0)

        if(it == 2)
            BA[2] = (BA[2] + v[1])/(1.0+n_p[2])
        end
        if(it > 2)
            BA[it] = ((1.0+r[it-1])*BA[it-1] + v[it-1])/(1.0+n_p[it])
        end
    end

end



# initialize variables and government parameters
initialize()

# compute initial long-run equilibrium
get_SteadyState()

utility(0)

#=
# write output
output_file =  open("prog07_03.txt","w")
close(output_file)

output_file =  open("prog07_03.txt","a")

output(0, output_file)

# output(0, 20)
=#

# calculate transition path
get_Transition()



#=
for it in 1:TT
    output(it, output_file)
end
close(output_file)
=#

# calculate transition with LSRA payments
global lsra_on = true

get_Transition()

for it in 1:TT
    HEV = ((util[1,it]/util[1,0])^(1.0/egam)-1.0)*100
    println(round(HEV, digits =2))
end 