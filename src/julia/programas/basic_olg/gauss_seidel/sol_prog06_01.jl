#=#############################################################################
! PROGRAM OLG_GAUSS_SEIDEL
!
! ## The OLG model with Gauss-Seidel method for initial steady state
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
=##############################################################################

using OffsetArrays

# model parameters
global TT = 25
global gamma = 0.5
global egam = 1.0 - 1.0/gamma
global beta = 0.9
global alpha = 0.3
global delta = 0.0
global tol = 1e-5
global damp = 0.25
global itermax = 1000

# model variables
global lsra_on = false

for param = [:by, :kappa, :n_p, :tax, :w, :r, :wn, :Rn, :p, :tauw, :taur, :tauc, :taup, :pen, :KK, :LL, :YY, :AA, :CC, :II, :BB, :GG, :BA]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT) ;
end

for param = [:a, :c, :util]
    @eval global $param = OffsetArray(zeros(3, TT+1), 1:3, 0:TT);
end

global v = OffsetArray(zeros(TT + 2), -1:TT)

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


# initialize variables and government parameters
function initialize()

    # set baseline parameters
    global g = [0.12, 0.12, 0.0]
    by[0:TT]    .= 0.0
    kappa[0:TT] .= 0.0
    n_p[0:TT]   .= 0.2
    tax[0:TT]   .= 1
    #lsra_on     = .false.

    # set reform values (uncomment respective line for different tables)
    #tax[1:TT] .= 2                               # Table 6.2
    tax[1:TT] .= 4                               # Table 6.3
    #tax[1:TT] .= 3                               # Table 6.4
    #tax[1:TT] .= 3 ; by[1:TT] = -0.058857      # Table 6.5
    #kappa[1:TT] .= 0.5                         # Table 6.6
    #by[1:TT] .= 0.0986                       # Table 6.7
    #n_p[1:TT] .= 0.0                           # Table 6.8
    #n_p[1:TT] .= 0.0 ; kappa = 0.5 ; g = 0.0 # Table 6.9

    # get labor supply and pension payments
    LL .= (2.0.+n_p)./(1.0.+n_p)
    taup .= kappa./((2.0.+n_p).*(1.0.+n_p))

    # initialize tax rates
    tauc .= 0.0
    tauw .= 0.0
    taur .= 0.0

    # initialize assets, LSRA payments and debt holdings
    a .= 0.0
    v .= 0.0
    BA .= 0.0
end 

# subroutine for calculating factor prices in a certain year
function factor_prices(it)

    # factor prices and pension payments in year t
    r[it] = alpha*(KK[it]/LL[it])^(alpha-1.0)-delta
    w[it] = (1.0-alpha)*(KK[it]/LL[it])^alpha
    wn[it] = w[it]*(1.0-tauw[it]-taup[it])
    Rn[it] = 1.0+r[it]*(1.0-taur[it])
    p[it] = 1.0 + tauc[it]
    pen[it] = kappa[it]*w[max(it-1, 0)]

end 

# subroutine for calculating individual decisions in a certain year
function decisions(it)

    # calculate future and past years
    it1 = year(it, 1, 2)
    it2 = year(it, 1, 3)
    itm = year(it, 2, 1)

    # individual decisions
    PVI = wn[it] + wn[it1]/Rn[it1] + pen[it2]/(Rn[it1]*Rn[it2]) + v[it]
    PSI = p[it]*(1.0 + beta^gamma*(p[it1]/p[it]/Rn[it1])^(1.0-gamma) + beta^(2*gamma)*(p[it2]/p[it]/Rn[it1]/Rn[it2])^(1.0-gamma))
    c[1, it] = PVI/PSI

    if(it == 1)
        PVI = Rn[it]*a[2, 0] + wn[it] + pen[it1]/Rn[it1] + v[0]
        PSI = p[it]*(1.0 + beta^gamma*(p[it1]/p[it]/Rn[it1])^(1.0-gamma))
        c[2,it] = PVI/PSI
        c[3,it] = (pen[it] + Rn[it]*a[3, itm] + v[-1])/p[it]
        a[2,it] = wn[itm] - p[itm]*c[1, itm]
    else
        c[2,it] = (beta*Rn[it]*p[itm]/p[it])^gamma*c[1, itm]
        c[3,it] = (beta*Rn[it]*p[itm]/p[it])^gamma*c[2, itm]
        a[2,it] = wn[itm] + v[itm] - p[itm]*c[1, itm]
    end

    if(it == 2)
        a[3,it] = wn[itm] + Rn[itm]*a[2, itm] + v[0] - p[itm]*c[2, itm]
    else
        a[3,it] = wn[itm] + Rn[itm]*a[2, itm] - p[itm]*c[2, itm]
    end

end 

# subroutine for calculating quantities in a certain year
function quantities(it)

    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    # individual decisions
    YY[it] = KK[it]^alpha * LL[it]^(1.0-alpha)
    CC[it] = c[1, it] + c[2, it]/(1.0+n_p[it]) + c[3, it]/((1.0+n_p[it])*(1.0+n_p[itm]))
    GG[it] = g[1] + g[2]/(1.0+n_p[it]) + g[3]/((1.0+n_p[it])*(1.0+n_p[itm]))
    AA[it] = a[2, it]/(1.0+n_p[it]) + a[3, it]/((1.0+n_p[it])*(1.0+n_p[itm]))
    BB[it] = by[itm]*YY[it]

    # damping for Gauss-Seidel procedure
    KK[it] = damp*(AA[it] - BB[it] - BA[it])+(1.0-damp)*KK[it]
    II[it] = (1.0+n_p[itp])*KK[itp]-(1.0-delta)*KK[it]

end 

# subroutine for calculating government parameters
function government(it)

    itm = year(it, 2, 1)

    # get budget balancing tax rate
    if(tax[it] == 1)
        tauc[it] = ((1.0+r[it])*BB[it] + GG[it] - (tauw[it]*w[it]*LL[it] + taur[it]*r[it]*AA[it] + (1.0+n_p[it])*BB[min(it+1,TT)] ))/CC[it]
    elseif(tax[it] == 2)
        tauw[it] = ((1.0+r[it])*BB[it] + GG[it] - (tauc[it]*CC[it] + (1.0+n_p[it])*BB[min(it+1,TT)]))/(w[it]*LL[it] + r[it]*AA[it])
        taur[it] = tauw[it]
    elseif(tax[it] == 3)
        tauw[it] = ((1.0+r[it])*BB[it] + GG[it] - (tauc[it]*CC[it] + taur[it]*r[it]*AA[it] + (1.0+n_p[it])*BB[in(it+1,TT)] ))/(w[it]*LL[it])
    else
        taur[it] = ((1.0+r[it])*BB[it] + GG[it] - (tauc[it]*CC[it] + tauw[it]*w[it]*LL[it] + (1.0+n_p[it])*BB[min(it+1,TT)] ))/(r[it]*AA[it])
    end

    # get budget balancing social security contribution
    taup[it] = (pen[it]/((2.0+n_p[it])*(1.0+n_p[itm])))/w[it]

end 

# subroutine to compute household utility
function utility(it)

    # get future years
    if(it == 0)
        it1 = 0
        it2 = 0
    else
        it1 = min(it+1, TT)
        it2 = min(it+2, TT)
    end

    # oldest cohort
    util[3, it] = c[3, it]^egam/egam

    # middle cohort
    util[2, it] = c[2, it]^egam/egam + beta*c[3, it1]^egam/egam

    # youngest cohort
    util[1, it] = c[1, it]^egam/egam + beta*c[2, it1]^egam/egam + beta^2*c[3, it2]^egam/egam

end 

# solves for steady state using Gauss-Seidel
function get_SteadyState()

    println("INITIAL EQUILIBRIUM")

    # initial guess for capital stock
    KK[0] = 1.0

    for iter in 1:itermax

        # get prices, decisions and quantities
        factor_prices(0)
        decisions(0)
        quantities(0)
        government(0)

        # check whether goods market is in equilibrium
        if(abs(YY[0] - CC[0] - II[0] - GG[0])/YY[0] < tol) 
                break
        end
    end

end 



# subroutine for calculating lsra transfers
function lsra()

    # calculate utility for each generation
    for it in 1:TT
        utility(it)
    end

    global PVI = OffsetArray(zeros(TT + 2), -1:TT)

    # transfers to old generations
    PVI[-1] = Rn[1]*a[3,1] + pen[1] + v[-1]
    PVI[0]  = Rn[1]*a[2,1] + wn[1] + pen[2]/Rn[2] + v[0]
    v[-1] = v[-1] + PVI[-1]*((util[3,0]/util[3,1])^(1.0/egam)-1.0)
    v[0]  = v[0] + PVI[0] *((util[2,0]/util[2,1])^(1.0/egam)-1.0)
    BA[2] = v[-1]/((1.0+n_p[0])*(1.0+n_p[1])) + v[0]/(1.0+n_p[1])

    # long run equilibrium
    PVI[TT] = wn[TT] + wn[TT]/Rn[TT] + pen[TT]/Rn[TT]^2 + v[TT]
    PVV  = v[TT]*(1.0+r[TT])/(r[TT]-n_p[TT])
    sum1 = PVI[TT]*(1.0+r[TT])/(r[TT]-n_p[TT])
    sum2 = PVI[TT]*(util[1,TT]*egam)^(-1.0/egam)*(1.0+r[TT])/ (r[TT]-n_p[TT])

    # transition path
    for it in TT-1:-1:1 
        it1 = year(it, 1, 2)
        it2 = year(it, 1, 3)
        PVI[it] = wn[it] + wn[it1]/Rn[it1] + pen[it2]/(Rn[it1]*Rn[it2]) + v[it]
        PVV = PVV*(1.0+n_p[it1])/(1.0+r[it1]) + v[it]
        sum1 = sum1*(1.0+n_p[it1])/(1.0+r[it1]) + PVI[it]
        sum2 = sum2*(1.0+n_p[it1])/(1.0+r[it1]) + PVI[it]* (util[1,it]*egam)^(-1.0/egam)
    end

    # calculate ustar for future generations
    global  ustar = ((sum1-BA[2]-PVV)/sum2)^egam/egam

    # calculate transfers to future generations and debt of LSRA
    for it in 1:TT
        v[it] = v[it] + PVI[it]*((ustar/util[1, it])^(1.0/egam)-1.0)
        if(it == 2)
            BA[2] = (BA[2] + v[1])/(1.0+n_p[2])
        end
        if(it > 2)
            BA[it] = ((1.0+r[it-1])*BA[it-1] + v[it-1])/(1.0+n_p[it])
        end
    end

end 

# solves for transition path using Gauss-Seidel
function get_Transition()

    if(lsra_on)
        println("TRANSITION PATH WITH LSRA")
    else
        println("TRANSITION PATH WITHOUT LSRA")
    end

    # initialize values from initial equilibrium
    KK[:] .= KK[0]
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
            if(abs(YY[it] - CC[it] - II[it] - GG[it] )/YY[it] < tol*10.0) 
                nmarket = nmarket + 1
            end
        end

        if(nmarket == TT)
            break
        end
    end

end 


# initialize variables and government parameters
initialize()

# compute initial long-run equilibrium
get_SteadyState()

utility(0)

# calculate transition path
get_Transition()


# calculate transition with LSRA payments
global lsra_on = true

get_Transition()

for it in 1:TT
    HEV = ((util[1,it]/util[1,0])^(1.0/egam)-1.0)*100
    println(round(HEV, digits =2))
end 