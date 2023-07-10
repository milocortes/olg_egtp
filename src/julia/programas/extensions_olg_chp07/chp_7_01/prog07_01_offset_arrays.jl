#=##############################################################################
! PROGRAM TRVL_OLG
!
! ## The OLG model with variable labor supply
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
=##############################################################################

using OffsetArrays
using Printf
using DataFrames
using StatsPlots

# model parameters
global TT      = 24
global JJ      = 3
global JR      = 3
global gamma   = 0.5
global egam    = 1.0 - 1.0/gamma
global beta    = 0.9
global nu      = 1.5
global rho     = 0.6
global erho    = 1.0 - 1.0/rho
global alpha   = 0.3
global delta   = 0.0
global tol     = 1e-5
global damp    = 0.30
global itermax = 200
global gy = 0.0
global lsra_on = false
global smopec = false


# model variables

for param = [:w, :r, :Rn, :p, :tauw, :taur, :tauc, :taup, :by, :kappa, :n_p, :tauk, :KK, :LL, :YY, :AA, :CC, :II, :BB, :GG, :BA, :BF, :TB, :Tpen, :TXR, :tax, :eps]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT) ;
end

for param = [:mu, :wn, :pen, :m, :a, :c, :util, :l]
    @eval global $param = OffsetArray(zeros(TT+1, 3), 0:TT, 1:JJ);
end

global v = OffsetArray(zeros(TT+2) , -JJ+2:TT);
#global h = zeros(TT);


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
    gy = 0.0
    by .= 0.0
    kappa .= 0.0
    tauk .= 0.0
    eps .= 0
    tax .= 1
    #lsra_on = false
    #smopec = false


    # human capital profiles for Table 7.1 (uncomment respective line for different rows)
    #global h = [1.0, 1.0, 1.0]
    #h = (/0.5d0, 1.0d0, 1.5d0/)
    #h = (/1.5d0, 1.0d0, 0.5d0/)
    #h = (/0.5d0, 2.0d0, 0.5d0/)


    # setup for Table 7.2
    #h = (/2.0d0, 2.0d0, 0.0d0/)
    #gy = 0.195d0
    #tax(1:TT) = 3
    global h = [2.0, 2.0, 0.0]
    global gy = 0.195
    tax[1:TT] .= 3

    # setup for Table 7.3
    #h = (/2.0d0, 2.0d0, 0.0d0/)
    #gy = 0.195d0
    #kappa(1:TT) = 0.5d0


    # for everything else you might wanna do
    #by(1:TT) = 0d0
    #eps(1:TT) = 0d0
    #n_p(1:TT) = 0.2d0

    # initialize tax rates shadow wages and pensions
    tauc .= 0.0
    tauw .= 0.0
    taur .= 0.0
    taup .= 0.0
    pen .= 0.0
    mu .= 0.0

    if(nu > 0.0)
        #mu(JR:JJ,:) = 0.5d0
        mu[:, JR:JJ] .= 0.5
    end 

    # initialize assets, LSRA payments and debt holdings
    a = 0.0
    c = 0.0
    l = 0.0
    v = 0.0
    YY = 0.0
    BA = 0.0
    BF = 0.0
    TB = 0.0
    TXR = 0.0

    # size of cohorts in specific year
    for it in 0:TT
        #m(1, it) = 1d0
        m[it, 1] = 1.0
        itm = year(it, 2, 1)

        for ij in 2:JJ
            #m(ij, it) = m(ij-1, itm)/(1d0 + n_p(it))
            m[it, ij] = m[itm, ij-1]/(1.0 + n_p[it])
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

    # derive after tax prices (incl. shadow wage)
    for ij in 1:JJ
        #wn(ij, it) = (h(ij)*w[it] + mu(ij, it))*(1d0-tauw[it]-taup[it])
        wn[it,ij] = (h[ij]*w[it] + mu[it,ij])*(1.0-tauw[it]-taup[it])
    end

    Rn[it] = 1.0 + r[it]*(1.0-taur[it])
    p[it] = 1.0 + tauc[it]

end 



# present value of resources for household aged ij in year it
function get_W(ij, it)

    # get current value of resources
    get_W = wn[it, ij] + pen[it, ij]

    if(it == 1 && ij > 1)
        get_W = get_W + Rn[it]*a[it, ij] + v[-ij+2]
    end

    if(it >= 1 && ij == 1)
        get_W = get_W + v[it]
    end

    # iterate over remainder of life span
    PRn = 1.0

    for ijp in ij+1:JJ
        itp = year(it, ij, ijp)
        PRn = PRn*Rn[itp]
        get_W = get_W + (wn[itp, ijp]  + pen[itp, ijp])/PRn
    end

    return get_W
end 


# marginal consumption for household aged ij in year it
function get_Psi(ij, it)

    vv = zeros(JJ)

    get_Psi = 0.0
    PRn = 1.0

    for ijp in ij:JJ
        itp = year(it, ij, ijp)
        itp1 = year(it, ij, ijp+1)

        vv[ijp] = (1.0 + nu^rho*(wn[itp, ijp ]/p[itp])^(1.0-rho))^((rho-gamma)/(1.0-rho))

        get_Psi = get_Psi + beta^((ijp-ij)*gamma)*(p[itp]/PRn/p[it])^(1.0-gamma)* vv[ijp]^((1.0-gamma)/(rho-gamma))
        
        PRn = PRn*Rn[itp1]
    end
    get_Psi = vv[ij]/p[it]/get_Psi

    return get_Psi
end 


# subroutine for calculating the optimal consumption path
function get_path(ij, it)

    vv = zeros(JJ)
    # get leisure at age ij
    PRn = 1.0
    vv[ij] = (1.0 + nu^rho*(wn[it, ij]/p[it])^(1.0-rho))^((rho-gamma)/(1.0-rho))
    
    if(nu > 0.0)
        l[it, ij] = (wn[it, ij]/nu/p[it])^(-rho)*c[it, ij]
    end

    # determine consumption and leisure path for remainder of lifetime
    for ijp in ij+1:JJ

        # get future and previous year as well as interest factor
        itp = year(it, ij, ijp)
        itm = year(it, ij, ijp-1)
        PRn = PRn*Rn[itp]

        # get consumption and leisure
        vv[ijp] = (1.0 + nu^rho*(wn[itp, ijp ]/p[itp])^(1.0-rho))^((rho-gamma)/(1.0-rho))
        c[itp, ijp ] = vv[ijp]/vv[ij]*(beta^(ijp-ij)*PRn*p[it]/p[itp])^gamma*c[it, ij]
        
        if(nu > 0.0)
            l[itp, ijp ] = (wn[itp, ijp ]/nu/p[itp])^(-rho)*c[itp, ijp ]
        end

        # calculate assets
        a[itp, ijp ] = wn[itm, ijp-1]*(1.0-l[itm, ijp-1] ) + pen[itm, ijp-1] +  Rn[itm]*a[itm, ijp-1] - p[itm]*c[itm, ijp-1]

        if(itp == 2)
            a[itp, ijp ] = a[itp, ijp ] + v[-ijp+3]
        end

        if(itp > 2 && ijp == 2)
            a[itp, ijp ] = a[itp, ijp ] + v[itm]
        end
    end

end 

# compute shadow wages
function shadw(ij, it)

    mu_new = 0
    do_all = false

    # iterate over all ages
    for ijp in ij:JJ
        itp = year(it, ij, ijp)

        # check whether leisure satisfies constraint
        if(l[itp, ijp ] > 1.0 || do_all)
            mu_new = c[itp, ijp ]^(1.0/rho)*nu*p[itp]/ (1.0-tauw[itp]-taup[itp]) - h[ijp]*w[itp]
            mu[itp, ijp ] = max((1.0-damp)*mu[itp, ijp ] + damp*mu_new, 0.0)
            # once leisure doesn't satisfy constraint, adjust all future shadow wages
            do_all = true
        else
            mu[itp, ijp ] = (1.0-damp)*mu[itp, ijp ]
        end
    end
end 

# subroutine for calculating individual decisions in a certain year
function decisions(it)

    # solve cohort that just entered the economy
    c[it, 1] = get_Psi(1, it)*get_W(1, it)
    
    get_path(1, it)
    
    if(nu > 0.0)
        shadw(1, it)
    end 

    # derive behavior for all other cohorts in year 1 of transition
    if(it == 1)
        for ij in 2:JJ
            c[it, ij] = get_Psi(ij, it)*get_W(ij, it)
            get_path(ij, it)
            if(nu > 0.0)
                shadw(ij, it)
            end
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

    for ij in 1:JJ
        CC[it] = CC[it] + c[it, ij]*m[it, ij]
        AA[it] = AA[it] + a[it, ij]*m[it, ij]
        LL[it] = LL[it] + h[ij]*(1.0-l[it, ij])*m[it, ij]
    end

    YY[it] = KK[it]^alpha*LL[it]^(1.0-alpha)
    BB[it] = by[itm]*YY[it]

    # derive capital in small open or closed economy
    if(smopec && it > 0)
        KK[it] = LL[it]*((r[it]*(1.0-eps[it]*tauk[it])/(1.0-tauk[it]) + delta)/alpha)^(1.0/(alpha-1.0))
        BF[it] = AA[it] - KK[it] - BA[it] - BB[it]
        TB[it] = (1.0+n_p[itp])*BF[itp] - (1.0+r[it])*BF[it]
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
        tauw[it] = ((1.0+r[it])*BB[it] + GG[it] - (taxrev[1] + taxrev[3] + taxrev[4] + (1.0+n_p[itp])*BB[itp]))/(w[it]*LL[it])
    else
        taur[it] = ((1.0+r[it])*BB[it] + GG[it] - (taxrev[1] + taxrev[2] + taxrev[4] + (1.0+n_p[itp])*BB[itp]))/(r[it]*AA[it])
    end

    TXR[it] = sum(taxrev)

    # get budget balancing social security contribution
    pen[it, JR:JJ] .= kappa[it]*w[it]
    
    Tpen[it] = 0.0

    for ij in JR:JJ
        Tpen[it] = Tpen[it] + pen[it, ij]*m[it, ij]
    end

    taup[it] = Tpen[it]/w[it]/LL[it]

end 


# subroutine to compute household utility
function utility(it)

    # for first generation
    util[it, 1] = 0.0
    for ij in 1:JJ
        itp = year(it, 1, ij)
        util[it, 1] = util[it, 1] + beta^(ij-1)* (c[itp,ij]^erho + nu*l[itp,ij]^erho)^(egam/erho)/egam
    end

    # for current total popilation if year = 0 or 1
    if(it < 2)
        for ij in 2:JJ
            util[it,ij] = 0.0
            for ijp in ij:JJ
                itp = year(it, ij, ijp)
                util[it,ij] = util[it,ij] + beta^(ijp-ij)* (c[itp,ijp]^erho + nu*l[itp,ijp]^erho)^(egam/erho)/egam
            end
        end
    end

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
    a[1,:] = a[0,:]
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

        #if(lsra_on)call lsra()

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

# subroutine for writing output to a file
function output(it, fname)


    # household utility
    utility(it)

    itp = year(it, 1, 2)

    # difference on goods market
    diff = YY[it] - CC[it] - GG[it] - II[it] - TB[it]


    io_buf = IOBuffer()

    @printf(io_buf, "Equilibrium: Year %i\n\n", it)
    @printf(io_buf, " Goods Market  \n")

    @printf(io_buf, "      Y      C      G      I      TB      DIFF\n")
    @printf(io_buf, "   %0.2f    %0.2f    %0.2f    %0.2f   %0.2f   %0.2f \n", YY[it], CC[it], GG[it], II[it], TB[it], diff)
    @printf(io_buf, " %0.2f   %0.2f   %0.2f    %0.2f   %0.2f\n\n", YY[it]/YY[it]*100, CC[it]/YY[it]*100, GG[it]/YY[it]*100, II[it]/YY[it]*100, TB[it]/YY[it]*100)

    @printf(io_buf, " Capital Market  \n")
    @printf(io_buf, "      A      K      BB      BA      BF      r\n")
    @printf(io_buf, "   %0.2f    %0.2f    %0.2f    %0.2f   %0.2f   %0.2f \n", AA[it], KK[it], BB[it], BA[it], BF[it], r[it])
    @printf(io_buf, "  %0.2f   %0.2f    %0.2f    %0.2f   %0.2f\n\n", AA[it]/YY[it]*100, KK[it]/YY[it]*100, BB[it]/YY[it]*100, BA[it]/YY[it]*100, BF[it]/YY[it]*100)

    #repeat("   %.2f",3)
    @printf(io_buf, " Labor Market  \n")
    @printf(io_buf,"     LL      w   util\n")
    @printf(io_buf, "   %.2f   %.2f   %.2f\n\n",LL[it], w[it], util[it, 1])

    #repeat("   %.2f",8)
    @printf(io_buf, " GOVERMENT  \n")
    @printf(io_buf,"   tauc   tauw   taur   taup   tauk    TXR     DD     rB  \n")
    @printf(io_buf, "   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f\n\n", tauc[it], tauw[it], taur[it], taup[it], tauk[it], TXR[it]/YY[it]*100.0, (1.0+n_p[itp])*BB[itp]-BB[it]/YY[it]*100.0, r[it]*BB[it]/YY[it]*100.0)

    # write individual output
    #repeat("   %.2f",8)
    @printf(io_buf," Age    cons     leis     wn       mu       pen      asset    Diff      \n")

    for ij in 1:JJ
        itp = year(it, 1, 2)

        if(ij < JJ)
            diff = Rn[it]*a[it, ij] + wn[it, ij]*(1.0-l[it, ij]) + pen[it, ij] - a[itp, ij+1] - p[it]*c[it, ij]
        else
            diff = Rn[it]*a[it, ij] + wn[it, ij]*(1.0-l[it, ij]) + pen[it, ij] - p[it]*c[it, ij]
        end

        @printf(io_buf, "  %i     %.2f     %.2f     %.2f     %.2f     %.2f     %.2f     %.2f\n", ij, c[it, ij], l[it, ij], wn[it, ij], mu[it, ij], pen[it, ij], a[it, ij], diff)
    end

    @printf(io_buf, "\n\n")

    write(fname, take!(io_buf))
    close(io_buf)

end 

# initialize variables and government parameters
initialize()

# compute initial long-run equilibrium
get_SteadyState()


# write output
output_file =  open("prog07_01.txt","w")
close(output_file)

output_file =  open("prog07_01.txt","a")

output(0, output_file)

# output(0, 20)

# calculate transition path
get_Transition()

for it in 1:TT
    output(it, output_file)
end

### save results
results_1D_df = []

for param = [:w, :r, :Rn, :p, :tauw, :taur, :tauc, :taup, :by, :kappa, :n_p, :tauk, :KK, :LL, :YY, :AA, :CC, :II, :BB, :GG, :BA, :BF, :TB, :Tpen, :TXR, :tax, :eps]
    push!(results_1D_df,
        @eval DataFrame($param = [$param...])  
    )
end

results_1D_df = hcat(results_1D_df...)

macro Name(arg)
    string(arg)
end

results_2D_df = []

for param = [:mu, :wn, :pen, :m, :a, :c, :util, :l]
    @eval col_names = [string(@Name($param),i) for i in 1:size($param.parent)[2]]
    push!(results_2D_df,
        @eval DataFrame(OffsetArrays.no_offset_view($param), col_names)
    )
end

results_2D_df = hcat(results_2D_df...)

results_all = hcat([results_1D_df, results_2D_df]...)


### Plots
## Aggregate variables effects
## Consumption effects
@df results_all plot([:YY :CC :GG :II], title =  "From consumption to labour-income taxation\n(with variable labour supply)")
plot!(legend=:outerbottom, legendcolumns=4)
xlabel!("time")


## Consumption effects
@df results_all plot([:c1 :c2 :c3], title =  "From consumption to labour-income taxation\n(with variable labour supply)")
plot!(legend=:outerbottom, legendcolumns=3)
xlabel!("time")

## Taxes effects
@df results_all plot([:tauw :tauc], title =  "From consumption to labour-income taxation\n(with variable labour supply)")
plot!(legend=:outerbottom, legendcolumns=2)
xlabel!("time")

