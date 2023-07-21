#=##############################################################################
! PROGRAM TRHCEG_OLG
!
! ## The OLG model with endogenous growth
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
global xi      = 2.0
global upsi    = 0.5
global epsi    = 0.0
global mu      = 0.8
global tol     = 1e-6
global damp    = 0.30
global itermax = 200
global gy = 0.0
global smopec = false

# model variables

for param = [:w, :r, :wn, :Rn, :p, :tauw, :taur, :tauc, :taup, :by, :kappa, :n_p, :tauk, :KK, :LL, :YY, :AA, :CC, :II, :BB, :GG, :BA, :BF, :TB, 
            :Tpen, :TXR, :taus, :e, :HH, :n_e, :tax, :eps]
    @eval global $param = OffsetArray(zeros(TT+1), 0:TT) ;
end

for param = [:winc, :pen, :m, :ma, :a, :c]
    @eval global $param = OffsetArray(zeros(TT+1, 3), 0:TT, 1:JJ);
end


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
    n_p .= 0.0
    by .= 0.0
    kappa .= 0.0
    eps .= 0
    tauk .= 0.0
    taus .= 0.0
    tax .= 1
    smopec = false


    # impact of mu Table 7.8
    # change mu in the parameter section at beginning of program
    global gy = 0.00


    # setup for Table 7.9
    # set mu = 0.8d0 at beginnig of program section
    #global gy = 0.22
    #tax[1:TT] .= 3


    # setup for Table 7.9
    # set mu = 0.8d0 at beginnig of program section
    #gy = 0.22d0
    #taus(1:TT) = 0.25d0


    # initialize tax rates and pensions
    tauc = 0.0
    tauw = 0.0
    taur = 0.0
    taup = 0.0
    pen = 0.0

    # initialize assets, LSRA payments and debt holdings
    a .= 0.0
    e .= 0.0
    winc .= 0.0
    n_e .= 0.0
    YY .= 0.0
    HH .= 1.0
    BF .= 0.0
    TB .= 0.0
    TXR .= 0.0
    ma .= 0.0

    # size of cohorts in specific year
    for it in 0:TT
        m[it,1] = 1.0
        itm = year(it, 2, 1)
        for ij in 2:JJ
            m[it,ij] = m[itm,ij-1]/(1.0 + n_p[it])
        end
    end

end 



# subroutine for calculating factor prices in a certain year
function factor_prices(it)

    # factor prices and pension payments in year t
    if(smopec && it > 0)
        r[it] = r[0]
    else
        r[it] = (1.0-tauk[it])/(1.0-eps[it]*tauk[it])*(alpha*(KK[it]/LL[it])^(alpha-1.0)*HH[it]^epsi-delta)
    end

    w[it] = (1.0-alpha)*(KK[it]/LL[it])^alpha*HH[it]^epsi

    # derive after tax prices
    wn[it] = w[it]*(1.0-tauw[it]-taup[it])
    winc[it, 1] = (1.0 - (1.0-taus[it])*e[it])*wn[it]
    
    for ij in 2:JR-1
        itp = year(it, ij, 2)
        winc[it, ij] = wn[it]*(1.0+n_e[itp])/mu
    end

    Rn[it] = 1.0 + r[it]*(1.0-taur[it])
    p[it] = 1.0 + tauc[it]

end 



# present value of resources for household aged ij in year it
function get_W(ij, it)

    # get current value of resources
    get_W = winc[it, ij] + pen[it, ij]

    if(it == 1 && ij > 1)
        get_W = get_W + Rn[it]*a[it, ij]  
    end

    # iterate over remainder of life span
    PRn = 1.0

    for ijp in ij+1:JJ
        itp = year(it, ij, ijp)
        PRn = PRn*Rn[itp]
        get_W = get_W + (winc[itp, ijp]  + pen[itp, ijp])/PRn
    end

    return get_W
end 


# marginal consumption for household aged ij in year it
function get_Psi(ij, it)

    get_Psi = 1.0
    PRn = 1.0

    for ijp in ij+1:JJ
        itp = year(it, ij, ijp)
        PRn = PRn*Rn[itp]

        get_Psi = get_Psi + beta^((ijp-ij)*gamma)*(p[itp]/PRn/p[it])^(1.0-gamma)

    end
    get_Psi = 1.0/p[it]/get_Psi

    return get_Psi
end 


# subroutine for calculating the optimal consumption path
function get_path(ij, it)

    PRn = 1.0

    # determine consumption and leisure path for remainder of lifetime
    for ijp in ij+1:JJ

        # get future and previous year as well as interest factor
        itp = year(it, ij, ijp)
        itm = year(it, ij, ijp-1)
        PRn = PRn*Rn[itp]

        # get consumption and assets
        #c(ijp, itp) = (beta**(ijp-ij)*PRn*p(it)/p(itp))**gamma*c(ij, it)
        c[itp, ijp] = (beta^(ijp-ij)*PRn*p[it]/p[itp])^gamma*c[it, ij]
        
        #a(ijp, itp) = winc(ijp-1, itm) + pen(ijp-1, itm) + Rn(itm)*a(ijp-1, itm) - p(itm)*c(ijp-1, itm)
        a[itp, ijp] = winc[itm, ijp-1] + pen[itm, ijp-1] + Rn[itm]*a[itm, ijp-1] - p[itm]*c[itm, ijp-1]
    end

end 


# subroutine for calculating individual decisions in a certain year
function decisions(it)

    # derive education decision and corresponding growth rate
    itp = year(it, 1, 2)
    #e(it) = (xi*upsi*wn(itp)/(wn(it)*(1d0-taus(it))*Rn(itp)))**(1d0/(1d0-upsi))
    e[it] = (xi*upsi*wn[itp]/(wn[it]*(1.0-taus[it])*Rn[itp]))^(1.0/(1.0-upsi))

    #n_e(itp) = mu*(1d0+xi*e(it)**upsi) - 1d0
    n_e[itp] = mu*(1.0+xi*e[it]^upsi) - 1.0

    # consumption path for cohort that just entered the economy
    
    #c(1, it) = get_Psi(1, it)*get_W(1, it)
    c[it, 1] = get_Psi(1, it)*get_W(1, it)
    get_path(1, it)

    # consumption for all other cohorts in year 1 of transition
    if(it == 1)
        for ij in 2:JJ
            #c(ij, it) = get_Psi(ij, it)*get_W(ij, it)
            c[it, ij] = get_Psi(ij, it)*get_W(ij, it)
            get_path(ij, it)
        end
    end

end 

# subroutine for calculating quantities in a certain year
function quantities(it)

    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    # government consumption
    if(it == 0)
        GG[it] = gy*YY[it]
    else
        GG[it] = GG[0]
    end

    # calculate adjusted aggregator
    #ma(1, it) = 1.0
    ma[it,1] = 1.0

    for ij in 2:JJ
        #ma(ij, it) = ma(ij-1, itm)/((1.0+n_e[it])*(1.0+n_p[it]))
        ma[it,ij] = ma[itm,ij-1] /((1.0+n_e[it])*(1.0+n_p[it]))
    end

    # calculate economy wide aggregates
    CC[it] = 0.0
    AA[it] = 0.0
    LL[it] = 0.0
    HH[it] = 0.0
    
    for ij in 1:JJ
        #CC[it] = CC[it] + c(ij, it)*ma(ij, it)
        CC[it] = CC[it] + c[it, ij]*ma[it,ij]

        #AA[it] = AA[it] + a(ij, it)*ma(ij, it)
        AA[it] = AA[it] + a[it,ij]*ma[it,ij]

        if(ij < JR)
            #LL[it] = LL[it] + mu^dble(1-ij)*m(ij, it)
            LL[it] = LL[it] + mu^Float64(1-ij)*m[it,ij]
        end 

        if(ij < JR)
            #HH[it] = HH[it] + m(ij, it)
            HH[it] = HH[it] + m[it,ij]
        end
    end

    LL[it] = LL[it] - e[it]
    HH[it] = LL[it]/(HH[it] - e[it])

    # derive output and government debt
    YY[it] = KK[it]^alpha*LL[it]^(1.0-alpha)*HH[it]^epsi
    BB[it] = by[itm]*YY[it]

    # capital stock (in closed or open economy)
    if(smopec && it > 0)
        KK[it] = LL[it]*((r[it]*(1.0-eps[it]*tauk[it])/(1.0-tauk[it])+ delta)/alpha)^(1.0/(alpha-1.0))
        BF[it] = AA[it] - KK[it] - BB[it]
        TB[it] = (1.0+n_p[itp])*(1.0+n_e[itp])*BF[itp] - (1.0+r[it])*BF[it]
    else
        KK[it] = damp*(AA[it]-BB[it]) + (1.0-damp)*KK[it]
    end
    II[it] = (1.0+n_p[itp])*(1.0+n_e[itp])*KK[itp]-(1.0-delta)*KK[it]


end 

# subroutine for calculating government parameters
function government(it)

    itp = year(it, 1, 2)

    taxrev = zeros(5)

    taxrev[1] = tauc[it]*CC[it]
    taxrev[2] = tauw[it]*w[it]*LL[it]
    taxrev[3] = taur[it]*r[it]*AA[it]
    taxrev[4] = tauk[it]*(YY[it] - w[it]*LL[it] - (delta+eps[it]*r[it])*KK[it])
    taxrev[5] = taus[it]*e[it]*wn[it]

    # total growth rate
    global n_g = (1.0+n_p[itp])*(1.0+n_e[itp]) - 1.0

    # get budget balancing tax rate
    if(tax[it] == 1)
        tauc[it] = (taxrev[5] + (1.0+r[it])*BB[it] + GG[it] - (taxrev[2] + taxrev[3] + taxrev[4] + (1.0+n_g)*BB[itp]))/CC[it]
    elseif(tax[it] == 2)
        tauw[it] = (taxrev[5] + (1.0+r[it])*BB[it] + GG[it] - (taxrev[1] + taxrev[4] + (1.0+n_g)*BB[itp]))/(w[it]*LL[it] + r[it]*AA[it])
        taur[it] = tauw[it]
    elseif(tax[it] == 3)
        tauw[it] = (taxrev[5] + (1.0+r[it])*BB[it] + GG[it] - (taxrev[1] + taxrev[3] + taxrev[4] + (1.0+n_g)*BB[itp]))/(w[it]*LL[it])
    else
        taur[it] = (taxrev[5] + (1.0+r[it])*BB[it] + GG[it] - (taxrev[1] + taxrev[2] + taxrev[4] + (1.0+n_g)*BB[itp]))/(r[it]*AA[it])
    end

    TXR[it] = sum(taxrev[1:4]) - taxrev[5]

    # get budget balancing social security contribution
    #pen(JR:JJ, it) = kappa[it]*w[it]
    pen[it, JR:JJ] .= kappa[it]*w[it]

    Tpen[it] = 0.0
    for ij in JR:JJ
        #Tpen[it] = Tpen[it] + pen(ij, it)*ma(ij, it)
        Tpen[it] = Tpen[it] + pen[it,ij]*ma[it,ij]
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


# subroutine for writing output to a file
function output(it, fname)


    itp = year(it, 1, 2)

    # total growth rate
    n_g = (1.0+n_p[itp])*(1.0+n_e[itp]) - 1.0

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
    @printf(io_buf,"     LL      e   n_e   w\n")
    @printf(io_buf, "   %.2f   %.2f   %.2f   %.2f\n\n",LL[it], e[it], n_e[it], w[it])

    #repeat("   %.2f",8)
    @printf(io_buf, " GOVERMENT  \n")
    @printf(io_buf,"   tauc   tauw   taur   taup   tauk    TXR     DD     rB  \n")
    @printf(io_buf, "   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f\n\n", tauc[it], tauw[it], taur[it], taup[it], tauk[it], TXR[it]/YY[it]*100.0, (1.0+n_g)*BB[itp]-BB[it]/YY[it]*100.0, r[it]*BB[it]/YY[it]*100.0)

    # write individual output
    #repeat("   %.2f",8)
    @printf(io_buf," Age    cons     wn       pen      asset    Diff      \n")

    for ij in 1:JJ
        itp = year(it, 1, 2)

        if(ij < JJ)
            diff = Rn[it]*a[it,ij] + wn[it]*((1.0+n_e[it])/mu)^(ij-1) + pen[it,ij] - a[itp,ij+1] - p[it]*c[it,ij]

            if(ij == 1)
                diff = diff - (1.0-taus[it])*e[it]*wn[it]
            end

            @printf(io_buf, "  %i     %.2f     %.2f     %.2f     %.2f     %.2f\n", ij, c[it, ij],  wn[it]*((1.0+n_e[it])/mu)^(ij-1), pen[it, ij], a[it, ij], diff)
        else
            diff = Rn[it]*a[it,ij] + pen[it,ij] - p[it]*c[it,ij]
            @printf(io_buf, "  %i     %.2f     %.2f     %.2f     %.2f     %.2f\n", ij, c[it, ij],  0.0, pen[it, ij], a[it, ij], diff)

        end

    end

    @printf(io_buf, "\n\n")

    write(fname, take!(io_buf))
    close(io_buf)

end 

# solves for transition path using Gauss-Seidel
function get_Transition()


    # initialize values from initial equilibrium
    a[1,:] = a[0,:]
    n_e[1] = n_e[0]
    KK[:] .= KK[0]
    LL[:] .= LL[0]
    HH[:] .= HH[0]

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



# initialize variables and government parameters
initialize()

# compute initial long-run equilibrium
get_SteadyState()


# write output
output_file =  open("prog07_03.txt","w")
close(output_file)

output_file =  open("prog07_03.txt","a")

output(0, output_file)

# output(0, 20)

# calculate transition path
get_Transition()

for it in 1:TT
    output(it, output_file)
end
close(output_file)

### save results
using DataFrames
using StatsPlots

results_1D_df = []

for param = [:w, :r, :wn, :Rn, :p, :tauw, :taur, :tauc, :taup, :by, :kappa, :n_p, :tauk, :KK, :LL, :YY, :AA, :CC, :II, :BB, :GG, :BA, :BF, :TB, 
            :Tpen, :TXR, :taus, :e, :HH, :n_e, :tax, :eps]
    push!(results_1D_df,
        @eval DataFrame($param = [$param...])  
    )
end

results_1D_df = hcat(results_1D_df...)

macro Name(arg)
    string(arg)
end

results_2D_df = []

for param = [:winc, :pen, :m, :ma, :a, :c]
    @eval col_names = [string(@Name($param),i) for i in 1:size($param.parent)[2]]
    push!(results_2D_df,
        @eval DataFrame(OffsetArrays.no_offset_view($param), col_names)
    )
end

results_2D_df = hcat(results_2D_df...)

results_all = hcat([results_1D_df, results_2D_df]...)