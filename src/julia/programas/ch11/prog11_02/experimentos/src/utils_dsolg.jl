using CSV 
using DataFrames 
using Query
using OffsetArrays
using Roots
using HTTP.WebSockets

#=#############################################################################
# SUBROUTINE discretize_AR
#
# Discretizes an AR(1) process of the form z_j = \rho*z_{j-1} + eps using
#     the Rouwenhorst method.
#
# REFERENCE: Kopecky, K.A., Suen, R.M.H., Finite state Markov-chain 
#            approximations to highly persistent processes, Review of Economic
#            Dynamics, Vol. 13, No. 3, 2010, 701-714.
=##############################################################################

function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real=0.0)
    σ_y = σ / sqrt(1-ρ^2)
    p  = (1+ρ)/2
    ψ = sqrt(N-1) * σ_y
    m = μ / (1 - ρ)

    state_values, p = _rouwenhorst(p, p, m, ψ, N)
    #return p, state_values

    sigma_eta = σ /(1.0-ρ^2)
    psi_val = sqrt(N-1)*sqrt(sigma_eta)

    w_est = zeros(N)
    w_est = 1.0/float(N)
    for in in 1:10000
        w_est = *(transpose(p), w_est)
    end
 
    # [-psi_val + 2.0*psi_val*float(i-1)/float(n-1) for i in 1:N]
    return p , [-psi_val + (2.0*psi_val* ((i-1)/(N-1))) for i in 1:N], w_est
end

function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
    if n == 2
        return [m-Δ, m+Δ],  [p 1-p; 1-q q]
    else
        _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n-1)
        θN = p    *[θ_nm1 zeros(n-1, 1); zeros(1, n)] +
             (1-p)*[zeros(n-1, 1) θ_nm1; zeros(1, n)] +
             q    *[zeros(1, n); zeros(n-1, 1) θ_nm1] +
             (1-q)*[zeros(1, n); θ_nm1 zeros(n-1, 1)]

        θN[2:end-1, :] ./= 2

        return range(m-Δ, stop=m+Δ, length=n), θN
    end
end


###########################################################################
# FUNCTION get_tomorrow
#
# Calculates value of function that should be integrated for pis.
###########################################################################
function get_tomorrow(pi_model)

    ###### INPUT/OUTPUT VARIABLES #########################################

    # transition probabilities
    #real*8, intent(in) :: pi(:)

    # tomorrows shock
    #integer :: get_tomorrow


    ###### OTHER VARIABLES ################################################

    #real*8 :: rand
    #integer :: i1


    ###### ROUTINE CODE ###################################################

    # get random number
    #call random_number(rand)

    random_generated = rand()
    # get tomorrows value
    for i1 in 1:size(pi_model, 1)-1

        if(random_generated <= sum(pi_model[1:i1]))
            return i1
        end
    end

    # else choose last value
    return size(pi_model, 1)-1

end 


###############################################################################
# SUBROUTINE simulate_AR
#
# Simulates a discrete AR(1) process.
###############################################################################
function simulate_AR(pi_model, shocks)

        
    ###### INPUT/OUTPUT VARIABLES #############################################
    
    # transition matrix
    #real*8, intent(in) :: pi(:, :)
    
    # simulated schocks
    #integer, intent(out) :: shocks(:)
    
    # should the random seed be initialized at a fixed values
    #logical, optional :: fixed
    
    
    ###### OTHER VARIABLES ####################################################
    
    #integer :: T, n, j
    
    
    ###### ROUTINE CODE #######################################################
    
    # assert size equality and get number of simulated schocks
    #n = assert_eq(size(pi,1), size(pi_model,2), 'tauchen')
    n = size(pi_model,1)
    T = size(shocks,1)
    
    #= initialize the random seed
    if(tbox_seed)then
        if(present(fixed))then
            call init_random_seed(fixed)
        else
            call init_random_seed()
        endif
        tbox_seed = .false.
    endif
    =#

    # get first entry
    shocks[1] = round(Int, n ÷ (2+1))  
    
    # now calculate other shocks
    for j in 2:T
        shocks[j] = get_tomorrow( pi_model[round(Int,shocks[j-1]) , :] )
    end


    return shocks
end 

###############################################################################
# SUBROUTINE grid_Cons_Grow
#
# Constructs a growing grid on [left, right].
###############################################################################
function grid_Cons_Grow(a, n, left, right, growth)
    ccall((:grid_Cons_Grow, "./fortran_wrappers/mod_julfort.so"),
                Cvoid,
                (Ptr{Float64},Ref{Int64},Ref{Float64}, Ref{Float64}, Ref{Float64}),
                a, n,  left, right,growth)
end

###############################################################################
# subroutine linint_Grow
#
# Calculates linear interpolant on a growing grid.
###############################################################################
function linint_Grow(x, left, right, growth, n, il, ir, phi)

    ccall((:linint_Grow, "./fortran_wrappers/mod_julfort.so"),
            Cvoid,
            (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64} ),
            x, left, right, growth, n, il, ir, phi)

    return il[1], ir[1], phi[1]
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

# function which computes the year in which the household lives
function year2(it, addit)


    year2 = it + addit

    if(year2 > TT)
        year2 = TT
    end 

    if(year2 < 0)
        year2 = 0
    end 

    if(it == 0)
        year2 = 0
    end
    if(it == TT)
        year2 = TT
    end

    return year2
end 

# the first order condition
function foc(x_in)
    global ij_com
    global cons_com
    global lab_com

    # calculate tomorrows assets
    a_plus  = x_in

    # get tomorrows year
    itp = year(it_com, ij_com, ij_com+1)

    # get lsra transfer payment
    v_ind = v[ij_com, ia_com, ip_com, is_com, it_com]

    # calculate the wage rate
    wage = wn[it_com]*eff[ij_com]*theta[ip_com]*eta[is_com]

    # calculate available resources
    available = (1.0+rn[it_com])*a[ia_com] + beq[ij_com, it_com] + pen[ij_com, it_com] + v_ind

    # determine labor
    if(ij_com < JR)
        global lab_com = min( max( nu + (1.0-nu)*(a_plus-available)/wage, 0.0) , 1.0-1e-10)
    else
        global lab_com = 0.0
    end

    # calculate consumption
    global cons_com = max( (available + wage*lab_com - a_plus)/p[it_com] , 1e-10)

    # calculate linear interpolation for future part of first order condition
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = max(varphi*RHS[ij_com+1, ial, ip_com, is_com, itp] + (1.0-varphi)*RHS[ij_com+1, iar, ip_com, is_com, itp], 0.0)

    # calculate first order condition for consumption
    return margu(cons_com, lab_com, it_com)^(-gamma) - tomorrow

end 

# calculates marginal utility of consumption
function margu(cons, lab, it)
    return nu*(cons^nu*(1.0-lab)^(1.0-nu))^egam/(p[it]*cons)
end 

# calculates the value function
function valuefunc(a_plus, cons, lab, ij, ip, is, it)

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0),1.0-1e-10)

    # get tomorrows year
    itp = year(it, ij, ij+1)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0
    if(ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ial, ip, is, itp] + (1.0-varphi)*EV[ij+1, iar, ip, is, itp], 1e-10)^egam/egam
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help)^(1.0-nu))^egam/egam + beta*psi[ij+1, itp]*valuefunc

end 


# computes the initial steady state of the economy
function get_SteadyState()

    # initialize remaining variables
    #call initialize()
    println("INITIAL EQUILIBRIUM")
    println("ITER     H     K/Y     C/Y     I/Y       r       w        DIFF")

    # iterate until value function converges
    for iter in 1:itermax
        # derive prices
        prices(0)

        # solve the household problem
        solve_household(1, 0)

        # calculate the distribution of households over state space
        get_distribution(0)

        # aggregate individual decisions
        aggregation(0)

        # determine the government parameters
        government(0)

        

        println(iter,"     ",round(digits = 2, HH[0]),"   ", round(digits = 2, 5.0*KK[0]/YY[0]*100.0), "   ", round(digits = 2, CC[0]/YY[0]*100.0), "   ", round(digits = 2, II[0]/YY[0]*100.0), "   ", round(digits = 2, ((1.0+r[0])^0.2-1.0)*100.0), "   ", round(digits = 2, w[0]), "   ", round(digits = 6, DIFF[0]/YY[0]*100.0))

        if(abs(DIFF[0]/YY[0])*100.0 < sig)
            break
        end

    end

end 


# subroutine for calculating prices
function prices(it)

    r[it] = Omega2*alpha*(KK[it]/LL[it])^(alpha-1.0)-delta
    w[it] = Omega2*(1.0-alpha)*(KK[it]/LL[it])^alpha
    rn[it] = r[it]*(1.0-taur[it])
    wn[it] = w[it]*(1.0-tauw[it]-taup[it])
    p[it] = 1.0 + tauc[it]

end 


# determines the solution to the household optimization problem
function solve_household(ij_in, it_in)

    global cons_com
    global lab_com
    # get decision in the last period of life
    it = year(it_in, ij_in, JJ)

    #bequest for the olds
    beq[JJ, it] = damp*GAM[JJ, it]*BQ[it] + (1.0-damp)*beq[JJ, it]

    for ia in 0:NA
        aplus[JJ, ia, :, :, it] .= 0.0
        #c[JJ, ia, :, :, it] .= ((1.0+rn[it] )*a[ia] + pen[JJ, it] .+ v[JJ, ia, :, :, it])/p[it]
        c[JJ, ia, :, :, it] .= ((1.0+rn[it])*a[ia]+ beq[JJ, it] .+ pen[JJ, it] .+ v[JJ, ia, :, :, it])/p[it]
        l[JJ, ia, :, :, it] .= 0.0
        VV[JJ, ia, :, :, it] .= valuefunc(0.0, c[JJ, ia, 1, 1, it],l[JJ, ia, 1, 1, it], JJ, 1, 1, it)
    end

    # interpolate individual RHS
    interpolate(JJ, it)
    
    for ij in JJ-1:-1:ij_in

        it = year(it_in, ij_in, ij)

        #bequest for the olds
        beq[ij, it] = damp*GAM[ij, it]*BQ[it] + (1.0-damp)*beq[ij, it]

        # check about how many is to iterate
        if(ij >= JR)
            ip_max = 1
            is_max = 1
        else
            ip_max = NP
            is_max = NS
        end

        for ia in 0:NA

            # determine decision for zero assets at retirement without pension
            if(ij >= JR && ia == 0 && kappa[it] <= 1e-10)
                aplus[ij, ia, :, :, it] .= 0.0
                c[ij, ia, :, :, it] .= 0.0
                l[ij, ia, :, :, it] .= 0.0
                VV[ij, ia, :, :, it] .= valuefunc(0.0, 0.0, 0.0, ij, 1, 1, it)
                continue
            end

            for ip in 1:ip_max
                for is in 1:is_max

                    # get initial guess for the individual choices
                    x_in = aplus[ij, ia, ip, is, it]

                    # set up communication variables
                    global ij_com = ij
                    global ia_com = ia
                    global ip_com = ip
                    global is_com = is
                    global it_com = it

                    # solve the household problem using rootfinding
                    x_root = fzero(foc, x_in)
                    foc(x_root)
                    # write screen output in case of a problem
                    #if(check)write(*,'(a, 5i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is, it

                    # check for borrowing constraint
                    if(x_root < 0.0)
                        x_root = 0.0
                        wage = wn[it]*eff[ij]*theta[ip]*eta[is]
                        v_ind = v[ij, ia, ip, is, it]
                        #available = (1.0+rn[it])*a[ia] + pen[ij, it] + v_ind
                        available = (1.0+rn[it])*a[ia] + beq[ij, it] + pen[ij, it] + v_ind
                        if(ij < JR)
                            global lab_com = min( max(nu-(1.0-nu)*available/wage , 0.0) , 1.0-1e-10)
                        else
                            global lab_com = 0.0
                        end
                        cons_com = max( (available + wage*lab_com)/p[it] , 1e-10)
                    end

                    # copy decisions
                    aplus[ij, ia, ip, is, it] = x_root
                    c[ij, ia, ip, is, it] = cons_com
                    l[ij, ia, ip, is, it] = lab_com
                    VV[ij, ia, ip, is, it] = valuefunc(x_root, cons_com, lab_com, ij, ip, is, it)

                end

                # copy decision in retirement age
                if(ij >= JR)
                    aplus[ij, ia, :, :, it] .= aplus[ij, ia, 1, 1, it]
                    c[ij, ia, :, :, it] .= c[ij, ia, 1, 1, it]
                    l[ij, ia, :, :, it] .= l[ij, ia, 1, 1, it]
                    VV[ij, ia, :, :, it] .= VV[ij, ia, 1, 1, it]
                end
            end
        end

        # interpolate individual RHS
        interpolate(ij, it)
    end

end 


# for calculating the rhs of the first order condition at age ij
function interpolate(ij, it)


    for ia in 0:NA
        for ip in 1:NP
            for is in 1:NS

                # calculate the RHS of the first order condition
                RHS[ij, ia, ip, is, it] = 0.0
                EV[ij, ia, ip, is, it] = 0.0
                for is_p in 1:NS
                    chelp = max(c[ij, ia, ip, is_p, it],1e-10)
                    lhelp = max(l[ij, ia, ip, is_p, it],1e-10)
                    RHS[ij, ia, ip, is, it] = RHS[ij, ia, ip, is, it] + pi_model[is, is_p]*margu(chelp, lhelp, it)
                    EV[ij, ia, ip, is, it]  = EV[ij, ia, ip, is, it] + pi_model[is, is_p]*VV[ij, ia, ip, is_p, it]
                end

                RHS[ij, ia, ip, is, it] = max((1.0+rn[it])*beta*psi[ij,it]*RHS[ij, ia, ip, is, it], 1e-10)^(-gamma)
                #EV[ij, ia, ip, is, it] = ((1.0-1.0/gamma)*EV[ij, ia, ip, is, it])^(1.0/(1.0-1.0/gamma))
                EV[ij, ia, ip, is, it] = (egam*EV[ij, ia, ip, is, it])^(1.0/egam)
            end
        end
    end
end 


# determines the invariant distribution of households
function get_distribution(it)

    # get yesterdays year
    itm = year(it, 2, 1)

    # set distribution to zero
    phi[:, :, :, :, it] .= 0.0

    # get initial distribution in age 1
    for ip in 1:NP
        phi[1, 0, ip, is_initial, it] = dist_theta[ip]
    end

    # successively compute distribution over ages
    for ij in 2:JJ

        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS

                    # interpolate yesterday's savings decision
                    ial, iar, varphi = linint_Grow(aplus[ij-1, ia, ip, is, itm], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                    # restrict values to grid just in case
                    #ial = max(min(ial, NA-1), 0)
                    #iar = max(min(iar, NA), 1)
                    #varphi = max(min(varphi, 1.0), 0.0)
                    ial = min(ial, NA)
                    iar = min(iar, NA)
                    varphi = min(varphi, 1.0)

                    # redistribute households
                    for is_p in 1:NS
                        phi[ij, ial, ip, is_p, it] = phi[ij, ial, ip, is_p, it] + pi_model[is, is_p]*varphi*phi[ij-1, ia, ip, is, itm]
                        phi[ij, iar, ip, is_p, it] = phi[ij, iar, ip, is_p, it] + pi_model[is,is_p]*(1.0-varphi)*phi[ij-1,ia,ip,is,itm]
                    end
                end
            end
        end
    end

end 


# subroutine for calculating quantities in a certain
function aggregation(it)

    m_coh = zeros(JJ)

    # get tomorrow's year
    itp = year(it, 1, 2)
    LL_old = LL[it]

    # calculate cohort aggregates
    c_coh[:, it]  .= 0.0
    l_coh[:, it]  .= 0.0
    y_coh[:, it]  .= 0.0
    a_coh[:, it]  .= 0.0
    VV_coh[:, it] .= 0.0
    m_coh[:]      .= 0.0
    FLC[:,it]     .= 0.0
    beq_coh[:,it] .= 0.0

    for ij in 1:JJ
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS

                    c_coh[ij, it] = c_coh[ij, it] + c[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]
                    l_coh[ij, it] = l_coh[ij, it] + l[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]
                    y_coh[ij, it] = y_coh[ij, it] + eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]
                    a_coh[ij, it] = a_coh[ij, it] + a[ia]*phi[ij, ia, ip, is, it]

                    # exclude households who die
                    if(ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                        continue
                    end

                    if(aplus[ij, ia, ip, is, it] < 1e-4)
                        FLC[ij, it] = FLC[ij, it] + phi[ij, ia, ip, is, it]
                    end 
 
                    VV_coh[ij, it] = VV_coh[ij, it] + VV[ij, ia, ip, is, it]*phi[ij, ia, ip, is, it]
                    m_coh[ij]      = m_coh[ij] + phi[ij, ia, ip, is, it]
                    beq_coh[ij, it] = beq_coh[ij, it] + a[ia]*(1.0 + rn[it])* (1.0 - psi[ij, it])*phi[ij, ia, ip, is, it]/psi[ij, it] 
                end
            end
        end
    end

    # normalize VV_coh (because hh excluded)
    VV_coh[:, it] = VV_coh[:, it]./m_coh[:]
    FLC[:, it] = FLC[:, it]./m_coh[:]

    # calculate aggregate quantities
    CC[it] = 0.0
    LL[it] = 0.0
    HH[it] = 0.0
    AA[it] = 0.0
    workpop = 0.0
    BQ[it] = 0.0
    for ij in 1:JJ
        CC[it] = CC[it] + c_coh[ij, it]*m[ij, it]
        LL[it] = LL[it] + y_coh[ij, it]*m[ij, it]
        HH[it] = HH[it] + l_coh[ij, it]*m[ij, it]
        #AA[it] = AA[it] + a_coh[ij, it]*m[ij, it]
        AA[it] = AA[it] + a_coh[ij, it]*m[ij, it]/psi[ij, it]

        BQ[it] = BQ[it] + beq_coh[ij, it]*m[ij, it]   

        if(ij < JR)
            workpop = workpop + m[ij, it]
        end
    end

    # damping and other quantities
    KK[it] = damp*(AA[it]-BB[it]-BA[it])+(1.0-damp)*KK[it]
    LL[it] = damp*LL[it] + (1.0-damp)*LL_old
    II[it] = (1.0+n_p)*KK[itp] - (1.0-delta)*KK[it]
    YY[it] = Omega2 * KK[it]^alpha * LL[it]^(1.0-alpha)

    # get average income and average working hours
    INC[it] = w[it]*LL[it]/workpop
    HH[it]  = HH[it]/workpop    # damping and other quantities


    # get difference on goods market
    DIFF[it] = YY[it]-CC[it]-II[it]-GG[it]

end 



# subroutine for calculating government parameters
function government(it)

    # last year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    # set government quantities and pension payments
    GG[it] = gy*YY[0]
    BB[it] = by*YY[0]
    pen[JR:JJ, it] .= kappa[it]*INC[itm]
    PP[it] = 0.0

    for ij in JR:JJ
        PP[it] = PP[it] + pen[ij, it]*m[ij, it]
    end

    # calculate government expenditure
    expend = GG[it] + (1.0+r[it])*BB[it] - (1.0+n_p)*BB[itp]

    # get budget balancing tax rate
    if(tax[it] == 1)
        tauc[it] = (expend - (tauw[it]*w[it]*LL[it] + taur[it]*r[it]*AA[it]))/CC[it]
        p[it]    = 1.0 + tauc[it]
    elseif(tax[it] == 2)
        tauw[it] = (expend - tauc[it]*CC[it])/(w[it]*LL[it] + r[it]*AA[it])
        taur[it] = tauw[it]
    elseif(tax[it] == 3)
        tauw[it] = (expend - (tauc[it]*CC[it] + taur[it]*r[it]*AA[it]))/(w[it]*LL[it])
    else
        taur[it] = (expend - (tauc[it]*CC[it] + tauw[it]*w[it]*LL[it]))/(r[it]*AA[it])
    end

    taxrev[1, it] = tauc[it]*CC[it]
    taxrev[2, it] = tauw[it]*w[it]*LL[it]
    taxrev[3, it] = taur[it]*r[it]*AA[it]
    taxrev[4, it] = sum(taxrev[1:3, it])

    # get budget balancing social security contribution
    taup[it] = PP[it]/(w[it]*LL[it])

end 



# computes the transition path of the economy
function get_transition()
    global lsra_comp
    global lsra_all
    global Lstar
    global lsra_on

    # initialize remaining variables
    if(!lsra_on)
        initialize_trn()
    else
        println("ITER    COMP_OLD  EFFICIENCY        DIFF")
    end

    # start timer
    #call tic()

    # iterate until value function converges
    for iter in 1:itermax

        # derive prices
        for it in 1:TT
            prices(it)
        end

        # solve the household problem
        for ij in JJ:-1:2
            solve_household(ij, 1)
        end
        for it in 1:TT
            solve_household(1, it)
        end

        # calculate the distribution of households over state space
        for it in 1:TT
            get_distribution(it)
        end

        # calculate lsra transfers if needed
        if(lsra_on)
            LSRA()
        end 

        # aggregate individual decisions
        for it in 1:TT
            aggregation(it)
        end

        # determine the government parameters
        for it in 1:TT
            government(it)
        end

        # get differences on goods markets
        itcheck = 0
        for it in 1:TT
            if(abs(DIFF[it]/YY[it])*100.0 < sig)
                itcheck = itcheck + 1
            end
        end

        #itmax = maxloc(abs(DIFF[1:TT]/YY[1:TT]), 1)
        itmax = argmax(abs.(DIFF[1:TT]./YY[1:TT]))

        # check for convergence and write screen output
        if(!lsra_on)

            check = iter > 1 && itcheck == TT && abs(DIFF[itmax]/YY[itmax])*100.0 < sig*100.0
            println(iter,"     ",round(digits = 2, HH[TT]),"   ", round(digits = 2, 5.0*KK[TT]/YY[TT]*100.0), "   ", round(digits = 2, CC[TT]/YY[TT]*100.0), "   ", round(digits = 2, II[TT]/YY[TT]*100.0), "   ", round(digits = 2, ((1.0+r[TT])^0.2-1.0)*100.0), "   ", round(digits = 2, w[TT]), "   ", round(digits = 5, DIFF[itmax]/YY[itmax]*100.0))

        else
            check = iter > 1 && itcheck == TT && lsra_comp/lsra_all > 0.99999 && abs(DIFF[itmax]/YY[itmax])*100.0 < sig*100.0
            println(iter,"     ",round(digits = 5, lsra_comp/lsra_all*100.0), "    ", round(digits = 5, (Lstar^(1.0/(1.0-1.0/gamma))-1.0)*100.0), "     ", round(digits = 5, DIFF[itmax]/YY[itmax]*100.0))
        end

        # check for convergence
        if(check)
            for it in 1:TT
                if(!lsra_on)
                    output(it)
                end
            end
            #call output_summary()
            break
        end
    end

    #call toc
    for it in 1:TT
        if(!lsra_on)
            output(it)
        end
    end
    #call output_summary()


end 


# initializes transitional variables
function initialize_trn()


    println("TRANSITION PATH")

    println("ITER       H     K/Y     C/Y     I/Y       r       w        DIFF")


    for it in 1:TT

        taup[it] = taup[0]
        tauc[it] = tauc[0]
        tauw[it] = tauw[0]
        taur[it] = taur[0]
        r[it] = r[0]
        rn[it] = r[it]*(1.0-taur[it])
        w[it] = w[0]
        wn[it] = w[it]*(1.0-tauw[it]-taup[it])
        p[it] = 1.0 + tauc[it]
        KK[it] = KK[0]
        AA[it] = AA[0]
        BB[it] = BB[0]
        LL[it] = LL[0]
        HH[it] = HH[0]
        YY[it] = YY[0]
        CC[it] = CC[0]
        II[it] = II[0]
        GG[it] = GG[0]
        INC[it] = INC[0]
        pen[:,it] = pen[:, 0]
        PP[it] = PP[0]
        taxrev[:,it] = taxrev[:, 0]
        c_coh[:, it] = c_coh[:, 0]
        l_coh[:, it] = l_coh[:, 0]
        y_coh[:, it] = y_coh[:, 0]
        a_coh[:, it] = a_coh[:, 0]
        aplus[:, :, :, :, it] = aplus[:, :, :, :, 0]
        c[:, :, :, :, it] = c[:, :, :, :, 0]
        l[:, :, :, :, it] = l[:, :, :, :, 0]
        phi[:, :, :, :, it] = phi[:, :, :, :, 0]
        VV[:, :, :, :, it] = VV[:, :, :, :, 0]
        RHS[:, :, :, :, it] = RHS[:, :, :, :, 0]

        beq[:,it] = beq[:,0]
        BQ[it] = BQ[0]
        
    end

end 

# subroutine for calculating lsra payments
function LSRA()
    global lsra_comp
    global lsra_all
    global Lstar
    global lsra_on

    # initialize variables
    SV[:] .= 0.0
    v_coh[:, :] .= 0.0

    # initialize counters
    lsra_comp     = 0.0
    lsra_all      = 0.0

    for ij in 2:JJ
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS

                    # do not do anything for an agent at retirement without pension and savings
                    if(ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                        v[ij, ia, ip, is, 1] = 0.0
                        continue
                    end

                    # get today's utility
                    VV_1 = VV[ij, ia, ip, is, 1]

                    # get target utility
                    VV_0 = VV[ij, ia, ip, is, 0]

                    # get derivative of the value function
                    dVV_dv = margu(c[ij, ia, ip, is, 1],l[ij, ia, ip, is, 1], 1)

                    # calculate change in transfers
                    v_tilde = (VV_0-VV_1)/dVV_dv

                    # restrict z_tilde to income maximum
                    v_tilde = max(v_tilde, -((1.0+rn[1])*a[ia] + pen[ij, 1] + wn[1]*eff[ij]*theta[ip]*eta[is]*0.99 + v[ij, ia, ip, is, 1]))

                    # check whether individual is already compensated
                    lsra_all = lsra_all + phi[ij, ia, ip, is, 1]*m[ij, 1]

                    if(abs((VV_1-VV_0)/VV_0)*100.0 < sig) 
                        lsra_comp = lsra_comp + phi[ij, ia, ip, is, 1]*m[ij, 1]
                    end
                    # calculate total transfer
                    v[ij, ia, ip, is, 1] = v[ij, ia, ip, is, 1] + damp*v_tilde

                    # aggregate transfers by cohort
                    v_coh[ij, 1] = v_coh[ij, 1] + v[ij, ia, ip, is, 1]*phi[ij, ia, ip, is, 1]

                end
            end
        end
    end

    # aggregate transfers in year 1
    for ij in 2:JJ
        SV[1] = SV[1] + v_coh[ij, 1]*m[ij, 1]
    end

    # initialize present value variables
    PV_t = 0.0
    PV_0 = 0.0
    PV_trans = 0.0

    # calculate present value of utility changes (in monetary values)
    for it in TT:-1:1

        # get today's ex ante utility
        EVV_t = damp*VV_coh[1, it]

        # get damped target utility
        EVV_0 = damp*VV_coh[1, 0]

        # get derivative of expected utility function
        dEVV_dv = 0.0

        for ip in 1:NP
            for is in 1:NS
                dEVV_dv = dEVV_dv + margu(c[1, 0, ip, is, it], l[1, 0, ip, is, it], it)*phi[1, 0, ip, is, it]
            end
        end

        # calculate present values
        if(it == TT)
            PV_t     = EVV_t/dEVV_dv    *(1.0+r[it])/(r[it]-n_p)
            PV_0     = EVV_0/dEVV_dv    *(1.0+r[it])/(r[it]-n_p)
            PV_trans = v[1, 0, 1, 1, it]*(1.0+r[it])/(r[it]-n_p)
        else
            PV_t     = PV_t    *(1.0+n_p)/(1.0+r[it+1]) + EVV_t/dEVV_dv
            PV_0     = PV_0    *(1.0+n_p)/(1.0+r[it+1]) + EVV_0/dEVV_dv
            PV_trans = PV_trans*(1.0+n_p)/(1.0+r[it+1]) + v[1, 0, 1, 1, it]
        end
    end

    # calculate the constant utility gain/loss for future generations
    Lstar = (PV_t-PV_trans-SV[1])/PV_0

    # calculate compensation payments for future cohorts
    for it in TT:-1:1

        # get today's ex ante utility
        EVV_t = damp*VV_coh[1, it]

        # get target utility
        EVV_0 = damp*VV_coh[1, 0]*Lstar

        # get derivative of expected utility function
        dEVV_dv = 0.0

        for ip in 1:NP
            for is in 1:NS
                dEVV_dv = dEVV_dv + margu(c[1, 0, ip, is, it], l[1, 0, ip, is, it], it)*phi[1, 0, ip, is, it]
            end
        end

        # compute change in transfers (restricted)
        v_tilde = (EVV_0-EVV_t)/dEVV_dv

        # calculate cohort transfer level
        v[1, 0, :, :, it] .= v[1, 0, :, :, it] .+ v_tilde

        # aggregate transfers
        v_coh[1, it] = v[1, 0, 1, 1, it]
        SV[it] = SV[it] + v_coh[1, it]*m[1, it]

    end

    # determine sequence of LSRA debt/savings
    BA[2] = SV[1]/(1.0+n_p)
    for it in 3:TT
        BA[it] = ((1.0+r[it-1])*BA[it-1] + SV[it-1])/(1.0+n_p)
    end

end 


# subroutine for writing output
function output(it)

    iamax = zeros(JJ)
    # check for the maximium grid point used
    check_grid(iamax, it)

end 

# subroutine that checks for the maximum gridpoint used
function check_grid(iamax, it)

    iamax .= 0
    
    for ij in 1:JJ
        # check for the maximum asset grid point used at a certain age
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS
                    if(phi[ij, ia, ip, is, it] > 1e-8)
                        iamax[ij] = ia
                    end
                end
            end
        end
    end

end 


### Build parameters
function build_parameters(country, data_file_name)
    df_param = DataFrame(CSV.File(joinpath(pwd(), "datos",data_file_name)))

    df_pais = @from i in df_param begin
        @where i.countrycode == country
        @select i
        @collect DataFrame
   end


    ### Calcula alpha --> capital share in production
    # Labour-income share --> (labsh : Share of labour compensation in GDP at current national prices)
    alpha = 1 - df_pais.labsh[1] 

    ### Calcula Omega (normaliza w = 1)
    ## w = (1 - alpha) * Omega * (K/L)**alpha
    ## Omega = (w/(1 - alpha))/((K/L)**alpha)

    K_Y = (df_pais.cn/df_pais.cgdpo)[1]
    w = 1 
    Omega = (w/(1 - alpha))/((K_Y)^alpha)


    ### Calcula tasa de depreciación
    ## I/Y = (np + delta)(K/Y) ---> delta = (I/Y)/(K/Y) - np  
    #delta = ((df_pais.csh_i[1]/(K_Y/5)) - df_pais.np[1])^(1/5)
    delta = df_pais.delta[1]

    ### Calcula nu --> is the parameter which governs the preference of the household for consumption as compared to leissure
    #=
    The autors set nu equal to 0.335, which leads labour hours to average at approximately one third of the time endowment

    This share is derived from assuming a maximum weekly working-time endowment of 110 hours as well as 50 working weeks per year. 
    We relate this average annual working hours per employee of around 1800.

    nu = 1800/(110*50)
    =#

    nu = (df_pais.avh[1]/(110*50))
    
    params_names = ["alpha", "Omega", "delta", "nu", "np", "gy" , "tauc", "tauw_min", "tauw_max", "tauw_mean", "kappa"]
    params_values = [alpha, Omega, delta, nu, df_pais.fertility_rate[1], df_pais.csh_g[1], df_pais.vat_tax[1], df_pais.tax_w_ocde_min[1], df_pais.tax_w_ocde_max[1],df_pais.tax_w_ocde_medio[1], df_pais.kappa[1]]
    
    return Dict(zip(params_names, params_values))
end 


macro Name(arg)
    string(arg)
end

function DSOLG(calib_params, stress_params)
    try
        # number of transition periods
        global TT = 140

        # number of years the household lives
        global JJ = 16

        # number of years the household retires
        global JR = 10

        # number of persistent shock process values
        global NP = 2

        # number of transitory shock process values
        global NS = 5

        # number of points on the asset grid
        global NA = 100

        # household preference parameters
        global gamma = calib_params["gamma"]
        global egam = 1.0 - 1.0/gamma
        global nu    = calib_params["nu"]
        global beta  = 0.998^5

        # household risk process
        global sigma_theta = 0.23
        global sigma_eps   = 0.05
        global rho         = 0.98

        # production parameters
        global alpha = calib_params["alpha"]
        global delta = 1.0-(1.0-calib_params["delta"])^5
        global Omega2 = calib_params["Omega"]

        # size of the asset grid
        global a_l    = 0.0
        global a_u    = 50.0
        global a_grow = 0.05

        # demographic parameters
        global n_p   = (1.0+calib_params["np"])^5-1.0

        # simulation parameters
        global damp    = 0.40
        global sig     = 1e-8
        global itermax = 50

        # counter variables
        #integer :: iter

        # macroeconomic variables
        # prices 
        for param = [:r, :rn, :w, :wn, :p]
            @eval global $param = OffsetArray(zeros(TT+1), 0:TT);
        end

        # capital market
        for param = [:KK, :AA, :BB, :LL, :HH]
            @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
        end

        # good market
        for param = [:YY, :CC, :II, :GG, :INC, :BQ]
            @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
        end

        # government variables
        for param = [:tauc, :tauw, :taur, :taup, :kappa, :PP]
            @eval global $param = OffsetArray(zeros(TT+1), 0:TT)
        end

        global gy 
        global by 

        global pen = OffsetArray(zeros(JJ,TT+1), 1:JJ,0:TT)

        global taxrev = OffsetArray(zeros(4,TT+1), 1:4,0:TT)
        global tax = OffsetArray(Int.(zeros(TT+1)), 0:TT)

        # LSRA variables
        global BA = OffsetArray(zeros(TT+1), 0:TT)
        global SV = OffsetArray(zeros(TT+1), 0:TT) 

        global lsra_comp
        global lsra_all
        global Lstar
        global lsra_on
        global Vstar

        # cohort aggregate variables
        for param = [:c_coh, :y_coh, :l_coh, :a_coh, :v_coh, :VV_coh]
            @eval global $param = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT) ;
        end

        for param = [:GAM, :beq, :beq_coh]
            @eval global $param = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT) ;
        end

        global omega = zeros(JJ)

        global psi = OffsetArray(zeros(JJ+1, TT+1),1:JJ+1, 0:TT)

        # the shock process
        global dist_theta = zeros(NP)
        global theta = zeros(NP)
        #pi(NS, NS), eta(NS)
        global is_initial = 3

        # demographic and other model parameters
        global eff = zeros(JJ)

        for param = [:m, :pop]
            @eval global $param = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT) ;
        end

        # individual variables

        global a = OffsetArray(zeros(NA+1),0:NA)

        global aplus = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT)

        for param = [:c, :l, :phi, :VV, :v]
            @eval global $param = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT)
        end 

        global FLC = OffsetArray(zeros(JJ,TT+1), 1:JJ,0:TT)

        # numerical variables
        global RHS = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT) 
        global EV = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT) 

        for param = [:ij_com, :ia_com, :ip_com, :is_com, :it_com, :cons_com, :lab_com]
            @eval global $param
        end

        global DIFF = OffsetArray(zeros(TT+1),0:TT)


        ##### initializes the remaining model parameters and variables

        # survival probabilities
        psi[1:6,0:TT] .= 1.00000000
        psi[7,0:TT] .= 0.98972953
        psi[8,0:TT] .= 0.98185396
        psi[9,0:TT] .= 0.97070373
        psi[10,0:TT] .= 0.95530594
        psi[11,0:TT] .= 0.93417914
        psi[12,0:TT] .= 0.90238714
        psi[13,0:TT] .= 0.83653436
        psi[14,0:TT] .= 0.71048182
        psi[15,0:TT] .= 0.52669353
        psi[16,0:TT] .= 0.31179803
        psi[17,0:TT] .= 0.00000000

        # set bequest distribution
        omega[1] = 1.0/6.0
        omega[2] = 1.0/6.0
        omega[3] = 1.0/6.0
        omega[4] = 1.0/6.0
        omega[5] = 1.0/6.0
        omega[6] = 1.0/6.0
        #        omega(7) = 1d0/9d0
        #        omega(8) = 1d0/9d0
        #        omega(9) = 1d0/9d0
        omega[7:16] .= 0.0

        # set up population structure
        for it in 0:TT
            m[1,it] = 1.0
            GAM[1,it] = omega[1]
            itm = year2(it, -1)
            for ij in 2:JJ
                m[ij,it] = m[ij-1, itm]*psi[ij, it]/(1.0+n_p)
                GAM[1,it] = GAM[1, it] + omega[ij]*m[ij, it]
            end
            for ij in JJ:-1:1
                GAM[ij, it] = omega[ij]/GAM[1, it]
            end
        end

        # initialize asset grid
        grid_Cons_Grow(a, NA+1, a_l, a_u, a_grow)


        # get initial guess for savings decision
        for ij in 1:JJ
            for ip in 1:NP
                for is in 1:NS
                    aplus[ij, :, ip, is, 0] .= max.(a/2, a[1]/2) 
                end
            end
        end

        # initialize age earnings process
        eff[1] = 1.0000
        eff[2] = 1.3527
        eff[3] = 1.6952
        eff[4] = 1.8279
        eff[5] = 1.9606
        eff[6] = 1.9692
        eff[7] = 1.9692
        eff[8] = 1.9392
        eff[9] = 1.9007
        eff[JR:JJ] .= 0.0

        # initialize fixed effect
        dist_theta .= 1.0/float(NP)
        theta[1]   = -sqrt(sigma_theta)
        theta[2]   = sqrt(sigma_theta)
        theta .= exp.(theta)

        # calculate the shock process
        global pi_model
        global eta

        pi_model, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
        eta = exp.(eta)

        # tax and transfers
        tax   .= 3
        tauc  .= calib_params["tauc"]
        tauw  .= 0.0
        taur  .= 0.0
        taup  .= 0.1
        kappa .= calib_params["kappa"]
        gy    = calib_params["gy"]
        by    = calib_params["by"]/5.0

        beq[:,0] .= 0.0
        BQ[0] = 0.0

        # initial guesses for macro variables
        KK .= 1.0
        LL .= 1.0
        YY .= 1.0
        II .= (n_p+delta)*KK

        GG .= gy*YY[0]
        BB .= by*YY[0]

        pen .= 0.0
        pen[JR:JJ, 0] .= kappa[0]


        global ial_v = Array{Int64}(undef, 1)
        global iar_v = Array{Int64}(undef, 1)
        global varphi_v = zeros(1)



        # calculate initial equilibrium
        get_SteadyState()

        #### ESTRESA PARÁMETROS

        # household preference parameters
        global gamma = stress_params["gamma"]
        global egam = 1.0 - 1.0/gamma
        global nu    = stress_params["nu"]

        # production parameters
        global alpha = stress_params["alpha"]
        global delta = 1.0-(1.0-stress_params["delta"])^5
        global Omega2 = stress_params["Omega"]

        # demographic parameters
        global n_p   = (1.0+stress_params["np"])^5-1.0

        # tax and transfers
        tauc  .= stress_params["tauc"]
        kappa[1:TT] .= stress_params["kappa"]
        gy    = stress_params["gy"]
        by    = stress_params["by"]/5.0

        psi[7:17,1:TT] .= psi[7:17,1:TT]*stress_params["psi_scalar"]
        psi[psi.>1.00] .= 1.00

        # Define la política 
        #=
        if stress_params["policy"] == 1
            println("BAU")
        elseif stress_params["policy"] == 2
            kappa[1:TT] .= 1.10*kappa[0];
        elseif stress_params["policy"] == 3
            kappa[1:TT] .= 1.20*kappa[0];
        end 
        =#
        global sig     = 1e-4

        # calculate transition path without lsra
        lsra_on = false;

        get_transition()

        # calculate transition path with lsra
        #lsra_on = true;
        #get_transition()

        #### GUARDAMOS LOS RESULTADOS

        ## Variables agregadas

        results_1D_df = []

        for param = [:r, :rn, :w, :wn, :p, :KK, :AA, :BB, :LL, :HH, :YY, :CC, :II, :GG, :INC, :BQ, :tauc, :tauw, :taur, :taup, :kappa, :PP]
            push!(results_1D_df,
                @eval DataFrame($param = [$param...])  
            )
        end

        results_1D_df = hcat(results_1D_df...)

        ## Variables por cohorte

        results_2D_df = []

        for param = [:pen, :c_coh, :y_coh, :l_coh, :a_coh, :v_coh, :VV_coh]
            @eval col_names = [string(@Name($param),"_cohort_",i) for i in 1:size($param.parent)[1]]
            push!(results_2D_df,
                @eval DataFrame(OffsetArrays.no_offset_view($param)', col_names)
            )
        end

        results_2D_df = hcat(results_2D_df...)

        ## Recaudación de impuestos
        df_tax_rev = DataFrame(OffsetArrays.no_offset_view(taxrev)', ["tauc_rev", "tauw_rev", "taur_rev", "total_tax_rev"])

        results_all = hcat([ DataFrame(time = 0:TT, run_id = [stress_params["iteracion"] for i in 0:TT]), results_1D_df, results_2D_df, df_tax_rev]...)

        # Agrega ganancia de eficiencia
        
        #results_all[:, :hicksian] .= (Lstar^(1.0/egam)-1.0)*100.0
        
        file_name = "DSOLG_output_"*stress_params["pais"]*"_run_id_"*string(stress_params["iteracion"])*".csv"
        
        full_save_path = joinpath(abspath(pwd(), "..", "output"), stress_params["pais"], file_name)
        
        #return results_all
        CSV.write(full_save_path, results_all)

    catch e
        println("ERROR")
    end

end


function get_run_id()
    global id
    WebSockets.open("ws://127.0.0.1:8081") do ws
           send(ws, "Hello")
           global id = parse(Int,receive(ws))
           
       end;
    return id
end