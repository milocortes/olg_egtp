using CSV 
using DataFrames 
using Query

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
function get_tomorrow(pi)

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
    for i1 in 1:size(pi, 1)-1

        if(random_generated <= sum(pi[1:i1]))
            return i1
        end
    end

    # else choose last value
    return size(pi, 1)-1

end 


###############################################################################
# SUBROUTINE simulate_AR
#
# Simulates a discrete AR(1) process.
###############################################################################
function simulate_AR(pi, shocks)

        
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
    #n = assert_eq(size(pi,1), size(pi,2), 'tauchen')
    n = size(pi,1)
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
        shocks[j] = get_tomorrow( pi[round(Int,shocks[j-1]) , :] )
    end


    return shocks
end 

###############################################################################
# SUBROUTINE grid_Cons_Grow
#
# Constructs a growing grid on [left, right].
###############################################################################
function grid_Cons_Grow(a, n, left, right, growth)
    ccall((:grid_Cons_Grow, "fortran_wrappers/mod_julfort.so"),
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

    ccall((:linint_Grow, "fortran_wrappers/mod_julfort.so"),
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
    global ik_com 

    # calculate tomorrows assets
    a_plus  = x_in

    # get tomorrows year
    itp = year(it_com, ij_com, ij_com+1)

    # get lsra transfer payment
    v_ind = v[ij_com, ia_com, ip_com, is_com, it_com, ik_com]

    # calculate the wage rate
    wage = wn[it_com]*eff[ij_com, ik_com]*theta[ip_com]*eta[is_com]

    # calculate available resources
    available = (1.0+rn[it_com])*a[ia_com] + beq[ij_com, it_com, ik_com] + pen[ij_com, it_com, ik_com] + v_ind

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

    tomorrow = max(varphi*RHS[ij_com+1, ial, ip_com, is_com, itp, ik_com] + (1.0-varphi)*RHS[ij_com+1, iar, ip_com, is_com, itp, ik_com], 0.0)

    # calculate first order condition for consumption
    return margu(cons_com, lab_com, it_com)^(-gamma) - tomorrow

end 

# calculates marginal utility of consumption
function margu(cons, lab, it)
    return nu*(cons^nu*(1.0-lab)^(1.0-nu))^egam/(p[it]*cons)
end 

# calculates the value function
function valuefunc(a_plus, cons, lab, ij, ip, is, it, ik)

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
        valuefunc = max(varphi*EV[ij+1, ial, ip, is, itp, ik] + (1.0-varphi)*EV[ij+1, iar, ip, is, itp, ik], 1e-10)^egam/egam
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
        for ik in 1:SS
            solve_household(1, 0, ik)
        end 

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
function solve_household(ij_in, it_in, ik_in)

    global cons_com
    global lab_com
    # get decision in the last period of life
    it = year(it_in, ij_in, JJ)

    #bequest for the olds
    beq[JJ, it, ik_in] = damp*GAM[JJ, it, ik_in]*BQ[it] + (1.0-damp)*beq[JJ, it, ik_in]

    for ia in 0:NA
        aplus[JJ, ia, :, :, it, ik_in] .= 0.0
        #c[JJ, ia, :, :, it] .= ((1.0+rn[it] )*a[ia] + pen[JJ, it] .+ v[JJ, ia, :, :, it])/p[it]
        c[JJ, ia, :, :, it, ik_in] .= ((1.0+rn[it])*a[ia]+ beq[JJ, it, ik_in] .+ pen[JJ, it, ik_in] .+ v[JJ, ia, :, :, it, ik_in])/p[it]
        l[JJ, ia, :, :, it, ik_in] .= 0.0
        VV[JJ, ia, :, :, it, ik_in] .= valuefunc(0.0, c[JJ, ia, 1, 1, it, ik_in],l[JJ, ia, 1, 1, it, ik_in], JJ, 1, 1, it, ik_in)
    end

    # interpolate individual RHS
    interpolate(JJ, it, ik_in)

    for ij in JJ-1:-1:ij_in

        it = year(it_in, ij_in, ij)

        #bequest for the olds
        beq[ij, it, ik_in] = damp*GAM[ij, it, ik_in]*BQ[it] + (1.0-damp)*beq[ij, it, ik_in]

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
                aplus[ij, ia, :, :, it, ik_in] .= 0.0
                c[ij, ia, :, :, it, ik_in] .= 0.0
                l[ij, ia, :, :, it, ik_in] .= 0.0
                VV[ij, ia, :, :, it, ik_in] .= valuefunc(0.0, 0.0, 0.0, ij, 1, 1, it, ik_in)
                continue
            end

            for ip in 1:ip_max
                for is in 1:is_max

                    # get initial guess for the individual choices
                    x_in = aplus[ij, ia, ip, is, it, ik_in]

                    # set up communication variables
                    global ij_com = ij
                    global ia_com = ia
                    global ip_com = ip
                    global is_com = is
                    global it_com = it
                    global ik_com = ik_in

                    # solve the household problem using rootfinding
                    x_root = fzero(foc, x_in, xatol=1e-12)
                    foc(x_root)
                    # write screen output in case of a problem
                    #if(check)write(*,'(a, 5i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is, it

                    # check for borrowing constraint
                    if(x_root < 0.0)
                        x_root = 0.0
                        wage = wn[it]*eff[ij, ik_in]*theta[ip]*eta[is]
                        v_ind = v[ij, ia, ip, is, it, ik_in]
                        #available = (1.0+rn[it])*a[ia] + pen[ij, it] + v_ind
                        available = (1.0+rn[it])*a[ia] + beq[ij, it, ik_in] + pen[ij, it, ik_in] + v_ind
                        if(ij < JR)
                            global lab_com = min( max(nu-(1.0-nu)*available/wage , 0.0) , 1.0-1e-10)
                        else
                            global lab_com = 0.0
                        end
                        cons_com = max( (available + wage*lab_com)/p[it] , 1e-10)
                    end

                    # copy decisions
                    aplus[ij, ia, ip, is, it, ik_in] = x_root
                    c[ij, ia, ip, is, it, ik_in] = cons_com
                    l[ij, ia, ip, is, it, ik_in] = lab_com
                    VV[ij, ia, ip, is, it, ik_in] = valuefunc(x_root, cons_com, lab_com, ij, ip, is, it, ik_in)

                end

                # copy decision in retirement age
                if(ij >= JR)
                    aplus[ij, ia, :, :, it, ik_in] .= aplus[ij, ia, 1, 1, it, ik_in]
                    c[ij, ia, :, :, it, ik_in] .= c[ij, ia, 1, 1, it, ik_in]
                    l[ij, ia, :, :, it, ik_in] .= l[ij, ia, 1, 1, it, ik_in]
                    VV[ij, ia, :, :, it, ik_in] .= VV[ij, ia, 1, 1, it, ik_in]
                end
            end
        end

        # interpolate individual RHS
        interpolate(ij, it, ik_in)
    end

end 


# for calculating the rhs of the first order condition at age ij
function interpolate(ij, it, ik)
    for ia in 0:NA
        for ip in 1:NP
            for is in 1:NS
                # calculate the RHS of the first order condition
                RHS[ij, ia, ip, is, it, ik] = 0.0
                EV[ij, ia, ip, is, it, ik] = 0.0
                for is_p in 1:NS
                    chelp = max(c[ij, ia, ip, is_p, it, ik],1e-10)
                    lhelp = max(l[ij, ia, ip, is_p, it, ik],1e-10)
                    RHS[ij, ia, ip, is, it, ik] = RHS[ij, ia, ip, is, it, ik] + pi[is, is_p]*margu(chelp, lhelp, it)
                    EV[ij, ia, ip, is, it, ik]  = EV[ij, ia, ip, is, it, ik] + pi[is, is_p]*VV[ij, ia, ip, is_p, it, ik]
                end
                RHS[ij, ia, ip, is, it, ik] = max((1.0+rn[it])*beta*psi[ij,it]*RHS[ij, ia, ip, is, it, ik], 1e-10)^(-gamma)
                #EV[ij, ia, ip, is, it] = ((1.0-1.0/gamma)*EV[ij, ia, ip, is, it])^(1.0/(1.0-1.0/gamma))
                EV[ij, ia, ip, is, it, ik] = (egam*EV[ij, ia, ip, is, it, ik])^(1.0/egam)
            end
        end
    end
end 


# determines the invariant distribution of households
function get_distribution(it)

    # get yesterdays year
    itm = year(it, 2, 1)

    # set distribution to zero
    phi[:, :, :, :, it, :] .= 0.0

    # get initial distribution in age 1
    for ip in 1:NP
        phi[1, 0, ip, is_initial, it, :] .= dist_theta[ip]
    end

    # successively compute distribution over ages
    for ij in 2:JJ

        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS

                    for ik in 1:SS
                        # interpolate yesterday's savings decision
                        ial, iar, varphi = linint_Grow(aplus[ij-1, ia, ip, is, itm, ik], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

                        # restrict values to grid just in case
                        #ial = max(min(ial, NA-1), 0)
                        #iar = max(min(iar, NA), 1)
                        #varphi = max(min(varphi, 1.0), 0.0)
                        ial = min(ial, NA)
                        iar = min(iar, NA)
                        varphi = min(varphi, 1.0)

                        # redistribute households
                        for is_p in 1:NS
                            phi[ij, ial, ip, is_p, it, ik] = phi[ij, ial, ip, is_p, it, ik] + pi[is, is_p]*varphi*phi[ij-1, ia, ip, is, itm, ik]
                            phi[ij, iar, ip, is_p, it, ik] = phi[ij, iar, ip, is_p, it, ik] + pi[is,is_p]*(1.0-varphi)*phi[ij-1,ia,ip,is,itm, ik]
                        end
                    end
                end
            end
        end
    end

end 


# subroutine for calculating quantities in a certain
function aggregation(it)

    # calculate fraction of low and high households skills
    for ij in 1:JJ
        for ik in 1:SS
            frac_phi[ij, it, ik] = sum(phi[ij, :, :, :, it, ik])
        end
    end

    global SS
    m_coh = zeros(JJ,SS)

    # get tomorrow's year
    itp = year(it, 1, 2)
    LL_old = LL[it]

    # calculate cohort aggregates
    for ik in 1:SS
        c_coh[:, it, ik]  .= 0.0
        l_coh[:, it, ik]  .= 0.0
        y_coh[:, it, ik]  .= 0.0
        a_coh[:, it, ik]  .= 0.0
        VV_coh[:, it, ik] .= 0.0
        m_coh[:, ik]      .= 0.0
        FLC[:,it, ik]     .= 0.0
        beq_coh[:,it, ik] .= 0.0
    end

    for ij in 1:JJ
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS
                    for ik in 1:SS
                        c_coh[ij, it, ik] = c_coh[ij, it, ik] + c[ij, ia, ip, is, it, ik]*phi[ij, ia, ip, is, it, ik]/frac_phi[ij, it, ik]
                        l_coh[ij, it, ik] = l_coh[ij, it, ik] + l[ij, ia, ip, is, it, ik]*phi[ij, ia, ip, is, it, ik]/frac_phi[ij, it, ik]
                        y_coh[ij, it, ik] = y_coh[ij, it, ik] + eff[ij, ik]*theta[ip]*eta[is]*l[ij, ia, ip, is, it, ik]*phi[ij, ia, ip, is, it, ik]/frac_phi[ij, it, ik]
                        a_coh[ij, it, ik] = a_coh[ij, it, ik] + a[ia]*phi[ij, ia, ip, is, it, ik]/frac_phi[ij, it, ik]

                        # exclude households who die
                        if(ij >= JR && ia == 0 && (kappa[0] <= 1e-10 || kappa[1] <= 1e-10))
                            continue
                        end

                        if(aplus[ij, ia, ip, is, it, ik] < 1e-4)
                            FLC[ij, it, ik] = FLC[ij, it, ik] + phi[ij, ia, ip, is, it, ik]/frac_phi[ij, it, ik]
                        end 
    
                        VV_coh[ij, it, ik] = VV_coh[ij, it, ik] + VV[ij, ia, ip, is, it, ik]*phi[ij, ia, ip, is, it, ik]/frac_phi[ij, it, ik]
                        m_coh[ij, ik]      = m_coh[ij, ik] + phi[ij, ia, ip, is, it, ik]/frac_phi[ij, it, ik]
                        beq_coh[ij, it, ik] = beq_coh[ij, it, ik] + a[ia]*(1.0 + rn[it])* (1.0 - psi[ij, it])*phi[ij, ia, ip, is, it, ik]/frac_phi[ij, it, ik]/psi[ij, it] 
                    end
                end
            end
        end
    end

    # normalize VV_coh (because hh excluded)

    for ik in 1:SS
        VV_coh[:, it, ik] = VV_coh[:, it, ik]./m_coh[:, ik]
        FLC[:, it, ik] = FLC[:, it, ik]./m_coh[:, ik]
    end

    # calculate aggregate quantities
    CC[it] = 0.0
    LL[it] = 0.0
    HH[it] = 0.0
    AA[it] = 0.0
    workpop = 0.0
    BQ[it] = 0.0
    for ij in 1:JJ
        for ik in 1:SS
            CC[it] = CC[it] + c_coh[ij, it, ik]*m[ij, it, ik]
            LL[it] = LL[it] + y_coh[ij, it, ik]*m[ij, it, ik]
            HH[it] = HH[it] + l_coh[ij, it, ik]*m[ij, it, ik]
            #AA[it] = AA[it] + a_coh[ij, it]*m[ij, it]
            AA[it] = AA[it] + a_coh[ij, it, ik]*m[ij, it, ik]/psi[ij, it]

            BQ[it] = BQ[it] + beq_coh[ij, it, ik]*m[ij, it, ik]   

            if(ij < JR)
                workpop = workpop + m[ij, it, ik]
            end
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

    global universal

    # last year
    itm = year(it, 2, 1)
    itp = year(it, 1, 2)

    # set government quantities and pension payments
    GG[it] = gy*YY[0]
    BB[it] = by*YY[0]

    #for ik in 1:SS
    #    pen[JR:JJ, it, ik] .= kappa[it]*INC[itm]
    #end

    
    if universal
        pen[JR:JJ, it, 1] .= (kappa[it]*0.5)*INC[itm]
        pen[JR:JJ, it, 2] .= (kappa[it]*0.5)*INC[itm]
    else
        #pen[JR:JJ, it, 1] .= kappa[it]*(0.75*INC[itm])
        #pen[JR:JJ, it, 2] .= (kappa[it]*0.25)*(0.25*INC[itm])
        pen[JR:JJ, it, 1] .= kappa[it]*INC[itm]  
    end
    
    PP[it] = 0.0

    for ij in JR:JJ
        for ik in 1:SS
            PP[it] = PP[it] + pen[ij, it, ik]*m[ij, it, ik]
        end
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
            for ik in 1:SS
                solve_household(ij, 1, ik)
            end
        end
        for it in 1:TT
            for ik in 1:SS
                solve_household(1, it, ik)
            end
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
        PP[it] = PP[0]
        taxrev[:,it] = taxrev[:, 0]
        BQ[it] = BQ[0]

        for ik in 1:SS
            pen[:, it, ik] = pen[:, 0, ik]
            c_coh[:, it, ik] .= c_coh[:, 0, ik]
            l_coh[:, it, ik] = l_coh[:, 0, ik]
            y_coh[:, it, ik] = y_coh[:, 0, ik]
            a_coh[:, it, ik] = a_coh[:, 0, ik]
            aplus[:, :, :, :, it, ik] = aplus[:, :, :, :, 0, ik]
            c[:, :, :, :, it, ik] = c[:, :, :, :, 0, ik]
            l[:, :, :, :, it, ik] = l[:, :, :, :, 0, ik]
            phi[:, :, :, :, it, ik] = phi[:, :, :, :, 0, ik]
            VV[:, :, :, :, it, ik] = VV[:, :, :, :, 0, ik]
            RHS[:, :, :, :, it, ik] = RHS[:, :, :, :, 0, ik]
            beq[:,it, ik] = beq[:,0, ik]
        end
        
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
                    for ik in 1:SS
                        if(phi[ij, ia, ip, is, it, ik] > 1e-8)
                            iamax[ij] = ia
                        end
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