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
    psi = sqrt(N-1)*sqrt(sigma_eta)

    w_est = zeros(N)
    w_est = 1.0/float(N)
    for in in 1:10000
        w_est = *(transpose(p), w_est)
    end
 
    # [-psi + 2.0*psi*float(i-1)/float(n-1) for i in 1:N]
    return p , [-psi + (2.0*psi* ((i-1)/(N-1))) for i in 1:N], w_est
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

# the first order condition
function foc(x_in)
    
    global ij_com
    global cons_com
    global lab_com
    
    # calculate tomorrows assets
    a_plus  = x_in

    # calculate the wage rate
    wage = wn*eff[ij_com]*theta[ip_com]*eta[is_com]

    # calculate available resources
    available = (1.0+rn)*a[ia_com] + pen[ij_com]

    # determine labor
    if(ij_com < JR)
        global lab_com = min( max( nu + (1.0-nu)*(a_plus - available)/wage, 0.0), 1.0-1e-10)
    else
        global lab_com = 0.0
    end

    # calculate consumption
    cons_com = max( (available + wage*lab_com - a_plus)/p , 1e-10)
    # calculate linear interpolation for future part of first order condition
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    tomorrow = varphi*RHS[ij_com+1, ial, ip_com, is_com] + (1.0-varphi)*RHS[ij_com+1, iar, ip_com, is_com]

    # calculate first order condition for consumption
    return margu(cons_com, lab_com)^(-gamma) - tomorrow

end 


# calculates marginal utility of consumption
function margu(cons, lab)

    return nu/p*(cons^nu*(1.0-lab)^(1.0-nu))^egam/cons

end 


# calculates the value function
function valuefunc(a_plus, cons, lab, ij, ip, is)


    global ial_v 
    global iar_v 
    global varphi_v

    # check whether consumption or leisure are too small
    c_help = max(cons, 1e-10)
    l_help = min(max(lab, 0.0),1.0-1e-10)

    # get tomorrows utility
    ial, iar, varphi = linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)

    # calculate tomorrow's part of the value function
    valuefunc = 0.0
    if(ij < JJ)
        valuefunc = max(varphi*EV[ij+1, ial, ip, is] + (1.0-varphi)*EV[ij+1, iar, ip, is], 1e-10)^egam/egam
    end

    # add todays part and discount
    return (c_help^nu*(1.0-l_help)^(1.0-nu))^egam/egam + beta*valuefunc

end 

# for calculating the rhs of the first order condition at age ij
function interpolate(ij)

    for ia in 0:NA
        for ip in 1:NP
            for is in 1:NS

                # calculate the RHS of the first order condition
                RHS[ij, ia, ip, is] = 0.0
                EV[ij, ia, ip, is] = 0.0
                for is_p in 1:NS
                    chelp = max(c[ij, ia, ip, is_p],1e-10)
                    lhelp = max(l[ij, ia, ip, is_p],1e-10)
                    RHS[ij, ia, ip, is] = RHS[ij, ia, ip, is] + pi[is, is_p]*margu(chelp, lhelp)
                    EV[ij, ia, ip, is]  = EV[ij, ia, ip, is] + pi[is, is_p]*V[ij, ia, ip, is_p]
                end
                RHS[ij, ia, ip, is] = ((1.0+rn)*beta*RHS[ij, ia, ip, is])^(-gamma)
                EV[ij, ia, ip, is] = (egam*EV[ij, ia, ip, is])^(1.0/egam)
            end
        end
    end

end 


function grid_Cons_Grow(a, n, left, right, growth)
    ccall((:grid_Cons_Grow, "./mod_julfort.so"),
                Cvoid,
                (Ptr{Float64},Ref{Int64},Ref{Float64}, Ref{Float64}, Ref{Float64}),
                a, n,  left, right,growth)
end


function linint_Grow(x, left, right, growth, n, il, ir, phi)

    ccall((:linint_Grow, "./mod_julfort.so"),
            Cvoid,
            (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64} ),
            x, left, right, growth, n, il, ir, phi)

    return il[1], ir[1], phi[1]
end