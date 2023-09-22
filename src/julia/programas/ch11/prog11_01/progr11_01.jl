#=##############################################################################
# PROGRAM SOLG_LR
#
# ## Long-run equilibria in the stochastic OLG model
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr and Fabian Kindermann
#          contact@ce-fortran.com
#
=##############################################################################
include("utils.jl")

using OffsetArrays
using Roots

global ial_v = Array{Int64}(undef, 1)
global iar_v = Array{Int64}(undef, 1)
global varphi_v = zeros(1)

# number of years the household lives
global JJ = 12

# number of years the household retires
global JR = 10

# number of persistent shock process values
global NP = 2

# number of transitory shock process values
global NS = 5

# number of points on the asset grid
global NA = 100

# household preference parameters
global gamma = 0.50
global egam = 1.0 - 1.0/gamma
global nu    = 0.335
global beta  = 0.998^5

# household risk process
global sigma_theta = 0.23
global sigma_eps   = 0.05
global rho         = 0.98

# production parameters
global alpha = 0.36
global delta = 1.0-(1.0-0.0823)^5
global Omega = 1.60

# size of the asset grid
global a_l    = 0.0
global a_u    = 35.0
global a_grow = 0.05

# demographic parameters
global n_p   = (1.0+0.01)^5-1.0

# simulation parameters
global damp    = 0.30
global sig     = 1e-4
global itermax = 50


# macroeconomic variables
# prices
for param = [:r, :rn, :w, :wn, :p]
    @eval global $param ;
end

# capital market
for param = [:KK, :AA, :BB, :LL, :HH]
    @eval global $param ;
end

# goods market
for param = [:YY, :CC, :II, :GG, :INC]
    @eval global $param;
end


# government variables
for param = [:tauc, :tauw, :taur, :taup, :kappa]
    @eval global $param;
end

global gy
global by
global pen = zeros(JJ)
global PP
global taxrev = zeros(4)
global tax
global reform_on

# cohort aggregate variables
for param = [:c_coh, :y_coh, :l_coh, :a_coh, :v_coh]
    @eval global $param = zeros(JJ);
end

# the shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)
global pi 
global eta 
global is_initial = 3

# demographic and other model parameters
global m = zeros(JJ)
global eff = zeros(JJ)

# individual variables
global a = OffsetArray(zeros(NA+1), 0:NA)
global aplus =  OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS);
global c =  OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS);
global l =  OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS);
global phi =  OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS);
global V =  OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS);

# numerical variables
global RHS =  OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS);
global EV =  OffsetArray(zeros(JJ, NA+1, NP, NS), 1:JJ, 0:NA, 1:NP, 1:NS);

for param = [:c_coh, :y_coh, :l_coh, :a_coh, :v_coh]
    @eval global $param = zeros(JJ);
end

for param = [:ij_com, :ia_com, :ip_com, :is_com, :it_com]
    @eval global $param 
end

for param = [:cons_com, :lab_com, :DIFF, :INC_init]
    @eval global $param 
end


# initialize variables

# set up population structure
for ij in 1:JJ
    m[ij] = (1.0+n_p)^(1.0-ij)
end 

# initialize asset grid
grid_Cons_Grow(a, NA+1, a_l, a_u, a_grow)

# get initial guess for savings decision
for ij in 1:JJ
    for ip in 1:NP
        for is in 1:NS
            aplus[ij, :, ip, is] .= max.(a/2, a[1]/2)
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
theta = exp.(theta)

# calculate the shock process
#call discretize_AR(rho, 0d0, sigma_eps, eta, pi)
#eta = exp(eta)
pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
eta = exp.(eta)

# tax and transfers
tax   = 2
tauc  = 0.075
tauw  = 0.0
taur  = 0.0
taup  = 0.1
kappa = 0.5
gy    = 0.19
by    = 0.60/5.0

# initial guesses for macro variables
KK = 1.0
LL = 1.0
YY = 1.0
II = (n_p+delta)*KK

GG = gy*YY
BB = by*YY

pen .= 0.0
pen[JR:JJ] .= kappa

# calculate initial equilibrium
reform_on = false

# computes the initial steady state of the economy
# subroutine for calculating prices
function prices()

    global r = Omega*alpha*(KK/LL)^(alpha-1.0)-delta
    global w = Omega*(1.0-alpha)*(KK/LL)^alpha
    global rn = r*(1.0-taur)
    global wn = w*(1.0-tauw-taup)
    global p = 1.0 + tauc

end 


# solve the household problem

# determines the solution to the household optimization problem
function solve_household()

    global cons_com
    global lab_com
    # get decision in the last period of life
    for ia in 0:NA
        aplus[JJ, ia, :, :] .= 0.0
        c[JJ, ia, :, :] .= ((1.0+rn)*a[ia] + pen[JJ])/p
        l[JJ, ia, :, :] .= 0.0
        V[JJ, ia, :, :] .= valuefunc(0.0, c[JJ, ia, 1, 1], l[JJ, ia, 1, 1], JJ, 1, 1)
    end

    # interpolate individual RHS
    interpolate(JJ)

    for ij in JJ-1:-1:1

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
            if(ij >= JR && ia == 0 && kappa <= 1e-10)
                aplus[ij, ia, :, :] .= 0.0
                c[ij, ia, :, :] .= 0.0
                l[ij, ia, :, :] .= 0.0
                V[ij, ia, :, :] .= valuefunc(0.0, 0.0, 0.0, ij, 1, 1)
                continue
            end

            for ip in 1:ip_max
                for is in 1:is_max

                    # get initial guess for the individual choices
                    x_in = aplus[ij, ia, ip, is]

                    # set up communication variables
                    global ij_com = ij
                    global ia_com = ia
                    global ip_com = ip
                    global is_com = is

                    # solve the household problem using rootfinding
                    x_root = fzero(foc, x_in)

                    #println(ij, " ", ia, " ", ip, " ", is, " ", foc(x_root), x_root)

                    # write screen output in case of a problem
                    #if(check)write(*,'(a, 4i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is

                    # check for borrowing constraint
                    if(x_root < 0.0)
                        x_root = 0.0
                        wage = wn*eff[ij]*theta[ip]*eta[is]
                        available = (1.0+rn)*a[ia] + pen[ij]

                        if(ij < JR)
                            lab_com = min( max(nu-(1.0-nu)*available/wage , 0.0) , 1.0-1e-10)
                        else
                            lab_com = 0.0
                        end

                        cons_com = max( (available + wage*lab_com)/p , 1e-10)
                    end
                    # copy decisions
                    aplus[ij, ia, ip, is] = x_root
                    c[ij, ia, ip, is] = cons_com
                    l[ij, ia, ip, is] = lab_com
                    V[ij, ia, ip, is] = valuefunc(x_root, cons_com, lab_com, ij, ip, is)

                end

                # copy decision in retirement age
                if(ij >= JR)
                    aplus[ij, ia, :, :] .= aplus[ij, ia, 1, 1]
                    c[ij, ia, :, :] .= c[ij, ia, 1, 1]
                    l[ij, ia, :, :] .= l[ij, ia, 1, 1]
                    V[ij, ia, :, :] .= V[ij, ia, 1, 1]
                end
            end
        end
        # interpolate individual RHS
        interpolate(ij)
    end

end 

# calculate the distribution of households over state space
# determines the invariant distribution of households
function get_distribution()

    # set distribution to zero
    phi[:, :, :, :] .= 0.0

    # get initial distribution in age 1
    for ip in 1:NP
        phi[1, 0, ip, is_initial] = dist_theta[ip]
    end

    # successively compute distribution over ages
    for ij in 2:JJ

        # iterate over yesterdays gridpoints
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS

                    # interpolate yesterday's savings decision
                    ial, iar, varphi = linint_Grow(aplus[ij-1, ia, ip, is], a_l, a_u, a_grow, NA, ial_v, iar_v, varphi_v)
                    # restrict values to grid just in case
                    ial = max(min(ial, NA-1), 0)
                    iar = max(min(iar, NA), 1)
                    varphi = max(min(varphi, 1.0), 0.0)

                    # redistribute households
                    for is_p in 1:NS
                        phi[ij, ial, ip, is_p] = phi[ij, ial, ip, is_p] + pi[is, is_p]*varphi*phi[ij-1, ia, ip, is]
                        phi[ij, iar, ip, is_p] = phi[ij, iar, ip, is_p] + pi[is, is_p]*(1.0-varphi)*phi[ij-1, ia, ip, is]
                    end
                end
            end
        end
    end

end 

# aggregate individual decisions over cohorts
# subroutine for calculating quantities
function aggregation()

    global LL
    global KK

    LL_old = LL

    # calculate cohort aggregates
    c_coh[:] .= 0.0
    l_coh[:] .= 0.0
    y_coh[:] .= 0.0
    a_coh[:] .= 0.0
    v_coh[:] .= 0.0

    for ij in 1:JJ
        for ia in 0:NA
            for ip in 1:NP
                for is in 1:NS
                    c_coh[ij] = c_coh[ij] + c[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    l_coh[ij] = l_coh[ij] + l[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    y_coh[ij] = y_coh[ij] + eff[ij]*theta[ip]*eta[is]*l[ij, ia, ip, is]*phi[ij, ia, ip, is]
                    a_coh[ij] = a_coh[ij] + a[ia]*phi[ij, ia, ip, is]
                    v_coh[ij] = v_coh[ij] + V[ij, ia, ip, is]*phi[ij, ia, ip, is]
                end
            end
        end
    end

    # calculate aggregate quantities
    global CC = 0.0
    global LL = 0.0
    global HH = 0.0
    global AA = 0.0
    global workpop = 0.0

    for ij in 1:JJ
        global CC = CC + c_coh[ij]*m[ij]
        global LL = LL + y_coh[ij]*m[ij]
        global HH = HH + l_coh[ij]*m[ij]
        global AA = AA + a_coh[ij]*m[ij]
        if(ij < JR)
            global workpop = workpop + m[ij]
        end
    end

    # damping and other quantities
    global KK = damp*(AA-BB) + (1.0-damp)*KK
    global LL = damp*LL + (1.0-damp)*LL_old
    global II = (n_p+delta)*KK
    global YY = Omega * KK^alpha * LL^(1.0-alpha)

    # get average income and average working hours
    global INC = w*LL/workpop
    global HH  = HH/workpop

    # get difference on goods market
    global DIFF = YY-CC-II-GG

end 
# determine the government parameters
# subroutine for calculating government parameters
function government()

    global tauc

    # set government quantities and pension payments
    if(!reform_on)
        global GG = gy*YY
        global BB = by*YY
    end

    # calculate government expenditure
    global expend = GG + (1.0+r)*BB - (1.0+n_p)*BB

    # get budget balancing tax rate
    if(tax == 1)
        global tauc = (expend - (tauw*w*LL + taur*r*AA))/CC
        global p    = 1.0 + tauc
    elseif(tax == 2)
        global tauw = (expend - tauc*CC)/(w*LL + r*AA)
        global taur = tauw
    elseif(tax == 3)
        global tauw = (expend - (tauc*CC + taur*r*AA))/(w*LL)
    else
        global taur = (expend - (tauc*CC + tauw*w*LL))/(r*AA)
    end

    taxrev[1] = tauc*CC
    taxrev[2] = tauw*w*LL
    taxrev[3] = taur*r*AA
    taxrev[4] = sum(taxrev[1:3])

    # get budget balancing social security contribution
    pen[JR:JJ] .= kappa*INC
    global PP = 0.0
    for ij in JR:JJ
        global PP = PP + pen[ij]*m[ij]
    end

    taup = PP/(w*LL)

end 


# computes the initial steady state of the economy
function get_SteadyState()

    # iterate until value function converges
    for iter in 1:itermax

        # derive prices
        prices()

        # solve the household problem
        solve_household()

        # calculate the distribution of households over state space
        get_distribution()

        # aggregate individual decisions over cohorts
        aggregation()

        # determine the government parameters
        government()

        println(iter,"     ", 5.0*KK/YY*100.0, "   ", CC/YY*100.0, "   ", II/YY*100.0, "   ", r, "   ", w, "   ", DIFF/YY*100.0)

        if(abs(DIFF/YY)*100.0 < sig)
            break
        end

    end

end 

println("INITIAL EQUILIBRIUM")
println("ITER     K/Y     C/Y     I/Y       r       w        DIFF")
get_SteadyState()

using Plots 

plot([i for i in 20:5:75], c_coh, label = "Consumption", title = "Average life-cycle", xlabel = "Age j", ylabel = "Mean")
plot!([i for i in 20:5:75], l_coh, label = "Hours Worked")
plot!([i for i in 20:5:75], y_coh + pen, label = "Labour-related Income")


# set reform variables
reform_on = true

# set reform values (Table 11.3, different rows)
tax = 1    ;    tauw = 0.0    ;    taur = 0.0  # Simulation (1)
#tax = 3    ;    taur = 0d0                     ! Simulation (2)
#kappa = 0d0                                    ! Simulation (3)

# calculate final equilibrium
println("REFORM")
println("ITER     K/Y     C/Y     I/Y       r       w        DIFF")

get_SteadyState()


