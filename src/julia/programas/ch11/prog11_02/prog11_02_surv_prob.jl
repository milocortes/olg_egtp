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

include("utils_surv_prob.jl")
using OffsetArrays
using Roots

# number of transition periods
global TT = 40

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
global Omega2 = 1.60

# size of the asset grid
global a_l    = 0.0
global a_u    = 50.0
global a_grow = 0.05

# demographic parameters
global n_p   = (1.0+0.02)^5-1.0

# simulation parameters
global damp    = 0.30
global sig     = 1e-4
global itermax = 70

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
pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
eta = exp.(eta)

# tax and transfers
tax   .= 2
tauc  .= 0.075
tauw  .= 0.0
taur  .= 0.0
taup  .= 0.1
kappa .= 0.5
gy    = 0.19
by    = 0.60/5.0

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


# set reform parameter (adjsust accordingly for Figure 11.7)
#kappa[1:TT] .= 0.0
kappa[1:TT] .= 0.5;


# calculate transition path without lsra
lsra_on = false;

get_transition()



# calculate transition path with lsra
lsra_on = true;
get_transition()




plot([i for i in 20:5:75],  l_coh[:,0], title = "Average life-cycle", label = "Hours Worked - Pre-Reforma")
plot!([i for i in 20:5:75],  l_coh[:,40], label = "Hours Worked - Post-Reforma")


plot([i for i in 20:5:75],  c_coh[:,0], title = "Average life-cycle", label = "Consumption - Pre-Reforma")
plot!([i for i in 20:5:75],  c_coh[:,40], label = "Consumption- Post-Reforma")

plot([i for i in 20:5:75], y_coh[:,0] + pen[:,0], title = "Average life-cycle", label = "Labour-related Income - Pre-Reforma")
plot!([i for i in 20:5:75], y_coh[:,40] + pen[:,40], label = "Labour-related Income - Post-Reforma")


plot([i for i in 20:5:75], c_coh, label = "Consumption", title = "Average life-cycle", xlabel = "Age j", ylabel = "Mean")

plot!([i for i in 20:5:75], y_coh + pen, label = "Labour-related Income")
