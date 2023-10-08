###############################################################################
# PROGRAM SOLG_TR
#
# ## The stochastic OLG model with transitional dynamics
#
# This code is published under the GNU General Public License v3
#                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Authors: Hans Fehr and Fabian Kindermann
#          contact@ce-fortran.com
#
###############################################################################
include("utils.jl")
using OffsetArrays
using Roots

# Get parameters
alpha_param, Omega_param, delta_param, nu_param, np_param, gy_param , tauc_param = build_parameters("MEX","parametros_olg.csv")
by_param = 0.44
kappa_param = 0.3

# number of transition periods
global TT = 40

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
global gamma = 0.22
global egam = 1.0 - 1.0/gamma
global nu    = nu_param
global beta  = 0.998^5

# household risk process
global sigma_theta = 0.23
global sigma_eps   = 0.05
global rho         = 0.98

# production parameters
global alpha = alpha_param
global delta = 1.0-(1.0-delta_param)^5
global Omega = Omega_param

# size of the asset grid
global a_l    = 0.0
global a_u    = 35.0
global a_grow = 0.05

# demographic parameters
global n_p   = (1.0+np_param)^5-1.0

# simulation parameters
global damp    = 0.30
global sig     = 1e-4
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
for param = [:YY, :CC, :II, :GG, :INC]
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

# cohort aggregate variables
for param = [:c_coh, :y_coh, :l_coh, :a_coh, :v_coh, :VV_coh]
    @eval global $param = OffsetArray(zeros(JJ, TT+1), 1:JJ, 0:TT) ;
end

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


# Calculate initial equilibrium
# set up population structure
for ij in 1:JJ
    pop[ij, 0] = 1.0/(1.0+n_p)^(ij-1)
end

for ij in 1:JJ
    m[ij, 0] = pop[ij, 0]/pop[1, 0]
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
tax   .= 1
tauc  .= 0.0
tauw  .= 0.3
taur  .= 0.0
taup  .= 0.1
kappa .= kappa_param
gy    = gy_param
by    = by_param/5.0

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

capital_market = DataFrame( 
          etiqueta = ["valor", "(in %)"],
          K = [KK[0], KK[0]/YY[0]*500],
          A = [AA[0], AA[0]/YY[0]*500],  
          B = [BB[0], BB[0]/YY[0]*500], 
          BA = [BA[0], BA[0]/YY[0]*500],   
          r = [r[0], ""],
          pa = [((1.0+r[0])^(1.0/5.0)-1.0)*100.0,  ""]) 

labour_market = DataFrame( 
    etiqueta = ["valor"],
    L = [LL[0]],
    HH = [HH[0]*100],
    INC = [INC[0]],
    w = [w[0]]
) 

good_market = DataFrame( 
          etiqueta = ["valor", "(in %)"],
          Y = [YY[0], YY[0]/YY[0]*100],
          C = [CC[0], CC[0]/YY[0]*100],  
          I = [II[0], II[0]/YY[0]*100], 
          G = [GG[0], GG[0]/YY[0]*100]) 

gov_accounts = DataFrame(
    etiqueta = ["valor", "(in %)", "(rate)"],
    TAUC = [taxrev[1,0], taxrev[1,0]/YY[0]*100, tauc[0]*100],   
    TAUW = [taxrev[2,0], taxrev[2,0]/YY[0]*100, tauw[0]*100],
    TAUR = [taxrev[3,0], taxrev[3,0]/YY[0]*100, taur[0]*100],
    TOTAL = [taxrev[4,0], taxrev[4,0]/YY[0]*100, ""], 
    G = [GG[0], GG[0]/YY[0]*100, ""],
    B = [BB[0], (BB[0]*5.0)/YY[0]*100, ""]
)

pension_system = DataFrame(
    etiqueta = ["valor", "(in %)"],
    TAUP = [taup[0]*w[0]*LL[0], taup[0]*100], 
    PEN = [pen[JR, 0], kappa[0]],
    PP = [PP[0], PP[0]/YY[0]*100]
)
capital_market
labour_market
good_market
gov_accounts
pension_system
# set reform parameter (adjsust accordingly for Figure 11.7)
#kappa[1:TT] .= 0.0
kappa[1:TT] .= 0.33


# calculate transition path without lsra
lsra_on = false

get_transition()



# calculate transition path with lsra
lsra_on = true
get_transition()




plot([i for i in 20:5:75],  l_coh[:,0], title = "Average life-cycle", label = "Hours Worked - Pre-Reforma")
plot!([i for i in 20:5:75],  l_coh[:,40], label = "Hours Worked - Post-Reforma")


plot([i for i in 20:5:75],  c_coh[:,0], title = "Average life-cycle", label = "Consumption - Pre-Reforma")
plot!([i for i in 20:5:75],  c_coh[:,40], label = "Consumption- Post-Reforma")

plot([i for i in 20:5:75], y_coh[:,0] + pen[:,0], title = "Average life-cycle", label = "Labour-related Income - Pre-Reforma")
plot!([i for i in 20:5:75], y_coh[:,40] + pen[:,40], label = "Labour-related Income - Post-Reforma")


plot([i for i in 20:5:75], c_coh, label = "Consumption", title = "Average life-cycle", xlabel = "Age j", ylabel = "Mean")

plot!([i for i in 20:5:75], y_coh + pen, label = "Labour-related Income")
