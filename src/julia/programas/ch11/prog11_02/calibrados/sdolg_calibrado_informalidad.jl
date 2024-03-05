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

include("utils_informalidad.jl")
using OffsetArrays
using Roots

# Get parameters
pais = "MEX"
OLG_params = build_parameters(pais, "parametros_olg.csv")

if pais == "MEX"
    OLG_params["by"] = 0.442
    #OLG_params["kappa"] = 0.3
elseif pais == "CHL"
    OLG_params["by"] = 0.303
elseif pais == "CRI"
    OLG_params["by"] = 0.478
end

#OLG_params["kappa"] = 0.3

# number of transition periods
global TT = 100

# number of years the household lives
global JJ = 16

# number of years the household retires
global JR = 10

# number of persistent shock process values
global NP = 2

# number of transitory shock process values
global NS = 5

# number of points on the asset grid
global NA = 200

# number of skill group
global SS = 2

# household preference parameters
global gamma = 0.18
global egam = 1.0 - 1.0/gamma
global nu    = OLG_params["nu"]
global beta  = 0.998^5

# household risk process
global sigma_theta = 0.23
global sigma_eps   = 0.05
global rho         = 0.98

# production parameters
global alpha = OLG_params["alpha"]
global delta = 1.0-(1.0-OLG_params["delta"])^5
#global delta = 1.0-(1.0-0.08)^5
#global Omega2 = OLG_params["Omega"]
OLG_params["Omega"] = 1.45
global Omega2 = OLG_params["Omega"]

# size of the asset grid
global a_l    = 0.0
global a_u    = 50.0
global a_grow = 0.05

# demographic parameters
global n_p   = (1.0+OLG_params["np"])^5-1.0

# simulation parameters
global damp    = 0.30
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

## Agregamos una seguna dimensión al arreglo pen (JJ, TT, SS)
global pen = OffsetArray(zeros(JJ,TT+1, SS), 1:JJ,0:TT, 1:SS)

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
for param = [:c_coh, :y_coh, :l_coh, :a_coh, :v_coh, :VV_coh, :frac_phi]
    @eval global $param = OffsetArray(zeros(JJ, TT+1, SS), 1:JJ, 0:TT, 1:SS) ;
end

for param = [:GAM, :beq, :beq_coh]
    @eval global $param = OffsetArray(zeros(JJ, TT+1, SS), 1:JJ, 0:TT, 1:SS) ;
end

global omega = zeros(JJ)

global psi = OffsetArray(zeros(JJ+1, TT+1),1:JJ+1, 0:TT)

# the shock process
global dist_theta = zeros(NP)
global theta = zeros(NP)
#pi(NS, NS), eta(NS)
global is_initial = 3

# demographic and other model parameters
# Agregamos una dimensión adicional
global eff = zeros(JJ,SS)

for param = [:m, :pop]
    @eval global $param = OffsetArray(zeros(JJ, TT+1, SS), 1:JJ, 0:TT, 1:SS) ;
end

# individual variables

global a = OffsetArray(zeros(NA+1),0:NA)

# Agregamos una dimensión adicional
global aplus = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1, SS), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT, 1:SS)

for param = [:c, :l, :phi, :VV, :v]
    @eval global $param = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1, SS), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT, 1:SS)
end 

global FLC = OffsetArray(zeros(JJ,TT+1,SS), 1:JJ,0:TT, 1:SS)

# numerical variables
# Agregamos dimensión adicional
global RHS = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1, SS), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT, 1:SS) 
global EV = OffsetArray(zeros(JJ, NA+1, NP, NS, TT+1, SS), 1:JJ, 0:NA, 1:NP, 1:NS, 0:TT, 1:SS) 

# Agregamos una variable adicional para comunicar el grupo de skill
for param = [:ij_com, :ia_com, :ip_com, :is_com, :it_com, :cons_com, :lab_com, :ik_com]
    @eval global $param
end

global DIFF = OffsetArray(zeros(TT+1),0:TT)

global universal = false

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
    for ik in 1:SS
        m[1,it,ik] = 0.5
        GAM[1,it,ik] = omega[1]
        itm = year2(it, -1)
        for ij in 2:JJ
            m[ij,it,ik] = m[ij-1, itm,ik]*psi[ij, it]/(1.0+n_p)
            GAM[1,it, ik] = GAM[1, it, ik] + omega[ij]*m[ij, it,ik]
        end
        for ij in JJ:-1:1
            GAM[ij, it, ik] = omega[ij]/GAM[1, it, ik]
        end
    end
end

# initialize asset grid
grid_Cons_Grow(a, NA+1, a_l, a_u, a_grow)


# get initial guess for savings decision
for ij in 1:JJ
    for ip in 1:NP
        for is in 1:NS
            for ik in 1:SS
                aplus[ij, :, ip, is, 0, ik] .= max.(a/2, a[1]/2) 
            end
        end
    end
end

# initialize age earnings process
eff[1,1] = 1.0000
eff[1,2] = 1.0000
eff[2,1] = 1.3527
eff[3,1] = 1.6952
eff[4,1] = 1.8279
eff[5,1] = 1.9606
eff[6,1] = 1.9692
eff[7,1] = 1.9692
eff[8,1] = 1.9392
eff[9,1] = 1.9007
eff[JR:JJ,1] .= 0.0

eff[2:end,2] = 0.7.*eff[2:end,1]

# initialize fixed effect
dist_theta .= 1.0/float(NP)
theta[1]   = -sqrt(sigma_theta)
theta[2]   = sqrt(sigma_theta)
theta .= exp.(theta)

# calculate the shock process
pi, eta = rouwenhorst(NS, rho, sigma_eps, 0.0);
eta = exp.(eta)

# tax and transfers
tax   .= 3
tauc  .= OLG_params["tauc"]
tauw  .= 0.0
taur  .= 0.0
taup  .= 0.1
kappa .= OLG_params["kappa"]
gy    = OLG_params["gy"]
by    = OLG_params["by"]/5.0

beq[:,0,1] .= 0.0
beq[:,0,2] .= 0.0

BQ[0] = 0.0

# initial guesses for macro variables
KK .= 1.0
LL .= 1.0
YY .= 1.0
II .= (n_p+delta)*KK

GG .= gy*YY[0]
BB .= by*YY[0]

pen .= 0.0
pen[JR:JJ, 0, 1] .= kappa[0]
#pen[JR:JJ, 0, 2] .= kappa[0]
pen[JR:JJ, 0, 2] .= 0

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
          pa = [((1.0+r[0])^(1.0/5.0)-1.0)*100.0,  ""]
);

labour_market = DataFrame( 
    etiqueta = ["valor"],
    L = [LL[0]],
    HH = [HH[0]*100],
    INC = [INC[0]],
    w = [w[0]]
);

good_market = DataFrame( 
          etiqueta = ["valor", "(in %)"],
          Y = [YY[0], YY[0]/YY[0]*100],
          C = [CC[0], CC[0]/YY[0]*100],  
          I = [II[0], II[0]/YY[0]*100], 
          G = [GG[0], GG[0]/YY[0]*100]
); 

gov_accounts = DataFrame(
    etiqueta = ["valor", "(in %)", "(rate)"],
    TAUC = [taxrev[1,0], taxrev[1,0]/YY[0]*100, tauc[0]*100],   
    TAUW = [taxrev[2,0], taxrev[2,0]/YY[0]*100, tauw[0]*100],
    TAUR = [taxrev[3,0], taxrev[3,0]/YY[0]*100, taur[0]*100],
    TOTAL = [taxrev[4,0], taxrev[4,0]/YY[0]*100, ""], 
    G = [GG[0], GG[0]/YY[0]*100, ""],
    B = [BB[0], (BB[0]*5.0)/YY[0]*100, ""]
);

pension_system = DataFrame(
    etiqueta = ["valor", "(in %)"],
    TAUP = [taup[0]*w[0]*LL[0], taup[0]*100], 
    PEN = [pen[JR, 0], kappa[0]],
    PP = [PP[0], PP[0]/YY[0]*100]
);
capital_market
labour_market
good_market
gov_accounts
pension_system


for param = ["capital_market", "labour_market", "good_market", "gov_accounts", "pension_system"]
    file_name = "DSOLG_probs_"*param*"_"*pais*".csv" 
    full_save_path = joinpath(pwd(), "output","csvs", pais, file_name)

    df_to_save = Symbol(param)

    CSV.write(full_save_path, @eval $df_to_save)
end

df_OLG_params = DataFrame(OLG_params)

### Agregamos el parámetro gamma
df_OLG_params.gamma .= gamma

### Agregamos el gasto público que, en estado de equilibrio, refleja el costo necesitado para mantener el nivel de deuda constante.
### Se calcula como : (r - n_p)*BB
### Lo representamos como %PIB 
df_OLG_params.deuda_necesaria .= ((r[0] - n_p)*BB[0]/YY[0])*100


file_name = "DSOLG_probs_OLG_params_"*pais*".csv"
full_save_path = joinpath(pwd(), "output","csvs", pais, file_name)

CSV.write(full_save_path, df_OLG_params)

#### Guardamos datos a nivel de cohorte
## Cohort relative size
m_to_df = OffsetArrays.no_offset_view(m[:,0])
## assets by cohort
a_coh_to_df = OffsetArrays.no_offset_view(a_coh[:,0])
## income by cohort
y_coh_to_df = OffsetArrays.no_offset_view(y_coh[:,0])
## labour hours by cohort
l_coh_to_df = OffsetArrays.no_offset_view(l_coh[:,0])
## consumption by cohort
c_coh_to_df = OffsetArrays.no_offset_view(c_coh[:,0])

## join data in data frame
df_by_cohort = DataFrame(m = m_to_df, a_coh = a_coh_to_df, y_coh = y_coh_to_df, l_coh = l_coh_to_df, c_coh = c_coh_to_df)

## save dataframe
file_name = "DSOLG_probs_OLG_params_by_cohort_"*pais*".csv"
full_save_path = joinpath(pwd(), "output","csvs", pais, file_name)
CSV.write(full_save_path, df_by_cohort)

# set reform parameter (adjsust accordingly for Figure 11.7)
kappa[1:TT] .= 0.5;

global sig     = 1e-4
# calculate transition path without lsra
lsra_on = false;
global universal = false

get_transition()




plot([i for i in 20:5:75],  l_coh[1:12,0,1], title = "Average life-cycle", label = "Hours Worked - Pre-Reforma")
plot!([i for i in 20:5:75],  l_coh[1:12,40,1], label = "Hours Worked - Post-Reforma")


plot([i for i in 20:5:75], y_coh[1:12,0,1] + pen[1:12,0,1], ylimits=(0.2,1), title = "Average life-cycle", label = "Labour-related Income - High Skill (Pre-Reforma)", 
                                                            linestyle=:dot, 
                                                            marker = :circle, 
                                                            markersize = 4,
                                                            seriescolor = :blue)

plot!([i for i in 20:5:75], y_coh[1:12,0,2] + pen[1:12,0,2], label = "Labour-related Income - Low Skill (Pre-Reforma)", 
                                                            linestyle=:dot, 
                                                            marker = :circle, 
                                                            markersize = 4,
                                                            seriescolor = :red)

plot!([i for i in 20:5:75], y_coh[1:12,0,2] , label = "Labour-related Income - Low Skill-Sin Pensión(Pre-Reforma)", 
                                                            linestyle=:dot, 
                                                            marker = :circle, 
                                                            markersize = 4,
                                                            seriescolor = :orange)

plot!([i for i in 20:5:75], y_coh[1:12,40,1] + pen[1:12,40,1], label = "Labour-related Income - High Skill (Post-Reforma)", 
                                                            linestyle=:dash, 
                                                            marker = :sticks, 
                                                            markersize = 4,
                                                            seriescolor = :blue)

plot!([i for i in 20:5:75], y_coh[1:12,40,2] + pen[1:12,40,2], label = "Labour-related Income - Low Skill (Post-Reforma)", 
                                                            linestyle=:dash, 
                                                            marker = :sticks, 
                                                            markersize = 4,
                                                            seriescolor = :red)

