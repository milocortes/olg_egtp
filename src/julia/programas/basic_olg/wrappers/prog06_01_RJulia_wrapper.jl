#=#############################################################################
! PROGRAM LR_OLG
!
! ## Long-run equilibria in the overlapping generations model
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
=##############################################################################

using NLsolve
using DataFrames

function OLG(γ_p, β_p, α_p, δ_p, a1, a2, a3, tax_p , by_p , κ_p , nₚ_p )
	# model parameters
	global γ = γ_p
	global egam = 1-(1/γ)
	global β = β_p
	global α = α_p
	global δ = δ_p
	global by = by_p
	global κ = κ_p
	global nₚ = nₚ_p
	global tax = tax_p
	
	global g = [a1, a2, a3]
	
	# solve the steady state equation system
	
	function eqns(F, x)
	    # model variables
	    global w = 0.0
	    global r = 0.0
	    global wn = 0.0
	    global Rₙ = 0.0
	    global p = 0.0
	    global τw = 0.0
	    global τr = 0.0
	    global τc = 0.0
	    global τp = 0.0
	    global pen = 0.0
	
	    global K = 0.0
	    global L = 0.0 
	    global Y = 0.0
	    global A = 0.0
	    global C = 0.0
	    global B = 0.0
	    global G = 0.0
	    global I = 0.0
	
	    global a = [0.0, 0.0, 0.0]
	    global c = [0.0, 0.0, 0.0] 
	
	    #initialize labor supply, pension payments and tax rates
	
	    global L = (2.0+nₚ)/(1.0+nₚ)
	    global τp = κ/((2.0+nₚ)*(1.0+nₚ))
	    global τc = 0.0
	    global τw = 0.0
	    global τr = 0.0
	
	    K = x[1]
	    
	    if tax == 1
	        τc = x[2]
	    elseif tax==2
	        τw = x[2]
	        τr = x[2]
	    elseif tax == 3
	        τw = x[2]
	    else tax > 3
	        τr = x[2]
	    end
	    
	    # factor prices and pension payments
	    r = α*(K/L)^(α-1.0)-δ
	    w = (1.0-α)*(K/L)^α
	    wn = w*(1.0-τw-τp)
	    Rₙ = 1.0+r*(1.0-τr)
	    p = 1.0+τc
	    pen = κ*w
	
	    # individual decisions
	    PVI = wn + wn/Rₙ + pen/Rₙ^2
	    PSI = p*(1.0 + (β*Rₙ)^γ/Rₙ + (β*Rₙ)^(2*γ)/Rₙ^2)
	
	    c[1] = PVI/PSI
	    c[2] = (β*Rₙ)^γ*c[1]
	    c[3] = (β*Rₙ)^γ*c[2]
	    a[1] = 0.0
	    a[2] = wn - p*c[1]
	    a[3] = wn + Rₙ*a[2] - p*c[2]
	
	    # quantities
	    Y = K^α * L^(1.0-α)
	    C = c[1] + c[2]/(1.0+nₚ) + c[3]/(1.0+nₚ)^2
	    G = g[1] + g[2]/(1.0+nₚ) + g[3]/(1.0+nₚ)^2
	    A = a[2]/(1.0+nₚ) + a[3]/(1.0+nₚ)^2
	    I = (nₚ+δ)*K
	    B = by*Y
	
	    # get equations defining general equilibrium
	    F[1] = K + B - A
	    F[2] = τc*C + τw*w*L + τr*r*A - (r-nₚ)*B - G
	end
	
	x = nlsolve(eqns, [ 0.1; 0.1])


	function calculate_household_utility()
		# Calculate household utility
		global util = 0
		
		for i in 1:3
		    util +=  β^(i-1)*c[i]^egam/egam
		end
		return util
	end

	util = calculate_household_utility()

	results = ( by = round(by; digits=3), kappa = round(κ; digits=2), np = round(nₚ; digits=2), tauw = round(τw; digits=2), 
	tawr = round(τr; digits=2), tauc = round(τc; digits=2), c1 = round(c[1]; digits=2), 
	c2 = round(c[2]; digits=2), c3 = round(c[3]; digits=2), K = round(K; digits=2) , U = round(util; digits=2))

	return	results
     
end
