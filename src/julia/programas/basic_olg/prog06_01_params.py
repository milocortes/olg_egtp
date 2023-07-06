"""
python3 -m pip install --user julia

import julia
julia.install()
"""

from julia.api import Julia
jl = Julia(compiled_modules=False)

import os 
PATH = os.getcwd()

from julia import Main
Main.include(os.path.join(PATH, "wrappers", "prog06_01_PyJulia_wrapper.jl")  )

from typing import Dict
import numpy as np 

def OLG(γ : float, β : float, α : float, δ : float, a1 : float, a2 : float, a3 : float, tax : int, by : float, κ :float, nₚ : float) -> Dict[str, np.array] :
    """
    *****************************************************************************
    ## Long-run equilibria in the overlapping generations model
    
     This code is published under the GNU General Public License v3
                             (https://www.gnu.org/licenses/gpl-3.0.en.html)
    
     Authors: Hans Fehr and Fabian Kindermann
              contact@ce-fortran.com
    *****************************************************************************

    ------------------------
            Arguments
    ------------------------
    - γ     :   Coefficient of relative risk aversion (or the reciprocal of the elasticity of intertemporal substitution)
    - β     :   Time discount factor
    - α     :   Capital intensity elasticity
    - δ     :   Depreciation rate
    - a1    :   Per capita coeficient of cohort 1 for public good provision
    - a2    :   Per capita coeficient of cohort 2 for public good provision
    - a3    :   Per capita coeficient of cohort 3 for public good provision
    - tax   :   Tax regimen
    - by    :   Deficit path
    - κ     :   Replacement rate of the pension system.
    - nₚ    :   Population growth rate

    ------------------------
            Output 
    ------------------------

    - np    :   Population growth rate
    - c3    :   Consumption cohort 3
    - tauw  :   Cabor-income tax rate
    - by    :   Deficit path
    - tauc  :   Consumption tax rate
    - c1    :   Consumption cohort 1
    - c2    :   Consumption cohort 2
    - U     :   Utility
    - kappa :   Replacement rate of the pension system.
    - tawr  :   Capital-income tax rate
    - K     :   Capital
    """

    
    return Main.OLG(γ, β, α, δ, a1, a2, a3, tax, by, κ, nₚ)

model_parameters = {
    "γ" : 0.5,   # Coefficient of relative risk aversion (or the reciprocal of the elasticity of intertemporal substitution)
    "β" : 0.9,   # Time discount factor
    "α" : 0.3,   # Capital intensity elasticity
    "δ" : 0.0,   # Depreciation rate
    "a1" : 0.12, # Per capita coeficient of cohort 1 for public good provision
    "a2" : 0.12, # Per capita coeficient of cohort 2 for public good provision
    "a3" : 0.0   # Per capita coeficient of cohort 3 for public good provision
}

## 0) Consumer taxes (by=kappa=0.0, np = 0.20, tax = 1)
OLG(**model_parameters, tax = 1, by = 0.0, κ = 0.0, nₚ = 0.2)

## 1) Wage and capital taxes (by=kappa=0.0, np = 0.20, tax = 2)
OLG(**model_parameters, tax = 2, by = 0.0, κ = 0.0, nₚ = 0.2)

## 2) Capital taxes (by=kappa=0.0, np = 0.20, tax = 4)
OLG(**model_parameters, tax = 4, by = 0.0, κ = 0.0, nₚ = 0.2)

## 3) Wage taxes (by=kappa=0.0, np = 0.20, tax = 3)
OLG(**model_parameters, tax = 3, by = 0.0, κ = 0.0, nₚ = 0.2)

## 4) Wage taxes (by=-0.059, kappa=0.0, np = 0.20, tax = 1)
OLG(**model_parameters, tax = 3, by = -0.059, κ = 0.0, nₚ = 0.2)

## 5) Consumer taxes (by=0.0, kappa=0.50, np = 0.20, tax = 4)
OLG(**model_parameters, tax = 1, by = 0.0, κ = 0.50, nₚ = 0.2)

## 6) Consumer taxes (by=0.099, kappa=0.0, np = 0.20, tax = 4)
OLG(**model_parameters, tax = 1, by = 0.099, κ = 0.0, nₚ = 0.2)

## 7) Consumer taxes (by=0.0, kappa=0.0, np = 0.0, tax = 4)
OLG(**model_parameters, tax = 1, by = 0.0, κ = 0.0, nₚ = 0.0)

