library("JuliaConnectoR")

PATH = "/home/milo/Documents/egap/BID/OLG/olg_egtp/src/julia/programas/basic_olg"

juliaCall("include", file.path(PATH, "wrappers", "prog06_01_RJulia_wrapper.jl"))

OLG <- function(γ, β, α, δ, a1, a2, a3, tax, by, κ, nₚ){
    return (data.frame(juliaGet(juliaFun("OLG")(γ, β, α, δ, a1, a2, a3, tax, by, κ, nₚ))))
}

## 0) Consumer taxes (by=kappa=0.0, np = 0.20, tax = 1)
OLG(γ = 0.5, β = 0.9, α = 0.3, δ = 0.0, a1 = 0.12, a2 = 0.12, a3 = 0.0, tax = 1, by = 0.0, κ = 0.0, nₚ = 0.2)

## 1) Wage and capital taxes (by=kappa=0.0, np = 0.20, tax = 2)
OLG(γ = 0.5, β = 0.9, α = 0.3, δ = 0.0, a1 = 0.12, a2 = 0.12, a3 = 0.0, tax = 2, by = 0.0, κ = 0.0, nₚ = 0.2)

## 2) Capital taxes (by=kappa=0.0, np = 0.20, tax = 4)
OLG(γ = 0.5, β = 0.9, α = 0.3, δ = 0.0, a1 = 0.12, a2 = 0.12, a3 = 0.0, tax = 4, by = 0.0, κ = 0.0, nₚ = 0.2)

## 3) Wage taxes (by=kappa=0.0, np = 0.20, tax = 3)
OLG(γ = 0.5, β = 0.9, α = 0.3, δ = 0.0, a1 = 0.12, a2 = 0.12, a3 = 0.0, tax = 3, by = 0.0, κ = 0.0, nₚ = 0.2)

## 4) Wage taxes (by=-0.059, kappa=0.0, np = 0.20, tax = 1)
OLG(γ = 0.5, β = 0.9, α = 0.3, δ = 0.0, a1 = 0.12, a2 = 0.12, a3 = 0.0, tax = 3, by = -0.059, κ = 0.0, nₚ = 0.2)

## 5) Consumer taxes (by=0.0, kappa=0.50, np = 0.20, tax = 4)
OLG(γ = 0.5, β = 0.9, α = 0.3, δ = 0.0, a1 = 0.12, a2 = 0.12, a3 = 0.0, tax = 1, by = 0.0, κ = 0.50, nₚ = 0.2)

## 6) Consumer taxes (by=0.099, kappa=0.0, np = 0.20, tax = 4)
OLG(γ = 0.5, β = 0.9, α = 0.3, δ = 0.0, a1 = 0.12, a2 = 0.12, a3 = 0.0, tax = 1, by = 0.099, κ = 0.0, nₚ = 0.2)

## 7) Consumer taxes (by=0.0, kappa=0.0, np = 0.0, tax = 4)
OLG(γ = 0.5, β = 0.9, α = 0.3, δ = 0.0, a1 = 0.12, a2 = 0.12, a3 = 0.0, tax = 1, by = 0.0, κ = 0.0, nₚ = 0.0)
