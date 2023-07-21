### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ e574ec5e-21a5-11ee-16f0-353b47d26bd3
md"""
# Introduction to Computational Economics Using Fortran
"""

# ╔═╡ 610bfad9-97f6-4734-81b9-4199949364cd
md"""
# Cap 7: Extending the OLG model
* Assumption of a certain lifespan in the previous OLG models can be justified by arguing that households buy a perfect annuity insurance (annuity puzzle).
* ¿Implies for public policy? Consequences of the introduction of a perfect annuity market.
* In the previous programs of this chapter, agents knew with certainty when their life was over.
* In a setup with uncertain survival, households may die prior to their maximum lifespan (J).
* Therefore, households leave unintended bequests $b_{j,t}$ unintended bequest an agent at age j receives at date t.
"""

# ╔═╡ 3dda2831-d9d9-4810-b2ca-b91192eeeccc
md"""

## 7.3 Longevity risk and annuitization

* Survival probabilities, before: $\psi_{j,t}=1$

* A household that enters the labor market in period t reaches adulthood with certainty:

```math
\psi_{1,t}=1-q_{1,t}
```

```math
\psi_{2,t+1}=1-q_{2,t+1}
```

```math
\psi_{3,t+2}=1-q_{3,t+2}
```

* Lifespan is restricted to three periods 
```math
\psi_{4,t+3}=0
```

* Cohorts weights:
```math
m_{1,t}=1
```

```math
m_{j,t}= \dfrac{\psi_{j,t}}{1+n_{p,t}}*m_{j-1,t-1}
```
"""

# ╔═╡ 0f862bf8-9958-4528-a253-4ab00ad5dfdb
md"""

### The Household decision problem

* Households born in period t maximize the time-separable utility function.

```math
U_{t}=U(c_{1,t},c_{2,t+1},c_{3,t+2})=\sum_{j=1}^{3}\beta^{j-1}(\Pi_{i=1}^{j}\psi_{i,j})u(c_{j,s})
```

* Now future consumption is further discounted by the survival probabilities
* Intertemporal budget constraint:

```math
\sum_{j=1}^{3}{\dfrac{p_{s}c_{j,s}}{\Pi_{k=t+1}^{s}R_{k}^{n}}}=\sum_{j=1}^{3}{\dfrac{w_{s}^{n}+b_{j,s}+pen_{j,s}}{\Pi_{k=t+1}^{s}R_{k}^{n}}}=:W_{1,t}
```

```math
pen_{j,s}=0 for j<3
```

Total assets are distributed over the surviving cohort members with a flexible distribution scheme:

```math
\Gamma_{j,s}=\dfrac{\omega_{b,j}}{\sum_{j=1}^{3}\omega_{b,i}m_{i,s}}
```

Cohort amount of unintended bequests:
```math
b_{j,t}=\Gamma_{j,t}BQ_{t}
```


* First order conditions:

```math
p_{t+1}u'(c_{1,t})=\beta\psi_{2,t+1}R_{t+1}^{n}p_{t}u'{c_{2,t+1}}
```

```math
p_{t+2}u'(c_{2,t+1})=\beta\psi_{3,t+2}R_{t+2}^{n}p_{t}u'{c_{3,t+2}}
```

The CRRA utility function given by

```math
	 u(c) = \dfrac{c^{1-\frac{1}{\gamma}}}{1-\frac{1}{\gamma}} \quad 
```

* Consumption path

```math
	 c_{2,t+1}= 
```



"""

# ╔═╡ Cell order:
# ╟─e574ec5e-21a5-11ee-16f0-353b47d26bd3
# ╟─610bfad9-97f6-4734-81b9-4199949364cd
# ╟─3dda2831-d9d9-4810-b2ca-b91192eeeccc
# ╠═0f862bf8-9958-4528-a253-4ab00ad5dfdb
