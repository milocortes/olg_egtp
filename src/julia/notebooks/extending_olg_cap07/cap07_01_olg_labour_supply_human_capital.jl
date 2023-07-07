### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 29cebbb7-6986-4d47-96ad-0ccebb993dfc
html"<button onclick='present()'>present</button>"

# ╔═╡ 63da7e09-da7d-441a-b967-892d7dcf1f78
md"""

# Introduction to Computational Economics Using Fortran
"""

# ╔═╡ 75ba6402-0628-11ee-02ad-bbd3ec9278e3
md"""
# Cap 7: Extending the OLG model
* Households not only decide on their saving, but also on their time use. 
* Given a specific time endowment, agents can either work in the market (and earn income), go to school (and acquire human capital for future income generation), or consume leisure.
* Public policy may distort all of these decisions.

"""

# ╔═╡ 4278b868-655e-4cf8-9719-4e3b072e3764
md"""

## 7.1 Accounting for variable labour supply
* The households can decide how many hours to work in each period.
* Leisure demand in each period of the life cycle strongly depends on the respective value of human capital $h_j$, which measure the value of the time endowment in terms of labour market productivity.
* Hence agents may work the same number of hours, but they may be differently productive, so that they earn a different wage per time unit. 
* In order to guarantee that the time endowment is met, we calculate a so-called *shadow wage* $\mu_{j,s}$ (whenever the wage a household earns in the labour market is very small, the household might want to consume more leisure than the actual time endowment).
* The shadow wage is added to the regular wage of the household and calculated such that the household's optimal decision consists in consuming the households's total endowment of time as leisure.
"""

# ╔═╡ 692c88e3-b132-4b79-9e6f-57eaa2600cd0
md"""
### The household decision problem

* Households now have to decide how much to consume and save in each period and how to split up their total time endowment of one unit between working time and leisure  $l_{j,s}$.
* The dynamic budget contraints of the individual entering the labour market in period $t$ now change to:

```math
	p_{t}c_{1,t} = (1-l_{1,t})w_{1,t}^{n} - a_{2,t+1}
```
```math
	p_{t+1}c_{2,t+1} = R_{t+1}^{n} a_{2,t+1} + (1-l_{2,t+1})w_{2,t+1}^{n} - a_{3,t+2}
```

```math
	p_{t+2}c_{3,t+2} = R_{t+2}^{n} a_{3,t+2} + (1-l_{3,t+2})w_{3,t+2}^{n}+pen_{t+2}
```

Where

```math
	w_{j,s}^{n} = [w_s h_j + \mu_{j,s}](1-\tau_s^w-\tau_s^p)
```
defines the net wage of a $j$-year-old worker in period $s=t+j-1$. $w_s$ denotes the gross wage per unit of human capital. Households of different ages may earn more or less *per hour* because of differences in human-capital endowments.

Thus $w_s h_j$ may be interpreted as the individual's gross wage rate.

The human-capital profile $h_1$, $h_2$, $h_3$ is set exogenously and is separate from the general level of wages $w_s$.

If, for a given wage level, leisure demand $l_{j,s}$ exceeds the total time endowment of unity, we compute a shadow wage rate $\mu_{j,s}$ that reduce leisure demand exactly to the time endowment. If $l_{j,s}$ is less than the time endowment, then $\mu_{j,s}$ is zero, *i.e.*

```math
	l_{j,s} ≤ 1, \quad \mu_{j,s} ≥ 0 \quad \text{and} \quad \mu_{j,s}(1-l_{j,s}) = 0
```

* The dynamic budget constraints in the basic model was given by

```math
	p_{t}c_{1,t} = w_{t}^{n} - a_{2,t+1}
```
```math
	p_{t+1}c_{2,t+1} = R_{t+1}^{n} a_{2,t+1} +w_{t+1}^{n} - a_{3,t+2}
```

```math
	p_{t+2}c_{3,t+2} = R_{t+2}^{n} a_{3,t+2} + pen_{t+2}
```

* Finally, aggregate labour supply in period $t$-now measured in units of human capital-is computed from

```math
	L_t = \dfrac{1}{N_t} \Big[ h_1(1-l_{1,t})N_t + h_2(1-l_{2,t})N_{t-1} + h_3(1-l_{3,t})N_{t-2} \Big]
```
"""

# ╔═╡ 3ed1ace5-ce3b-4419-9e66-ccb716d0c3cb
md"""
### Funtional forms and numerical implementation

The utility function in a specific period $s$ is given by

```math
u(c_{j,s}, l_{j,s}) = \dfrac{1}{1-\frac{1}{\gamma}} \Big[ c_{j,s}^{1-\frac{1}{\rho}} + \nu l_{j,s}^{1-\frac{1}{\rho}}\Big]^{\dfrac{1-\frac{1}{\gamma}} { 1- \frac{1}{\rho}}}
```

where $\gamma$ again denotes the intertemporal elasticity of substitution between consumption in different years and $\rho$ defines the intratemporal elasticity of substitution between consumption and leisure. The leisure preference parameter $\nu$ is used to calibrate a realistic fraction of working time. 

From the first-order conditions we derive

```math
	l_{j,s} = \Big( \dfrac{w_{j,s}^n}{\nu p_s}\Big)^{-\rho} c_{j,s} \quad \text{for} \quad j=1,2,3
```

Substituting en the others first order conditions, yields
```math
c_{2,t+1} = \dfrac{\nu_{2,t+1}}{\nu_{1,t}} \Big[ \beta R_{t+1}^n \dfrac{p_t}{p_{t+1}}\Big]^\gamma c_{1,t}
```
```math
c_{3,t+2} = \dfrac{\nu_{3,t+2}}{\nu_{1,t}} \Big[ \beta^2 R_{t+1}^nR_{t+2}^n \dfrac{p_t}{p_{t+2}}\Big]^\gamma c_{1,t}
```

with $\nu_{j,1} = \Big[ 1 + \nu^\rho \Big(  \dfrac{w_{j,s}^n}{p_s} \Big)^{1-\rho} \Big]^{\frac{\rho-\gamma}{1-\rho}}$

FALTA SUSTITUIR EN LA RESTRICCIÓN PRESUPUESTARIA
"""

# ╔═╡ c47014fe-2fd4-4ac7-a141-b802f5ff3c8e
md"""
## 7.2 Human capital and the growth process
* Return to the model with exogenous labour supply.
* Assume that in the first working period households can use a share of their time endowment for education, which increases their human capital in the second period of life. 
* Education can be seen as a private investment, where households invest time they could otherwise have used to generate labour income in order to yield higher labour productivity and therefore income in the future.
* On top of this private return to education, we assume that there may also be social externality associated with the average stock of human capital per worker which increases output and factor returns.
* Since the initial endowment of human capital be constant and equal to unity, the long-run growth rate of the economy is constant and cannot be influenced by public policy. 
* In next step, we model human-capital spillovers from parent cohort to the generation children. Succesive cohorts then start with different initial human-capital endowments. As a results, the growth rate of the economy is endogenous.
"""

# ╔═╡ 96b0b8c0-801d-4b7e-adae-34fcb86aea54
md"""
### 7.2.1 Education investment and externalities

Households have to decide how to split their time endowment in the first period between working time and education investment $e_t$ and how much consume and save in each period. Their productivity in the second period $h_{2,t+1}$ is related to their training decision by

```math
	h_{2,t+1} = h_1\big[ 1 +  \Phi(e_t) \big]
```

where $\Phi(e_t)$ is a decreasing returns to scale function measuring the agent's productivity growth for a given education investment. As already mentioned, we assume $h_1$ to be constant for each new cohort. The budget constraint and household decisions are very similar as in the basic model covered in Chapter 6. The first-period budget constraint now changes to

```math
	p_t c_{1,t} = \big[ 1-(1-\tau_t^s)e_t \big] w_t^n h_1 - a_{2,t+1}
```

where $\tau_t^s$ denotes the subsidy rate of the goverment for education costs which are forgone net wages.

Given the time-separable utility function from Chapter 6 as well as the fact theres is no leisure consumption, we can solve the household optimization problem in two steps.

Firs, we maximize the household's lifetime income by choosing the optimal education investment

```math
	\max_{e_t} \quad W_{1,t} = \big[ 1 -(1-\tau_t^s)e_t \big] w_t^n h_1 + \dfrac{w_{t+1}^n h_1 \big[ 1+ \Phi(e_t) \big] }{R_{t+1}^n} + \dfrac{pen_{t+2}}{R_{t+1}^n R_{t+2}^n}
```

The resulting first-order condition:

```math
	\dfrac{w_{t+1}^n \Phi'(e_t)}{R_{t+1}^n} = (1-\tau_t^s)w_t^n
```

shows that optimal education investment balances the benefits from higher wages in the future against the current cost from forgone wages (net of subsidiers). Note that in steady state, the optimal education decision is independent of net wages. If we specify the education function as 

```math
	\Phi(e_t) = \xi e_{t}^v \quad \xi > 0 \quad \text{and} \quad 0<v<1
```

the latest first order condition yields an optimal education time of

```math
	e_t = \Big( \dfrac{\xi v w_{t+1}^n}{w_t^n(1-\tau_t^s)R_{t+1}^n} \Big)^{\frac{1}{1-v}}
```

The education investment increases with the discounted level of future net wages and decreases with the current net wages (including the subsidy rate) which represents an opportunity cost.

With the optimal education time $e_t$ and the implied lifetime resources $W_{1,t}$, the second step of solving the household problem consists in deriving the optimal path of consumption and saving over the life cycle.
"""

# ╔═╡ 6f0fb77a-f5a8-4b06-a7b3-ca17bf74c42a
md"""
Human capital may affect aggregate output via two distinct channels. On the one hand, it increases aggregare labour supply $L_t$ which is measured in units of human capital. On the other hand, we assume that there may be an additional externality associated with human capital. More specifically, we let the average stock of human capital per worker $H_t$ feature directly in the production funtion, so that

```math
	Y_t = F(K_t, L_t, H_t) = K^\alpha L_t^{1-\alpha}H_t^\epsilon
```

where $\epsilon$ denotes the strength of the externality. Labour input and the average stock of human capital are computed from

```math
	L_t = h_1(1-e_t) + \dfrac{h_{2,t}}{1+n_{p,t}} \quad \text{and} \quad H_t = \dfrac{L_t}{1-e_t+\frac{1}{1+n_{p,t}}}
```

The wage rate and the interest rate then read

```math
	w_t = (1-\alpha) k_t^\alpha H_t^\epsilon \quad \text{and} \quad r_t = \alpha k_l^{\alpha - 1}H_t^\epsilon - \delta
```

Finally, the education subsidy has to be included in the budget constraint of the goverment and equation (6.9) changes to

```math
	T_t + (1+n_{p, t+1})B_{t+1}-B_t = G_t + r_t B_t + \tau_t^s e_t w_t^n h_{1,t}
```
"""

# ╔═╡ eea81ee5-a24a-494a-b874-d0882b6859c6
md"""
### 7.2.3 Human-Capital spillovers and endogenous growth

In the before extension the initial human-capital stock $h_1$ was constant and human capital was only endogenous in the second working-period. 

In this subsection we develop a model variant in which the growth rate of the economy is driven by human-capital accumulation and the respective spillovers from the parent cohort to their children. 

In contrast to the equation 

```math
	h_{2,t+1} = h_1\big[ 1 +  \Phi(e_t) \big]
```

We assume that the second human capital $h_{2, t+1}$ increases by

```math
	h_{2,t+1} = h_{1,t}\big( 1 +  \xi e_{t}^u \big)
```

so that the initial human-capital endowment $h_{1,t}$ is cohort-specific. The most simple way to model the spillover process of human capital from parents to children is to assume a linear relationship

```math
	h_{1,t} = \mu h_{2,t}
```

$\mu >0$ defines the fraction of parental human capital that is acquired by the children generation. Substituting the last two equations we can derive the (endogenous) human capital growth rate $n_{e, t+1}$ as

```math
	h_{1,t+1} =  h_{1,t} \mu \big( 1 +  \xi e_{t}^u \big)= h_{1,t}(1+n_{e,t+1})
```

The way households make their optimal life-cycle choices is very similar to Section 7.2.1. Ther first choose their optimal education investment by maximizing lifetime resources

```math
	\max_{e_t} \quad \widetilde{W}_{1,t} = \big[ 1 -(1-\tau_t^s)e_t \big] w_t^n h_{1,t} + \dfrac{w_{t+1}^n h_{2,t+1} }{R_{t+1}^n} + \dfrac{\widetilde{pen}_{t+2}}{R_{t+1}^n R_{t+2}^n}
```

For every cohort total resources $\widetilde{W}_{1,t}$ grow owing to the rising human-capital stock. 

"""

# ╔═╡ b4d74fff-aa46-4d89-af63-f54432c7f46f
md"""
As before, available resources are spent for private consumption, so that we have

```math
	\widetilde{W}_{1,t} = p_t \widetilde{c}_{1,t} + \dfrac{p_{t+1} \widetilde{c}_{2,t+1}}{R_{t+1}^n} + \dfrac{p_{t+2} \widetilde{c}_{3,t+2}}{R_{t+1}^nR_{t+2}^n}
```

In order to compute a long-run equilibrium of the economy, we normalize all individual variables in thers of the human-capital endowment of a cohort, meaning we calculate

```math
\begin{gather}
	W_{1,t} = \dfrac{\widetilde{W}_{1,t}}{h_{1,t}} = [1-(1-\tau_t^s)e_t]w_t^n +  \dfrac{w_{t+1}^n \frac{(1+n_{e,t+1})}{\mu} }{R_{t+1}^n} + \dfrac{pen_{t+2}}{R_{t+1}^n R_{t+2}^n}\notag\\
	= p_t c_{1,t} + \dfrac{p_{t+1} c_{2,t+1}}{R_{t+1}^n} + \dfrac{p_{t+2} c_{3,t+2}}{R_{t+1}^nR_{t+2}^n}
\end{gather}
```


Note that $W_{1,t}$, $c_{j,1} = \dfrac{\widetilde{c}_{j,s}}{h_{1,t}}$, and $pen_{t+1} = \dfrac{\widetilde{pen}_{t+2}}{h_{1,t}}$ then define total resources, consumption as well as pensions measured in units of initial human capital, respectively. Note that the resources of a household aged $j=2$ need to be computed from

```math
	W_{2,t} = w_t^n \dfrac{(1+n_{e,t})}{\mu} + \dfrac{pen_{t+1}}{R_{t+1}^n}
```

so that one has to account for the increase in resources compared to the previous cohort. 

Of course, we then also have to normalize aggregate variables in a similar fashion, therefore distinguishing between total labour input

```math
	\widetilde{L}_t = (1-e_{t}) h_{1,t} N_t + h_{2,t} N_{t-1}
```

and aggregate labour input in units of a initial human capital and normalized by the size of the young cohort

```math
	L_t = \dfrac{\widetilde{L}_t}{h_{1,t} N_t} = (1-e_t) + \dfrac{1}{\mu (1-n_{p,t})}
```
"""

# ╔═╡ 82933d61-8af8-496c-95a3-36288bdbfa24


# ╔═╡ 6bbf0754-994f-4855-ab5d-a9396b646fd4


# ╔═╡ Cell order:
# ╟─29cebbb7-6986-4d47-96ad-0ccebb993dfc
# ╟─63da7e09-da7d-441a-b967-892d7dcf1f78
# ╟─75ba6402-0628-11ee-02ad-bbd3ec9278e3
# ╟─4278b868-655e-4cf8-9719-4e3b072e3764
# ╟─692c88e3-b132-4b79-9e6f-57eaa2600cd0
# ╟─3ed1ace5-ce3b-4419-9e66-ccb716d0c3cb
# ╟─c47014fe-2fd4-4ac7-a141-b802f5ff3c8e
# ╟─96b0b8c0-801d-4b7e-adae-34fcb86aea54
# ╟─6f0fb77a-f5a8-4b06-a7b3-ca17bf74c42a
# ╟─eea81ee5-a24a-494a-b874-d0882b6859c6
# ╟─b4d74fff-aa46-4d89-af63-f54432c7f46f
# ╠═82933d61-8af8-496c-95a3-36288bdbfa24
# ╠═6bbf0754-994f-4855-ab5d-a9396b646fd4
