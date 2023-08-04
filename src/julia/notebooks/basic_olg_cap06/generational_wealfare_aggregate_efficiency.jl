### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 952ec228-325a-11ee-28c3-2584f8e12945
md"""
# Generational Wealfare and Aggregate Efficiency
"""

# ╔═╡ ba9f4fdd-f15c-4ef8-b804-c0181f0e3b03
md"""
Es útil introducir una medida de bienestar que nos permita comparar e interpretar las consecuencias en el bienestar de ciertas políticas de reforma para diferentes cohortes viviendo a lo largo de la trayectoria de transición y en el nuevo equilibrio de largo plazo. Los puntos de referencia son los niveles de utilidad de los tres cohortes que estaban viviendo en el estado estacionario inicial antes de la reforma (*i.e.* en el periodo 0) $U_{j0,}j = 1,2,3$. Llevando a cabo una reforma de política tiene un impacto en diferentes tipos de cohortes. Nos referimos a las actuales generaciones como los que ya estaban haciendo decisiones económicas en el equilibrio inicial y que son "sorprendidas" por el cambio de política durante su ciclo de vida. Sus respectivos nivels de utilidad son $U_{21}$ y $U_{31}$ después de la reforma. Generaciones futuras, por otra parte, son aquellas quienes entran durante la fase de transición. Cómo la utilidad no está escalada cardialmente, la comparación directa de los niveles de utilidad sólo indicarían si un cohorte específico se beneficia o empeora de una reforma específica. Con el objetivo de derivar una medida cuantitativa, se utiliza la propiedad de homogeneidad de la función de utilidad y se deriva el incremento (decremento) relativo en el ingreso que es necesario para alcanzar el nivel de utilidad post-reforma dado los niveles de precios pre-reforma. La medida indica si la reforma ha incrementado el bienestar del hogar considerado en la misma forma que el aumento de los recursos iniciales en $\Delta_t$ por ciento.
De manera que $\Delta_t$ es llamada la variación equivalente en el ingreso (income equivalent variation) o -dado que fue introducido por Hicks- *variación hicksiana equivalente* (*Hicksian equivalent variation*, HEV). La HEV está dada por 

```math
\Delta_t = \Bigg(\dfrac{U_t}{U_0}\Bigg)^{\dfrac{1}{1-\frac{1}{\gamma}}} -1
```

El cambio en bienestar expresado como un porcentaje de los recursos iniciales de un cohorte es el principal indicador de las consecuencias de una reforma en el bienestar intergeneracional.

"""

# ╔═╡ 2c59d525-c51a-4614-b8a8-5a7a07fd8e16
md"""
## Efectos de eficiencia

Nos gustaría saber si las consecuencias en bienestar son debidas unicamente por una redistribución intergeneracional o si la reforma ha mejorado (o empeorado) la asignación de recursos. 

Con el objetivo de evaluar esos efectos de eficiencia, se procede de la siguiente forma:

1. Primero, se calculan las  *lump-sum transfers* (transferencias de suma fija) $\nu_1$ y $\nu_0$  que son necesarias otorgar a las generaciones actuales para que estuvieran tan bien como en el equilibrio inicial.
2. Se calculan las transferencias $\nu_t$, $t=1,\dots,T$ a las generaciones futuras para que todos estén en igualdad de condiciones, *i.e.* todas experimenten el mismo nivel de utilidad $U^*$. El nivel de utilidad $U^*$ está determinado por el requerimiento que los valores presentes de todas las transferencias sumen cero. Los efectos de este *esquema de compensación* cueden ser vistas en la figura.

En esta situación, todas las generaciones actuales recibirían las transferencias de suma fija mientras que las cohortes futuras habrían de pagar impuestos de suma fija, de manera que 

```math
	\Delta^* = \Bigg(\dfrac{U^*}{U_0}\Bigg)^{\dfrac{1}{1-\frac{1}{\gamma}}} -1 > 0.
```

$\Delta^*$ puede ser interpretada como una medida de la ganancia de eficiencia agregada de la reforma de política como porcentaje de los recursos iniciales. Decimos que una reforma es *efficiency improving* o *Pareto superior después de compensación* si $\Delta^* >0$.

"""

# ╔═╡ 2052005b-81a1-4302-af20-b5210ba14834
md"""
Derivar $U^*$ y por lo tanto $\Delta*$ no es tarea fácil. Se calcula incorporando un nuevo agente denominado *Lump-Sum Redistribution Authority* (LSRA). Se permite a esta autoridad pagar transferencias de suma fija a todas las generaciones, actuales y futuras. El LSRA sigue los siguientes tres pasos:

1. Calcula las transferencias requeridas que llevan a las cohortes de mayor edad en el año de la reforma desde su nivel de utilidad actual a sus utilidades de equilibrio iniciales $U_{i0,}i = 2,3$. 

```math
\nu_{-1} = \Bigg[ \Bigg( \dfrac{U_{30}}{U_{31}} \Bigg)^{\frac{1}{1-\frac{1}{\gamma}}} -1\Bigg] W_{31} \qquad \text{y} \qquad \nu_{0} = \Bigg[ \Bigg( \dfrac{U_{20}}{U_{21}} \Bigg)^{\frac{1}{1-\frac{1}{\gamma}}} -1\Bigg] W_{21} 
```

donde $W_{21}$ y $W_{31}$ define los recursos de los hogares de mayor edad *después* de la reforma en el periodo 1.

"""

# ╔═╡ 5bc33213-3b84-427e-9224-40ba9115f7d3
md"""

2. Se calculan las transferencias a las cohortes futuras mediante el presupuesto intertemporal de la agencia de redistribución. Este presupuesto garantiza que la suma descontada de las transferencias es exactamente iguale a cero, *i.e.* el LSRA no obtiene ninguna ganancia. En consecuencia, a partir del periodo 1

3. Dado $U^*$ las compensaciones transferencias de suma fija para las cohortes futuras pueden ser calculadas

```math
\nu_t = \Bigg[ \Bigg( \dfrac{U^*}{U_t} \Bigg)^{\frac{1}{1-\frac{1}{\gamma}}} -1 \Bigg]W_t \qquad \text{para} \qquad t=1,\dots,T
```

Los pagos de compensación $\nu_t$ tienen que ser incluidos en las respectivas restricciones presupuestales individuales cuando se calculan los valores de ahorro y consumo. Además, se asume que el LSRA es un actor regular del modelo, *i.e.* tenemos que tomar en cuenta que la suma inicial de transferencias $\nu_{-1}N_{-1}+\nu_0N_0+\nu_1N_1$ necesita ser financiado en el mercado de capital. Consecuentemente, la deuda (o activos) de la agencia, per cápita a los trabajadores en el periodo 2, está dada por

```math
B_2^a = \dfrac{1}{1 + n_{p,2}} \Bigg[ \dfrac{\nu_{-1}}{(1+n_{p,0}) (1+n_{p,1})} + \dfrac{\nu_0}{1+n_{p,1}} + \nu_1  \Bigg]
```

En todos los periodos futuros $t= 2,\dots, T$ la agencia tiene que pagar (o recibir) intereses y financiar transferencias, de manera que la restricción presupuestaria periódica es

```math
(1+n_{p,t+1})B_{t+1}^a = (1+r_t)B_t^a + \nu_t
```

Con el objetivo de contabilizar los activos y pasivos del LSRA, cambiamos la condición de equilibrio del mercado de capitales 

```math
A_t = K_t + B_t + B_t^a \qquad \text{para} \qquad t=2,\dots,T
```

"""

# ╔═╡ Cell order:
# ╟─952ec228-325a-11ee-28c3-2584f8e12945
# ╟─ba9f4fdd-f15c-4ef8-b804-c0181f0e3b03
# ╟─2c59d525-c51a-4614-b8a8-5a7a07fd8e16
# ╟─2052005b-81a1-4302-af20-b5210ba14834
# ╟─5bc33213-3b84-427e-9224-40ba9115f7d3
