import pandas as pd 

## labsh ---> Share of labour compensation in GDP at current national prices
## avh ---> Average annual hours worked by persons engaged
## cn ----> Capital stock at current PPPs (in mil. 2017US$)
## ctfp -->  TFP level at current PPPs (USA=1)
## rtfpna --->	TFP at constant national prices (2017=1)
## csh_c ---> Share of household consumption at current PPPs
## csh_g ---> Share of government consumption at current PPPs
## cgdpo ---> Output-side real GDP at current PPPs (in mil. 2017US$)
## delta ---> Average depreciation rate of the capital stock
## csh_c ---> Share of household consumption at current PPPs
## csh_i ---> Share of gross capital formation at current PPPs
## csh_g ---> Share of government consumption at current PPPs

df = pd.read_excel("pwt1001.xlsx", sheet_name = "Data")

pwt_ids = ["countrycode", "country", "currency_unit", "year"]
pwt_vars = ["labsh", "avh", "cn", "ctfp", "rtfpna", "csh_c", "csh_g", "cgdpo", "cn", "delta", "csh_c", "csh_i", "csh_g"]
df = df[pwt_ids + pwt_vars].query("year==2019")
df = df[df.countrycode.isin(["CHL", "CRI", "MEX"])]

df["tauc"] = [0.19, 0.13, 0.16]
df["np"] = [0.019,0.018,0.022]      

df.to_csv("parametros_olg.csv", index = False)

### Calcula alpha --> capital share in production
# Labour-income share --> (labsh : Share of labour compensation in GDP at current national prices)
alpha = 1 - df["labsh"].values[0] 

### Calcula Omega (normaliza w = 1)
## w = (1 - alpha) * Omega * (K/L)**alpha
## Omega = (w/(1 - alpha))/((K/L)**alpha)

K_Y = (df["cn"]/df["cgdpo"]).values[0]
w = 1 
Omega = (w/(1 - alpha))/((K_Y)**alpha)

### Calcula tasa de depreciación
## I/Y = (np + delta)(K/Y) ---> delta = (I/Y)/(K/Y) - np  
delta = ((df["csh_i"]/(K_Y/5)) - df["np"])**(1/5)


### Calcula nu --> is the parameter which governs the preference of the household for consumption as compared to leissure
"""
The autors set nu equal to 0.335, which leads labour hours to average at approximately one third of the time endowment

This share is derived from assuming a maximum weekly working-time endowment of 110 hours as well as 50 working weeks per year. 
We relate this average annual working hours per employee of around 1800.

nu = 1800/(110*50)
"""

nu = (df["avh"]/(110*50)).values[0]