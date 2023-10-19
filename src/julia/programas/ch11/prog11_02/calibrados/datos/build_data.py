import pandas as pd 

year_analysis = 2015

iso_map = {"Chile" : "CHL", "Mexico" : "MEX", "Costa Rica" : "CRI"}

"""
+++++++++++++++++++++++++++++
    Data to calibrate
+++++++++++++++++++++++++++++

(In % of GDP if not indicated otherwise)

Capital market:
    * Interest rate (in % p.a)

Goods market:
    * Private Consumption (csh_c variable in Penn World Table)
    * Public Consumption (csh_g variable in Penn World Table)
    * Investment (csh_i in Penn World Table)

Tax rates:
    * Consumption (in %) (Source : OECD Tax Database-https://www.oecd.org/tax/tax-policy/tax-database/)
    * Labour earnings (in %) (Source : OECD Tax Database-https://www.oecd.org/tax/tax-policy/tax-database/)
    * Capital income (in %) (Source : OECD Tax Database-https://www.oecd.org/tax/tax-policy/tax-database/)

Labour market:
    * Average hours worked (in % of time endowment) (avh in Penn World Table)

Pension system:
    * Total pension payments:
        - (OECD-Founded Pension Indicators-Contributions as a % of GDP)
        - (OECD-Pensions at Glance-Public expenditure on pensions as a % of GDP)
Tax revenue:
    * Consumption tax (OECD-Revenue Statistics-5000 Taxes on good and services as % GDP)
    * Labour-earnings tax (OECD-Revenue Statistics-1000 Taxes on income, profits and capital gains as % GDP)    
    * Capital-income tax (OECD-Revenue Statistics-1000 Taxes on income, profits and capital gains as % GDP)    

+++++++++++++++++++++++++++++
    Parameters
+++++++++++++++++++++++++++++

Calibrated:
    * alpha : Capital intensity elasticity
    * Omega : Being the technology level which is constant over time (Yₜ = ΩKₜᵅLₜ⁽¹⁻ᵅ))
    * delta : Depreciation rate
    * nu : Parameter which governs the preference of the household for consumption as compared to leissure

Historical:
    * np : 
        + Description : Population growth rate
        + Source : 
    * gy : 
        + Public spending (equals a fixed fraction of production in the initial steady state and remains constant throughout the transition).
        + Source : csh_g variable in Penn World Table
    * by : 
        + Description : Public debt rate
        + Source : 
    * tauc : 
        + Description : Consumption tax rate
        + Source : 
    * tauw : 
        + Description : Labour income tax rate
        + Source :
    * taur : 
        + Description : Capital tax rate
        + Source : 
    * kappa :
        + Description : Pension replacement rates
        + Source : (OECD-Pension at a Glance-Pension replacement rates)
"""


#### PEN WORLD TABLE SOURCE

"""
Indicators in PWT:
    + labsh ---> Share of labour compensation in GDP at current national prices
    + avh ---> Average annual hours worked by persons engaged
    + cn ----> Capital stock at current PPPs (in mil. 2017US$)
    + ctfp -->  TFP level at current PPPs (USA=1)
    + rtfpna --->	TFP at constant national prices (2017=1)
    + csh_c ---> Share of household consumption at current PPPs
    + csh_g ---> Share of government consumption at current PPPs
    + cgdpo ---> Output-side real GDP at current PPPs (in mil. 2017US$)
    + delta ---> Average depreciation rate of the capital stock
    + csh_c ---> Share of household consumption at current PPPs
    + csh_i ---> Share of gross capital formation at current PPPs
    + csh_g ---> Share of government consumption at current PPPs
"""

pwt = pd.read_excel("pwt/pwt1001.xlsx", sheet_name = "Data")

pwt_ids = ["countrycode", "country", "currency_unit", "year"]
pwt_vars = ["labsh", "avh", "cn", "ctfp", "rtfpna", "csh_c", "csh_g", "cgdpo", "delta", "csh_c", "csh_i", "csh_g"]
pwt = pwt[pwt_ids + pwt_vars].query(f"year=={year_analysis}")
pwt = pwt[pwt.countrycode.isin(["CHL", "CRI", "MEX"])]

#### OECD Statistics

"""
OECD-Founded Pension Indicators:
    + Pension system-Total pension payments
        - (OECD-Founded Pension Indicators-Contributions as a % of GDP)
        - (OECD-Pensions at Glance-Public expenditure on pensions as a % of GDP)
"""

oecd_fpi = pd.read_csv("oecd/founded_pensions_indicators/contribution_as_percentage_gdp/PNNI_NEW_08102023231430986.csv") 
oecd_fpi = oecd_fpi[oecd_fpi.Country.isin(['Chile', "Costa Rica", "Mexico"])].query(f"Year == {year_analysis}")[["Country","Year","Value"]]
oecd_fpi["Country"] = oecd_fpi["Country"].replace(iso_map)
oecd_fpi = oecd_fpi[["Country","Value"]].rename(columns = {"Country" : "countrycode", "Value" : "cont_pension_perc_gdp"})
oecd_fpi["cont_pension_perc_gdp"] /=100 

oecd_pep = pd.read_csv("oecd/founded_pensions_indicators/public_expenditure_on_pensions/public_expenditure_pensions.csv")
oecd_pep = oecd_pep[oecd_pep.Country.isin(['Chile', "Costa Rica", "Mexico"])][["Country", "2017"]].rename(columns = {"Country" : "countrycode", "2017" : "public_expenditure_pensions_perc_gdp"})
oecd_pep["public_expenditure_pensions_perc_gdp"] /= 100
oecd_pep["countrycode"] = oecd_pep["countrycode"].replace(iso_map)

"""
OECD-Revenue Statistics:
    * Consumption tax (OECD-Revenue Statistics-5000 Taxes on good and services as % GDP)
    * Labour-earnings tax (OECD-Revenue Statistics-1000 Taxes on income, profits and capital gains as % GDP)    
    * Capital-income tax (OECD-Revenue Statistics-1000 Taxes on income, profits and capital gains as % GDP)    
"""
oecd_rs = pd.read_csv("oecd/taxes_as_percentage_gdp/REV_08102023224830389.csv")
oecd_rs = oecd_rs[oecd_rs.Country.isin(['Chile', "Costa Rica", "Mexico"])]

oecd_rs_gpd = oecd_rs[(oecd_rs["Tax revenue"]=='1000 Taxes on income, profits and capital gains') | (oecd_rs["Tax revenue"]=='5000 Taxes on goods and services')].query(f"GOV == 'FED' and VAR =='TAXGDP' and Year =={year_analysis}")
oecd_rs_total_tax = oecd_rs[(oecd_rs["Tax revenue"]=='1000 Taxes on income, profits and capital gains') | (oecd_rs["Tax revenue"]=='5000 Taxes on goods and services')].query(f"GOV == 'FED' and VAR =='TAXPER' and Year =={year_analysis}")

oecd_rs_gpd = oecd_rs_gpd[["COU", "TAX","Value"]]
oecd_rs_total_tax = oecd_rs_total_tax[["COU", "TAX", "Value"]]

oecd_rs_gpd["TAX"] = oecd_rs_gpd["TAX"].replace({"1000" : "consump_tax_rev_per_gdp" , 5000 : "income_tax_rev_per_gdp"})
oecd_rs_total_tax["TAX"] = oecd_rs_total_tax["TAX"].replace({"1000" : "consump_tax_rev_per_total_tax" , 5000 : "income_tax_rev_per_total_tax"})

oecd_rs_gpd["Value"] /= 100
oecd_rs_total_tax["Value"] /= 100

oecd_rs_gpd = oecd_rs_gpd.pivot(index='COU', columns='TAX', values='Value').reset_index()
oecd_rs_total_tax = oecd_rs_total_tax.pivot(index='COU', columns='TAX', values='Value').reset_index()

"""
OECD-Pension replacement rates:
    + Pension replacement rates : Pension replacement rates (OECD-Pension at a Glance-Pension replacement rates)
    NOTE: Average Worker wage (AW)
"""

oecd_prr = pd.read_csv("oecd/pensions_at_a_glance/PAG_08102023230102857.csv")
oecd_prr = oecd_prr[oecd_prr.Country.isin(['Chile', "Costa Rica", "Mexico"])].query(f"Indicator == 'Net pension replacement rate, Male, 1.50 of AW'")
oecd_prr = oecd_prr[["COU","Value"]]
oecd_prr = oecd_prr.rename(columns = {"Value" : "kappa"})
oecd_prr["kappa"] /= 100

"""
OECD-Tax Database:
    + pop_growth_rate
    + fertility_rate
    + vat_tax
    + w_tax

Impuestos sobre los salarios :
https://www.oecd.org/tax/tax-policy/Impuestos-sobre-los-salarios-costa-rica.pdf
https://www.oecd.org/tax/tax-policy/Impuestos-sobre-los-salarios-chile.pdf
https://www.oecd.org/tax/tax-policy/Impuestos-sobre-los-salarios-mexico.pdf

"""

oecd_tax_db = pd.read_csv("tax_data/output/tax_data.csv")

#### Imputamos la tasa de impuesto sobre el VA de Costa Rica constante en todo el periodo
oecd_tax_db.loc[oecd_tax_db.countrycode=='CRI',"vat_tax"] = 13

"""
#### Ajustamos datos de impuestos laborales
    + MEX aplica las siguientes tasas de impuestos:
        - [1.92, 6.4, 10.88, 16.0, 17.92, 21.36, 23.52, 30.0, 32.0, 34.0, 35.0]
    + CHI aplica las siguientes tasas de impuestos:
        - [0.0, 4.0, 8.0, 13.5, 23.0, 30.4, 35.0, 40.0]
    + CRI aplica una tasa impositiva neta media de los trabajadores es una medida del impuesto neto sobre la renta salarial pagado directamente por el trabajador.
    En Costa Rica, en 2021, el trabajador soltero promedio estuvo sujeto a una tasa impositiva neta media de 10.5%, en comparación con la media de la OCDE de 24.6%.

Se construyen tres tasas de impuestos sobre el ingreso: 
    1) Tasa mínima de impuesto (corresponde al threshold más bajo de los datos de OCDE)
    2) Tasa media de impuesto (corresponde al dato de las infografías de impuestos sobre los salarios de OCDE)
    3) Tasa máxima impuesto (corresponde al threshold más alto de los datos de OCDE).

NOTA: Dado que CRI no tiene dato de thresholds, usamos el valor máximo y mínimo de impuesto sobre el ingreso de OCDE 

"""
oecd_tax_db["tax_w_ocde_min"] = 0.0
oecd_tax_db["tax_w_ocde_max"] = 0.0
oecd_tax_db["tax_w_ocde_medio"] = 0.0

oecd_tax_db.loc[oecd_tax_db.countrycode=='MEX', "tax_w_ocde_min"] = 0.0192
oecd_tax_db.loc[oecd_tax_db.countrycode=='CHL', "tax_w_ocde_min"] = 0.0
oecd_tax_db.loc[oecd_tax_db.countrycode=='CRI', "tax_w_ocde_min"] = 0.0

oecd_tax_db.loc[oecd_tax_db.countrycode=='MEX', "tax_w_ocde_max"] = 0.35
oecd_tax_db.loc[oecd_tax_db.countrycode=='CHL', "tax_w_ocde_max"] = 0.40
oecd_tax_db.loc[oecd_tax_db.countrycode=='CRI', "tax_w_ocde_max"] = 0.246

oecd_tax_db.loc[oecd_tax_db.countrycode=='MEX', "tax_w_ocde_medio"] = 0.102
oecd_tax_db.loc[oecd_tax_db.countrycode=='CHL', "tax_w_ocde_medio"] = 0.070
oecd_tax_db.loc[oecd_tax_db.countrycode=='CRI', "tax_w_ocde_medio"] = 0.105

oecd_tax_db = oecd_tax_db.query(f"year=={year_analysis}")
oecd_tax_db = oecd_tax_db[["countrycode", "pop_growth_rate", "fertility_rate", "vat_tax", "tax_w_ocde_min", "tax_w_ocde_max", "tax_w_ocde_medio"]]
oecd_tax_db.pop_growth_rate /= 100.0
oecd_tax_db.fertility_rate /= 100.0
oecd_tax_db.vat_tax = oecd_tax_db.vat_tax.astype(float)/100

#### Join data
df = pd.concat([pwt.set_index("countrycode"),
                oecd_fpi.set_index("countrycode"),
                oecd_rs_gpd.set_index("COU"), 
                oecd_rs_total_tax.set_index("COU"), 
                oecd_prr.set_index("COU"), 
                oecd_pep.set_index("countrycode"),
                oecd_tax_db.set_index("countrycode")], axis = 1)

df.index.name = 'countrycode'       
df = df.reset_index()         
   
df.to_csv("parametros_olg.csv", index = False)

