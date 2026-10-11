rm(list=ls())

# 0. INPUT DATA

# Utilisations [-]
U_glc <- 0.67          # glucose [-]
U_NH3 <- 0.89          # ammonia [-]
U_O2  <- 0.12          # oxygen [-]

# Prices [EUR/kg]
p_glc <- 0.19          # glucose [EUR/kg]
p_AS  <- 1.09          # ammonium sulphate [EUR/kg]
p_CA  <- 9.54          # coumaric acid [EUR/kg]

# Yields on glucose consumed [g product / g glucose]  (mg/g / 1000)
Y_X    <- 74  / 1000   # biomass [g/g]
Y_EtOH <- 71  / 1000   # ethanol [g/g]
Y_Gly  <- 122 / 1000   # glycerol [g/g]
Y_CA   <- 69  / 1000   # coumaric acid [g/g]

# Biomass composition: CH1.78 O0.66 N0.158 
h_X <- 1.78; o_X <- 0.66; n_X <- 0.158 #; p_X <- 0.009; s_X <- 0.0035; k_X <- 0.0015;   # H, O, N, P, S, K per C-mol [mol/C-mol]

# Target production [kg/yr]
P_CA <- 4950           # coumaric acid [kg/yr]

# Molecular weights [g/mol]
MW_glc  <- 180.16    # glucose [g/mol]
MW_O2   <- 31.998    # oxygen [g/mol]
MW_X    <- 27.027    # biomass unit (C-mol) [g/mol]
MW_NH3  <- 17.03     # ammonia [g/mol]
MW_EtOH <- 46.07     # ethanol [g/mol]
MW_Gly  <- 92.09     # glycerol [g/mol]
MW_CA   <- 164.16    # coumaric acid [g/mol]
MW_AS   <- 132.14    # ammonium sulphate (NH4)2SO4 [g/mol] (not given; standard value)

# Cooling data
Ws_kW   <- 10.7        # impeller power [kW] = [kJ/s]
dh_vap  <- 2430.7      # latent heat of water [kJ/kg]
m_vap   <- 127000      # water evaporated [kg/yr]
t_year  <- 8760 * 3600 # seconds per year [s/yr]

# Standard heats of combustion [kJ/mol] (negative by convention)
dhc <- c(glc = -2805, X = -582.16158, NH3 = -382.6,
         EtOH = -1366.8, Gly = -1655.4, CA = -4152)




# 1. STOICHIOMETRY PER MOL OF GLUCOSE CONSUMED  (basis: 1 mol glc)

# Moles of each product per mol glucose [mol/mol glc]
c_X  <- Y_X    * MW_glc / MW_X      # biomass [C-mol/mol glc]
f1   <- Y_EtOH * MW_glc / MW_EtOH   # ethanol [mol/mol glc]
f2   <- Y_Gly  * MW_glc / MW_Gly    # glycerol [mol/mol glc]
f3   <- Y_CA   * MW_glc / MW_CA     # coumaric acid [mol/mol glc]

# Elemental balances  (glucose C6H12O6; ethanol C2H6O; glycerol C3H8O3; coumaric acid C9H8O3; biomass CH1.78O0.66N0.158)
b <- n_X * c_X                                         # N balance -> NH3 [mol/mol glc]
d <- 6 - c_X - 2*f1 - 3*f2 - 9*f3                      # C balance -> CO2 [mol/mol glc]
e <- (12 + 3*b - h_X*c_X - 6*f1 - 8*f2 - 8*f3) / 2     # H balance -> H2O [mol/mol glc]
a <- (o_X*c_X + 2*d + e + f1 + 3*f2 + 3*f3 - 6) / 2    # O balance -> O2 [mol/mol glc]

stoich <- data.frame(
  species = c("O2 (a)", "NH3 (b)", "Biomass (c)", "CO2 (d)", "H2O (e)",
              "Ethanol (f1)", "Glycerol (f2)", "Coumaric acid (f3)"),
  mol_per_mol_glc = c(a, b, c_X, d, e, f1, f2, f3))
cat("\n Stoichiometric coefficients [mol / mol glucose consumed] \n")
print(stoich, digits = 5, row.names = FALSE)

# 2. GLUCOSE (Q3) - needed first because everything scales with it

m_glc_cons <- P_CA / Y_CA                       # glucose consumed [kg/yr]
n_glc_cons <- m_glc_cons * 1000 / MW_glc        # = SMR, glucose reacted [mol/yr]
m_glc_feed <- m_glc_cons / U_glc                # glucose fed [kg/yr]
n_glc_feed <- m_glc_feed * 1000 / MW_glc        # glucose fed [mol/yr]

# 3. OXYGEN (Q1)

n_O2_stoich <- a * n_glc_cons                   # O2 consumed [mol/yr]
m_O2_stoich <- n_O2_stoich * MW_O2 / 1000       # O2 consumed [kg/yr]
m_O2_feed   <- m_O2_stoich / U_O2               # O2 supplied [kg/yr]
n_O2_feed   <- m_O2_feed * 1000 / MW_O2         # O2 supplied [mol/yr]


# 4. AMMONIA (Q2)
  
n_NH3_stoich <- b * n_glc_cons                  # NH3 consumed [mol/yr]
m_NH3_stoich <- n_NH3_stoich * MW_NH3 / 1000    # NH3 consumed [kg/yr]
m_NH3_feed   <- m_NH3_stoich / U_NH3            # NH3 fed [kg/yr]
n_NH3_feed   <- m_NH3_feed * 1000 / MW_NH3      # NH3 fed [mol/yr]

# Equivalent ammonium sulphate: 1 mol (NH4)2SO4 supplies 2 mol N
n_AS <- n_NH3_feed / 2                          # [mol/yr]
m_AS <- n_AS * MW_AS / 1000                     # [kg/yr]


# 5. OTHER STREAMS (for the process diagram)  [mol/yr] and [kg/yr]

streams <- data.frame(
  stream = c("Biomass", "Ethanol", "Glycerol", "Coumaric acid", "CO2", "H2O produced",
             "Unreacted glucose", "Unreacted NH3", "Unused O2"),
  mol_per_yr = c(c_X*n_glc_cons, f1*n_glc_cons, f2*n_glc_cons, f3*n_glc_cons,
                 d*n_glc_cons, e*n_glc_cons,
                 n_glc_feed - n_glc_cons, n_NH3_feed - n_NH3_stoich,
                 n_O2_feed - n_O2_stoich))
MWs <- c(MW_X, MW_EtOH, MW_Gly, MW_CA, 44.01, 18.015, MW_glc, MW_NH3, MW_O2)
streams$kg_per_yr <- streams$mol_per_yr * MWs / 1000
streams$kg_per_h  <- streams$kg_per_yr / 8760
cat("\nOutlet streams\n")
print(streams, digits = 6, row.names = FALSE)

# 6. ECONOMICS (Q4, Q5)

economics <- function(U_glc = 0.67, Y_CA = 0.069, U_NH3 = 0.89) {
  m_g  <- P_CA / Y_CA                           # glucose consumed [kg/yr]
  n_g  <- m_g * 1000 / MW_glc                   # [mol/yr]
  m_gf <- m_g / U_glc                           # glucose fed [kg/yr]
  m_nf <- b * n_g * MW_NH3 / 1000 / U_NH3       # NH3 fed [kg/yr]
  m_as <- (m_nf * 1000 / MW_NH3) / 2 * MW_AS / 1000   # (NH4)2SO4 [kg/yr]
  cost_g  <- m_gf * p_glc                       # [EUR/yr]
  cost_as <- m_as * p_AS                        # [EUR/yr]
  rev_CA  <- P_CA * p_CA                        # [EUR/yr]
  c(glucose_cost = cost_g, AS_cost = cost_as, revenue = rev_CA,
    margin = rev_CA - cost_g - cost_as)         # [EUR/yr]
}
base <- economics()
cat("\nEconomics, base case [EUR/yr]\n"); print(round(base, 0))
cat(sprintf("Break-even CA price: %.2f EUR/kg\n",
            (base["glucose_cost"] + base["AS_cost"]) / P_CA))

# Sensitivity (Q5)
sens <- data.frame(
  case = c("Base", "Glucose util. 50%", "Glucose util. 60%", "Glucose util. 80%",
           "CA yield -20%", "CA yield +20%", "NH3 util. 70%"),
  margin_EUR_per_yr = c(
    base["margin"],
    economics(U_glc = 0.50)["margin"], economics(U_glc = 0.60)["margin"],
    economics(U_glc = 0.80)["margin"],
    economics(Y_CA = 0.069*0.8)["margin"], economics(Y_CA = 0.069*1.2)["margin"],
    economics(U_NH3 = 0.70)["margin"]))
cat("\nSensitivity of gross margin\n")
print(sens, digits = 6, row.names = FALSE)

# 7. COOLING DUTY (Q6-Q10)

# Q6: dH_rxn = sum(reactants * dhc) - sum(products * dhc)   [kJ/mol glc]
# (O2, CO2 and H2O have dhc = 0)
dH_rxn <- (dhc["glc"]*1 + dhc["NH3"]*b) -
  (dhc["X"]*c_X + dhc["EtOH"]*f1 + dhc["Gly"]*f2 + dhc["CA"]*f3)
cat(sprintf("\nQ6: dH_rxn = %.1f kJ/mol glucose (negative = exothermic)\n", dH_rxn))

# Q7: annual reaction heat [MJ/yr] = dH_rxn [kJ/mol] * SMR [mol/yr] / 1000
H_rxn <- dH_rxn * n_glc_cons / 1000
cat(sprintf("Q7: dH_rxn,annual = %.0f MJ/yr (heat released)\n", H_rxn))
# Alternative if the marker uses glucose FED instead of glucose reacted:
H_rxn_fed <- dH_rxn * n_glc_feed / 1000
cat(sprintf("(fed-glucose basis would give %.0f MJ/yr)\n", H_rxn_fed))

# Q8: shaft work [MJ/yr] = kW * s/yr / 1000
Ws <- Ws_kW * t_year / 1000
cat(sprintf("Q8: Ws = %.0f MJ/yr (heat added)\n", Ws))

# Q9: evaporative heat removal [MJ/yr] = kg/yr * kJ/kg / 1000
H_evap <- m_vap * dh_vap / 1000
cat(sprintf("Q9: m_vap*dh_vap = %.0f MJ/yr (heat removed)\n", H_evap))

# Q10: net cooling duty [MJ/yr]  Q = -dH_rxn + Ws - m_vap*dh_vap
Q <- -H_rxn + Ws - H_evap
Q_noevap <- -H_rxn + Ws
cat(sprintf("Q10: Q = %.0f MJ/yr (positive = cooling required)\n", Q))
cat(sprintf("Without evaporation: %.0f MJ/yr; evaporation reduces duty by %.1f %%\n",
            Q_noevap, 100 * H_evap / Q_noevap))
cat(sprintf("Share of heat load: reaction %.1f %%, shaft work %.1f %%\n",
            100 * (-H_rxn) / Q_noevap, 100 * Ws / Q_noevap))

# 8. SUMMARY TABLE

summary_tbl <- data.frame(
  item = c("O2 supplied", "NH3 fed", "Glucose fed",
           "Glucose reacted", "Glucose cost", "(NH4)2SO4 cost", "CA revenue",
           "Gross margin", "Reaction heat", "Shaft work", "Evaporation", "Net cooling Q"),
  value = c(m_O2_feed, m_NH3_feed, m_glc_feed, m_glc_cons,
            base["glucose_cost"], base["AS_cost"], base["revenue"], base["margin"],
            H_rxn, Ws, -H_evap, Q),
  unit = c("kg/yr", "kg/yr", "kg/yr", "kg/yr", "EUR/yr", "EUR/yr", "EUR/yr",
           "EUR/yr", "MJ/yr", "MJ/yr", "MJ/yr", "MJ/yr"))
cat("\nSUMMARY\n"); print(summary_tbl, digits = 7, row.names = FALSE)
