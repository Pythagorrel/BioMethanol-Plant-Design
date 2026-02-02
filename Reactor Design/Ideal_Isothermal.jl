using DifferentialEquations
using PyPlot
using Printf

# ==============================================================================
# 1. IMPORT THERMODYNAMICS
# ==============================================================================
include("BioMethEQK.jl")
using .BioMethEQK

# ==============================================================================
# 2. OPERATING PARAMETERS (Base Case - Stream 5)
# ==============================================================================
const T_op = 523.15      # 250 C
const P_op = 50.0        # 50 bar
const R_gas = 8.314

# FEED (Stream 5)
const F_CO2_in = 514.0   # mol/s
const F_H2_in = 1542.0  # mol/s 
const F_CO_in = 0.0     # Minimal CO in fresh feed
const F_H2O_in = 0.0
const F_MeOH_in = 0.0

# ==============================================================================
# 3. KINETICS (Hybrid: Van-Dal Rate + BioMethEQK Thermo)
# ==============================================================================
function get_kinetics_hybrid(T)
    # Rate Constants (Van-Dal 2013)
    k1 = 1.07 * exp(40000.0 / (R_gas * T))
    k_rwgs = 1.22e10 * exp(-98084.0 / (R_gas * T))

    # Adsorption (Van-Dal 2013)
    K_water_h2 = 3453.38
    K_h2_sqrt = 0.499 * exp(17197.0 / (R_gas * T))
    K_water = 6.62e-11 * exp(124119.0 / (R_gas * T))

    # Equilibrium (Rigorous Gibbs)
    K_eq1 = EQConst1(T)
    K_eq3 = EQConst3(T)

    return k1, k_rwgs, K_water_h2, K_h2_sqrt, K_water, K_eq1, K_eq3
end

# ==============================================================================
# 4. SOLVER
# ==============================================================================
function pbr_integrated!(dF, F, p, W)
    F_CO2, F_H2, F_CO, F_MeOH, F_H2O = F[1], F[2], F[3], F[4], F[5]
    F_total = sum(F)

    p_CO2 = P_op * (F_CO2 / F_total)
    p_H2 = P_op * (F_H2 / F_total)
    p_CO = P_op * (F_CO / F_total)
    p_MeOH = P_op * (F_MeOH / F_total)
    p_H2O = P_op * (F_H2O / F_total)

    k1, k_rwgs, K_wh, K_h, K_w, K_eq1, K_eq3 = get_kinetics_hybrid(T_op)

    p_H2_safe = max(p_H2, 1e-5)
    DEN_base = 1.0 + (K_wh * p_H2O / p_H2_safe) + (K_h * sqrt(p_H2_safe)) + (K_w * p_H2O)

    # Rates
    drive_1 = (p_CO2 * p_H2_safe) * (1.0 - (1.0 / K_eq1) * (p_MeOH * p_H2O) / (p_CO2 * p_H2_safe^3))
    r1 = (k1 * drive_1) / (DEN_base^3)

    drive_2 = p_CO2 * (1.0 - (1.0 / K_eq3) * (p_CO * p_H2O) / (p_CO2 * p_H2_safe))
    r2 = (k_rwgs * drive_2) / (DEN_base)

    dF[1] = -r1 - r2
    dF[2] = -3 * r1 - r2
    dF[3] = r2
    dF[4] = r1
    dF[5] = r1 + r2
end

# ==============================================================================
# 5. EXECUTION & YIELD ANALYSIS
# ==============================================================================
println("----------------------------------------------------------")
println("Biomethanol Reactor Design (Yield vs Conversion Analysis)")
println("----------------------------------------------------------")

u0 = [F_CO2_in, F_H2_in, F_CO_in, F_MeOH_in, F_H2O_in]
W_span = (0.0, 50000.0)
prob = ODEProblem(pbr_integrated!, u0, W_span)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

# Extract Data
W_plot = sol.t
# 1. Conversion X = (In - Out) / In
X_CO2 = [(F_CO2_in - u[1]) / F_CO2_in for u in sol.u]

# 2. Yield Y = Methanol_Out / CO2_In
Yield_MeOH = [u[4] / F_CO2_in for u in sol.u]

# 3. Production Rate (kmol/hr)
MeOH_flow_kmol = [u[4] * 3600 / 1000 for u in sol.u]

# Find Point for Target PRODUCTION (Which implies Target Yield)
target_prod = 394.1
idx = findfirst(x -> x >= target_prod, MeOH_flow_kmol)

if !isnothing(idx)
    W_req = W_plot[idx]
    Actual_Prod = MeOH_flow_kmol[idx]

    # Calculate Tube Count (Van-Dal Table 2 Data)
    rho_bulk = 1775.0 * (1.0 - 0.4)
    cat_per_tube = (pi * (0.02)^2 * 6.0) * rho_bulk
    n_tubes = ceil(Int, W_req / cat_per_tube)

    println("[SUCCESS] Target Production Reached!")
    @printf("Methanol Production:    %.2f kmol/hr\n", Actual_Prod)
    @printf("Required Catalyst Mass: %.2f kg\n", W_req)
    @printf("Tube Count:             %d\n", n_tubes)
    println("----------------------------------------")
    @printf("METHANOL YIELD:         %.2f %% (Target ~21.3%%)\n", Yield_MeOH[idx] * 100)
    @printf("CO2 CONVERSION:         %.2f %%\n", X_CO2[idx] * 100)
    @printf("Selectivity:            %.2f %%\n", (Yield_MeOH[idx] / X_CO2[idx]) * 100)
    println("----------------------------------------")
else
    println("[WARNING] Target production not reached (Equilibrium Limit).")
    max_prod = maximum(MeOH_flow_kmol)
    max_yield = maximum(Yield_MeOH)
    max_conv = maximum(X_CO2)

    @printf("Max Production: %.2f kmol/hr\n", max_prod)
    @printf("Max Yield:      %.2f %%\n", max_yield * 100)
    @printf("Max Conversion: %.2f %%\n", max_conv * 100)
end

# Plot Yield vs Conversion
figure(figsize=(8, 6))
plot(W_plot, X_CO2, label="Total CO2 Conversion", color="blue", linewidth=2, linestyle="--")
plot(W_plot, Yield_MeOH, label="Methanol Yield", color="green", linewidth=2)
axhline(y=target_prod / (F_CO2_in * 3.6 / 1000), color="red", linestyle=":", label="Target Yield")

title("Reactor Performance: Yield vs Conversion", fontsize=14)
xlabel("Catalyst Weight (kg)", fontsize=12)
ylabel("Fraction (-)", fontsize=12)
grid(true, linestyle="--", alpha=0.6)
legend(loc="right")
tight_layout()
display(gcf())