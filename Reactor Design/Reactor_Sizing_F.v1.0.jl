using DifferentialEquations
using Printf
using Statistics

# ==============================================================================
# 1. IMPORT MODULES
# ==============================================================================
include("BioMethEQK.jl")
include("FugacityCo.jl") 
using .BioMethEQK
using .FugacityCo

# ==============================================================================
# 2. CALIBRATED CATALYST FACTORS (Determined via Reconciliation)
# ==============================================================================
# These account for pore diffusion and specific catalyst selectivity
const ETA_MEOH = 0.3143  
const ETA_RWGS = 0.1960  

# ==============================================================================
# 3. OPERATING PARAMETERS (SI UNITS)
# ==============================================================================
const T_op = 523.15      # 250 C
const P_in = 50.0        # 50 bar Inlet
const R_gas = 8.314      

# FEED (Stream 5) - Specified in mol/s
const F_CO2_in = 744.303 
const F_H2_in  = 2232.909
const F_CO_in  = 0.0      
const F_H2O_in = 0.0
const F_MeOH_in= 0.0

# ==============================================================================
# 4. GEOMETRY & PHYSICS
# ==============================================================================
const Tube_Length = 6.0         
const Tube_ID     = 0.04        
const Ac_Tube     = (pi * Tube_ID^2) / 4.0

const D_p_nominal = 0.0055      
const Void_Frac   = 0.4         
const Psi         = 1.0         
const D_eff       = D_p_nominal * Psi 

const Particle_Density = 1775.0 
const Cat_Bulk_Density = Particle_Density * (1.0 - Void_Frac) 
const Viscosity = 1.8e-5        

# ==============================================================================
# 5. KINETICS (Calibrated Van-Dal 2013)
# ==============================================================================
function get_kinetics_real(T, P)
    # 1. Intrinsic Rates
    k1_intr = 1.07 * exp(40000.0 / (R_gas * T))
    k_rwgs_intr = 1.22e10 * exp(-98084.0 / (R_gas * T))

    # 2. Apply Calibrated Effectiveness Factors
    k1 = k1_intr * ETA_MEOH
    k_rwgs = k_rwgs_intr * ETA_RWGS

    # 3. Adsorption & Equilibrium
    K_water_h2 = 3453.38
    K_h2_sqrt = 0.499 * exp(17197.0 / (R_gas * T))
    K_water = 6.62e-11 * exp(124119.0 / (R_gas * T))

    K_eq1_ideal = EQConst1(T) 
    K_eq3_ideal = EQConst3(T) 
    
    safe_P = clamp(P, 1.0, 200.0)
    phis = fugacity(T, safe_P) 
    
    K_phi_1 = (phis["CH3OH"] * phis["H2O"]) / (phis["CO2"] * phis["H2"]^3)
    K_phi_3 = (phis["CO"] * phis["H2O"]) / (phis["CO2"] * phis["H2"])

    K_eq1_real = K_eq1_ideal / K_phi_1
    K_eq3_real = K_eq3_ideal / K_phi_3

    return k1, k_rwgs, K_water_h2, K_h2_sqrt, K_water, K_eq1_real, K_eq3_real
end

# ==============================================================================
# 6. SOLVER (SI UNITS: mol/s, kg)
# ==============================================================================
function pbr_integrated!(dF, F, p, W)
    Ac_Total_Current = p[1]
    
    F_CO2, F_H2, F_CO, F_MeOH, F_H2O, P_bar = F[1], F[2], F[3], F[4], F[5], F[6]
    F_total = sum(F[1:5])
    
    if P_bar <= 2.0; dF .= 0.0; return; end

    k1, k_rwgs, K_wh, K_h, K_w, K_eq1, K_eq3 = get_kinetics_real(T_op, P_bar)
    
    p_CO2  = P_bar * (F_CO2 / F_total)
    p_H2   = P_bar * (F_H2  / F_total)
    p_CO   = P_bar * (F_CO  / F_total)
    p_MeOH = P_bar * (F_MeOH/ F_total)
    p_H2O  = P_bar * (F_H2O / F_total)
    
    DEN_base = 1.0 + (K_wh * p_H2O / p_H2) + (K_h * sqrt(p_H2)) + (K_w * p_H2O)
    
    drive_1 = (p_CO2 * p_H2) * (1.0 - (1.0/K_eq1) * (p_MeOH * p_H2O) / (p_CO2 * p_H2^3))
    r1 = (k1 * drive_1) / (DEN_base^3)

    drive_2 = p_CO2 * (1.0 - (1.0/K_eq3) * (p_CO * p_H2O) / (p_CO2 * p_H2))
    r2 = (k_rwgs * drive_2) / (DEN_base)

    dF[1] = -r1 - r2
    dF[2] = -3*r1 - r2
    dF[3] = r2
    dF[4] = r1
    dF[5] = r1 + r2

    # Momentum Balance (Ergun)
    MW_mix = (F_CO2*0.04401 + F_H2*0.002016 + F_CO*0.02801 + F_MeOH*0.03204 + F_H2O*0.018015) / F_total
    rho_g = (P_bar * 1e5 * MW_mix) / (R_gas * T_op)
    G = (F_total * MW_mix) / Ac_Total_Current
    
    term1 = (150 * (1 - Void_Frac) * Viscosity) / D_eff
    term2 = 1.75 * G
    constant_group = - (G / (rho_g * D_eff)) * ((1 - Void_Frac) / Void_Frac^3)
    
    dF[6] = (constant_group * (term1 + term2) * 1e-5) / (Cat_Bulk_Density * Ac_Total_Current)
end

# ==============================================================================
# 7. DIAGNOSTIC SEARCH
# ==============================================================================
println("----------------------------------------------------------")
println("DIGITAL TWIN: Simulating Calibrated Performance")
println("Constraints: Eta_MeOH=$(ETA_MEOH), Eta_RWGS=$(ETA_RWGS)")
println("----------------------------------------------------------")

const Target_Flow_kmol_hr = 394.1
const N_min = 1000
const N_max = 5000

# Tracking "Best" Configuration
best_N = 0
max_MeOH_seen = 0.0
best_sol = nothing
target_hit = false

for N in N_min:10:N_max
    Current_Ac_Total = Ac_Tube * N
    W_Total_Fixed = Current_Ac_Total * Tube_Length * Cat_Bulk_Density
    
    u0 = [F_CO2_in, F_H2_in, F_CO_in, F_MeOH_in, F_H2O_in, P_in]
    p = [Current_Ac_Total]
    
    prob = ODEProblem(pbr_integrated!, u0, (0.0, W_Total_Fixed), p)
    sol = solve(prob, Rodas4(autodiff=false), reltol=1e-6, abstol=1e-6)
    
    Final_MeOH_kmol_hr = sol.u[end][4] * 3.6
    
    if Final_MeOH_kmol_hr > max_MeOH_seen
        max_MeOH_seen = Final_MeOH_kmol_hr
        best_N = N
        best_sol = sol
    end
    
    if Final_MeOH_kmol_hr >= Target_Flow_kmol_hr
        target_hit = true
        break
    end
end

# ==============================================================================
# 8. GENERATE REPORT (For Best N Found)
# ==============================================================================
if best_sol !== nothing
    # Metrics
    W_Total_Best = Ac_Tube * best_N * Tube_Length * Cat_Bulk_Density
    
    # Peak Analysis
    MeOH_Profile = [u[4] * 3.6 for u in best_sol.u]
    Peak_MeOH, Peak_Idx = findmax(MeOH_Profile)
    
    Cat_Used_At_Peak = best_sol.t[Peak_Idx]
    
    Final_MeOH = best_sol.u[end][4] * 3.6
    Final_H2O  = best_sol.u[end][5] * 3.6
    Final_CO   = best_sol.u[end][3] * 3.6
    
    Peak_H2O   = best_sol.u[Peak_Idx][5] * 3.6
    Peak_CO    = best_sol.u[Peak_Idx][3] * 3.6
    
    Conversion_At_Peak = (F_CO2_in - best_sol.u[Peak_Idx][1]) / F_CO2_in * 100.0
    Yield_At_Peak      = (best_sol.u[Peak_Idx][4] / F_CO2_in) * 100.0
    dP_Total           = P_in - best_sol.u[end][6]
    
    if target_hit
        println("\n>>> TARGET REACHED <<<")
    else
        println("\n>>> TARGET MISSED - SHOWING BEST ATTEMPT <<<")
        println("(Expected due to calibrated efficiency factors)")
        Diff = Target_Flow_kmol_hr - max_MeOH_seen
        @printf("Shortfall: %.2f kmol/hr\n", Diff)
    end
    
    println("\nConfiguration:")
    @printf("  Number of Tubes:        %d\n", best_N)
    @printf("  Tube Length:            %.2f m\n", Tube_Length)
    @printf("  Total Catalyst Mass:    %.2f kg\n", W_Total_Best)
    
    # --- RESTORED SECTIONS ---
    println("\nProduction Statistics:")
    @printf("  Peak CH3OH Production:  %.2f kmol/hr\n", Peak_MeOH)
    @printf("  Catalyst used at Peak:  %.2f kg\n", Cat_Used_At_Peak)
    
    println("\nSpecies at Peak (kmol/hr):")
    @printf("  Water Production:       %.2f\n", Peak_H2O)
    @printf("  CO Production:          %.2f\n", Peak_CO)
    # -------------------------

    println("\nReactor Outlet:")
    @printf("  Final CH3OH Production: %.2f kmol/hr\n", Final_MeOH)
    @printf("  Final CO Production:    %.2f kmol/hr\n", Final_CO)
    @printf("  Final Water Production: %.2f kmol/hr\n", Final_H2O)

    println("\nPerformance:")
    @printf("  CO2 Conversion (Peak):  %.2f %%\n", Conversion_At_Peak)
    @printf("  Methanol Yield (Peak):  %.2f %%\n", Yield_At_Peak)
    @printf("  Total Pressure Drop:    %.2f bar\n", dP_Total)

else
    println("Simulation Failed to run.")
end
println("----------------------------------------------------------")