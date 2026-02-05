using DifferentialEquations
using Printf

# ==============================================================================
# 1. IMPORT MODULES
# ==============================================================================
include("BioMethEQK.jl")
include("FugacityCo.jl") 

using .BioMethEQK
using .FugacityCo

# ==============================================================================
# 2. OPERATING PARAMETERS
# ==============================================================================
const T_op = 523.15      # 250 C
const P_in = 50.0        # 50 bar Inlet
const R_gas = 8.314      

# FEED (Stream 5)
const F_CO2_in = 514.0   
const F_H2_in  = 1542.0   
const F_CO_in  = 0.0     
const F_H2O_in = 0.0
const F_MeOH_in= 0.0

# ==============================================================================
# 3. INDUSTRIAL CATALYST & GEOMETRY (Table 2 Values)
# ==============================================================================
# Catalyst Properties (Industrial Pellet)
const D_p = 0.0055          # 5.5 mm 
const Void_Frac = 0.4       # 0.4 
const Particle_Density = 1775.0 # kg/m3

# Bulk Density = 1775 * (1 - 0.4) = 1065 kg/m3
const Cat_Bulk_Density = Particle_Density * (1.0 - Void_Frac) 

# Reactor Geometry
const N_Tubes = 2348.0      
const Tube_ID = 0.04        # 40mm
const Ac_Tube = (pi * Tube_ID^2) / 4.0
const Ac_Total = Ac_Tube * N_Tubes 

# Fluid Properties
const Viscosity = 1.8e-5    # Pa.s 

# ==============================================================================
# 4. KINETICS WRAPPER (Real Gas)
# ==============================================================================
function get_kinetics_real(T, P)
    k1 = 1.07 * exp(40000.0 / (R_gas * T))
    k_rwgs = 1.22e10 * exp(-98084.0 / (R_gas * T))

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
# 5. SOLVER (Real Gas + Ergun Pressure Drop)
# ==============================================================================
function pbr_integrated!(dF, F, p, W)
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

    # --- PRESSURE DROP (FIXED) ---
    MW_mix = (F_CO2*44.01 + F_H2*2.016 + F_CO*28.01 + F_MeOH*32.04 + F_H2O*18.015) / F_total / 1000.0
    rho_g = (P_bar * 1e5 * MW_mix) / (R_gas * T_op)

    # FIXED LINE: Removed the extra * 1000.0 multiplication
    m_dot_total = F_total * MW_mix # This is now correctly in kg/s
    G = m_dot_total / Ac_Total

    term1 = (150 * (1 - Void_Frac) * Viscosity) / D_p
    term2 = 1.75 * G
    
    constant_group = - (G / (rho_g * D_p)) * ((1 - Void_Frac) / Void_Frac^3)
    dP_dz_Pa = constant_group * (term1 + term2)
    
    dF[6] = (dP_dz_Pa * 1e-5) / (Cat_Bulk_Density * Ac_Total)
end

# ==============================================================================
# 6. EXECUTION
# ==============================================================================
println("----------------------------------------------------------")
println("Biomethanol Reactor: Final Industrial Check")
println("----------------------------------------------------------")

u0 = [F_CO2_in, F_H2_in, F_CO_in, F_MeOH_in, F_H2O_in, P_in]
W_span = (0.0, 100000.0) 
prob = ODEProblem(pbr_integrated!, u0, W_span)

sol = solve(prob, Rodas4(autodiff=false), reltol=1e-8, abstol=1e-8, maxiters=1e7)

W_plot = sol.t
MeOH_flow_kmol = [u[4] * 3.6 for u in sol.u] 
P_profile = [u[6] for u in sol.u]

max_prod, max_idx = findmax(MeOH_flow_kmol)
W_at_max = W_plot[max_idx]
P_at_max = P_profile[max_idx]

# Calculated Tube Length
Vol_Used = W_at_max / Cat_Bulk_Density
Calc_Length = Vol_Used / Ac_Total

@printf("Peak Production:      %.2f kmol/hr\n", max_prod)
@printf("Catalyst at Peak:     %.2f kg\n", W_at_max)
@printf("Calculated Length:    %.2f m\n", Calc_Length)
@printf("Pressure at Peak:     %.2f bar\n", P_at_max)
@printf("Total Pressure Drop:  %.2f bar\n", P_in - P_at_max)
println("----------------------------------------------------------")