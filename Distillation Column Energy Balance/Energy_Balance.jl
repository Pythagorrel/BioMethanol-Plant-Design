struct Watson_parameters
   ΔH_n::Float64  # Latent Heat at Normal Boiling Point (J/mol)
    T_n::Float64        # Normal Boiling Point (K)
    T_c::Float64        # Critical Temperature (K)
end

methanol_watson = Watson_parameters(35210.0, 337.8, 512.5)
water_watson = Watson_parameters(40660.0, 373.2, 647.1)

function ΔH_LV_pure(comp::Watson_parameters, T)
    # 1. Calculate Reduced Temperatures
    Tr_actual = T / comp.T_c
    Tr_normal = comp.T_n / comp.T_c

    # 2. Apply Watson's Equation
    # We prevent negative roots if T > Tc (unlikely in distillation but good safety)
    if T >= comp.T_c
        return 0.0
    end

    factor = ((1 - Tr_actual) / (1 - Tr_normal))^0.38
    ΔH = comp.ΔH_n * factor

    return ΔH
end

function ΔH_LV_mix(ΔH_m,ΔH_w,x_m)
    ΔH_LV_mix = (x_m*ΔH_m)+((1-x_m)*ΔH_w)
    return ΔH_LV_mix
end

struct L_cp_parameters 
    A::Float64
    B::Float64
    C::Float64
end

methanol_cp = L_cp_parameters(13.431,-51.28e-3,131.13e-6)
water_cp = L_cp_parameters(8.712,1.25e-3,-0.18e-6)

const x_D1 = 0.952; const x_f = 0.425; const x_B1 = 0.407; 
const R = 8.314 # Jmol⁻¹
const T_ref = 298.15 #K
const T_f = 360.75
const T_D = 342.94
const T_B = 352.65
const D = 31.12 #kmol/hr
const B = 895.51 #kmol/hr
const Reflux = 0.754
const F = 926.63
const V = 54.58 

function cp_pure(sys::L_cp_parameters,T)
    cp_pure = R*(sys.A+(sys.B*T)+(sys.C*T^2))
    return cp_pure
end

function cp_mix(cp_m,cp_w,x_m)
    cp_mix = (x_m*cp_m)+((1-x_m)*cp_w)
    return cp_mix
end

#-------------------Cₚ calculations for the various mixtures-----------------
Fcp_m = cp_pure(methanol_cp,T_f)
Fcp_w = cp_pure(water_cp, T_f)
Fcp_mix = cp_mix(Fcp_m, Fcp_w,x_f) #Cp of the feed mixture 

Dcp_m = cp_pure(methanol_cp,T_D)
Dcp_w = cp_pure(water_cp, T_D)
Dcp_mix = cp_mix(Dcp_m, Dcp_w,x_D1) #Cp of the feed mixture 

Bcp_m = cp_pure(methanol_cp,T_B)
Bcp_w = cp_pure(water_cp, T_B)
Bcp_mix = cp_mix(Bcp_m, Bcp_w, x_B1) #Cp of the feed mixture 
#------------------------Enthalpy calculations-------------------------------
h_F = Fcp_mix * (T_f - T_ref)
h_D = Dcp_mix * (T_D - T_ref)
h_B = Bcp_mix * (T_B - T_ref)
#------------------------Enthalpy calculations-------------------------------

DΔH_m = ΔH_LV_pure(methanol_watson, T_D)
DΔH_w = ΔH_LV_pure(water_watson, T_D)
DΔH_mix = ΔH_LV_mix(DΔH_m, DΔH_w, x_D1) 

Q_C = V *DΔH_mix #kJ
Q_B = (D*h_D+B*h_B+Q_C)-F*h_F
# QB​=[D⋅hD​+B⋅hB​+QC​]−F⋅hF​ (final equation)

Q_C_W = Q_C/3600
Q_B_W = Q_B/3600

println("$Q_C_W KW is the condenser duty.")
println("$Q_B_W KW is the reboiler duty.")
