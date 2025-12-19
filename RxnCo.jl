"""
Updated solver that accepts a dictionary of fugacity coefficients (phi).
"""
module RxnCo
export rxnco
include("BioMethEQK.jl")
using .BioMethEQK
using NLsolve
function rxnco(T, P, phi)
    # --- 1. Define Constants ---
    P0  = 1.0
    
    # Calculate Equilibrium constants for current T
    K1 = EQConst1(T)
    K2 = EQConst2(T)
    K3 = EQConst3(T)
    
    # --- 2. Define the System of Equations ---
    function f!(F, x)
        e1, e2, e3 = x[1], x[2], x[3]

        # Denominator (Total Moles)
        D = 4 - 2*e1 - 2*e2
        if D ≈ 0; D = 1e-9; end

        # Mole Fractions
        y_CO2   = (1 - e1 - e3) / D
        y_H2    = (3 - 3*e1 - 2*e2 - e3) / D
        y_CO    = (-e2 + e3) / D
        y_CH3OH = (e1 + e2) / D
        y_H2O   = (e1 + e3) / D

        # Reaction Quotients using the PASSED `phi` values
        # Eq 1
        A_num = (phi["CH3OH"] * y_CH3OH) * (phi["H2O"] * y_H2O)
        A_den = (phi["CO2"]   * y_CO2)   * (phi["H2"]  * y_H2)^3
        A = A_num / A_den
        
        # Eq 2
        B_num = (phi["CH3OH"] * y_CH3OH)
        B_den = (phi["CO"]    * y_CO)    * (phi["H2"]  * y_H2)^2
        B = B_num / B_den
        
        # Eq 3
        C_num = (phi["CO"]    * y_CO)    * (phi["H2O"] * y_H2O)
        C_den = (phi["CO2"]   * y_CO2)   * (phi["H2"]  * y_H2)
        C = C_num / C_den

        # Residuals
        F[1] = ((P0/P)^2 * A) - K1
        F[2] = ((P0/P)^2 * B) - K2
        F[3] = ((P0/P)^0 * C) - K3 
    end

    # --- 3. Solve ---
    x0 = [0.1, 0.05, 0.1]
    result = nlsolve(f!, x0)
    return result.zero
end
end