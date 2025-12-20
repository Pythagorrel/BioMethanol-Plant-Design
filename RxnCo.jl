module RxnCo

include("BioMethEQK.jl")
using .BioMethEQK
using NLsolve

export rxnco

"""
    rxnco(T, P, phi, init_guess, n_CO2_0, n_H2_0)

Solves the 2-reaction methanol synthesis system for ANY feed ratio.

Arguments:
- T, P: Temperature (K), Pressure (bar)
- phi: Fugacity coefficient dictionary
- init_guess: [epsilon_2, epsilon_3]
- n_CO2_0: Initial moles of CO2 (Feed)
- n_H2_0:  Initial moles of H2  (Feed)
"""
function rxnco(T, P, phi, init_guess, n_CO2_0, n_H2_0)
    P0 = 1.0
    
    # Calculate Equilibrium Constants (from BioMethEQK.jl)
    K2 = EQConst2(T)
    K3 = EQConst3(T)

    # Calculate Total Initial Moles for the Denominator
    n_total_0 = n_CO2_0 + n_H2_0

    function f!(F, x)
        # Solve for e2 and e3. Force e1 to 0 (Implied).
        e2, e3 = x[1], x[2]
        
        # --- Stability Guards ---
        # e3 must be > e2 for CO to exist.
        if e2 < -0.1 || e3 < -0.1 || e2 > 1.2 * n_CO2_0 || e3 > 1.2 * n_CO2_0
             F[1] = 1e6 * (e2 + 10)
             F[2] = 1e6 * (e3 + 10)
             return
        end

        # Denominator (Total Moles)
        # Formula: n_total_0 + (-2)*e2 + (0)*e3
        D = max(n_total_0 - 2*e2, 1e-6)
        
        # Mole Fractions (Dynamic Feed)
        # CO2: Consumed by Rxn 3
        y_CO2   = max(1e-12, (n_CO2_0 - e3) / D)       
        
        # H2: Consumed by Rxn 2 (2x) and Rxn 3 (1x)
        y_H2    = max(1e-12, (n_H2_0 - 2*e2 - e3) / D) 
        
        # CO: Produced by e3, consumed by e2. (Intermediate)
        y_CO    = max(1e-12, (e3 - e2) / D)       
        
        # MeOH: Produced by Rxn 2
        y_CH3OH = max(1e-12, (e2) / D)            
        
        # H2O: Produced by Rxn 3
        y_H2O   = max(1e-12, (e3) / D)            
        
        # Reaction Quotients
        Q2 = (phi["CH3OH"]*y_CH3OH) / (phi["CO"]*y_CO * (phi["H2"]*y_H2)^2)
        Q3 = (phi["CO"]*y_CO * phi["H2O"]*y_H2O) / (phi["CO2"]*y_CO2 * phi["H2"]*y_H2)
            
        # Residuals
        term_P = (P/P0)^(-2)
        F[1] = (Q2 * term_P) / K2 - 1.0
        F[2] = Q3 / K3 - 1.0
    end
    
    return nlsolve(f!, init_guess, method = :trust_region, autoscale=true)
end

end # module