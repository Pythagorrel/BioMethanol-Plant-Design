# ==========================================
# 1. LOAD MODULES
# ==========================================
# Assuming these are in the same folder. 
# If they are in the same file, you can just remove the 'include' lines.
include("BioMethEQK.jl")
include("FugacityCo.jl")
include("Yield.jl")
include("Conversion.jl") 
include("RxnCo.jl")
using .BioMethEQK
using .FugacityCo
using .Yield
using .Conversion
using .RxnCo
using NLsolve
using Printf

# 2. FEED PARAMETERS 
FEED_CO2 = 1.0
FEED_H2  = 12.0  

# ==========================================
# 3. MAIN EXECUTION LOOP
# ==========================================

T_range = 473.15:20.0:573.15 
P_range = 50.0:10.0:100.0   
current_guess = [0.45, 0.50]

for T in T_range
    println("\n==================================================================================")
    # UPDATED HEADER: Now includes the Feed Ratio
    println("Results for T = $(T) K  |  Feed Ratio CO2:H2 = 1:12")
    println("==================================================================================")
    
    @printf("%-10s %-12s %-12s %-15s %-15s %-15s\n", 
            "P (bar)", "epsilon_2", "epsilon_3", "CH3OH Yield", "S_CH3OH/CO (%)", "CO2 Conv")
    println("-"^85)
    
    for P in P_range
        current_phi = fugacity(T, P)
        sol = rxnco(T, P, current_phi, current_guess, FEED_CO2, FEED_H2)
        
        if converged(sol)
            roots = sol.zero
            
            # Calculate metrics
            yield_val = yieldF(roots)
            sel_val   = selectivityF(roots)
            conv_val  = conv(roots)
            
            @printf("%-10.1f %-12.4f %-12.4f %-15.4f %-15.2f %-15.4f\n", 
                    P, roots[1], roots[2], yield_val, sel_val, conv_val)
            
            global current_guess = roots
        else
            @printf("%-10.1f %s\n", P, "FAILED TO CONVERGE")
            global current_guess = [0.45, 0.50]
        end
    end
    println("\n")
end
