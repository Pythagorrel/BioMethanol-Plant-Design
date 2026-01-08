module BubblePointSolver

    # Dependencies for Saturation Pressure

    using Reexport 

    @reexport using AntoineMethanol
    @reexport using AntoineWater

    export solve_bubble_point

    """
        solve_bubble_point(x1, P_target, A12, A21)
    
    Calculates T_bubble and y1 for a specific liquid composition x1.
    Optimized: Calculates Activity Coefficients once (outside the T-loop).
    """
    function solve_bubble_point(x1, P_target, A12, A21)
        
        # --- 1. Pre-Calculate Activity Coefficients (Model) ---
        # We do this once because Gamma is assumed independent of T.
        x2 = 1.0 - x1
        
        # Margules Equations
        ln_gamma1 = (x2^2) * (A12 + 2 * (A21 - A12) * x1)
        ln_gamma2 = (x1^2) * (A21 + 2 * (A12 - A21) * x2)
        
        gamma1 = exp(ln_gamma1)
        gamma2 = exp(ln_gamma2)

        # --- 2. Initialize Solver ---
        # Search Range for Temperature (°C)
        T_low = 40.0
        T_high = 160.0
        epsilon = 1e-5 
        
        T_solution = 0.0
        P1_sat_final = 0.0

        # --- 3. Iteration Loop (Find T) ---
        for iter in 1:100
            
            T_guess = (T_low + T_high) / 2.0
            
            # Only Psat needs to be updated in the loop
            P1_sat = AntoineMethanol.PsatM(T_guess)
            P2_sat = AntoineWater.PsatW(T_guess)
            
            # Check Pressure
            P_calc = (x1 * gamma1 * P1_sat) + (x2 * gamma2 * P2_sat)
            
            if abs(P_calc - P_target) < epsilon
                T_solution = T_guess
                P1_sat_final = P1_sat
                break
            elseif P_calc < P_target
                T_low = T_guess
            else
                T_high = T_guess
            end
        end
        
        # --- 4. Final Result ---
        y1 = (x1 * gamma1 * P1_sat_final) / P_target
        
        return T_solution, y1
    end

end