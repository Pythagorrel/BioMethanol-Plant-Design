module LnGamma

    # Import the necessary modules from the parent scope
	using Reexport

	@reexport using AntoineMethanol
	@reexport using AntoineWater
	@reexport using VLE_MeH2O 	  

    export get_ln_gammas

    """
        get_ln_gammas()

    Calculates the natural logarithm of activity coefficients (ln gamma) for the 
    Methanol-Water system.
    
    Returns:
        (ln_gamma1_list, ln_gamma2_list) :: Tuple{Vector{Float64}, Vector{Float64}}
    """
    function get_ln_gammas()
        
        # [cite_start]Access data from the VLE module [cite: 7, 8]
        T_data   = VLE_MeH2O.T_degC
        x1_data  = VLE_MeH2O.x_methanol
        y1_data  = VLE_MeH2O.y_methanol
        P_system = VLE_MeH2O.P_kPa
        
        # Define epsilon for endpoint correction (prevents log(0) errors)
        epsilon = 1e-7
        
        # Initialize arrays to store ln(gamma) results
        ln_gamma1_list = Float64[]
        ln_gamma2_list = Float64[]

        for i in 1:length(T_data)
            t = T_data[i]
            
            # Read raw data
            x1_raw = x1_data[i]
            y1_raw = y1_data[i]
            
            # Apply epsilon correction for calculation safety
            # Necessary because ln(0) is undefined
            x1 = max(x1_raw, epsilon)
            y1 = max(y1_raw, epsilon)
            x2 = max(1.0 - x1_raw, epsilon)
            y2 = max(1.0 - y1_raw, epsilon)
            
            # [cite_start]Calculate Saturation Pressures (Antoine Eqn) [cite: 20, 29]
            P1_sat = AntoineMethanol.PsatM(t)
            P2_sat = AntoineWater.PsatW(t)
            
            # Calculate Activity Coefficients (Gamma)
            # [cite_start]Formula: gamma_i = (y_i * P) / (x_i * P_sat_i) [cite: 2]
            gamma1 = (y1 * P_system) / (x1 * P1_sat)
            gamma2 = (y2 * P_system) / (x2 * P2_sat)
            
            # Calculate Natural Logarithm and Store
            push!(ln_gamma1_list, log(gamma1))
            push!(ln_gamma2_list, log(gamma2))
        end
        
        # Return the lists of natural logs
        return ln_gamma1_list, ln_gamma2_list
    end

end