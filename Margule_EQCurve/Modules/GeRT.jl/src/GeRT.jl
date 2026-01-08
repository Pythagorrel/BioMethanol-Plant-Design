module GeRT

    # Import necessary modules
    
    using Reexport 


    @reexport using VLE_MeH2O
    @reexport using LnGamma


    export get_GeRT

    """
        get_GeRT()

    Calculates the Excess Gibbs Energy (Ge/RT) for the Methanol-Water system.
    Formula: Ge/RT = x1*ln(gamma1) + x2*ln(gamma2)
    
    Returns:
        ge_rt_list :: Vector{Float64}
    """
    function get_GeRT()
        
        # Access composition data
        x1_data = VLE_MeH2O.x_methanol
        
        # Get natural log of gammas from previous module
        ln_gamma1, ln_gamma2 = LnGamma.get_ln_gammas()
        
        # Define epsilon for endpoint correction
        epsilon = 1e-7
        
        # Initialize array to store Ge/RT results
        ge_rt_list = Float64[]

        # Loop through data points
        for i in 1:length(x1_data)
            
            # Read raw composition
            x1_raw = x1_data[i]
            
            # Apply epsilon correction to match the logic used in LnGamma
            # strictly: x2 = 1 - x1
            x1 = max(x1_raw, epsilon)
            x2 = max(1.0 - x1_raw, epsilon)
            
            # Calculate Ge/RT
            # Formula: (x1 * ln_gamma1) + (x2 * ln_gamma2)
            val = (x1 * ln_gamma1[i]) + (x2 * ln_gamma2[i])
            
            push!(ge_rt_list, val)
        end
        
        return ge_rt_list
    end

end