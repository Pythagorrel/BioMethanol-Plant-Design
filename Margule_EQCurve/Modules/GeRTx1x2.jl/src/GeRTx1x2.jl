module GeRTx1x2

    # Import necessary modules

    using Reexport
    
    @reexport using VLE_MeH2O
    @reexport using GeRT 


    export get_GeRTx1x2

    """
        get_GeRTx1x2()

    Calculates the normalized Excess Gibbs Energy parameter (Ge / RT / x1 / x2)
    used for linearizing and fitting thermodynamic models.
    
    Returns:
        fitting_list :: Vector{Float64}
    """
    function get_GeRTx1x2()
        
        # Access composition data
        x1_data = VLE_MeH2O.x_methanol
        
        # Get Ge/RT values from previous module
        # Note: calling the function defined in your GeRT module
        ge_rt_data = get_GeRT()
        
        # Define epsilon for endpoint correction
        epsilon = 1e-7
        
        # Initialize array to store results
        fitting_list = Float64[]

        for i in 1:length(x1_data)
            
            # Read raw composition
            x1_raw = x1_data[i]
            ge_rt_val = ge_rt_data[i]
            
            # Apply epsilon correction to prevent division by zero
            # strictly: x2 = 1 - x1
            x1 = max(x1_raw, epsilon)
            x2 = max(1.0 - x1_raw, epsilon)
            
            # Calculate (Ge/RT) / (x1*x2)
            # This represents (A21*x1 + A12*x2) in the Margules model
            val = ge_rt_val / (x1 * x2)
            
            push!(fitting_list, val)
        end
        
        return fitting_list
    end

end