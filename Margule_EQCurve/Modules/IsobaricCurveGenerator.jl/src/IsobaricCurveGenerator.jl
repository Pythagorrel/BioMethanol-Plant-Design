module IsobaricCurveGenerator

    # 1. Import Dependencies
    
    using Reexport 
    @reexport using Final_Margules_Fit  # To get A12, A21
    @reexport using BubblePointSolver   # To solve T and y for each point

   

    export generate_isobaric_data

    """
        generate_isobaric_data(P_bar, num_points)

    Generates a full T-x-y equilibrium dataset at a specified pressure.
    
    Inputs:
        P_bar      : Target Pressure in bar (e.g., 1.2)
        num_points : Resolution of the curve (e.g., 1000)

    Returns:
        (x_list, y_list, T_list) :: Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    """
    function generate_isobaric_data(P_bar::Float64, num_points::Int)
        
        # A. Retrieve Margules Parameters from your regression module
        # This ensures we are always using the latest fitted values.
        A12, A21 = Final_Margules_Fit.perform_margules_fit()
        
        # B. Convert Pressure to kPa (Antoine calc requires kPa)
        P_target_kPa = P_bar * 100.0
        
        # C. Initialize Storage
        x_list = Float64[]
        y_list = Float64[]
        T_list = Float64[]
        
        println("\nGenerating Isobaric Data...")
        println("Pressure:  $P_bar bar")
        println("Points:    $num_points")
        println("...calculating...")

        # D. Generation Loop
        # Create a range from 0 to 1
        x_range = range(0.0, 1.0, length=num_points)

        for x1 in x_range
            
            # Delegate to the Solver
            # We pass the parameters we just retrieved
            T_sol, y_sol = BubblePointSolver.solve_bubble_point(x1, P_target_kPa, A12, A21)
            
            push!(x_list, x1)
            push!(y_list, y_sol)
            push!(T_list, T_sol)
        end
        
        println("Done.")
        return x_list, y_list, T_list
    end

end