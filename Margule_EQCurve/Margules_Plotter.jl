import Pkg
Pkg.activate(".")  # Activate the environment in the current folder

# Margules_Plotter

    # Import external libraries
    using PyPlot
    using Printf

    # Import your specific project modules

    using Margule_EQCurve

    """
        plot_margules_comparison()

    Generates a high-quality plot comparing experimental/calculated scatter data 
    against the smooth Margules model curves derived from fitted parameters A12 and A21.
    """
    function plot_margules_comparison()
        
        # ---------------------------------------------------------
        # 1. Retrieve Data (Scatter Points)
        # ---------------------------------------------------------
        x_exp      = VLE_MeH2O.x_methanol
        
        # Experimental/Calculated values
        y_fitting  = GeRTx1x2.get_GeRTx1x2() # Contains endpoint artifacts
        y_gert     = GeRT.get_GeRT()         # Excess Gibbs Energy
        lng1, lng2 = LnGamma.get_ln_gammas() # Natural log of activity coeffs

        # ---------------------------------------------------------
        # 2. Retrieve Fitted Parameters
        # ---------------------------------------------------------
        # We call the regression module to get the constants
        A12, A21 = Final_Margules_Fit.perform_margules_fit()
        
        println("\nPlotting using Margules Parameters:")
        @printf("A12: %.4f\n", A12)
        @printf("A21: %.4f\n", A21)

        # ---------------------------------------------------------
        # 3. Define Model Equations (Smooth Curves)
        # ---------------------------------------------------------
        # Generate 100 points for smooth lines
        x_model = range(0, 1, length=100)
        
        model_fit_param = Float64[]
        model_gert      = Float64[]
        model_lng1      = Float64[]
        model_lng2      = Float64[]

        for x1 in x_model
            x2 = 1.0 - x1
            
            # --- Explicit Margules Equations ---
            
            # Equation A: Fitting Parameter (Linearized form)
            # y = A21*x1 + A12*x2
            val_fit = (A21 * x1) + (A12 * x2)
            push!(model_fit_param, val_fit)
            
            # Equation B: Excess Gibbs Energy (Ge/RT)
            # Ge/RT = x1*x2 * (A21*x1 + A12*x2)
            val_gert = x1 * x2 * ((A21 * x1) + (A12 * x2))
            push!(model_gert, val_gert)
            
            # Equation C: ln(gamma1)
            # ln(g1) = x2^2 * [A12 + 2(A21 - A12)x1]
            val_lng1 = (x2^2) * (A12 + (2 * (A21 - A12) * x1))
            push!(model_lng1, val_lng1)
            
            # Equation D: ln(gamma2)
            # ln(g2) = x1^2 * [A21 + 2(A12 - A21)x2]
            val_lng2 = (x1^2) * (A21 + (2 * (A12 - A21) * x2))
            push!(model_lng2, val_lng2)
        end

        # ---------------------------------------------------------
        # 4. Plotting
        # ---------------------------------------------------------
        fig, ax1 = plt.subplots(figsize=(10, 7))

        # --- A. Scatter Plots (Data) ---
        # Note: We slice [2:end-1] for the fitting param to hide endpoint artifacts
        ax1.scatter(x_exp[2:end-1], y_fitting[2:end-1], color="black", marker="s", s=30, facecolors="none", label="_nolegend_")
        ax1.scatter(x_exp, y_gert,    color="green", marker="o", s=30, facecolors="none", label="_nolegend_") 
        ax1.scatter(x_exp, lng1,      color="blue",  marker="^", s=30, facecolors="none", label="_nolegend_")
        ax1.scatter(x_exp, lng2,      color="red",   marker="v", s=30, facecolors="none", label="_nolegend_")

        # --- B. Model Curves (Lines) ---
        ax1.plot(x_model, model_fit_param, "k--", linewidth=1)   # Fitting Line
        ax1.plot(x_model, model_gert,      "g-",  linewidth=1.5) # Ge/RT curve
        ax1.plot(x_model, model_lng1,      "b-",  linewidth=1.5) # ln(g1) curve
        ax1.plot(x_model, model_lng2,      "r-",  linewidth=1.5) # ln(g2) curve

        # --- C. Formatting ---
        # Minimalist Setup: No Title, No Legend
        ax1.set_xlabel(raw"$x_1$ (Mole Fraction Methanol)", fontsize=12)
        
        # X-Axis limits
        ax1.set_xlim(0.0, 1.0)
        
        # Y-Axis Setup (Start at 0, no negative quadrants)
        # Dynamic limit based on max value in curves + 10% buffer
        y_max = max(maximum(model_fit_param), maximum(model_lng1), maximum(model_lng2)) * 1.1
        ax1.set_ylim(0.0, y_max)
        
        # Grid
        ax1.grid(true, which="both", linestyle="--", linewidth=0.5)

        # Double Y-Axes (Mirror)
        ax2 = ax1.twinx()
        ax2.set_ylim(0.0, y_max)
        ax2.set_ylabel("") # Left blank as requested

        # Display
        plt.tight_layout()
        display(fig)
        
        # Optional: Return figure object if saving is needed later
        return fig
    end

plot_margules_comparison()