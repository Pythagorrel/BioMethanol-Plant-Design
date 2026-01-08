import Pkg
Pkg.activate(".")  # Activate the environment in the current folder

# EQCurve_1.2_BAr

    # 1. Import Dependencies
    using Margule_EQCurve
    using PyPlot
    
    
    """
        plot_final_curve()

    Orchestrates the full simulation pipeline:
    1. Calls IsobaricCurveGenerator to get data at 1.2 bar.
    2. Generates a publication-quality T-x-y equilibrium plot.
    """
    function plot_final_curve()
        
        # --- A. Define Conditions ---
        target_pressure_bar = 1.2
        resolution = 1000
        
        # --- B. Get Data (The Pipeline Trigger) ---
        # This single line triggers:
        # Master -> Generator -> Solver -> Antoine/Margules -> Final Data
        x_data, y_data, T_data = IsobaricCurveGenerator.generate_isobaric_data(target_pressure_bar, resolution)
        
        # --- C. Plotting ---
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # 1. Plot the Equilibrium Curve
        ax.plot(x_data, y_data, "k-", linewidth=2, label="Equilibrium Curve ($target_pressure_bar bar)")
        
        # 2. Plot Reference Diagonal (y=x)
        ax.plot([0, 1], [0, 1], "k--", linewidth=1.0, alpha=0.6)
        
        # 3. Formatting
        ax.set_xlabel(raw"$x_{MeOH}$", fontsize=14)
        ax.set_ylabel(raw"$y_{MeOH}$", fontsize=14)
        
        # Limits (Strict 0 to 1)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        
        # Grid
        ax.grid(true, which="major", linestyle="--", linewidth=0.5, alpha=0.5)
        
        # Minimalist Style (No title, simple legend if needed, or remove)
        # ax.legend(loc="upper left", frameon=false) 
        
        # Display
        plt.tight_layout()
        display(fig)
        
        # Return figure and data if you want to save them to CSV later
        #return fig, x_data, y_data
    end
plot_final_curve()