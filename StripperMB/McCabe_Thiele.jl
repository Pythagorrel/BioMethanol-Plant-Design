import Pkg
Pkg.activate(".")  # Activate the environment in the current folder

using StripperMB
using PyPlot

# ==========================================
# 1. SETUP
# ==========================================
target = 0.995  # <--- CHANGE THIS SINGLE VARIABLE TO UPDATE DESIGN

# Get the specific line coordinates for this target
op_line = calculate_op_line_points(target)

# Get the physics curve
X_eq, Y_eq = calculate_equilibrium_curve(100)

# ==========================================
# 2. PLOT
# ==========================================
figure(figsize=(10, 8))

# Plot Equilibrium
plot(X_eq, Y_eq, "k-", linewidth=2, label="Equilibrium (y=3.56x)")

# Plot Operating Line (Connecting the two coordinates)
# X-coords: [Bottom, Top]
# Y-coords: [Bottom, Top]
plot([op_line.X_btm, op_line.X_top], [op_line.Y_btm, op_line.Y_top], 
     "b-", linewidth=2, marker="o", label="OP Line ($(target*100)% Removal)")

# ==========================================
# 3. FORMATTING
# ==========================================
title("McCabe-Thiele Diagram (Target: $(target*100)%)")
xlabel("Liquid Mole Ratio (X)")
ylabel("Vapor Mole Ratio (Y)")
grid(true, which="both", linestyle="--")
minorticks_on()
legend()

# Display Values for verification
println("Drawing Line from:")
println("  Bottom: ($(op_line.X_btm), $(op_line.Y_btm))")
println("  Top:    ($(op_line.X_top), $(op_line.Y_top))")

display(gcf())