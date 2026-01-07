module StripperXYEQ

using Reexport
@reexport using StripperBasis

export calculate_op_line_points, calculate_equilibrium_curve

# A clean container for the four numbers needed to draw the line
struct OpLineCoords
    X_top::Float64    # Feed Liquid Ratio (X_0)
    Y_top::Float64    # Exit Vapor Ratio (Y_1)
    X_btm::Float64    # Exit Liquid Ratio (X_N)
    Y_btm::Float64    # Inlet Vapor Ratio (Y_N+1) - Usually 0.0
end

"""
calculate_op_line_points(target_removal::Float64)

Calculates the start and end coordinates for the Operating Line 
based on a single specific removal target.

Input: target_removal (e.g., 0.995)
Output: OpLineCoords struct
"""
function calculate_op_line_points(target_removal::Float64)
    
    # 1. GET CONSTANTS
    L0 = calculate_L0()
    L_prime = calculate_L_prime()
    V_prime = calculate_V_prime()
    
    # 2. CALCULATE TOP POINT (Feed Condition)
    # The Feed X is fixed by the Flash, independent of the stripper performance.
    xA0 = calculate_xA0()
    X_feed = xA0 / (1.0 - xA0)
    
    # 3. CALCULATE BOTTOM POINT (Lean Liquid Condition)
    CO2_inlet = L0 - L_prime
    CO2_removed = CO2_inlet * target_removal
    
    # Liquid Out
    LN = L0 - CO2_removed
    xaN = 1.0 - (L_prime / LN)
    
    if (1.0 - xaN) == 0
        X_bottom = 0.0
    else
        X_bottom = xaN / (1.0 - xaN)
    end
    
    # 4. CALCULATE Y_TOP (Rich Gas Condition)
    # Using the OP Line Slope equation: Y1 = (L'/V')(X0 - XN)
    slope = L_prime / V_prime
    Y_top = slope * (X_feed - X_bottom)
    
    # 5. RETURN COORDINATES
    # Line goes from (X_bottom, 0) -> (X_feed, Y_top)
    return OpLineCoords(X_feed, Y_top, X_bottom, 0.0)
end

# (Keep equilibrium curve function as-is, it's perfect)
function calculate_equilibrium_curve(num_points=100)
    # ... (Same code as before) ...
    m = calculate_equilibrium_m()
    xA0 = calculate_xA0()
    x_fracs = range(0, stop=xA0 * 1.1, length=num_points)
    X_ratios = Float64[]
    Y_ratios = Float64[]
    for x in x_fracs
        y = m * x
        if y >= 1.0; y = 0.999; end
        push!(X_ratios, x/(1-x))
        push!(Y_ratios, y/(1-y))
    end
    return X_ratios, Y_ratios
end

end # module