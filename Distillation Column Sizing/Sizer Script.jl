# --- Dependencies & Constants ---
const R_gas = 8314.0       # J/kmolK
const MW_MeOH = 32.04
const MW_Water = 18.02
const ρ_MeOH_std = 750.0   # kg/m3 
const ρ_Water_std = 978.0  # kg/m3 
const flooding_safety = 0.85 

# --- Structs ---
struct Section
    name::String
    T::Float64          # K
    P::Float64          # Pa
    L_flow::Float64     # kmol/hr
    V_flow::Float64     # kmol/hr
    x_MeOH::Float64     # Mole Frac
    y_MeOH::Float64     # Mole Frac
    K1::Float64         # User Input
    downcomer::Float64  # Fraction (0.05 or 0.12)
end

# --- Physics Functions ---
function get_mw_liquid(s::Section)
    return s.x_MeOH * MW_MeOH + (1 - s.x_MeOH) * MW_Water
end

function get_mw_vapor(s::Section)
    return s.y_MeOH * MW_MeOH + (1 - s.y_MeOH) * MW_Water
end

function get_rho_vapor(s::Section)
    return (s.P * get_mw_vapor(s)) / (R_gas * s.T)
end

function get_rho_liquid(s::Section)
    w_MeOH = (s.x_MeOH * MW_MeOH) / get_mw_liquid(s)
    w_Water = 1.0 - w_MeOH
    rho_mix = 1.0 / ( (w_MeOH / ρ_MeOH_std) + (w_Water / ρ_Water_std) )
    return rho_mix
end

function calculate_diameter(s::Section)
    # 1. Densities
    ρ_v = get_rho_vapor(s)
    ρ_l = get_rho_liquid(s)
    
    # 2. Velocities
    u_flood = s.K1 * sqrt((ρ_l - ρ_v) / ρ_v)
    u_design = u_flood * flooding_safety
    
    # 3. Volumetric Vapor Flow (m3/s)
    V_mass_kg_s = (s.V_flow * get_mw_vapor(s)) / 3600.0
    Q_v = V_mass_kg_s / ρ_v
    
    # 4. Areas
    Area_net = Q_v / u_design
    Area_total = Area_net / (1 - s.downcomer)
    
    # 5. Diameter
    D = sqrt((4 * Area_total) / π)
    
    return D, u_flood, u_design, Area_total
end

# --- DEFINE SECTIONS ---

# C1 Top: Low liquid, 5% DC
c1_top = Section("C1 Top", 342.94, 1.2e5, 23.46, 54.58, 0.952, 0.952, 
                 0.10, 0.01)

# C1 Bottom: High liquid (Flv=0.68), 12% DC required
c1_bot = Section("C1 Bottom", 352.65, 1.2e5, 950.09, 54.58, 0.407, 0.734, 
                 0.048, 0.12)

# C2 Top: Low liquid, 5% DC
c2_top = Section("C2 Top", 342.63, 1.2e5, 449.27, 823.66, 0.972, 0.972, 
                 0.102, 1e-5)

# C2 Bottom: Low/Med liquid, 5% DC ok (Flv=0.07)
c2_bot = Section("C2 Bottom", 377.77, 1.2e5, 1344.78, 823.66, 0.001, 0.007, 
                 0.101, 0.275)

# --- EXECUTION ---
sections = [c1_top, c1_bot, c2_top, c2_bot]

println("--- FINAL COLUMN SIZING ---")
for s in sections
    D, uf, ud, At = calculate_diameter(s)
    println("\n$(s.name):")
    println("  K1 Used:        $(s.K1)")
    println("  Design Velocity: $(round(ud, digits=2)) m/s")
    println("  Required Diam:   $(round(D, digits=3)) m")
end