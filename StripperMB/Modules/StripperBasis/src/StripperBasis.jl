module StripperBasis

# Export the specific calculation functions
export calculate_L0, calculate_VN_plus_1, calculate_xA0
export calculate_L_prime, calculate_V_prime
export calculate_equilibrium_m

# ==========================================
# RAW INPUTS: Update these numbers if PFD changes
# ==========================================
const FLASH_LIQ_MEOH = 394.1   # kmol/hr
const FLASH_LIQ_H2O  = 532.2   # kmol/hr
const FLASH_LIQ_CO2  = 66.25   # kmol/hr (Dissolved)
const FRESH_H2_FEED  = 1320.0  # kmol/hr (Stripping Gas)

# === PHYSICS CONSTANTS ===
const HENRYS_CONSTANT = 160.0 # bar
const STRIPPER_PRESSURE = 45.0 # bar


# ==========================================

"""
Calculates L0: The total liquid molar flow rate entering the top of the stripper.
"""
function calculate_L0()
    return FLASH_LIQ_MEOH + FLASH_LIQ_H2O + FLASH_LIQ_CO2
end

"""
Calculates VN+1: The total gas molar flow rate entering the bottom of the stripper.
"""
function calculate_VN_plus_1()
    return FRESH_H2_FEED
end

"""
Calculates xA0: The mole fraction of Solute A (CO2) in the liquid feed L0.
"""
function calculate_xA0()
    # We call our own function here to ensure consistency
    L_total = calculate_L0()
    return FLASH_LIQ_CO2 / L_total
end

"""
Calculates L': The flow rate of the liquid on a solute-free basis (Solvent only).
Formula: Moles Methanol + Moles Water.
"""
function calculate_L_prime()
    return FLASH_LIQ_MEOH + FLASH_LIQ_H2O
end

"""
Calculates V': The flow rate of the gas on a solute-free basis (Carrier only).
Formula: Equals Fresh H2 Feed (assuming 100% purity).
"""
function calculate_V_prime()
    return FRESH_H2_FEED
end

"""
Calculates the Equilibrium Constant (m-value) based on Henry's Law.
Formula: m = H / P
"""
function calculate_equilibrium_m()
    return HENRYS_CONSTANT / STRIPPER_PRESSURE
end

end # module