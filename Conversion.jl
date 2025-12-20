module Conversion

export conv

"""
   conv(roots)
   Calculates CO2 Conversion per pass.
   Formula: (Moles CO2 Consumed) / (Moles CO2 Feed)
          = e3 / 1.0
   
   Expects roots = [epsilon_3]
"""
function conv(roots)
    # roots[2] corresponds to epsilon_3 (the RWGS reaction extent)
    e3 = roots[2]

    # Since Feed is 1.0 and only Rxn 3 consumes CO2:
    X_CO2 = e3
    
    return X_CO2
end

end # module