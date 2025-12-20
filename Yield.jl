module Yield

export yieldF, selectivityF

"""
   yieldF(roots)
   Returns Methanol Yield (Moles MeOH / Moles CO2 feed).
   Since basis is 1.0 mol CO2, this is just epsilon_2.
"""
function yieldF(roots)
    e2 = roots[1]
    return e2
end

"""
   selectivityF(roots)
   Returns Carbon Selectivity to Methanol (%).
   Formula: (Moles MeOH / Total Carbon Converted) * 100
          = (e2 / e3) * 100
"""
function selectivityF(roots)
    e2 = roots[1]
    e3 = roots[2]
    
    # Avoid division by zero if e3 is 0 (at start)
    if e3 < 1e-9
        return 0.0
    end
    
    return (e2 / e3) * 100.0
end

end # module