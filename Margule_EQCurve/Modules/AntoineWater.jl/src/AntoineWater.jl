module AntoineWater

    export PsatW

    # [cite_start]Constants for Water from Table B.2 (Continued) [cite: 29]
    # [cite_start]Valid Temperature Range: 0 to 200 °C [cite: 29]
    const A = 16.3872
    const B = 3885.70
    const C = 230.170

    """
        PsatW(T_degC)

    Calculates the saturation pressure of Water in kPa given temperature in °C.
    Formula: ln(P_sat) = A - B / (T + C)
    """
    function PsatW(T_degC::Float64)
        # Antoine Equation
        ln_Psat = A - (B / (T_degC + C))
        
        # Return P_sat in kPa
        return exp(ln_Psat)
    end

end