module AntoineMethanol

    export PsatM

    # [cite_start]Constants for Methanol from Table B.2 [cite: 20]
    # Valid Temperature Range: -11 to 83 °C
    const A = 16.5785
    const B = 3638.27
    const C = 239.500

    """
        PsatM(T_degC)

    Calculates the saturation pressure of Methanol in kPa given temperature in °C.
    Formula: ln(P_sat) = A - B / (T + C)
    """
    function PsatM(T_degC::Float64)
        # Antoine Equation
        ln_Psat = A - (B / (T_degC + C))
        
        # Return P_sat in kPa
        return exp(ln_Psat)
    end

end