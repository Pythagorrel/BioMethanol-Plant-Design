module BioMethEQK

export EQConst1, EQConst2, EQConst3

const T_0 = 298.15
const R = 8.314

struct cp_parameters
    A::Float64
    B::Float64
    C::Float64
    D::Float64
end

const methanol_cp = cp_parameters(2.211, 12.216e-3, -3.45e-6, 0.0)
const co2_cp = cp_parameters(5.457, 1.045e-3, 0.0, -1.157e5)
const co_cp = cp_parameters(3.376, 0.557e-3, 0.0, -0.031e5)
const h2_cp = cp_parameters(3.249, 0.422e-3, 0.0, 0.083e5)
const h2o_cp = cp_parameters(3.470, 1.450e-3, 0.0, 0.121e5)

const formation_parameters = Dict("CH3OH" => (ΔG0=-161960.0, ΔH0=-200660.0), "CO2" => (ΔG0=-394359.0, ΔH0=-393509.0), "CO" => (ΔG0=-137169.0, ΔH0=-110525.0), "H2" => (ΔG0=0.0, ΔH0=0.0), "H2O" => (ΔG0=-228572.0, ΔH0=-241818.0))

const ΔA1 = (methanol_cp.A .+ h2o_cp.A) .- (co2_cp.A .+ 3 * h2_cp.A)
const ΔA2 = methanol_cp.A .- (co_cp.A .+ 2 * h2_cp.A)
const ΔA3 = (co_cp.A .+ h2o_cp.A) .- (co2_cp.A .+ h2_cp.A)

const ΔB1 = (methanol_cp.B .+ h2o_cp.B) .- (co2_cp.B .+ 3 * h2_cp.B)
const ΔB2 = methanol_cp.B .- (co_cp.B .+ 2 * h2_cp.B)
const ΔB3 = (co_cp.B .+ h2o_cp.B) .- (co2_cp.B .+ h2_cp.B)

const ΔC1 = (methanol_cp.C .+ h2o_cp.C) .- (co2_cp.C .+ 3 * h2_cp.C)
const ΔC2 = methanol_cp.C .- (co_cp.C .+ 2 * h2_cp.C)
const ΔC3 = (co_cp.C .+ h2o_cp.C) .- (co2_cp.C .+ h2_cp.C)

const ΔD1 = (methanol_cp.D .+ h2o_cp.D) .- (co2_cp.D .+ 3 * h2_cp.D)
const ΔD2 = methanol_cp.D .- (co_cp.D .+ 2 * h2_cp.D)
const ΔD3 = (co_cp.D .+ h2o_cp.D) .- (co2_cp.D .+ h2_cp.D)

const ΔG1 = (formation_parameters["CH3OH"].ΔG0 .+ formation_parameters["H2O"].ΔG0) .- (formation_parameters["CO2"].ΔG0 .+ 3 * formation_parameters["H2"].ΔG0)
const ΔG2 = formation_parameters["CH3OH"].ΔG0 .- (formation_parameters["CO"].ΔG0 .+ 2 * formation_parameters["H2"].ΔG0)
const ΔG3 = (formation_parameters["CO"].ΔG0 .+ formation_parameters["H2O"].ΔG0) .- (formation_parameters["CO2"].ΔG0 .+ formation_parameters["H2"].ΔG0)

const ΔH1 = (formation_parameters["CH3OH"].ΔH0 .+ formation_parameters["H2O"].ΔH0) .- (formation_parameters["CO2"].ΔH0 .+ 3 * formation_parameters["H2"].ΔH0)
const ΔH2 = formation_parameters["CH3OH"].ΔH0 .- (formation_parameters["CO"].ΔH0 .+ 2 * formation_parameters["H2"].ΔH0)
const ΔH3 = (formation_parameters["CO"].ΔH0 .+ formation_parameters["H2O"].ΔH0) .- (formation_parameters["CO2"].ΔH0 .+ formation_parameters["H2"].ΔH0)

struct rxn_params
    ΔA::Float64
    ΔB::Float64
    ΔC::Float64
    ΔD::Float64
    ΔG::Float64
    ΔH::Float64
end

rxn1 = rxn_params(ΔA1, ΔB1, ΔC1, ΔD1, ΔG1, ΔH1)
rxn2 = rxn_params(ΔA2, ΔB2, ΔC2, ΔD2, ΔG2, ΔH2)
rxn3 = rxn_params(ΔA3, ΔB3, ΔC3, ΔD3, ΔG3, ΔH3)


function EQConst(T::Float64, rxn::rxn_params)

    K₀ = exp.(-1.0 .* (rxn.ΔG) ./ (R .* T_0))
    K₁ = exp.((rxn.ΔH / (R .* T_0)) .* (1.0 .- T_0 ./ T))
    K₂ = exp.(rxn.ΔA .* ((log.(T ./ T_0) .- ((T .- T_0) ./ T))) .+ 0.5 * rxn.ΔB .* (((T .- T_0) .^ 2) ./ T) .+ (rxn.ΔC ./ 6) .* (((T .- T_0) .^ 2) .* ((T .+ 2 * T_0)) ./ T) .+ 0.5 * rxn.ΔD .* ((T .- T_0) .^ 2) ./ ((T .^ 2.) .* (T_0 .^ 2.)))

    K = K₀ .* K₁ .* K₂
    return K
end

function EQConst1(T)
    return EQConst(T, rxn1)
end


function EQConst2(T)
    return EQConst(T, rxn2)
end


function EQConst3(T)
    return EQConst(T, rxn3)
end

end