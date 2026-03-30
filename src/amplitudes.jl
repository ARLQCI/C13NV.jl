module Amplitudes


using ComponentArrays: ComponentVector, Axis
using QuantumControl.Controls: ParameterizedFunction
using QuantumPropagators.Shapes: flattop


struct LinearChirp <: ParameterizedFunction
    values::ComponentVector{Float64,Vector{Float64},Tuple{Axis{(t₀ = 1, α = 2)}}}
    parameters::Vector{Float64}
end


function LinearChirp(; t₀, α)
    return LinearChirp(ComponentVector(; t₀, α), Float64[1.0, 1.0])
end

function (control::LinearChirp)(t)
    t₀ = control.values.t₀ * control.parameters[1]
    α = control.values.α * control.parameters[2]
    return α * (t - t₀)
end


struct ConstantDrive <: ParameterizedFunction
    Ω::Float64
    parameters::Vector{Float64}
    function ConstantDrive(Ω)
        new(Ω, Float64[1.0])
    end
end


function (control::ConstantDrive)(t)
    return control.Ω
end


@kwdef struct SinglePump <: Function
    amplitude::Float64
    tlist::Vector{Float64}
    pump_width::Float64
    ramp_width::Float64
    pump_at::Float64
end


function SinglePump(amplitude; kwargs...)
    return SinglePump(; amplitude, kwargs...)
end


function (Λ::SinglePump)(t::Float64)
    T = Λ.tlist[end]
    w = Λ.pump_width
    r = Λ.ramp_width
    return Λ.amplitude *
           flattop(t; t₀ = Λ.pump_at * T, T = (Λ.pump_at + w) * T, t_rise = (r * w * T))
end

end
