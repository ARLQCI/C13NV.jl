using TestItems


# Test that we can reproduce the dynamics we've shown previously, e.g. on GRC
# poster. Derived from examples/Polarization_with_Pumping.jl


@testmodule PolarizationPumpingTest begin
    using C13NV.Units: kHz, ns
    using QuantumPropagators.Shapes: flattop

    const γcB = 128.4kHz
    const γcB_per_2π = 128.4kHz / 2π
    const γcB⁻¹ = 2π / 128.4kHz

    const Ω₀ = 17γcB_per_2π
    const t₀ = 10γcB⁻¹
    const α = 20γcB_per_2π / γcB⁻¹
    const Γ = (1 / (12ns))

    const tlist = collect(range(0, 4 * 10γcB⁻¹; length = 4001))

    @kwdef struct SinglePump <: Function
        amplitude::Float64
        tlist::Vector{Float64} = tlist
        pump_width::Float64 = 0.5
        ramp_width::Float64 = 0.1
        pump_at::Float64 = 0.5
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

    """Get matrix of population for `states`."""
    function get_pops(states)::Matrix{Float64}
        if ndims(states) == 1
            return get_pops(reshape(states, (length(states), 1)))
        end
        N², nt = size(states)
        N = isqrt(N²)
        pops = zeros(N, nt)
        for i = 1:nt
            for j = 0:(N-1)
                p = abs(states[(j*N)+j+1, i])
                pops[j+1, i] = p
            end
        end
        return pops
    end
end


@testitem "Polarization with optical pumping dynamics" setup = [PolarizationPumpingTest] begin
    using C13NV.Models: make_nv_system, ket
    using C13NV.Defaults: with_defaults
    using C13NV.Units: ns
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using QuantumPropagators: propagate, Newton
    using UnPack: @unpack

    @unpack Ω₀, t₀, α, Γ, tlist, SinglePump, get_pops = PolarizationPumpingTest

    Λ = SinglePump(Γ; pump_width = 0.4, ramp_width = 0.01)

    L, labels = make_nv_system(;
        with_defaults()...,
        Ω₋ = ConstantDrive(Ω₀),
        ω₋ = LinearChirp(; t₀, α),
        Γ = (1 / (12ns)),
        Λ,
        frame = :rwa,
    )

    N = length(labels)

    Ψ_up = ket(("G", "0", "↑"), labels)
    Ψ_down = ket(("G", "0", "↓"), labels)
    ρ₀ = reshape(0.5 * Ψ_up * Ψ_up' + 0.5 * Ψ_down * Ψ_down', N * N)

    states = propagate(ρ₀, L, tlist; method = Newton, check = true, storage = true)

    pops = get_pops(states)

    # Expected populations at t=2000 (during chirp, before pumping completes)
    expected_pops_2000 = [0.497, 0.000868, 0.497, 0.00531]

    # Check first 4 levels at t=2000 with 3 significant digits (rtol=1e-3)
    @test pops[1, 2000] ≈ expected_pops_2000[1] rtol = 1e-3
    @test pops[2, 2000] ≈ expected_pops_2000[2] rtol = 1e-3
    @test pops[3, 2000] ≈ expected_pops_2000[3] rtol = 1e-3
    @test pops[4, 2000] ≈ expected_pops_2000[4] rtol = 1e-3

    # Check that the first level has population > 0.98 at the end (successful polarization)
    @test pops[1, end] > 0.98
end

using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
