using TestItems

# Test that the `μ` parameter is correctly taken into account. We test this via
# simple Rabi cycling. We set up a system where we observer a π-pulse for the
# default μ = 1. The same system should then show a π/2 pulse if we set μ = 0.5


@testitem "Rabi cycling in |0⟩ ↔ |+1⟩ system" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Amplitudes: ConstantDrive
    using C13NV.Units: μs, kHz
    using C13NV.Defaults: DEFAULTS
    using QuantumPropagators: propagate, ExpProp

    T = 1μs
    tlist = collect(range(0, T; length = 101))

    # For a π pulse: Ω·T = π → Ω = π/T = 0.5 MHz = 500 kHz
    Ω_π = 500kHz

    H₊pi, labels₊ = make_nv_system(;
        DEFAULTS...,
        A_zz = 0.0,
        A_zx = 0.0,
        B = 0.0,
        Ω₊ = ConstantDrive(Ω_π),
        ω₊ = nothing,
        δ₊ = 0.0,
    )

    H₊pi_half, _ = make_nv_system(;
        DEFAULTS...,
        A_zz = 0.0,
        A_zx = 0.0,
        B = 0.0,
        Ω₊ = ConstantDrive(Ω_π),
        ω₊ = nothing,
        δ₊ = 0.0,
        μ = 0.5,
    )

    H₋pi, labels₋ = make_nv_system(;
        DEFAULTS...,
        A_zz = 0.0,
        A_zx = 0.0,
        B = 0.0,
        Ω₋ = ConstantDrive(Ω_π),
        ω₋ = nothing,
        δ₋ = 0.0,
    )

    H₋pi_half, _ = make_nv_system(;
        DEFAULTS...,
        A_zz = 0.0,
        A_zx = 0.0,
        B = 0.0,
        Ω₋ = ConstantDrive(Ω_π),
        ω₋ = nothing,
        δ₋ = 0.0,
        μ = 0.5,
    )

    for electronic_spin in ("+1", "-1")

        H_pi = (electronic_spin == "+1") ? H₊pi : H₋pi
        H_pi_half = (electronic_spin == "+1") ? H₊pi_half : H₋pi_half
        labels = (electronic_spin == "+1") ? labels₊ : labels₋

        for nuclear_spin in ("↑", "↓")

            Ψ₀ = ket(("G", electronic_spin, nuclear_spin), labels)
            Ψ = propagate(Ψ₀, H_pi, tlist; method = ExpProp)

            # After a π pulse: U = exp(-i(π/2)σₓ) = -iσₓ
            # So |±1⟩ → -i|0⟩
            Ψ_pi = -1im * ket(("G", "0", nuclear_spin), labels)

            @test Ψ ≈ Ψ_pi

            # After a π/2 pulse: U = exp(-i(π/4)σₓ) = (1/√2)(I - iσₓ)
            # So |±1⟩ → (1/√2)(|±1⟩ - i|0⟩)
            Ψ_pi_half =
                (
                    ket(("G", electronic_spin, nuclear_spin), labels) -
                    1im * ket(("G", "0", nuclear_spin), labels)
                ) / √2

            Ψ = propagate(Ψ₀, H_pi_half, tlist; method = ExpProp)

            @test Ψ ≈ Ψ_pi_half

        end

    end

end

using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
