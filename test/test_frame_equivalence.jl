using TestItems


@testitem "RWA and diagonal frame propagation equivalence" begin
    # Test that propagation in the RWA frame gives the same results as
    # propagation in the diagonal frame followed by transformation back to RWA.
    # This validates the frame transformation theory from notes/hamiltonian.qmd.

    using C13NV.Models: make_nv_system, ket
    using C13NV.Defaults: DEFAULTS
    using C13NV.Units: kHz, MHz, μs, ns
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.FrameTransformations: transform_frame
    using QuantumPropagators: propagate, ExpProp
    using LinearAlgebra: norm
    using UnPack: @unpack

    # Use the same setup as examples/Frame_Equivalence.jl
    γcB = 128.4kHz
    γcB_per_2π = 128.4kHz / 2π
    γcB⁻¹ = 2π / 128.4kHz

    Ω₀ = 17γcB_per_2π
    t₀ = 10γcB⁻¹
    α = 20γcB_per_2π / γcB⁻¹

    tlist = collect(range(0, 2 * 10γcB⁻¹; length = 501))

    # System in RWA frame
    H_rwa, labels = make_nv_system(;
        DEFAULTS...,
        Ω₋ = ConstantDrive(Ω₀),
        ω₋ = LinearChirp(; t₀, α),
        Γ = (1 / (12ns)),
        frame = :rwa,
    )

    # System in diagonal frame
    H_diag, _ = make_nv_system(;
        DEFAULTS...,
        Ω₋ = ConstantDrive(Ω₀),
        ω₋ = LinearChirp(; t₀, α),
        Γ = (1 / (12ns)),
        frame = :diag,
    )

    @unpack A_zz, A_zx = DEFAULTS

    # Test with spin-up initial state
    Ψ_up = ket(("G", "0", "↑"), labels)
    Ψ_up_diag = transform_frame(Ψ_up; A_zz, A_zx, to_frame = :diag)

    states_rwa_up = propagate(Ψ_up, H_rwa, tlist; method = ExpProp, storage = true)
    states_diag_up = propagate(Ψ_up_diag, H_diag, tlist; method = ExpProp, storage = true)

    # Transform diagonal frame results back to RWA
    states_diag_rwa_up = mapslices(
        Ψ -> transform_frame(Ψ; A_zz, A_zx, to_frame = :rwa),
        states_diag_up;
        dims = 1
    )

    @test norm(states_rwa_up - states_diag_rwa_up) < 1e-10

    # Test with spin-down initial state
    Ψ_down = ket(("G", "0", "↓"), labels)
    Ψ_down_diag = transform_frame(Ψ_down; A_zz, A_zx, to_frame = :diag)

    states_rwa_down = propagate(Ψ_down, H_rwa, tlist; method = ExpProp, storage = true)
    states_diag_down =
        propagate(Ψ_down_diag, H_diag, tlist; method = ExpProp, storage = true)

    states_diag_rwa_down = mapslices(
        Ψ -> transform_frame(Ψ; A_zz, A_zx, to_frame = :rwa),
        states_diag_down;
        dims = 1
    )

    @test norm(states_rwa_down - states_diag_rwa_down) < 1e-10
end


@testitem "RWA and diagonal frame equivalence with superposition" begin
    # Test frame equivalence with a superposition initial state

    using C13NV.Models: make_nv_system, ket
    using C13NV.Defaults: DEFAULTS
    using C13NV.Units: kHz, MHz, μs, ns
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.FrameTransformations: transform_frame
    using QuantumPropagators: propagate, ExpProp
    using LinearAlgebra: norm, normalize
    using UnPack: @unpack

    γcB = 128.4kHz
    γcB_per_2π = 128.4kHz / 2π
    γcB⁻¹ = 2π / 128.4kHz

    t₀ = 10γcB⁻¹
    α = 20γcB_per_2π / γcB⁻¹

    tlist = collect(range(0, 2 * 10γcB⁻¹; length = 501))

    H_rwa, labels = make_nv_system(;
        DEFAULTS...,
        Ω₋ = ConstantDrive(17γcB_per_2π),
        ω₋ = LinearChirp(; t₀, α),
        Γ = (1 / (12ns)),
        frame = :rwa,
    )

    H_diag, _ = make_nv_system(;
        DEFAULTS...,
        Ω₋ = ConstantDrive(17γcB_per_2π),
        ω₋ = LinearChirp(; t₀, α),
        Γ = (1 / (12ns)),
        frame = :diag,
    )

    @unpack A_zz, A_zx = DEFAULTS

    # Superposition of spin up and spin down
    Ψ_up = ket(("G", "0", "↑"), labels)
    Ψ_down = ket(("G", "0", "↓"), labels)
    Ψ_superpos = normalize(Ψ_up + im * Ψ_down)

    Ψ_superpos_diag = transform_frame(Ψ_superpos; A_zz, A_zx, to_frame = :diag)

    states_rwa = propagate(Ψ_superpos, H_rwa, tlist; method = ExpProp, storage = true)
    states_diag =
        propagate(Ψ_superpos_diag, H_diag, tlist; method = ExpProp, storage = true)

    states_diag_rwa = mapslices(
        Ψ -> transform_frame(Ψ; A_zz, A_zx, to_frame = :rwa),
        states_diag;
        dims = 1
    )

    @test norm(states_rwa - states_diag_rwa) < 1e-10
end


@testitem "RWA and diagonal frame equivalence without Zeeman" begin
    # Test frame equivalence when B=0 (no Zeeman term).
    # In this case, the only difference between frames is the hyperfine term,
    # which is correctly transformed. This helps isolate any errors in
    # transforming the B̂_I^{(n)} (which is less trivial than Â_I^{(n)})

    using C13NV.Models: make_nv_system, ket
    using C13NV.Units: kHz, MHz, μs, Gauss
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.FrameTransformations: transform_frame
    using QuantumPropagators: propagate, ExpProp
    using LinearAlgebra: norm

    A_zz = 0.8MHz
    A_zx = 0.4MHz
    B = 0.0  # No Zeeman term!
    γ_c = 1.07kHz / Gauss

    γcB_per_2π = 128.4kHz / 2π
    γcB⁻¹ = 2π / 128.4kHz

    t₀ = 8γcB⁻¹
    α = 15γcB_per_2π / γcB⁻¹

    tlist = collect(range(0, 2 * t₀; length = 401))

    H_rwa, labels = make_nv_system(;
        A_zz,
        A_zx,
        B,
        γ_c,
        Ω₋ = ConstantDrive(10γcB_per_2π),
        ω₋ = LinearChirp(; t₀, α),
        frame = :rwa,
    )

    H_diag, _ = make_nv_system(;
        A_zz,
        A_zx,
        B,
        γ_c,
        Ω₋ = ConstantDrive(10γcB_per_2π),
        ω₋ = LinearChirp(; t₀, α),
        frame = :diag,
    )

    Ψ₀ = ket(("G", "0", "↑"), labels)
    Ψ₀_diag = transform_frame(Ψ₀; A_zz, A_zx, to_frame = :diag)

    states_rwa = propagate(Ψ₀, H_rwa, tlist; method = ExpProp, storage = true)
    states_diag = propagate(Ψ₀_diag, H_diag, tlist; method = ExpProp, storage = true)

    states_diag_rwa = mapslices(
        Ψ -> transform_frame(Ψ; A_zz, A_zx, to_frame = :rwa),
        states_diag;
        dims = 1
    )

    # This should pass because B=0 means no Zeeman term to transform
    @test norm(states_rwa - states_diag_rwa) < 1e-10
end


@testitem "RWA and diagonal frame equivalence with diagonal hyperfine" begin
    # Test frame equivalence when A_zx = A_zy = 0 (hyperfine already diagonal).
    # In this case, R = I (identity), so both frames are identical.

    using C13NV.Models: make_nv_system, ket
    using C13NV.Units: kHz, MHz, μs, Gauss
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.FrameTransformations: transform_frame
    using QuantumPropagators: propagate, ExpProp
    using LinearAlgebra: norm

    A_zz = 1.0MHz
    A_zx = 0.0  # Already diagonal!
    B = 100.0Gauss
    γ_c = 1.07kHz / Gauss

    γcB_per_2π = 128.4kHz / 2π
    γcB⁻¹ = 2π / 128.4kHz

    t₀ = 8γcB⁻¹
    α = 15γcB_per_2π / γcB⁻¹

    tlist = collect(range(0, 2 * t₀; length = 401))

    H_rwa, labels = make_nv_system(;
        A_zz,
        A_zx,
        B,
        γ_c,
        Ω₋ = ConstantDrive(10γcB_per_2π),
        ω₋ = LinearChirp(; t₀, α),
        frame = :rwa,
    )

    H_diag, _ = make_nv_system(;
        A_zz,
        A_zx,
        B,
        γ_c,
        Ω₋ = ConstantDrive(10γcB_per_2π),
        ω₋ = LinearChirp(; t₀, α),
        frame = :diag,
    )

    Ψ₀ = ket(("G", "0", "↑"), labels)
    Ψ₀_diag = transform_frame(Ψ₀; A_zz, A_zx, to_frame = :diag)

    states_rwa = propagate(Ψ₀, H_rwa, tlist; method = ExpProp, storage = true)
    states_diag = propagate(Ψ₀_diag, H_diag, tlist; method = ExpProp, storage = true)

    states_diag_rwa = mapslices(
        Ψ -> transform_frame(Ψ; A_zz, A_zx, to_frame = :rwa),
        states_diag;
        dims = 1
    )

    # This should pass because A_zx=0 means R=I (identity transformation)
    @test norm(states_rwa - states_diag_rwa) < 1e-10
end


using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
