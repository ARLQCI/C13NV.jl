using TestItems

# Test that various options return the expected labels (quantum numbers) and
# the expected shape of dynamical generators


@testitem "Coherent |G⟩ manifold, no Ω, 1 spin" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Defaults: with_defaults
    using C13NV.Units: ms
    using QuantumPropagators.Interfaces: check_generator

    H, labels = make_nv_system(; with_defaults()...)

    @test labels == [("G", "0", "↑"), ("G", "0", "↓"),]
    @test size(H) == (2, 2)

    Ψ₀ = ket(("G", "0", "↑"), labels)
    tlist = [0.0, 1ms]
    @test check_generator(H; state = Ψ₀, tlist)

end

@testitem "Coherent |G⟩ manifold, no Ω, 2 spin" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Defaults: with_defaults
    using C13NV.Units: MHz, kHz, ms
    using QuantumPropagators.Interfaces: check_generator

    A₁ = [
        0        0     200kHz
        0        0     200kHz
        200kHz 200kHz   1MHz
    ]

    A₂ = [
        0        0     250kHz
        0        0     250kHz
        250kHz 250kHz  1.2MHz
    ]

    H, labels = make_nv_system(; with_defaults()..., hyperfine_tensors = [A₁, A₂])

    @test labels == [("G", "0", "↑↑"), ("G", "0", "↑↓"), ("G", "0", "↓↑"), ("G", "0", "↓↓")]
    @test size(H) == (4, 4)

    Ψ₀ = ket(("G", "0", "↑↑"), labels)
    @test size(Ψ₀) == (4,)
    tlist = [0.0, 1ms]
    @test check_generator(H; state = Ψ₀, tlist)

end


@testitem "Coherent |G⟩ manifold, no Ω, 3 spin" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Defaults: with_defaults
    using C13NV.Units: MHz, kHz, ms
    using QuantumPropagators.Interfaces: check_generator

    A₁ = [
        0        0     200kHz
        0        0     200kHz
        200kHz 200kHz   1MHz
    ]

    A₂ = [
        0        0     250kHz
        0        0     250kHz
        250kHz 250kHz  1.2MHz
    ]

    A₃ = [
        0        0     150kHz
        0        0     250kHz
        150kHz 250kHz  1.3MHz
    ]

    H, labels = make_nv_system(; with_defaults()..., hyperfine_tensors = [A₁, A₂, A₃])

    @test labels == [
        ("G", "0", "↑↑↑"),
        ("G", "0", "↑↑↓"),
        ("G", "0", "↑↓↑"),
        ("G", "0", "↑↓↓"),
        ("G", "0", "↓↑↑"),
        ("G", "0", "↓↑↓"),
        ("G", "0", "↓↓↑"),
        ("G", "0", "↓↓↓")
    ]
    @test size(H) == (8, 8)

    Ψ₀ = ket(("G", "0", "↑↑↑"), labels)
    @test size(Ψ₀) == (8,)
    tlist = [0.0, 1ms]
    @test check_generator(H; state = Ψ₀, tlist)

end


@testitem "Coherent |G⟩ manifold, Ω₊, 1 spin" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.Units: MHz, kHz, ms
    using C13NV.Defaults: with_defaults
    using QuantumPropagators.Interfaces: check_generator

    H, labels = make_nv_system(;
        Ω₊ = ConstantDrive(257kHz),
        ω₊ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        with_defaults()...
    )

    @test labels == [("G", "+1", "↑"), ("G", "+1", "↓"), ("G", "0", "↑"), ("G", "0", "↓"),]
    @test size(H.ops[1]) == (4, 4)

    Ψ₀ = ket(("G", "0", "↑"), labels)
    @test size(Ψ₀) == (4,)
    tlist = [0.0, 1ms]
    @test check_generator(H; state = Ψ₀, tlist)

end


@testitem "Coherent |G⟩ manifold, Ω₋, 1 spin" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.Units: MHz, kHz, ms
    using C13NV.Defaults: with_defaults
    using QuantumPropagators.Interfaces: check_generator

    H, labels = make_nv_system(;
        Ω₋ = ConstantDrive(257kHz),
        ω₋ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        with_defaults()...
    )

    @test labels == [("G", "0", "↑"), ("G", "0", "↓"), ("G", "-1", "↑"), ("G", "-1", "↓"),]
    @test size(H.ops[1]) == (4, 4)

    Ψ₀ = ket(("G", "0", "↑"), labels)
    @test size(Ψ₀) == (4,)
    tlist = [0.0, 1ms]
    @test check_generator(H; state = Ψ₀, tlist)

end

@testitem "Coherent |G⟩ manifold, Ω₊/Ω₋, 1 spin" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.Units: MHz, kHz, ms
    using C13NV.Defaults: with_defaults
    using QuantumPropagators.Interfaces: check_generator

    H, labels = make_nv_system(;
        Ω₋ = ConstantDrive(257kHz),
        ω₋ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        Ω₊ = ConstantDrive(257kHz),
        ω₊ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        with_defaults()...
    )

    @test labels == [
        ("G", "+1", "↑"),
        ("G", "+1", "↓"),
        ("G", "0", "↑"),
        ("G", "0", "↓"),
        ("G", "-1", "↑"),
        ("G", "-1", "↓"),
    ]
    @test size(H.ops[1]) == (6, 6)

    Ψ₀ = ket(("G", "0", "↑"), labels)
    @test size(Ψ₀) == (6,)
    tlist = [0.0, 1ms]
    @test check_generator(H; state = Ψ₀, tlist)

end


@testitem "Incoherent |G⟩ manifold, Ω₊/Ω₋, 2 spins" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.Units: MHz, kHz, ms, ns
    using C13NV.Defaults: with_defaults
    using QuantumPropagators.Interfaces: check_generator

    A₁ = [
        0        0     200kHz
        0        0     200kHz
        200kHz 200kHz   1MHz
    ]

    A₂ = [
        0        0     250kHz
        0        0     250kHz
        250kHz 250kHz  1.2MHz
    ]

    L, labels = make_nv_system(;
        Ω₋ = ConstantDrive(257kHz),
        ω₋ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        Ω₊ = ConstantDrive(257kHz),
        ω₊ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        with_defaults()...,
        hyperfine_tensors = [A₁, A₂],
        γ₊₁ = 1/400ns,
        γ₋₁ = 1/400ns,
    )

    @test labels == [
        ("G", "+1", "↑↑"),
        ("G", "+1", "↑↓"),
        ("G", "+1", "↓↑"),
        ("G", "+1", "↓↓"),
        ("G", "0", "↑↑"),
        ("G", "0", "↑↓"),
        ("G", "0", "↓↑"),
        ("G", "0", "↓↓"),
        ("G", "-1", "↑↑"),
        ("G", "-1", "↑↓"),
        ("G", "-1", "↓↑"),
        ("G", "-1", "↓↓"),
    ]
    @test size(L.ops[1]) == (12^2, 12^2)

    Ψ₀ = ket(("G", "0", "↑↑"), labels)
    @test size(Ψ₀) == (12,)
    tlist = [0.0, 1ms]
    ρ₀ = vec(Ψ₀ * Ψ₀')
    @test check_generator(L; state = ρ₀, tlist)

end


@testitem "Full manifold, no Ω, 1 spin, static Λ" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Defaults: with_defaults
    using C13NV.Units: ms
    using QuantumPropagators.Interfaces: check_generator

    L, labels = make_nv_system(; with_defaults()..., Λ = 1.0)

    @test labels == [
        ("G", "0", "↑"),
        ("G", "0", "↓"),
        ("E", "0", "↑"),
        ("E", "0", "↓"),
        ("M", "↑"),
        ("M", "↓"),
    ]
    @test size(L.ops[1]) == (6^2, 6^2)

    Ψ₀ = ket(("G", "0", "↑"), labels)
    @test size(Ψ₀) == (6,)
    tlist = [0.0, 1ms]
    ρ₀ = vec(Ψ₀ * Ψ₀')
    @test check_generator(L; state = ρ₀, tlist)

end


@testitem "Full manifold, Ω₊, 1 spin, static Λ" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.Defaults: with_defaults
    using C13NV.Units: MHz, kHz, ms
    using QuantumPropagators.Interfaces: check_generator

    L, labels = make_nv_system(;
        with_defaults()...,
        Ω₊ = ConstantDrive(257kHz),
        ω₊ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        Λ = 1.0,
    )

    @test labels == [
        ("G", "+1", "↑"),
        ("G", "+1", "↓"),
        ("G", "0", "↑"),
        ("G", "0", "↓"),
        ("E", "+1", "↑"),
        ("E", "+1", "↓"),
        ("E", "0", "↑"),
        ("E", "0", "↓"),
        ("M", "↑"),
        ("M", "↓"),
    ]

    @test size(L.ops[1]) == (10^2, 10^2)

    Ψ₀ = ket(("G", "0", "↑"), labels)
    @test size(Ψ₀) == (10,)
    tlist = [0.0, 1ms]
    ρ₀ = vec(Ψ₀ * Ψ₀')
    @test check_generator(L; state = ρ₀, tlist)

end


@testitem "Full manifold, Ω₊/Ω₋, 2 spin, time-dependent Λ" begin

    using C13NV.Models: make_nv_system, ket
    using C13NV.Amplitudes: ConstantDrive, LinearChirp
    using C13NV.Defaults: with_defaults
    using C13NV.Units: MHz, kHz, ms
    using QuantumPropagators.Shapes: flattop
    using QuantumPropagators.Interfaces: check_generator

    T = 1ms
    tlist = collect(range(0.0, T; length = 1001))

    Λ(t) = flattop(t; T, t_rise = 0.1T)

    A₁ = [
        0        0     200kHz
        0        0     200kHz
        200kHz 200kHz   1MHz
    ]

    A₂ = [
        0        0     250kHz
        0        0     250kHz
        250kHz 250kHz  1.2MHz
    ]

    L, labels = make_nv_system(;
        with_defaults()...,
        Ω₋ = ConstantDrive(257kHz),
        ω₋ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        Ω₊ = ConstantDrive(257kHz),
        ω₊ = LinearChirp(t₀ = 0.1568ms, α = 26.24MHz / ms),
        hyperfine_tensors = [A₁, A₂],
        Λ,
    )

    @test Λ ∈ L.amplitudes

    @test labels == [
        ("G", "+1", "↑↑"),
        ("G", "+1", "↑↓"),
        ("G", "+1", "↓↑"),
        ("G", "+1", "↓↓"),
        ("G", "0", "↑↑"),
        ("G", "0", "↑↓"),
        ("G", "0", "↓↑"),
        ("G", "0", "↓↓"),
        ("G", "-1", "↑↑"),
        ("G", "-1", "↑↓"),
        ("G", "-1", "↓↑"),
        ("G", "-1", "↓↓"),
        ("E", "+1", "↑↑"),
        ("E", "+1", "↑↓"),
        ("E", "+1", "↓↑"),
        ("E", "+1", "↓↓"),
        ("E", "0", "↑↑"),
        ("E", "0", "↑↓"),
        ("E", "0", "↓↑"),
        ("E", "0", "↓↓"),
        ("E", "-1", "↑↑"),
        ("E", "-1", "↑↓"),
        ("E", "-1", "↓↑"),
        ("E", "-1", "↓↓"),
        ("M", "↑↑"),
        ("M", "↑↓"),
        ("M", "↓↑"),
        ("M", "↓↓"),
    ]

    @test size(L.ops[1]) == (28^2, 28^2)

    Ψ₀ = ket(("G", "0", "↑↑"), labels)
    @test size(Ψ₀) == (28,)
    ρ₀ = vec(Ψ₀ * Ψ₀')
    @test check_generator(L; state = ρ₀, tlist)

end

using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
