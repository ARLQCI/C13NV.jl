using TestItems


@testmodule DissipatorTest begin
    using C13NV.Models: make_nv_system, ket
    using C13NV.Amplitudes: ConstantDrive
    using C13NV.Units: ns

    # Test parameters: τ = 100ns decay time
    const τ_test = 100.0ns
    const γ_test = 1.0 / τ_test
    const n_τ = 5       # simulate for 5τ
    const n_points = 501

    # Base parameters with all dissipation rates set to 0
    const base_params = Dict{Symbol,Float64}(
        :A_zz => 0.0,
        :A_zx => 0.0,
        :B => 0.0,
        :Γ => 0.0,
        :Γ₀ => 0.0,
        :Γ₊₁ => 0.0,
        :Γ₋₁ => 0.0,
        :Σ₀ => 0.0,
        :Σ₊₁ => 0.0,
        :Σ₋₁ => 0.0,
        :γ₊₁ => 0.0,
        :γ₋₁ => 0.0,
    )

    """Get matrix of populations from vectorized density matrices."""
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

    """Create uniform mixture density matrix over specified indices."""
    function uniform_mixture(labels, indices)
        N = length(labels)
        ρ = zeros(ComplexF64, N, N)
        p = 1.0 / length(indices)
        for i in indices
            ρ[i, i] = p
        end
        return reshape(ρ, N * N)
    end

    """Sum populations over given indices at a specific time index."""
    function combined_pop(pops, indices, t_idx)
        return sum(pops[i, t_idx] for i in indices)
    end

    """Calculate time index for a fraction of total time."""
    function time_index(t_frac)
        return max(1, Int(round(t_frac / n_τ * (n_points - 1))) + 1)
    end

end


# =============================================================================
# Test: Â_Γ (|E⟩ → |G⟩)
# =============================================================================
@testitem "Dissipator Â_Γ: |E⟩ → |G⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    # Create system with only Γ active; Λ=0.0 enables optical space without pumping
    L, labels = make_nv_system(; base_params..., Λ = 0.0, Γ = γ_test)

    # Find E and G state indices (all spin states)
    E_indices = findall(l -> l[1] == "E", labels)
    G_indices = findall(l -> l[1] == "G", labels)

    # Initial state: uniform mixture over all E states
    ρ₀ = uniform_mixture(labels, E_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_E(t_idx) = combined_pop(pops, E_indices, t_idx)
    P_G(t_idx) = combined_pop(pops, G_indices, t_idx)

    # Verify exponential decay: P_E(t) = exp(-γt), P_G(t) = 1 - exp(-γt)
    @test P_E(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_G(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_E(time_index(1)) ≈ exp(-1) rtol = 1e-2
    @test P_G(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2

    @test P_E(time_index(3)) ≈ exp(-3) rtol = 1e-2
    @test P_G(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2

    @test P_E(time_index(5)) ≈ exp(-5) rtol = 1e-2
    @test P_G(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation: ↑/↓ ratio preserved
    E_up = findall(l -> l[1] == "E" && l[end] == "↑", labels)
    E_down = findall(l -> l[1] == "E" && l[end] == "↓", labels)
    G_up = findall(l -> l[1] == "G" && l[end] == "↑", labels)
    G_down = findall(l -> l[1] == "G" && l[end] == "↓", labels)

    # At t=τ, check that ↑ and ↓ populations are equal (started from uniform mixture)
    P_up = combined_pop(pops, E_up, time_index(1)) + combined_pop(pops, G_up, time_index(1))
    P_down =
        combined_pop(pops, E_down, time_index(1)) +
        combined_pop(pops, G_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


# =============================================================================
# Test: Â_Γ₀ (|E,0⟩ → |M⟩)
# =============================================================================
@testitem "Dissipator Â_Γ₀: |E,0⟩ → |M⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    L, labels = make_nv_system(; base_params..., Λ = 0.0, Γ₀ = γ_test)

    # Find |E,0⟩ and |M⟩ state indices
    E0_indices = findall(l -> l[1] == "E" && l[2] == "0", labels)
    M_indices = findall(l -> l[1] == "M", labels)

    # Initial state: uniform mixture over |E,0⟩ states (both nuclear spins)
    ρ₀ = uniform_mixture(labels, E0_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_E0(t_idx) = combined_pop(pops, E0_indices, t_idx)
    P_M(t_idx) = combined_pop(pops, M_indices, t_idx)

    @test P_E0(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_M(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_E0(time_index(1)) ≈ exp(-1) rtol = 1e-2
    @test P_M(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2

    @test P_E0(time_index(3)) ≈ exp(-3) rtol = 1e-2
    @test P_M(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2

    @test P_E0(time_index(5)) ≈ exp(-5) rtol = 1e-2
    @test P_M(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation
    E0_up = findall(l -> l[1] == "E" && l[2] == "0" && l[end] == "↑", labels)
    E0_down = findall(l -> l[1] == "E" && l[2] == "0" && l[end] == "↓", labels)
    M_up = findall(l -> l[1] == "M" && l[end] == "↑", labels)
    M_down = findall(l -> l[1] == "M" && l[end] == "↓", labels)

    P_up =
        combined_pop(pops, E0_up, time_index(1)) + combined_pop(pops, M_up, time_index(1))
    P_down =
        combined_pop(pops, E0_down, time_index(1)) +
        combined_pop(pops, M_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


# =============================================================================
# Test: Â_Γ₋₁ (|E,-1⟩ → |M⟩)
# =============================================================================
@testitem "Dissipator Â_Γ₋₁: |E,-1⟩ → |M⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    # Need Ω₋ to include the -1 spin level
    L, labels = make_nv_system(;
        base_params...,
        Λ = 0.0,
        Ω₋ = ConstantDrive(0.0),  # Include level without coherent dynamics
        Γ₋₁ = γ_test,
    )

    E_m1_indices = findall(l -> l[1] == "E" && l[2] == "-1", labels)
    M_indices = findall(l -> l[1] == "M", labels)

    ρ₀ = uniform_mixture(labels, E_m1_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_Em1(t_idx) = combined_pop(pops, E_m1_indices, t_idx)
    P_M(t_idx) = combined_pop(pops, M_indices, t_idx)

    @test P_Em1(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_M(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_Em1(time_index(1)) ≈ exp(-1) rtol = 1e-2
    @test P_M(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2

    @test P_Em1(time_index(3)) ≈ exp(-3) rtol = 1e-2
    @test P_M(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2

    @test P_Em1(time_index(5)) ≈ exp(-5) rtol = 1e-2
    @test P_M(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation
    Em1_up = findall(l -> l[1] == "E" && l[2] == "-1" && l[end] == "↑", labels)
    Em1_down = findall(l -> l[1] == "E" && l[2] == "-1" && l[end] == "↓", labels)
    M_up = findall(l -> l[1] == "M" && l[end] == "↑", labels)
    M_down = findall(l -> l[1] == "M" && l[end] == "↓", labels)

    P_up =
        combined_pop(pops, Em1_up, time_index(1)) + combined_pop(pops, M_up, time_index(1))
    P_down =
        combined_pop(pops, Em1_down, time_index(1)) +
        combined_pop(pops, M_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


# =============================================================================
# Test: Â_Γ₊₁ (|E,+1⟩ → |M⟩)
# =============================================================================
@testitem "Dissipator Â_Γ₊₁: |E,+1⟩ → |M⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    # Need Ω₊ to include the +1 spin level
    L, labels =
        make_nv_system(; base_params..., Λ = 0.0, Ω₊ = ConstantDrive(0.0), Γ₊₁ = γ_test)

    E_p1_indices = findall(l -> l[1] == "E" && l[2] == "+1", labels)
    M_indices = findall(l -> l[1] == "M", labels)

    ρ₀ = uniform_mixture(labels, E_p1_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_Ep1(t_idx) = combined_pop(pops, E_p1_indices, t_idx)
    P_M(t_idx) = combined_pop(pops, M_indices, t_idx)

    @test P_Ep1(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_M(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_Ep1(time_index(1)) ≈ exp(-1) rtol = 1e-2
    @test P_M(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2

    @test P_Ep1(time_index(3)) ≈ exp(-3) rtol = 1e-2
    @test P_M(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2

    @test P_Ep1(time_index(5)) ≈ exp(-5) rtol = 1e-2
    @test P_M(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation
    Ep1_up = findall(l -> l[1] == "E" && l[2] == "+1" && l[end] == "↑", labels)
    Ep1_down = findall(l -> l[1] == "E" && l[2] == "+1" && l[end] == "↓", labels)
    M_up = findall(l -> l[1] == "M" && l[end] == "↑", labels)
    M_down = findall(l -> l[1] == "M" && l[end] == "↓", labels)

    P_up =
        combined_pop(pops, Ep1_up, time_index(1)) + combined_pop(pops, M_up, time_index(1))
    P_down =
        combined_pop(pops, Ep1_down, time_index(1)) +
        combined_pop(pops, M_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


# =============================================================================
# Test: Â_Σ₀ (|M⟩ → |G,0⟩)
# =============================================================================
@testitem "Dissipator Â_Σ₀: |M⟩ → |G,0⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    L, labels = make_nv_system(; base_params..., Λ = 0.0, Σ₀ = γ_test)

    M_indices = findall(l -> l[1] == "M", labels)
    G0_indices = findall(l -> l[1] == "G" && l[2] == "0", labels)

    ρ₀ = uniform_mixture(labels, M_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_M(t_idx) = combined_pop(pops, M_indices, t_idx)
    P_G0(t_idx) = combined_pop(pops, G0_indices, t_idx)

    @test P_M(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_G0(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_M(time_index(1)) ≈ exp(-1) rtol = 1e-2
    @test P_G0(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2

    @test P_M(time_index(3)) ≈ exp(-3) rtol = 1e-2
    @test P_G0(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2

    @test P_M(time_index(5)) ≈ exp(-5) rtol = 1e-2
    @test P_G0(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation
    M_up = findall(l -> l[1] == "M" && l[end] == "↑", labels)
    M_down = findall(l -> l[1] == "M" && l[end] == "↓", labels)
    G0_up = findall(l -> l[1] == "G" && l[2] == "0" && l[end] == "↑", labels)
    G0_down = findall(l -> l[1] == "G" && l[2] == "0" && l[end] == "↓", labels)

    P_up =
        combined_pop(pops, M_up, time_index(1)) + combined_pop(pops, G0_up, time_index(1))
    P_down =
        combined_pop(pops, M_down, time_index(1)) +
        combined_pop(pops, G0_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


# =============================================================================
# Test: Â_Σ₋₁ (|M⟩ → |G,-1⟩)
# =============================================================================
@testitem "Dissipator Â_Σ₋₁: |M⟩ → |G,-1⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    # Need Ω₋ to include the -1 spin level
    L, labels =
        make_nv_system(; base_params..., Λ = 0.0, Ω₋ = ConstantDrive(0.0), Σ₋₁ = γ_test)

    M_indices = findall(l -> l[1] == "M", labels)
    G_m1_indices = findall(l -> l[1] == "G" && l[2] == "-1", labels)

    ρ₀ = uniform_mixture(labels, M_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_M(t_idx) = combined_pop(pops, M_indices, t_idx)
    P_Gm1(t_idx) = combined_pop(pops, G_m1_indices, t_idx)

    @test P_M(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_Gm1(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_Gm1(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2
    @test P_M(time_index(1)) ≈ exp(-1) rtol = 1e-2

    @test P_Gm1(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2
    @test P_M(time_index(3)) ≈ exp(-3) rtol = 1e-2

    @test P_Gm1(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2
    @test P_M(time_index(5)) ≈ exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation
    M_up = findall(l -> l[1] == "M" && l[end] == "↑", labels)
    M_down = findall(l -> l[1] == "M" && l[end] == "↓", labels)
    Gm1_up = findall(l -> l[1] == "G" && l[2] == "-1" && l[end] == "↑", labels)
    Gm1_down = findall(l -> l[1] == "G" && l[2] == "-1" && l[end] == "↓", labels)

    P_up =
        combined_pop(pops, M_up, time_index(1)) + combined_pop(pops, Gm1_up, time_index(1))
    P_down =
        combined_pop(pops, M_down, time_index(1)) +
        combined_pop(pops, Gm1_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


# =============================================================================
# Test: Â_Σ₊₁ (|M⟩ → |G,+1⟩)
# =============================================================================
@testitem "Dissipator Â_Σ₊₁: |M⟩ → |G,+1⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    # Need Ω₊ to include the +1 spin level
    L, labels =
        make_nv_system(; base_params..., Λ = 0.0, Ω₊ = ConstantDrive(0.0), Σ₊₁ = γ_test)

    M_indices = findall(l -> l[1] == "M", labels)
    G_p1_indices = findall(l -> l[1] == "G" && l[2] == "+1", labels)

    ρ₀ = uniform_mixture(labels, M_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_M(t_idx) = combined_pop(pops, M_indices, t_idx)
    P_Gp1(t_idx) = combined_pop(pops, G_p1_indices, t_idx)

    @test P_M(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_Gp1(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_Gp1(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2
    @test P_M(time_index(1)) ≈ exp(-1) rtol = 1e-2

    @test P_Gp1(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2
    @test P_M(time_index(3)) ≈ exp(-3) rtol = 1e-2

    @test P_Gp1(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2
    @test P_M(time_index(5)) ≈ exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation
    M_up = findall(l -> l[1] == "M" && l[end] == "↑", labels)
    M_down = findall(l -> l[1] == "M" && l[end] == "↓", labels)
    Gp1_up = findall(l -> l[1] == "G" && l[2] == "+1" && l[end] == "↑", labels)
    Gp1_down = findall(l -> l[1] == "G" && l[2] == "+1" && l[end] == "↓", labels)

    P_up =
        combined_pop(pops, M_up, time_index(1)) + combined_pop(pops, Gp1_up, time_index(1))
    P_down =
        combined_pop(pops, M_down, time_index(1)) +
        combined_pop(pops, Gp1_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


# =============================================================================
# Test: Â_γ₊₁ (|G,+1⟩ → |G,0⟩)
# =============================================================================
@testitem "Dissipator Â_γ₊₁: |G,+1⟩ → |G,0⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    # Need Ω₊ to include the +1 spin level (no optical space needed)
    L, labels = make_nv_system(; base_params..., Ω₊ = ConstantDrive(0.0), γ₊₁ = γ_test)

    G_p1_indices = findall(l -> l[1] == "G" && l[2] == "+1", labels)
    G0_indices = findall(l -> l[1] == "G" && l[2] == "0", labels)

    ρ₀ = uniform_mixture(labels, G_p1_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_Gp1(t_idx) = combined_pop(pops, G_p1_indices, t_idx)
    P_G0(t_idx) = combined_pop(pops, G0_indices, t_idx)

    @test P_Gp1(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_G0(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_Gp1(time_index(1)) ≈ exp(-1) rtol = 1e-2
    @test P_G0(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2

    @test P_Gp1(time_index(3)) ≈ exp(-3) rtol = 1e-2
    @test P_G0(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2

    @test P_Gp1(time_index(5)) ≈ exp(-5) rtol = 1e-2
    @test P_G0(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation
    Gp1_up = findall(l -> l[1] == "G" && l[2] == "+1" && l[end] == "↑", labels)
    Gp1_down = findall(l -> l[1] == "G" && l[2] == "+1" && l[end] == "↓", labels)
    G0_up = findall(l -> l[1] == "G" && l[2] == "0" && l[end] == "↑", labels)
    G0_down = findall(l -> l[1] == "G" && l[2] == "0" && l[end] == "↓", labels)

    P_up =
        combined_pop(pops, Gp1_up, time_index(1)) + combined_pop(pops, G0_up, time_index(1))
    P_down =
        combined_pop(pops, Gp1_down, time_index(1)) +
        combined_pop(pops, G0_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


# =============================================================================
# Test: Â_γ₋₁ (|G,-1⟩ → |G,0⟩)
# =============================================================================
@testitem "Dissipator Â_γ₋₁: |G,-1⟩ → |G,0⟩" setup = [DissipatorTest] begin
    using QuantumPropagators: propagate, Newton
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using UnPack: @unpack

    @unpack base_params,
    γ_test,
    τ_test,
    n_τ,
    n_points,
    get_pops,
    uniform_mixture,
    combined_pop,
    time_index = DissipatorTest

    tlist = collect(range(0, n_τ * τ_test; length = n_points))

    # Need Ω₋ to include the -1 spin level (no optical space needed)
    L, labels = make_nv_system(; base_params..., Ω₋ = ConstantDrive(0.0), γ₋₁ = γ_test)

    G_m1_indices = findall(l -> l[1] == "G" && l[2] == "-1", labels)
    G0_indices = findall(l -> l[1] == "G" && l[2] == "0", labels)

    ρ₀ = uniform_mixture(labels, G_m1_indices)

    states = propagate(ρ₀, L, tlist; method = Newton, storage = true)
    pops = get_pops(states)

    P_Gm1(t_idx) = combined_pop(pops, G_m1_indices, t_idx)
    P_G0(t_idx) = combined_pop(pops, G0_indices, t_idx)

    @test P_Gm1(time_index(0)) ≈ 1.0 rtol = 1e-2
    @test P_G0(time_index(0)) ≈ 0.0 atol = 1e-2

    @test P_Gm1(time_index(1)) ≈ exp(-1) rtol = 1e-2
    @test P_G0(time_index(1)) ≈ 1 - exp(-1) rtol = 1e-2

    @test P_Gm1(time_index(3)) ≈ exp(-3) rtol = 1e-2
    @test P_G0(time_index(3)) ≈ 1 - exp(-3) rtol = 1e-2

    @test P_Gm1(time_index(5)) ≈ exp(-5) rtol = 1e-2
    @test P_G0(time_index(5)) ≈ 1 - exp(-5) rtol = 1e-2

    # Verify nuclear spin conservation
    Gm1_up = findall(l -> l[1] == "G" && l[2] == "-1" && l[end] == "↑", labels)
    Gm1_down = findall(l -> l[1] == "G" && l[2] == "-1" && l[end] == "↓", labels)
    G0_up = findall(l -> l[1] == "G" && l[2] == "0" && l[end] == "↑", labels)
    G0_down = findall(l -> l[1] == "G" && l[2] == "0" && l[end] == "↓", labels)

    P_up =
        combined_pop(pops, Gm1_up, time_index(1)) + combined_pop(pops, G0_up, time_index(1))
    P_down =
        combined_pop(pops, Gm1_down, time_index(1)) +
        combined_pop(pops, G0_down, time_index(1))
    @test P_up ≈ P_down rtol = 1e-2
end


using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
