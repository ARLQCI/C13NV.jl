using TestItems

@testmodule HyperfineCouplingTest begin
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using C13NV.Units: MHz, kHz, Gauss, ns

    # Base parameters for testing (from defaults but explicit)
    const test_params = (
        B = 120Gauss,
        γ_c = 1.07kHz / Gauss,
        δ₋ = 0.0,
        δ₊ = 0.0,
        Γ = (1 / (12ns)),
        Γ₀ = 0.0,
        Γ₊₁ = (1 / (24ns + 0.9ns)),
        Γ₋₁ = (1 / (24ns + 0.9ns)),
        Σ₀ = (1 / (219ns)),
        Σ₊₁ = (1 / (219ns)),
        Σ₋₁ = (1 / (219ns))
    )
end

@testitem "A_zz hyperfine coupling parameter" setup = [HyperfineCouplingTest] begin
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using C13NV.Units: MHz, kHz

    # Test different A_zz values while keeping A_zx constant
    A_zx_fixed = 0.3MHz
    Ω₊_test = ConstantDrive(1.0kHz)  # Small amplitude to focus on hyperfine structure

    # Test case 1: A_zz = 0 (no diagonal hyperfine coupling)
    H1, labels1 = make_nv_system(;
        A_zz = 0.0MHz,
        A_zx = A_zx_fixed,
        Ω₊ = Ω₊_test,
        HyperfineCouplingTest.test_params...
    )
    H1_matrix = H1.ops[1]  # Get static Hamiltonian part

    # Test case 2: Default A_zz = 1.0MHz
    H2, labels2 = make_nv_system(;
        A_zz = 1.0MHz,
        A_zx = A_zx_fixed,
        Ω₊ = Ω₊_test,
        HyperfineCouplingTest.test_params...
    )
    H2_matrix = H2.ops[1]

    # Test case 3: Large A_zz = 3.0MHz
    H3, labels3 = make_nv_system(;
        A_zz = 3.0MHz,
        A_zx = A_zx_fixed,
        Ω₊ = Ω₊_test,
        HyperfineCouplingTest.test_params...
    )
    H3_matrix = H3.ops[1]

    # Verify that A_zz affects diagonal elements as expected
    # From the Hamiltonian matrix structure:
    # H[1,1] = A_zz/2 - B*γ_c/2 + δ₊  (for |+1,↑⟩ state)
    # H[2,2] = -A_zz/2 + B*γ_c/2 + δ₊ (for |+1,↓⟩ state)
    # So H[1,1] - H[2,2] = A_zz - B*γ_c = A_zz - (120*1.07) = A_zz - 128.4kHz

    # The energy difference between spin-up and spin-down states should scale with A_zz
    energy_diff_1 = H1_matrix[1, 1] - H1_matrix[2, 2]
    energy_diff_2 = H2_matrix[1, 1] - H2_matrix[2, 2]
    energy_diff_3 = H3_matrix[1, 1] - H3_matrix[2, 2]

    # Expected differences: A_zz - B*γ_c (where B*γ_c ≈ 128.4kHz from defaults)
    B_gamma_term = 120 * 1.07  # From defaults: B=120Gauss, γ_c=1.07kHz/Gauss
    @test energy_diff_1 ≈ (0.0MHz - B_gamma_term * kHz) atol = 1e-6
    @test energy_diff_2 ≈ (1.0MHz - B_gamma_term * kHz) atol = 1e-6
    @test energy_diff_3 ≈ (3.0MHz - B_gamma_term * kHz) atol = 1e-6

    # Test that off-diagonal elements remain unchanged (only A_zx affects these)
    @test H1_matrix[1, 2] ≈ H2_matrix[1, 2] ≈ H3_matrix[1, 2]
    @test H1_matrix[2, 1] ≈ H2_matrix[2, 1] ≈ H3_matrix[2, 1]

    # All systems should have same labels and dimensions
    @test labels1 == labels2 == labels3
    @test size(H1_matrix) == size(H2_matrix) == size(H3_matrix)
end


@testitem "A_zx hyperfine coupling parameter" setup = [HyperfineCouplingTest] begin
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using C13NV.Units: MHz, kHz

    # Test different A_zx values while keeping A_zz constant
    A_zz_fixed = 1.0MHz
    Ω₊_test = ConstantDrive(1.0kHz)

    # Test case 1: A_zx = 0 (no off-diagonal hyperfine coupling)
    H1, labels1 = make_nv_system(;
        A_zz = A_zz_fixed,
        A_zx = 0.0MHz,
        Ω₊ = Ω₊_test,
        HyperfineCouplingTest.test_params...
    )
    H1_matrix = H1.ops[1]

    # Test case 2: Default A_zx = 0.3MHz
    H2, labels2 = make_nv_system(;
        A_zz = A_zz_fixed,
        A_zx = 0.3MHz,
        Ω₊ = Ω₊_test,
        HyperfineCouplingTest.test_params...
    )
    H2_matrix = H2.ops[1]

    # Test case 3: Large A_zx = 1.0MHz
    H3, labels3 = make_nv_system(;
        A_zz = A_zz_fixed,
        A_zx = 1.0MHz,
        Ω₊ = Ω₊_test,
        HyperfineCouplingTest.test_params...
    )
    H3_matrix = H3.ops[1]

    # Verify that A_zx affects off-diagonal elements as expected
    # From Hamiltonian: H[1,2] = H[2,1] = A_zx/2
    @test abs(H1_matrix[1, 2]) ≈ 0.0 atol = 1e-6        # Should be 0 when A_zx=0
    @test abs(H2_matrix[1, 2]) ≈ 0.15MHz atol = 1e-6     # Should be A_zx/2 = 0.3/2 = 0.15
    @test abs(H3_matrix[1, 2]) ≈ 0.5MHz atol = 1e-6      # Should be A_zx/2 = 1.0/2 = 0.5

    # Test Hermiticity: H[1,2] should equal conj(H[2,1])
    @test H1_matrix[1, 2] ≈ conj(H1_matrix[2, 1])
    @test H2_matrix[1, 2] ≈ conj(H2_matrix[2, 1])
    @test H3_matrix[1, 2] ≈ conj(H3_matrix[2, 1])

    # Test that diagonal elements remain unchanged (only A_zz affects these)
    @test H1_matrix[1, 1] ≈ H2_matrix[1, 1] ≈ H3_matrix[1, 1]
    @test H1_matrix[2, 2] ≈ H2_matrix[2, 2] ≈ H3_matrix[2, 2]

    # All systems should have same labels and dimensions
    @test labels1 == labels2 == labels3
    @test size(H1_matrix) == size(H2_matrix) == size(H3_matrix)
end


@testitem "Combined A_zz and A_zx effects" setup = [HyperfineCouplingTest] begin
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using C13NV.Units: MHz, kHz

    Ω₊_test = ConstantDrive(1.0kHz)

    # Test that A_zz and A_zx have independent effects
    H_ref, _ = make_nv_system(;
        A_zz = 1.0MHz,
        A_zx = 0.3MHz,
        Ω₊ = Ω₊_test,
        HyperfineCouplingTest.test_params...
    )
    H_ref_matrix = H_ref.ops[1]

    H_modified, _ = make_nv_system(;
        A_zz = 2.0MHz,     # Changed from 1.0 to 2.0
        A_zx = 0.6MHz,     # Changed from 0.3 to 0.6
        Ω₊ = Ω₊_test,
        HyperfineCouplingTest.test_params...
    )
    H_modified_matrix = H_modified.ops[1]

    # Check that the changes are as expected
    # Diagonal difference should reflect A_zz change (1.0MHz increase)
    diag_diff =
        (H_modified_matrix[1, 1] - H_modified_matrix[2, 2]) -
        (H_ref_matrix[1, 1] - H_ref_matrix[2, 2])
    @test abs(diag_diff) ≈ 1.0MHz atol = 1e-6

    # Off-diagonal difference should reflect A_zx change (0.3MHz increase → 0.15MHz in matrix element)
    offdiag_diff = abs(H_modified_matrix[1, 2]) - abs(H_ref_matrix[1, 2])
    @test abs(offdiag_diff) ≈ 0.15MHz atol = 1e-6
end


@testitem "Physical constraints on hyperfine parameters" setup = [HyperfineCouplingTest] begin
    using C13NV.Models: make_nv_system
    using C13NV.Amplitudes: ConstantDrive
    using C13NV.Units: MHz, kHz
    using LinearAlgebra: eigvals, Hermitian, norm

    Ω₊_test = ConstantDrive(1.0kHz)

    # Test that Hamiltonians are Hermitian for all parameter values
    test_params = [
        (A_zz = 0.5MHz, A_zx = 0.1MHz),
        (A_zz = 2.0MHz, A_zx = 0.8MHz),
        (A_zz = -0.5MHz, A_zx = 0.2MHz),  # Negative A_zz should work
        (A_zz = 1.0MHz, A_zx = -0.3MHz)   # Negative A_zx should work
    ]

    for params in test_params
        H, _ = make_nv_system(;
            A_zz = params.A_zz,
            A_zx = params.A_zx,
            Ω₊ = Ω₊_test,
            HyperfineCouplingTest.test_params...
        )
        H_matrix = H.ops[1]

        # Test Hermiticity
        @test norm(H_matrix - H_matrix') < 1e-12

        # Test that eigenvalues are real (consequence of Hermiticity)
        evals = eigvals(Hermitian(H_matrix))
        @test all(abs.(imag.(evals)) .< 1e-12)
    end

end

using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
