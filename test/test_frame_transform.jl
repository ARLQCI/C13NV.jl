using TestItems


@testitem "R unitarity" begin
    using C13NV.FrameTransformations: frame_transformation_matrix
    using LinearAlgebra: I, norm

    # General case
    R = frame_transformation_matrix(; A_zz = 0.5, A_zx = 0.3, A_zy = 0.1)
    @test norm(R' * R - I) < 1e-14
    @test norm(R * R' - I) < 1e-14

    # Another general case with different values
    R = frame_transformation_matrix(; A_zz = 1.0, A_zx = 0.5, A_zy = 0.2)
    @test norm(R' * R - I) < 1e-14
    @test norm(R * R' - I) < 1e-14
end


@testitem "R diagonalizes hyperfine tensor" begin
    using C13NV.FrameTransformations: frame_transformation_matrix
    using LinearAlgebra: norm, Diagonal

    A_zz, A_zx, A_zy = 0.5, 0.3, 0.1
    A = √(A_zz^2 + A_zx^2 + A_zy^2)

    # Construct Â_I from eq-A-I-matrix in notes/hamiltonian.qmd
    Â_I = 0.5 * [A_zz (A_zx-im*A_zy); (A_zx+im*A_zy) -A_zz]

    R = frame_transformation_matrix(; A_zz, A_zx, A_zy)

    # R† Â_I R should be diagonal with eigenvalues ±A/2
    diagonalized = R' * Â_I * R
    @test norm(diagonalized - Diagonal([A / 2, -A / 2])) < 1e-14
end


@testitem "R degenerate cases" begin
    using C13NV.FrameTransformations: frame_transformation_matrix
    using LinearAlgebra: I, norm

    # Zero hyperfine → identity
    R = frame_transformation_matrix(; A_zz = 0.0, A_zx = 0.0, A_zy = 0.0)
    @test norm(R - I) < 1e-14

    # Already diagonal (A_zz > 0) → identity
    R = frame_transformation_matrix(; A_zz = 1.0, A_zx = 0.0, A_zy = 0.0)
    @test norm(R - I) < 1e-14

    # Already diagonal (A_zz < 0) → swap matrix
    R = frame_transformation_matrix(; A_zz = -1.0, A_zx = 0.0, A_zy = 0.0)
    @test norm(R - [0 1; 1 0]) < 1e-14

    # Near-zero off-diagonal elements should also give identity
    R = frame_transformation_matrix(; A_zz = 1.0, A_zx = 1e-20, A_zy = 1e-20)
    @test norm(R - I) < 1e-12

    # Small but non-negligible off-diagonal should give proper transformation
    R = frame_transformation_matrix(; A_zz = 1.0, A_zx = 0.01, A_zy = 0.0)
    @test norm(R' * R - I) < 1e-10  # Still unitary (relaxed tolerance for small values)
end


@testitem "2-carbon R structure" begin
    using C13NV.FrameTransformations: frame_transformation_matrix
    using LinearAlgebra: kron, norm, I

    tensors =
        [[0.0 0.0 0.3; 0.0 0.0 0.1; 0.3 0.1 0.5], [0.0 0.0 0.2; 0.0 0.0 0.0; 0.2 0.0 0.4],]

    R = frame_transformation_matrix(tensors)

    # Should be R₁ ⊗ R₂
    R₁ = frame_transformation_matrix(; A_zz = 0.5, A_zx = 0.3, A_zy = 0.1)
    R₂ = frame_transformation_matrix(; A_zz = 0.4, A_zx = 0.2, A_zy = 0.0)
    @test norm(R - kron(R₁, R₂)) < 1e-14

    # Should be unitary
    @test norm(R' * R - I) < 1e-14
end


@testitem "3-carbon R structure" begin
    using C13NV.FrameTransformations: frame_transformation_matrix
    using LinearAlgebra: kron, norm, I

    tensors = [
        [0.0 0.0 0.3; 0.0 0.0 0.1; 0.3 0.1 0.5],
        [0.0 0.0 0.2; 0.0 0.0 0.0; 0.2 0.0 0.4],
        [0.0 0.0 0.1; 0.0 0.0 0.05; 0.1 0.05 0.6],
    ]

    R = frame_transformation_matrix(tensors)

    # Should be R₁ ⊗ R₂ ⊗ R₃
    R₁ = frame_transformation_matrix(; A_zz = 0.5, A_zx = 0.3, A_zy = 0.1)
    R₂ = frame_transformation_matrix(; A_zz = 0.4, A_zx = 0.2, A_zy = 0.0)
    R₃ = frame_transformation_matrix(; A_zz = 0.6, A_zx = 0.1, A_zy = 0.05)
    @test norm(R - kron(kron(R₁, R₂), R₃)) < 1e-14

    # Should be unitary (8x8 matrix for 3 carbons)
    @test size(R) == (8, 8)
    @test norm(R' * R - I) < 1e-14
end


@testitem "transform-frame round-trip single carbon" begin
    using C13NV.FrameTransformations: transform_frame
    using LinearAlgebra: norm, normalize

    A_zz, A_zx, A_zy = 0.5, 0.3, 0.0

    # Random state in ℋ_S ⊗ ℋ_I (dim_S=2, N=1 → 4-dimensional)
    Ψ_rwa = normalize(randn(ComplexF64, 4))

    # RWA → diag → RWA should recover original
    Ψ_diag = transform_frame(Ψ_rwa; A_zz, A_zx, A_zy, to_frame = :diag)
    Ψ_back = transform_frame(Ψ_diag; A_zz, A_zx, A_zy, to_frame = :rwa)

    @test norm(Ψ_back - Ψ_rwa) < 1e-14

    # Also test the other direction
    Ψ_diag2 = normalize(randn(ComplexF64, 4))
    Ψ_rwa2 = transform_frame(Ψ_diag2; A_zz, A_zx, A_zy, to_frame = :rwa)
    Ψ_back2 = transform_frame(Ψ_rwa2; A_zz, A_zx, A_zy, to_frame = :diag)

    @test norm(Ψ_back2 - Ψ_diag2) < 1e-14
end


@testitem "transform-frame round-trip multiple carbons" begin
    using C13NV.FrameTransformations: transform_frame
    using LinearAlgebra: norm, normalize

    tensors =
        [[0.0 0.0 0.3; 0.0 0.0 0.1; 0.3 0.1 0.5], [0.0 0.0 0.2; 0.0 0.0 0.0; 0.2 0.0 0.4],]

    # dim_S=2, N=2 → 2 * 4 = 8 dimensional
    Ψ_rwa = normalize(randn(ComplexF64, 8))

    Ψ_diag = transform_frame(Ψ_rwa, tensors; to_frame = :diag)
    Ψ_back = transform_frame(Ψ_diag, tensors; to_frame = :rwa)

    @test norm(Ψ_back - Ψ_rwa) < 1e-14
end


@testitem "transform-frame dimension error" begin
    using C13NV.FrameTransformations: transform_frame

    # Wrong dimension: 5 is not divisible by 2 (nuclear spin dimension for N=1)
    Ψ_wrong = zeros(ComplexF64, 5)

    @test_throws ErrorException transform_frame(
        Ψ_wrong;
        A_zz = 0.5,
        A_zx = 0.3,
        to_frame = :diag
    )
end


@testitem "transform_frame invalid frame error" begin
    using C13NV.FrameTransformations: transform_frame

    Ψ = zeros(ComplexF64, 4)

    @test_throws ErrorException transform_frame(
        Ψ;
        A_zz = 0.5,
        A_zx = 0.3,
        to_frame = :invalid
    )
end


@testitem "transform-frame consistent with direct matrix transformation" begin
    using C13NV.FrameTransformations: frame_transformation_matrix, transform_frame
    using LinearAlgebra: kron, norm, normalize, I

    # Using random test vectors

    A_zz, A_zx = 0.5, 0.3

    # Get the R matrix
    R_I = frame_transformation_matrix(; A_zz, A_zx)

    # Random state
    Ψ_rwa = normalize(randn(ComplexF64, 4))

    # Transform using the function
    Ψ_diag_func = transform_frame(Ψ_rwa; A_zz, A_zx, to_frame = :diag)

    # Transform using direct matrix multiplication: (𝟙_S ⊗ R†) Ψ
    𝟙_S = Matrix{ComplexF64}(I, 2, 2)
    R_full = kron(𝟙_S, R_I')
    Ψ_diag_direct = R_full * Ψ_rwa

    @test norm(Ψ_diag_func - Ψ_diag_direct) < 1e-14
end


@testitem "transform_frame no change when already diagonal" begin
    using C13NV.FrameTransformations: transform_frame
    using LinearAlgebra: norm, normalize

    # Using random test vectors

    # When A_zx = A_zy = 0, the hyperfine is already diagonal
    # so the transformation should be identity
    Ψ = normalize(randn(ComplexF64, 4))

    Ψ_diag = transform_frame(Ψ; A_zz = 1.0, A_zx = 0.0, A_zy = 0.0, to_frame = :diag)

    @test norm(Ψ_diag - Ψ) < 1e-14
end


@testitem "transform_frame changes state when hyperfine has off-diagonal" begin
    using C13NV.FrameTransformations: transform_frame
    using LinearAlgebra: norm, normalize

    # Using random test vectors

    # When A_zx ≠ 0, the transformation should change the state
    Ψ = normalize(randn(ComplexF64, 4))

    Ψ_diag = transform_frame(Ψ; A_zz = 0.5, A_zx = 0.3, to_frame = :diag)

    # State should be different (but same norm)
    @test norm(Ψ_diag - Ψ) > 0.01  # Should be noticeably different
    @test abs(norm(Ψ_diag) - 1.0) < 1e-14  # But still normalized
end


using TestItemRunner
@run_package_tests filter = ti -> ti.filename == @__FILE__
