module FrameTransformations

# Note: The `:diag` frame didn't work out as expected, since it doesn't
# actually diagonalize. This code is here for "legacy", and is
# private/undocumented

using ..Models: IdMatrix, ⊗, _frame_transformation_matrix

const 𝕚 = 1im


"""Construct the transformation matrix between the `:diag` and `:rwa` frames.

```julia
R⁽ⁿ⁾ = frame_transformation_matrix(; A_zz, A_zx, A_zy=0.0)
```

constructs the unitary matrix ``R`` that transforms states in the nuclear spin
subspace between the RWA frame and the diagonal frame for a single carbon. Or,
for multiple carbons:

```julia
R = frame_transformation_matrix(hyperfine_tensors)
```

The matrix ``R^{(n)}`` for a single carbon satisfies:

```math
Â_I^{(n)} = R^{(n)} \\begin{pmatrix} A/2 & 0 \\\\ 0 & -A/2 \\end{pmatrix} R^{(n)†}
```

where ``A = \\sqrt{A_{zz}^2 + A_{zx}^2 + A_{zy}^2}``.

For multiple nuclear spins, ``R = R^{(1)} ⊗ R^{(2)} ⊗ … ⊗ R^{(N)}``.

# Arguments

Either pass:
- `hyperfine_tensors`: Vector of 3×3 hyperfine tensor matrices

Or use keyword arguments for a single carbon:
- `A_zz`, `A_zx`, `A_zy`: Hyperfine coupling components

# Returns

A complex matrix ``R̂ ∈ ℋ_I`` of size ``2^N × 2^N`` where ``N`` is the number of
carbons (length of `hyperfine_tensors`).

# See also

* [`transform_frame`](@ref) to apply ``R̂`` to a state.
"""
function frame_transformation_matrix(hyperfine_tensors::Vector{Matrix{Float64}})
    N = length(hyperfine_tensors)
    R_list = [
        _frame_transformation_matrix(
            hyperfine_tensors[n][3, 3],
            hyperfine_tensors[n][3, 1],
            hyperfine_tensors[n][3, 2],
        ) for n = 1:N
    ]
    R = R_list[1]
    for n = 2:N
        R = R ⊗ R_list[n]
    end
    return R
end

function frame_transformation_matrix(; A_zz::Float64, A_zx::Float64, A_zy::Float64 = 0.0)
    hyperfine_tensors = [[
        0 0 A_zx
        0 0 A_zy
        A_zx A_zy A_zz
    ],]
    return frame_transformation_matrix(hyperfine_tensors)
end


"""Transform a state between the `:rwa` and `:diag` frames.

```julia
Ψ′ = transform_frame(Ψ, hyperfine_tensors; to_frame)
Ψ′ = transform_frame(Ψ; A_zz, A_zx, A_zy=0.0, to_frame)
```

transforms a Hilbert space vector between the rotating wave approximation (RWA)
frame and the diagonal frame where the hyperfine interaction is diagonalized.

# Arguments

- `Ψ`: The state vector to transform
- `hyperfine_tensors`: Vector of 3×3 hyperfine tensor matrices (or use keyword
  arguments `A_zz`, `A_zx`, `A_zy` for a single carbon)

# Keyword Arguments

- `to_frame`: Target frame, either `:rwa` or `:diag`

# Returns

The transformed state vector.

The transformation uses the unitary ``R = R^{(1)} ⊗ … ⊗ R^{(N)}`` that
diagonalizes the hyperfine interaction in the nuclear spin subspace:

- RWA → diag: ``|Ψ⟩_{\\text{diag}} = (𝟙 ⊗ R^†) |Ψ⟩_{\\text{RWA}}``
- diag → RWA: ``|Ψ⟩_{\\text{RWA}} = (𝟙 ⊗ R) |Ψ⟩_{\\text{diag}}``

# See also

* [`frame_transformation_matrix`](@ref) to obtain the transformation matrix ``R``.
"""
function transform_frame(
    Ψ::Vector{ComplexF64},
    hyperfine_tensors::Vector{Matrix{Float64}};
    to_frame::Symbol
)
    if to_frame ∉ (:rwa, :diag)
        error("`to_frame` must be :rwa or :diag, not $(repr(to_frame))")
    end

    R_I = frame_transformation_matrix(hyperfine_tensors)
    dim_I = size(R_I, 1)

    dim_Ψ = length(Ψ)
    dim_OS, remainder = divrem(dim_Ψ, dim_I)
    if remainder != 0
        error(
            "State vector dimension $dim_Ψ is not divisible by nuclear spin " *
            "dimension $dim_I (N_carbons=$(length(hyperfine_tensors)))"
        )
    end

    𝟙_OS = IdMatrix(ComplexF64, dim_OS)

    if to_frame == :diag
        # RWA → diag: Ψ_diag = (𝟙 ⊗ R†) Ψ_rwa
        R_full = 𝟙_OS ⊗ R_I'
    else
        # diag → RWA: Ψ_rwa = (𝟙 ⊗ R) Ψ_diag
        R_full = 𝟙_OS ⊗ R_I
    end

    return R_full * Ψ
end


function transform_frame(
    Ψ::Vector{ComplexF64};
    A_zz::Float64,
    A_zx::Float64,
    A_zy::Float64 = 0.0,
    to_frame::Symbol
)
    hyperfine_tensors = [[
        0 0 A_zx
        0 0 A_zy
        A_zx A_zy A_zz
    ],]
    return transform_frame(Ψ, hyperfine_tensors; to_frame)
end

end
