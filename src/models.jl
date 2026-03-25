module Models

import LinearAlgebra
using LinearAlgebra: norm
using LinearAlgebra: Diagonal
using QuantumPropagators: hamiltonian, liouvillian
using QuantumPropagators.Generators: dissipator, Generator
using ..Units: kHz, MHz, Gauss

const 𝕚 = 1im

function ⊗(A, B)
    LinearAlgebra.kron(A, B)
end


function ZeroMatrix(T, N)
    return Array{T}(zeros(T, N) * zeros(T, N)')
end

ZeroMatrix(N) = ZeroMatrix(Float64, N)

function IdMatrix(T, N)
    return Array{T}(LinearAlgebra.I(N))
end

IdMatrix(N) = IdMatrix(Float64, N)


function ⊕(A::T, B::T) where {T<:Array}
    N = size(A, 1)
    M = size(B, 1)
    @assert size(A, 2) == N
    @assert size(B, 2) == M
    C = ZeroMatrix(eltype(T), N + M)
    C[1:N, 1:N] .= A[1:N, 1:N]
    C[(N+1):(N+M), (N+1):(N+M)] .= B[1:M, 1:M]
    return C
end


function ⊕(A::TA, B::TB) where {TA<:AbstractMatrix,TB<:AbstractMatrix}
    ET = promote_type(eltype(TA), eltype(TB))
    if !(A isa Array)
        A = convert(Array{ET}, A)
    end
    if !(B isa Array)
        B = convert(Array{ET}, B)
    end
    return ⊕(A, B)
end


"""Construct an eigenstate

```julia
Ψ = ket(label, labels; strict = true)
```

constructs the canonical basis state identified by `label` in a Hilbert space
with a basis defined by `labels`. If `label` is not an element of `labels`,
return a zero-vector if `strict = false`, and error if `strict = true`
(default). The use of `strict = false` can be useful when working in a
truncated Hilbert space.
"""
function ket(label, labels; strict = true)
    N = length(labels)
    i = findfirst(isequal(label), labels)
    Ψ = zeros(ComplexF64, N)
    if isnothing(i) && !strict
        return Ψ
    else
        Ψ[i] = one(ComplexF64)
        return Ψ
    end
end

function bra(label, labels; strict = true)
    return adjoint(ket(label, labels; strict))
end


function ketbra(l1, l2, labels; strict = true)
    return ket(l1, labels; strict) * bra(l2, labels; strict)
end


# Tensor to ``𝟙 ⊗ … ⊗ Â ⊗ … ⊗ 𝟙`` where Â is in the `n`'th position of `N`
function _n(A, n, N)
    @assert size(A) == (size(A, 1), size(A, 2))
    @assert 1 ≤ n ≤ N
    𝟙 = Diagonal(ones(ComplexF64, size(A, 1)))
    result = (n == 1) ? A : 𝟙
    for n′ = 2:N
        result = result ⊗ ((n′ == n) ? A : 𝟙)
    end
    return result
end


# see `frame_transformation_matrix` in `FrameTransformations` for public
# function
function _frame_transformation_matrix(A_zz::Real, A_zx::Real, A_zy::Real)
    A = √(A_zz^2 + A_zx^2 + A_zy^2)
    if A ≈ 0
        return ComplexF64[1 0; 0 1]
    end
    # Handle degenerate case where off-diagonal hyperfine elements are negligible
    # relative to A_zz (i.e., matrix is effectively already diagonal)
    off_diag_norm = √(A_zx^2 + A_zy^2)
    if off_diag_norm < 1e-14 * abs(A_zz) || off_diag_norm < 1e-30
        if A_zz ≥ 0
            return ComplexF64[1 0; 0 1]
        else
            return ComplexF64[0 1; 1 0]
        end
    end
    # General case from eq-R in notes/hamiltonian.qmd
    denom_plus = √(A^2 + A * A_zz)
    denom_minus = √(A^2 - A * A_zz)
    R =
        (1 / √2) * ComplexF64[
            (A_zz+A)/denom_plus        (A_zz-A)/denom_minus
            (A_zx+𝕚*A_zy)/denom_plus   (A_zx+𝕚*A_zy)/denom_minus
        ]
    return R
end


"""Construct the system.

```julia
H_or_L, labels = make_nv_system(; kwargs...)
```

return a generator, and a list of labels, each label a tuple of strings.

# Keyword Arguments

* `hyperfine_tensors`: Required if `A_zz`, `A_zx`, `A_zy` are not given. A
  vector of 3×3  matrices. The number length of the vector, i.e., the number of
  matrices given determines the number of C-13 atoms in the model. Defaults to
  a single hyperfine tensor with elements `A_zz`, `A_zx`, `A_zy`.
* `A_zz`, `A_zx`, `A_zy`: Required if `hyperfine_tensors` is not given, the
  strength of the hyperfine coupling along the 3 spatial axes relevant in the
  RWA, for a single C-13 atom.
* `B`: Required. The magnitude of the magnetic field
* `ω₊ = nothing`: The control ``ω_{+}(t) ≡ ∂ϕ_{+}(t)/∂t``
* `ω₋ = nothing`: The control ``ω_{-}(t) ≡ ∂ϕ_{-}(t)/∂t``
* `Ω₊ = nothing`: The control ``Ω_{+}(t)``
* `Ω₋ = nothing`: The control ``Ω_{-}(t)``
* `μ = 1.0`: The reduction factor for ``Ω_{±}(t)``.
* `Λ = nothing`: The time-dependent optical drive. If given, implies the use of
  the full optical Hilbert space.
* `θ = 0.0`: The azimuthal angle of the magnetic field
* `ϕ = 0.0`: The polar angle of the magnetic field
* `γ_c = 1.07kHz/Gauss`: The C-13 nuclear gyromagnetic ratio
* `δ₋ = 0.0`: The detuning of the ``|0⟩ ↔ |-1⟩`` transition. Defined as
  ``δ_{-} = D - B γ_e cos(θ) - ω{-}``
* `δ₊ = 0.0`: The detuning of the ``|0⟩ ↔ |+1⟩`` transition. Defined as
  ``δ_{+} = D - B γ_e cos(θ) - ω{+}``
* `Γ = 0.0`: Rate for the ``|E⟩ → |G⟩`` spontaneous decay
* `Γ₀ = 0.0`: Rate for the ``|E,0⟩ → |M⟩`` spontaneous decay
* `Γ₊₁ = 0.0`: Rate for the ``|E,+1⟩ → |M⟩`` spontaneous decay
* `Γ₋₁ = 0.0`: Rate for the ``|E,-1⟩ → |M⟩`` spontaneous decay
* `Σ₀ = 0.0`: Rate for the ``|M⟩ → |G,0⟩`` spontaneous decay
* `Σ₊₁ = 0.0`: Rate for the ``|M⟩ → |G,+1⟩`` spontaneous decay
* `Σ₋₁ = 0.0`: Rate for the ``|M⟩ → |G,-1⟩`` spontaneous decay
* `γ₊₁ = 0.0`: Rate for the ``|G,+1⟩ → |G,0⟩`` spontaneous decay
* `γ₋₁ = 0.0`: Rate for the ``|G,-1⟩ → |G,0⟩`` spontaneous decay
"""
function make_nv_system(;
    A_zz::Float64,
    A_zx::Float64,
    A_zy::Float64 = 0.0,
    hyperfine_tensors::Vector{Matrix{Float64}} = [[
         0     0    A_zx
         0     0    A_zy
        A_zx  A_zy  A_zz
    ],],
    B::Float64,
    θ::Float64 = 0.0,
    ϕ::Float64 = 0.0,
    γ_c::Float64 = 1.07kHz/Gauss,
    δ₋::Float64 = 0.0,
    δ₊::Float64 = 0.0,
    ω₊ = nothing,
    ω₋ = nothing,
    Ω₊ = nothing,
    Ω₋ = nothing,
    μ = 1.0,
    Λ = nothing, # incoherent optical excitation (proportional to laser power)
    Γ::Float64 = 0.0,
    Γ₀::Float64 = 0.0,
    Γ₊₁::Float64 = 0.0,
    Γ₋₁::Float64 = 0.0,
    Σ₀::Float64 = 0.0,
    Σ₊₁::Float64 = 0.0,
    Σ₋₁::Float64 = 0.0,
    γ₊₁::Float64 = 0.0,
    γ₋₁::Float64 = 0.0,
    frame::Symbol = :rwa,  # or :diag
    # The `frame` argument is undocumented: this turns out to be a dead end,
    # and not really useful. It just felt bed to remove the code entirely, in
    # case we ever come back to it.
)

    for value in (Γ, Γ₀, Γ₊₁, Γ₋₁, Σ₀, Σ₊₁, Σ₋₁, γ₊₁, γ₋₁)
        @assert value ≥ 0.0 "Dissipation rates must be ≥ 0"
    end

    use_full_optical_space = !isnothing(Λ)
    has_dissipation = !isnothing(Λ) || (γ₋₁ > 0.0) || (γ₊₁ > 0.0)

    B⃗ = B .* [sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ)]
    B_x = B⃗[1]
    B_y = B⃗[2]
    B_z = B⃗[3]

    N = length(hyperfine_tensors)  # number of carbons

    labels_S = String["0"]
    if !isnothing(Ω₊)
        pushfirst!(labels_S, "+1")  #  before |0⟩
    end
    if !isnothing(Ω₋)
        push!(labels_S, "-1")  # after |0⟩
    end

    labels_OS = Tuple{String,Vararg{String}}[("G", label_S) for label_S in labels_S]
    if use_full_optical_space
        append!(labels_OS, [("E", label_S) for label_S in labels_S])
        push!(labels_OS, ("M",))
    end

    labels_I_n = ["↑", "↓"]
    labels_I = vec(
        collect(join(reverse(t)) for t in Iterators.product(ntuple(_ -> labels_I_n, N)...))
    )

    labels = [(label_OS..., label_I) for label_OS in labels_OS for label_I in labels_I]

    # Trivial identities
    𝟙_G = 1
    𝟙_E = 1
    𝟙_M = 1
    𝟘_M = 0

    # Operators for a single nuclear spin, ℋ_{I}^{(n)}

    Î_x = 0.5 * ComplexF64[0 1; 1 0]
    Î_y = 0.5 * ComplexF64[0 -1im; 1im 0]
    Î_z = 0.5 * ComplexF64[1 0; 0 -1]

    # Operators for electronic spin, ℋ_S

    𝟙_S = Diagonal(ones(ComplexF64, length(labels_S)))
    𝟘_S = Diagonal(zeros(ComplexF64, length(labels_S)))

    # Using `strict = false` to make "truncation" easy (depending on the
    # presence of Ω₊/Ω₋, see above)
    Ŝ_z =
        ketbra("+1", "+1", labels_S; strict = false) -
        ketbra("-1", "-1", labels_S; strict = false)

    # Operators for tensored nuclear spins, ℋ_I

    𝟙_I = _n(Diagonal(ones(ComplexF64, length(labels_I_n))), 1, N)
    @assert size(𝟙_I) == (length(labels_I), length(labels_I))
    𝟘_I = _n(Diagonal(zeros(ComplexF64, length(labels_I_n))), 1, N)
    @assert size(𝟘_I) == (length(labels_I), length(labels_I))
    @assert norm(𝟘_I) == 0.0

    # Drift Hamiltonian ∈ ℋ_S ⊗ ℋ_I

    δ̂ = (
        δ₊ * ketbra("+1", "+1", labels_S; strict = false) +
        δ₋ * ketbra("-1", "-1", labels_S; strict = false)
    )
    Ĥ_0 = δ̂ ⊗ 𝟙_I
    for n = 1:N
        B̂_I_n = B_x * Î_x + B_y * Î_y + B_z * Î_z
        @assert norm(
            # just to check that the math is consistent
            B̂_I_n -
            (B / 2) * [
                            cos(θ)                   (sin(θ)*cos(ϕ)-𝕚*sin(θ)*sin(ϕ))
                (sin(θ)*cos(ϕ)+𝕚*sin(θ)*sin(ϕ))                        -cos(θ)
            ]
        ) < 1e-14
        A_zz = hyperfine_tensors[n][3, 3]
        A_zx = hyperfine_tensors[n][3, 1]
        A_zy = hyperfine_tensors[n][3, 2]
        if frame == :rwa
            Ĥ_0 = Ĥ_0 - γ_c * (𝟙_S ⊗ _n(B̂_I_n, n, N))
            Â_I_n = 0.5 * [
                    (A_zz)    (A_zx-𝕚*A_zy)
                (A_zx+𝕚*A_zy)     (-A_zz)
            ]
            Ĥ_0 = Ĥ_0 + Ŝ_z ⊗ _n(Â_I_n, n, N)
        elseif frame == :diag
            R̂ = _frame_transformation_matrix(A_zz, A_zx, A_zy)
            A_n = √(A_zz^2 + A_zx^2 + A_zy^2)
            B̂_I_n_diag = R̂' * B̂_I_n * R̂
            if θ == ϕ == 0.0
                # just to check that the math is consistent
                B̂_I_n_diag_expected =
                    (B / (2*A_n)) * [
                                A_zz         -√(A_zx^2 + A_zy^2)
                        -√(A_zx^2 + A_zy^2)         -A_zz
                    ]
                diff = norm(B̂_I_n_diag - B̂_I_n_diag_expected)
                @assert diff < 1e-12 "Unexpected B in diagonal frame: error $diff"
            end
            Ĥ_0 = Ĥ_0 - γ_c * (𝟙_S ⊗ _n(B̂_I_n_diag, n, N))
            Ĥ_0 = Ĥ_0 + A_n * (Ŝ_z ⊗ _n(Î_z, n, N))
        else
            error("`frame` must be one of :rwa, :diag, not $repr(frame)")
        end
    end

    # Control Hamiltonians ∈ ℋ_S ⊗ ℋ_I (to match Hilbert space of Ĥ_0)

    Ĥ_ω₊ = -1.0 * ketbra("+1", "+1", labels_S; strict = false) ⊗ 𝟙_I
    Ĥ_ω₋ = -1.0 * ketbra("-1", "-1", labels_S; strict = false) ⊗ 𝟙_I
    Ĥ_Ω₊ =
        (μ / 2) * (
            ketbra("+1", "0", labels_S; strict = false) +
            ketbra("0", "+1", labels_S; strict = false)
        ) ⊗ 𝟙_I
    Ĥ_Ω₋ =
        (μ / 2) * (
            ketbra("-1", "0", labels_S; strict = false) +
            ketbra("0", "-1", labels_S; strict = false)
        ) ⊗ 𝟙_I


    if use_full_optical_space
        padding = (𝟙_E ⊗ 𝟘_S ⊗ 𝟘_I) ⊕ (𝟙_M ⊗ 𝟘_I)
        Ĥ_0 = (𝟙_G ⊗ Ĥ_0) ⊕ padding
        Ĥ_ω₋ = (𝟙_G ⊗ Ĥ_ω₋) ⊕ padding
        Ĥ_ω₊ = (𝟙_G ⊗ Ĥ_ω₊) ⊕ padding
        Ĥ_Ω₋ = (𝟙_G ⊗ Ĥ_Ω₋) ⊕ padding
        Ĥ_Ω₊ = (𝟙_G ⊗ Ĥ_Ω₊) ⊕ padding
    end

    parts = Any[Ĥ_0,]
    if !isnothing(ω₋)
        push!(parts, [Ĥ_ω₋, ω₋])
    end
    if !isnothing(ω₊)
        push!(parts, [Ĥ_ω₊, ω₊])
    end
    if !isnothing(Ω₋)
        push!(parts, [Ĥ_Ω₋, Ω₋])
    end
    if !isnothing(Ω₊)
        push!(parts, [Ĥ_Ω₊, Ω₊])
    end
    H = hamiltonian(parts...)

    if has_dissipation

        c_ops = Any[]
        if (Γ > 0) && use_full_optical_space
            Â_Γ = √Γ * ketbra("G", "E", ["G", "E"]) ⊗ 𝟙_S ⊗ 𝟙_I ⊕ 𝟘_M ⊗ 𝟙_I
            push!(c_ops, Â_Γ)
        end
        if (Γ₀ > 0) && use_full_optical_space
            Â_Γ₀ = √Γ₀ * ketbra(("M",), ("E", "0"), labels_OS) ⊗ 𝟙_I
            push!(c_ops, Â_Γ₀)
        end
        if (Γ₋₁ > 0) && use_full_optical_space && (!isnothing(Ω₋))
            Â_Γ₋₁ = √Γ₋₁ * ketbra(("M",), ("E", "-1"), labels_OS) ⊗ 𝟙_I
            push!(c_ops, Â_Γ₋₁)
        end
        if (Γ₊₁ > 0) && use_full_optical_space && (!isnothing(Ω₊))
            Â_Γ₊₁ = √Γ₊₁ * ketbra(("M",), ("E", "+1"), labels_OS) ⊗ 𝟙_I
            push!(c_ops, Â_Γ₊₁)
        end
        if (Σ₀ > 0) && use_full_optical_space
            Â_Σ₀ = √Σ₀ * ketbra(("G", "0"), ("M",), labels_OS) ⊗ 𝟙_I
            push!(c_ops, Â_Σ₀)
        end
        if (Σ₋₁ > 0) && use_full_optical_space && (!isnothing(Ω₋))
            Â_Σ₋₁ = √Σ₋₁ * ketbra(("G", "-1"), ("M",), labels_OS) ⊗ 𝟙_I
            push!(c_ops, Â_Σ₋₁)
        end
        if (Σ₊₁ > 0) && use_full_optical_space && (!isnothing(Ω₊))
            Â_Σ₊₁ = √Σ₊₁ * ketbra(("G", "+1"), ("M",), labels_OS) ⊗ 𝟙_I
            push!(c_ops, Â_Σ₊₁)
        end
        if (γ₊₁ > 0) && (!isnothing(Ω₊))
            if use_full_optical_space
                Â_γ₊₁ = √γ₊₁ * ketbra(("G", "0"), ("G", "+1"), labels_OS) ⊗ 𝟙_I
            else
                Â_γ₊₁ = √γ₊₁ * ketbra("0", "+1", labels_S) ⊗ 𝟙_I
            end
            push!(c_ops, Â_γ₊₁)
        end
        if (γ₋₁ > 0) && (!isnothing(Ω₋))
            if use_full_optical_space
                Â_γ₋₁ = √γ₋₁ * ketbra(("G", "0"), ("G", "-1"), labels_OS) ⊗ 𝟙_I
            else
                Â_γ₋₁ = √γ₋₁ * ketbra("0", "-1", labels_S) ⊗ 𝟙_I
            end
            push!(c_ops, Â_γ₋₁)
        end

        L = liouvillian(H, c_ops; convention = :TDSE)

        # Add time-dependent dissipative drive

        if !isnothing(Λ)
            L_Λ = liouvillian(
                nothing,
                [ketbra("E", "G", ["G", "E"]) ⊗ 𝟙_S ⊗ 𝟙_I ⊕ 𝟘_M ⊗ 𝟙_I];
                convention = :TDSE
            )
            if L isa Generator
                L = Generator([L.ops..., L_Λ], [L.amplitudes..., Λ])
            else
                # The original `L` may not have been time-dependent
                L = Generator([L, L_Λ], [Λ,])
            end
        end

        return L, labels

    else

        return H, labels

    end

end

end
