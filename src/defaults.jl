module Defaults

using ..Units: MHz, kHz, Gauss, ns

const A_zz::Float64 = 1.0MHz
const A_zx::Float64 = 0.3MHz
const B::Float64 = 120Gauss
const γ_c::Float64 = 1.07kHz / Gauss
const δ₋::Float64 = 0.0
const δ₊::Float64 = 0.0
const Γ::Float64 = (1 / (12ns))
const Γ₀::Float64 = 0.0
const Γ₊₁::Float64 = (1 / (24ns + 0.9ns))
const Γ₋₁::Float64 = (1 / (24ns + 0.9ns))
const Σ₀::Float64 = (1 / (219ns))
const Σ₊₁::Float64 = (1 / (219ns))
const Σ₋₁::Float64 = (1 / (219ns))

"""
A dict of keyword arguments with default values.

Calling `make_nv_system(; DEFAULTS...) is equivalent to

```
make_nv_system(;
    A_zz = 1.0MHz,
    A_zx = 0.3MHz,
    B = 120Gauss,
    γ_c = 1.07kHz/Gauss,
    δ₋ = 0.0,
    δ₊ = 0.0,
    Γ = (1 / (12ns)),
    Γ₀ = 0.0,
    Γ₊₁ = (1 / (24ns + 0.9ns)),
    Γ₋₁ = (1 / (24ns + 0.9ns)),
    Σ₀ = (1 / (219ns)),
    Σ₊₁ = (1 / (219ns)),
    Σ₋₁ = (1 / (219ns)),
)
```

Note that it is possible to override default values by passing it again _after_
the `DEFAULTS...`, e.g., `make_nv_system(; DEFAULTS..., A_zz = 1.5MHz)`
"""
DEFAULTS = Dict{Symbol,Float64}(
    :A_zz => A_zz,
    :A_zx => A_zx,
    :B => B,
    :γ_c => γ_c,
    :δ₋ => δ₋,
    :δ₊ => δ₊,
    :Γ => Γ,
    :Γ₀ => Γ₀,
    :Γ₊₁ => Γ₊₁,
    :Γ₋₁ => Γ₋₁,
    :Σ₀ => Σ₀,
    :Σ₊₁ => Σ₊₁,
    :Σ₋₁ => Σ₋₁,
)


# deprecated: from before I understood that you can override default keyword
# arguments with, e.g., `make_nv_system(; DEFAULTS..., A_zz = 1.5MHz)`
function with_defaults(; kwargs...)
    result = copy(DEFAULTS)
    for (k, v) in kwargs
        result[k] = v
    end
    return result
end

end
