# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: IJulia
#     language: julia
#     name: ijulia
# ---

# # Polarization with Optical Pumping

# We reproduce here the dynamics for a "known good" single chirp pulse. The dynamics include optical pumping after he chirp concludes. This produced the results shown in the poster "Microwave Pulse Trains for Quantum Metrology with C13 in Diamond" presented by M. Goerz at the GRC "Quantum Control of Light and Matter" 2025 in Newport, RI. The only difference is the relabeling of $|+1⟩ \leftrightarrow |-1⟩$, due to the non-standard use of quantum numbers in the poster, which was corrected here.

using C13NV.Models: make_nv_system, ket
using C13NV.Defaults: with_defaults
using C13NV.Units: kHz, MHz, ms, μs, ns, Gauss
using C13NV.Amplitudes: ConstantDrive, LinearChirp
using QuantumPropagators.Shapes: flattop

using Plots

const γcB = 128.4kHz
const γcB_per_2π = 128.4kHz / 2π # internal energy unit in Jabir's notebooks
const γcB⁻¹ = 2π / 128.4kHz;  # internal time unit in Jabir's notebooks)

const Ω₀ = 17γcB_per_2π
Ω₀ / MHz

const t₀ = 10γcB⁻¹
t₀ / μs

const α = 20γcB_per_2π / γcB⁻¹
α / (MHz / ms)

const Γ = (1 / (12ns))

tlist = collect(range(0, 4 * 10γcB⁻¹; length = 4001));

# +
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
# -

Λ = SinglePump(Γ; pump_width = 0.4, ramp_width = 0.01)

L, labels = make_nv_system(;
    with_defaults()...,
    Ω₋ = ConstantDrive(17γcB_per_2π),
    ω₋ = LinearChirp(; t₀, α),
    Γ = (1 / (12ns)),
    Λ,
    frame = :rwa,
);

labels

const N = length(labels)

Ψ_up = ket(("G", "0", "↑"), labels);
Ψ_down = ket(("G", "0", "↓"), labels);
ρ₀ = reshape(0.5 * Ψ_up * Ψ_up' + 0.5 * Ψ_down * Ψ_down', N * N);

using QuantumPropagators: propagate, Newton

states = propagate(
    ρ₀,
    L,
    tlist;
    method = Newton,
    check = true,
    storage = true,
    show_progress = true,
);

size(states)

"""Get matrix of population for `states`.

```
pops = get_pops(states)
```

takes `states` as returned by `propagate` with `storage=true`, assuming a
vectorized density matrix as the initial state. That is, `states` is assumed
to be a complex array with shape ``(N², nt)`` where ``N`` is the size of the
Hilbert space, and `nt` is the number of time grid points.

The returned `pops` is a real-valued matrix with shape ``(N, nt)`` containing
the population for each level.

If `states` is given as single vectorized density matrix, a time grid
dimension with `nt=1` is implied.
"""
function get_pops(states)::Matrix{Float64}
    if ndims(states) == 1
        # single state - add time dimension
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

pops = get_pops(states);

pop_E = sum([pops[i, :] for i in eachindex(labels) if labels[i][1] == "E"]);

plot(tlist, pops[findfirst(==(("G", "0", "↑")), labels), :]; label = "|0↑⟩")
plot!(tlist, pops[findfirst(==(("G", "0", "↓")), labels), :]; label = "|0↓⟩")
plot!(tlist, pops[findfirst(==(("G", "-1", "↑")), labels), :]; label = "|-1↑⟩")
plot!(tlist, pops[findfirst(==(("G", "-1", "↓")), labels), :]; label = "|-1↓⟩")
plot!(tlist, pop_E; label = "|E⟩", color = :red, ls = :dash)
plot!(tlist, t -> (Λ(t) / Λ.amplitude); label = "Λ(t)/|Λ|", color = :black)
plot!(; legend = :outertop, legend_column = -1, xlabel = "time (μs)", ylabel = "population")

pops[:, 2000]

pops[:, end]

# ## Diagonal Frame

# **TODO: why does this not reproduce the dynamics?**

# +
L_diag, labels = make_nv_system(;
    with_defaults()...,
    Ω₋ = ConstantDrive(17γcB_per_2π),
    ω₋ = LinearChirp(; t₀, α),
    Γ = (1 / (12ns)),
    Λ,
    frame = :diag,
);

states_diag = propagate(
    ρ₀,
    L_diag,
    tlist;
    method = Newton,
    check = true,
    storage = true,
    show_progress = true,
);

pops_diag = get_pops(states_diag);

plot(tlist, pops_diag[findfirst(==(("G", "0", "↑")), labels), :]; label = "|0↑⟩")
plot!(tlist, pops_diag[findfirst(==(("G", "0", "↓")), labels), :]; label = "|0↓⟩")
plot!(tlist, pops_diag[findfirst(==(("G", "-1", "↑")), labels), :]; label = "|-1↑⟩")
plot!(tlist, pops_diag[findfirst(==(("G", "-1", "↓")), labels), :]; label = "|-1↓⟩")
plot!(tlist, t -> Λ(t) / Λ.amplitude; label = "Λ(t)/|Λ|", color = :black)
plot!(;
    tite = "diagonal frame",
    legend = :outertop,
    legend_column = -1,
    xlabel = "time (μs)",
    ylabel = "population"
)

# -

pops_diag[:, 2000]
