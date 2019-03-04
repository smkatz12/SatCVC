"""
satCVC_const.jl
- define any constants needed by the simulation here
(e.g. discretized values of ψ and λ, tolerances, object geometry values, 
radius of sphere defining surface, etc.)
"""

ψdisc = collect(-90.0:1:90.0)
λdisc = collect(-180.0:1:180.0)

r = 100.0

l = 20.0
w = 20.0
h = 20.0

k = 1.0

p₀ = [-90.0 0.0; 90.0 0.0]