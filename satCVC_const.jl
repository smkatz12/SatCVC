"""
satCVC_const.jl
- define any constants needed by the simulation here
(e.g. discretized values of ψ and λ, tolerances, object geometry values, 
radius of sphere defining surface, etc.)
"""

ψdisc = collect(-90.0:1:90.0)
λdisc = collect(-180.0:1:179.0)

# ψdisc = collect(-90.0:30:60.0)
# λdisc = collect(-180.0:30:150.0)

λpartitions = length(λdisc)*ones(length(ψdisc))

r = 100.0

l = 20.0
w = 20.0
h = 20.0

k = 1.0

p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; 40.0 40.0; 50.0 50.0; 60.0 60.0; 70.0 70.0]