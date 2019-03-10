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

# Constants
r = 5.0
k = 4.0


# p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; 40.0 40.0; 50.0 50.0]


#fileName = "test" # DONT PUT .csv ANYMORE
#p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; -10.0 -10.0; -20.0 -20.0; -30.0 -30.0; 0.0 -10.0; 0.0 -20.0; 0.0 -30.0; 0.0 -40.0; -10.0 0.0;
#	-20.0 0.0; -30.0 0.0; -40.0 0.0; 50.0 50.0; 60.0 60.0; -50.0 -50.0; -45.0 45.0; 120.0 -120.0]

#p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; -10.0 -10.0; -20.0 -20.0; -30.0 -30.0; 0.0 -10.0; 0.0 -20.0; 0.0 -30.0; 0.0 -40.0; -10.0 0.0]

#fileName = "8cluster"
#p₀ = [1.0 0.2; -1.0 -0.2; 1.0 0.1; 0.0 -0.1; 0.0 0.3; 0.0 -0.3; 0.0 0.4; 0.0 -0.4]

# # fileName = "TwoRobotsSamePoint.csv"
# # p₀ = [0.0 0.0; 0.0 180.0]

fileName = "ThirtyRobotsAsteroid_ellipsoid"
#p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; 40.0 40.0; 50.0 50.0]
lats = 10.0*rand(20,1) .-5.0
longs = 10.0*rand(20,1) .-5.0
p₀ = hcat(lats,longs)
