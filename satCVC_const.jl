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

k = 4.0

# p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; 40.0 40.0; 50.0 50.0]

fileName = "20RobotsSamePoint.csv"
#p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; -10.0 -10.0; -20.0 -20.0; -30.0 -30.0; 0.0 -10.0; 0.0 -20.0; 0.0 -30.0; 0.0 -40.0; -10.0 0.0; 
#	-20.0 0.0; -30.0 0.0; -40.0 0.0; 50.0 50.0; 60.0 60.0; -50.0 -50.0; -45.0 45.0; 120.0 -120.0]

p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; -10.0 -10.0; -20.0 -20.0; -30.0 -30.0; 0.0 -10.0; 0.0 -20.0; 0.0 -30.0; 0.0 -40.0; -10.0 0.0]

# fileName = "TwoRobotsSamePoint.csv"
# p₀ = [0.0 0.0; 0.0 180.0]

fileName = "TwelveRobotsRandom.csv"
# p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; 40.0 40.0; 50.0 50.0]
lats = 180*rand(12,1) .-90
longs = 360*rand(12,1) .-180
p₀ = hcat(lats,longs)

