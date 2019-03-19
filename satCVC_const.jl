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

# Cubic target
r = 100.0
l = 20.0
w = 20.0
h = 20.0



# Flat Earth
#r = 100.0
#l = 20.0
#w = 20.0
#h = 1.0

k = 4.0


# p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; 40.0 40.0; 50.0 50.0]


#fileName = "test" # DONT PUT .csv ANYMORE
#p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; -10.0 -10.0; -20.0 -20.0; -30.0 -30.0; 0.0 -10.0; 0.0 -20.0; 0.0 -30.0; 0.0 -40.0; -10.0 0.0; 
#	-20.0 0.0; -30.0 0.0; -40.0 0.0; 50.0 50.0; 60.0 60.0; -50.0 -50.0; -45.0 45.0; 120.0 -120.0]

#p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; -10.0 -10.0; -20.0 -20.0; -30.0 -30.0; 0.0 -10.0; 0.0 -20.0; 0.0 -30.0; 0.0 -40.0; -10.0 0.0]

#fileName = "8cluster"
#p₀ = [1.0 0.2; -1.0 -0.2; 1.0 0.1; 0.0 -0.1; 0.0 0.3; 0.0 -0.3; 0.0 0.4; 0.0 -0.4]

# # fileName = "TwoRobotsSamePoint.csv"
#p₀ = [45.0 90.0; 45.0 -90.0]

fileName = "30_flatearth"
#p₀ = [0.0 0.0; 10.0 10.0; 20.0 20.0; 30.0 30.0; 40.0 40.0; 50.0 50.0]
lats = 10.0*rand(30,1) .-5.0
longs = 20.0*rand(30,1) .-10.0
p₀ = hcat(lats,longs)

