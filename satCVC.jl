"""
satCVC.jl
- Implements the functions necessary to perform centroidal Voronoi coverage
control to sense a space object

Things I am currently concerned about:
Angle wrapover in the control law (this is going to be super weird when it happens
and can probably cause some undesired behavior of the system)
"""

using Distances
using LinearAlgebra

include("create_environment.jl")
include("satCVC_const.jl")

"""
function simulate_cvc
	- main function to simulate satellite coverage
	- ******* A note on the compute_voronoi function:
		- If this were truly distributed, its input would be only the positions of its neighbors
		- We are not losing generality when we just input all positions because
		the positions any robots that are not neighbors will be irrelevant in the calculation
	- TODO: Add functionality to keep track of p on each iteration in a 3D array so that we can plot
	their entire trajectories instead of just final position
	- Anything after the semicolon is an optional input with a default value

	INPUTS:
	- p₀: initial positions (ψ, λ) for the robots
		- Each row is a different robot (so columns are ψ and λ)
	- ψdisc: discretized ψ points
	- λdisc: discretized λ points
	- k: controller gain (this will need tuning)
	- r: radius of the sphere
	-------------
	- max_iter: maximum number of time steps
	- tol: magnitude of change in position must be smaller than this to declare convergence
	- dt: time step

	OUTPUTS:
	- p: final position of the robots (set up the same as p0)
"""
function simulate_cvc(p₀::Array{Float64,2}, ψdisc::Array{Float64,1}, λdisc::Array{Float64,1}, k::Float64, r::Float64, areas;
	max_iter = 100, tol = 1, dt = 0.1)
	p = p₀
	for j = 1:max_iter # Using j so that I can use i for the robots (sorry this is backwards)
		# Initialize vector to fill with new positions
		pnew = zeros(size(p₀))

		# Compute the voronoi region of each robot
		V = compute_voronoi(p, ψdisc, λdisc, r)

		# Iterate though each robot
		for i = 1:size(p,1)
			# Next compute the centroid
			# if i == 1
			# 	println(findall(V[i][:,1].>0))
			# end
			CVᵢ = compute_centroid(V[i], ϕ, ψdisc, λdisc, areas)

			# Find ṗᵢ
			ṗᵢ = k*rel_vector(CVᵢ,p[i,:]) # k*(CVᵢ-p[i,:])

			# Update p
			pnew[i,:] = norm_p(p[i,:] + vec(ṗᵢ*dt))
		end

		# Check termination
		if norm(pnew - p) < tol
			return p
		end

		# Set up for next iteration
		println("norm:", norm(pnew-p))
		p = pnew

		println(p)

		println("position of robots:",p)
	end

	println("Hit maximum iterations before converging :(. Returning final position anyway ...")
	return p
end

"""
function compute_voronoi
	- return all of the discretized points in the voronoi region of Q
	(according to our geodesic distance function)
	- This should fully work once d is implemented (no need for further code)

	INPUTS:
	- p: current robot positions
	- ψdisc: discretized ψ points
	- λdisc: discretized λ points
	- r: radius of the sphere

	OUTPUTS:
	- V: Dictionary mapping robot index to array of points in Voronoi region
"""
function compute_voronoi(p::Array{Float64,2}, ψdisc::Array{Float64,1}, λdisc::Array{Float64,1}, r::Float64)
	V = Dict{Int64,Array{Float64,2}}()
	for i = 1:length(ψdisc)
		for j = 1:length(λdisc)
			min_dist = Inf
			closest_robot = 0
			for k = 1:size(p,1) # num robots
				dist = d([ψdisc[i],λdisc[j]], p[k,:], r) # Distance between discretized point and robot k
				if dist < min_dist
					min_dist = dist
					closest_robot = k
				end
			end
			if haskey(V, closest_robot)
				V[closest_robot] = vcat(V[closest_robot], [ψdisc[i] λdisc[j]]) # This may not work immediately (may have to mess with transposes/dimensions)
			else
				V[closest_robot] = [ψdisc[i] λdisc[j]]
			end
		end
	end
	return V
end

"""
function d - Keiko
	- computes distance between two points on the sphere
	- use haversine formula
		- https://en.wikipedia.org/wiki/Haversine_formula?fbclid=IwAR1opVdva9xhRW5JuWO2RMV5_uRwX2M31VJFbRNVZ7Yx81nRF2z6-Ngz254

	INPUTS:
	- p₁: first point of the form (ψ,λ)
	- p₂: second point of the form (ψ,λ)
	- r: radius of the sphere

	OUTPUTS:
	- dist: distance between the two points
"""
function d(p₁::Array{Float64,1}, p₂::Array{Float64,1}, r::Float64)
	x = [p₁[2], p₁[1]]
	y = [p₂[2], p₂[1]]
	d = abs(haversine(x,y,r))
	return d
end

"""
function compute_centroid - Keiko
	- compute the centroid of a Voronoi region based on its points

	INPUTS:
	- Vᵢ: Array of points (columns (ψ,λ)) in the Voronoi region	of robot i

	OUTPUTS:
	- CVᵢ: Centroid of the Voronoi region in the form (ψ,λ)
"""
function compute_centroid(Vᵢ::Array{Float64,2}, ϕ::Array{Float64,2}, ψdisc::Array{Float64,1}, λdisc::Array{Float64,1}, areas)
	#ϕ = [2 1 2; 1 2 1; 2 1 2] #test
	num_q_in_voronoi_cell = size(Vᵢ)[1]
	CVᵢ_num_XYZ = [0;0;0]
	CVᵢ_den = 0

	for k = 1:num_q_in_voronoi_cell
		q = [Vᵢ[k,1]; Vᵢ[k,2]]
		ψind = findfirst(ψdisc.==q[1])
		λind = findfirst(λdisc.==q[2])
		q_XYZ = convertToXYZ(ψdisc[ψind], λdisc[λind])
		CVᵢ_num_XYZ = CVᵢ_num_XYZ + q_XYZ*ϕ[ψind,λind]*areas[ψind]
		CVᵢ_den = CVᵢ_den + ϕ[ψind,λind]*areas[ψind]
	end

	CVᵢ_num = convertToLatLon(CVᵢ_num_XYZ[1], CVᵢ_num_XYZ[2], CVᵢ_num_XYZ[3])

	CVᵢ_XYZ = CVᵢ_num_XYZ./CVᵢ_den # Fill in!

	CVᵢ = convertToLatLon(CVᵢ_XYZ[1], CVᵢ_XYZ[2], CVᵢ_XYZ[3])
	println("final: $CVᵢ")
	return CVᵢ
end

function convertToXYZ(lat, lon)
	X = cosd(lat) * cosd(lon)
	Y = cosd(lat) * sind(lon)
	Z = sind(lat)
	return [X,Y,Z]
end

function convertToLatLon(x, y, z)
	Lon = atand(y, x)
	Hyp = sqrt(x * x + y * y)
	Lat = atand(z, Hyp)
	return [Lat,Lon]
end

"""
function rel_vector - Somrita
	- returns our version of CVᵢ - pᵢ that accounts for angle wrapover to be used in the control law

	INPUTS:
	- CVᵢ: centroid of the voronoi cell for robot i
	- pᵢ: position of robot i

	OUTPUTS:
	- rel_vector: vector from robot position to centroid
"""
function rel_vector(CVᵢ::Array{Float64,1}, pᵢ::Array{Float64,1})
	latDiff = CVᵢ[1] -  pᵢ[1]
	longDiff = CVᵢ[2] - pᵢ[2]
	longDiffCorr = longDiff
	if longDiff > 180
		longDiffCorr = longDiff - 360
	elseif longDiff < -180
		longDiffCorr = longDiff + 360
	end
	# e.g. to get from -175 degrees to 175 degrees, you can either go 175-(-175)=350 degrees or the equivalent -(360-(175-(-175)))=-10 degrees
	rel_vector = [latDiff longDiffCorr]
	return rel_vector
	# Tests:
	# rel_vector([89.0,175],[-89.0,-175])
end

"""
function norm_p
	- normalizes p to have a ψ between -90 and 90 and a
	λ between -180 and 180

	INPUTS:
	- p: position vector in the form of (ψ,λ)

	OUTPUTS:
	- normp: equivalent position vector with ψ and λ in the correct ranges
"""
function norm_p(p::Array{Float64,1})
	norm_p = p

	if norm_p[1] < -90
		norm_p[1] = -(norm_p[1]+180)
		norm_p[2] = norm_p[2] + 180
	end
	if norm_p[1] > 90
		norm_p[1] = -(norm_p[1]-180)
		norm_p[2] = norm_p[2] + 180
	end

	while norm_p[2] < -180
		norm_p[2] += 360
	end
	while norm_p[2] > 180
		norm_p[2] -= 360
	end

	return norm_p
end

# Run everything
ϕ = calcϕ(ψdisc, λdisc, l, w, h, r);
areas = get_areas(ψdisc, λpartitions, r);
p = simulate_cvc(p₀, ψdisc, λdisc, k, r, areas, max_iter = 100, tol = 0.4)
