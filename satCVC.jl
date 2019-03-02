"""
satCVC.jl
- Implements the functions necessary to perform centroidal Voronoi coverage
control to sense a space object

Things I am currently concerned about:
Angle wrapover in the control law (this is going to be super weird when it happens
and can probably cause some undesired behavior of the system)
"""

"""
function simulate_cvc
	- main function to simulate satellite coverage
	- ******* A note on the compute_voronoi function:
		- If this were truly distributed, its input would be only the positions of its neighbors
		- We are not losing generality when we just input all positions because
		the positions any robots that are not neighbors will be irrelavent in the calculation
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
function simulate_cvc(p₀::Array{Real,2}, ψdisc::Array{Real,1}, λdisc::Array{Real,1}, k::Real, r::Real;
	max_iter = 100, tol = 0.1, dt = 0.1)
	p = p₀
	for j = 1:max_iter # Using j so that I can use i for the robots (sorry this is backwards)
		# Initialize vector to fill with new positions
		pnew = zeros(size(p₀))

		# Compute the voronoi region of each robot
		V = compute_voronoi(p, ψdisc, λdisc, r)

		# Iterate though each robot
		for i = 1:size(p,1)
			# Next compute the centroid
			CVᵢ = compute_centroid(V[i])
			
			# Find ṗᵢ
			ṗᵢ = k*rel_vector(CVᵢ,p[i,:]) # k*(CVᵢ-p[i,:])

			# Update p
			pnew[i,:] = p[i,:] + ṗᵢ*dt
		end

		# Check termination
		if norm(pnew - p) < tol
			return p
		end

		# Set up for next iteration
		p = p_new
	end

	println("Hit maximum iterations before converging :(. Returning final position anyway ...")
	return p
end

"""
function rel_vector - Somrita
"""


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
function compute_voronoi(p::Array{Real,2}, ψdisc::Array{Real,1}, λdisc::Array{Real,1}, r::Real)
	V = Dict{Int64,Array{Real,2}}()
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
function d(p₁::Array{Real,1}, p₂::Array{Real,1}, r::Real)
	d = 0 # Fill in!
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
function compute_centroid(Vᵢ::Array{Real,2})
	CVᵢ = [0 0] # Fill in!
	return CVᵢ
end

