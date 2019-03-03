"""
create_environment.jl
- determines the environment for the sensor coverage problem by calculating
ϕ(q) based on the object geometry

General notes:
ψ defines elevation (defined like latitude ranging from -90 to 90 degrees)
λ defines azimuth (defined like longitude ranging from -180 to 180 degrees)
- defined like lat and lon so that we can use haversine formula as distance function

xyz coordinate system centered at center of sphere
- x is 0 when longitude is 0
- x pointing to 0°N,0°E, y pointing to 0°N,90°E, and z pointing to 90°N


Here is how we are going to define our object:
- It is a rectangular prism centered at (0,0)
- Front face = all x positive, back face = all x negative
- Right face = all y positive, left face = all y negative
- Top face = all z positive, bottom face = all z negative
- All vectors specifying parameters will then be in order of
[Front, back, right, left, top, bottom]
- To fully specify our geometry we just need the length (corresponds to y-direction),
width (corresponds to x-direction), and height (corresponds to z-direction)
"""

"""
function calcϕ
	- Calculate the sensory function at each point in the discretized
	space based on the object geometry (assuming some form of rectangular
	prism for now)

	INPUTS:
	- ψdisc: discretized ψ points
	- λdisc: discretized λ points
	- l: length of the sensed object
	- w: width of the sensed object
	- h: height of the sensed object
	- r: radius of the sphere to place the satellites on

	OUTOUTS:
	- ϕ: matrix that constains a value for each ϕ(ψ,λ)
	ψ varies in the rows while λ varies in the columns
	(e.g. ϕ[1,2] is ϕ of the first discretized value of ψ and the second discretized value of λ)
"""
function calcϕ(ψdisc::Array{Real}, λdisc::Array{Real}, l::Real, w::Real, h::Real, r::Real)
	# Initialize ϕ
	ϕ = zeros(length(ψdisc),length(θdisc))
	for i = 1:length(ψdisc)
		for j = 1:length(λdisc)
			# Get vector from point on sphere to center (the origin) in intertial cartesian coordinates
			vc = get_vector_to_center(ψdisc[i], λdisc[i], r)
			# Figure out which surface normal we care about based on the vector to the origin
			n = choose_surface_normal(vc, l, w, h)
			# Fill in ϕ using the dot product (normalize vector to center first)
			vc_norm = vc./r
			ϕ[i,j] = abs(dot(vc_norm, n))
		end
	end
	return ϕ
end

"""
function get_vector_to_center - Somrita
	- calculate the vector from the point on the sphere defined by
	ψ and λ to the center of the sphere (the origin)
	- This should just be a matter of getting the cartesian coordinates
	(careful not typical spherical coord conversion since using lat/lon coord system)
	- Check that norm is r to test function

	INPUTS:
	- ψ: latitude
	- λ: longitude
	- r: radius of the sphere

	OUTPUTS:
	- vc: vector to the origin
"""
function get_vector_to_center(ψ::Real, λ::Real, r::Real)
	# Formula from https://www.movable-type.co.uk/scripts/latlong-vectors.html
	vc = r*[cosd(ψ)*cosd(λ) cosd(ψ)*sind(λ) sind(ψ)]
	return vc
	# Tests:
	# get_vector_to_center(0,0,5)
	# get_vector_to_center(0,90,5)
	# get_vector_to_center(0,180,5)
	# get_vector_to_center(0,-179,5)
	# get_vector_to_center(45,90,5)
	# get_vector_to_center(45,0,5)
end

"""
function choose_surface_normal - Simon
	- returns the surface normal of the surface that is being observed
	from a particular position on the sphere
	- OUTWARD SURFACE NORMAL
	- should be the first surface that vc crosses

	INPUTS:
	- vc: vector from the point on the sphere to the center (origin)
	- l: length of the object
	- w: width of the object
	- h: height of the object

	OUTPUTS:
	- n: normal vector to the surface we care about (outward)
"""
function choose_surface_normal(vc::Array{Real,1}, l::Real, w::Real, h::Real)
	#flip vector to get outward vector from center

	r = -vc/norm(vc)
	dist2plane = [0 0 0 0 0 0] #distance to 6 planes +z, -z, +x, -x, +y, -y
	normals = [0 0 1;0 0 -1;1 0 0;-1 0 0;0 1 0;0 -1 0]
	# planes have equations of form ax+by+cz+d=0
	# z plane distances
	a = 0
	b = 0
	c = 1
	d = -c*h/2;
	dist = -d/(a*r[1]+b*r[2]+c*r[3])
	if dist>0
		dist2plane[1] = dist
		dist2plane[2] = dist + 2 #opposite plane distance adds penalty, not a minimum
	else
		dist2plane[2] = -dist
		dist2plane[1] = -dist + 2 #opposite plane distance adds penalty, not a minimum
	end

	# x plane distances
	a = 1
	b = 0
	c = 0
	d = -a*l/2
	dist = -d/(a*r[1]+b*r[2]+c*r[3])
	if dist>0
		dist2plane[3] = dist
		dist2plane[5] = dist + 2 #opposite plane distance adds penalty, not a minimum
	else
		dist2plane[5] = -dist
		dist2plane[3] = -dist + 2 #opposite plane distance adds penalty, not a minimum
	end

	# y plane distances
	a = 0
	b = 1
	c = 0
	d = -b*w/2
	dist = -d/(a*r[1]+b*r[2]+c*r[3])
	if dist>0
		dist2plane[4] = dist
		dist2plane[6] = dist + 2 #opposite plane distance adds penalty, not a minimum
	else
		dist2plane[6] = -dist
		dist2plane[4] = -dist + 2 #opposite plane distance adds penalty, not a minimum
	end

	min_dist = dist2plane[1]
	n = [0 0 1]

	for i in 1:1:6
		if dist2plane[i]<min_dist
			min_dist = dist2plane[i]
			n = normals[i:i,1:3]
		end
	return n
end
