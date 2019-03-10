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
sa- x pointing to 0°N,0°E, y pointing to 0°N,90°E, and z pointing to 90°N


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
function calcϕ(ψdisc::Array{Float64}, λdisc::Array{Float64}, l::Float64, w::Float64, h::Float64, r::Float64)
	# Initialize ϕ
	ϕ = zeros(length(ψdisc),length(λdisc))
	for i = 1:length(ψdisc)
		for j = 1:length(λdisc)
			# Get vector from point on sphere to center (the origin) in intertial cartesian coordinates
			vc = get_vector_to_center(ψdisc[i], λdisc[j], r)
			# Figure out which surface normal we care about based on the vector to the origin
			n = choose_surface_normal(vc, l, w, h)
			# Fill in ϕ using the dot product (normalize vector to center first)
			vc_norm = vc./r
			ϕ[i,j] = abs(dot(vc_norm, n))^4
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
function get_vector_to_center(ψ::Float64, λ::Float64, r::Float64)
	# Formula from https://www.movable-type.co.uk/scripts/latlong-vectors.html
	vc = -r*[cosd(ψ)*cosd(λ) cosd(ψ)*sind(λ) sind(ψ)]
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
function choose_surface_normal(vc::Array{Float64,2}, l::Float64, w::Float64, h::Float64)
	#flip vector to get outward vector from center

	r = -vc/norm(vc)
	dist2plane = [0.0 0.0 0.0 0.0 0.0 0.0] #distance to 6 planes +z, -z, +x, -x, +y, -y
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
		dist2plane[4] = dist + 2 #opposite plane distance adds penalty, not a minimum
	else
		dist2plane[4] = -dist
		dist2plane[3] = -dist + 2 #opposite plane distance adds penalty, not a minimum
	end

	# y plane distances
	a = 0
	b = 1
	c = 0
	d = -b*w/2
	dist = -d/(a*r[1]+b*r[2]+c*r[3])
	if dist>0
		dist2plane[5] = dist
		dist2plane[6] = dist + 2 #opposite plane distance adds penalty, not a minimum
	else
		dist2plane[6] = -dist
		dist2plane[5] = -dist + 2 #opposite plane distance adds penalty, not a minimum
	end

	min_dist = dist2plane[1]
	n = [0.0 0.0 1.0]

	for i in 1:1:6
		if dist2plane[i]<min_dist
			min_dist = dist2plane[i]
			n = normals[i:i,1:3]
		end
	end
	return n
end

"""
function get_areas
	- Calculate the area on the surface of the sphere associated with each point in the angle mesh

	INPUTS:
	- ψdisc: discretized ψ points
	- λpartitions: array of size equal to length(ψdisc), each index represents the number of equal longitude partitions
	at that latitutde
	- r: radius of the sphere to place the satellites on

	OUTOUTS:
	- areas: matrix that constains an area for each (ψ,λ)
	ψ varies in the rows while λ varies in the columns
	(e.g. areas[1,2] is surface area around the point on sphere with angles given by 
	first discretized value of ψ and the second discretized value of λ)
"""
function get_areas(ψdisc::Array{Float64}, λpartitions::Array{Float64}, r::Float64)
	# Initialize areas
	areas = zeros(length(ψdisc))
	for i = 1:length(ψdisc)
		
		# special consideration for first latitude (closest to -90)
		if i == 1
			#special consideration for first latitude, nearest to south pole
			if ψdisc[i] == -90
				#if first latitude is -90, get area of the circle surrounding south pole
				#find latitude change required to get half way to next latitude
				dlat = abs(ψdisc[i+1]-ψdisc[i])/2
				radius = 2*pi*r*dlat/360 	#radius is the arlength based on dlat
				areas[i] = pi*radius^2/λpartitions[i]		#divide area by number of times south pole will be added in integral
			else
				#if first latitude is not -90, get sector area of circle surrounding south pole
				#find latitude from -90 corresponding to half way to next latitude
				dlat = 90-abs(ψdisc[i+1]+ψdisc[i])/2
				radius = 2*pi*r*dlat/360	#radius is arclength based on dlat
				areas[i] = pi*radius^2/λpartitions[i] 	#divide into equal sectors based on partitions at this latitude
			end
		elseif i == length(ψdisc)
			#special consideration for last latitude nearest to north pole
			if ψdisc[i] == 90
				#if last latitude is 90, get area of the circle surrounding north pole
				#find latitude change required to get half way to previous latitude
				dlat = abs(ψdisc[i]-ψdisc[i-1])/2
				radius = 2*pi*r*dlat/360 	#radius is the arlength based on dlat
				areas[i] = pi*radius^2/λpartitions[i]		#divide area by number of times south pole will be added in integral
			else
				#if first latitude is not 90, get sector area of circle surrounding north pole
				#find latitude from 90 corresponding to half way to next latitude
				dlat = 90-abs(ψdisc[i]+ψdisc[i-1])/2
				radius = 2*pi*r*dlat/360	#radius is arclength based on dlat
				areas[i] = pi*radius^2/λpartitions[i] 	#divide into equal sectors based on partitions at this latitude
			end

		else
			#all other latitudes between first and last latitude
			#approximate each area by a trapezoid.
			#trapezoid top length is circumference of middle of next latitude and current latitude, divide by partitions
			lat_top = (ψdisc[i+1]+ψdisc[i])/2
			radius_top = r*cosd(abs(lat_top))	#get radius of level curve (lol) of top latitude
			length_top = 2*radius_top*pi/λpartitions[i] 	#length of top edge of trapezoid

			#trapezoid bottom length is circuference of middle of previous latitude and current latitude, divide by partitions
			lat_bot = (ψdisc[i]+ψdisc[i-1])/2
			radius_bot = r*cosd(abs(lat_bot))	#get radius of level curve (lol) of bot latitude
			length_bot = 2*radius_bot*pi/λpartitions[i] 	#length of bot edge of trapezoid

			#height of trapezoid is arclength between top and bottom latitude
			dlat = abs(lat_top-lat_bot)
			height = 2*pi*r*dlat/360

			areas[i] = (length_top+length_bot)/2*height

		end
		
		
	end
	return areas
end
