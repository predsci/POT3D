import math
# Taken from report generator's PredictiveScience.math python
# file. 

def deg_to_rad(angle):
	# designed to replicate ps_deg_rad.sh
	# Converts angle from degrees to radians.
	return math.radians(angle % 360)

def rad_to_deg(angle):
	# designed to replicate ps_rad_deg.sh
	# Converts angle from radians to degrees.
	return math.degrees(angle) % 360

class Math:
	@staticmethod
	def t_rad_lat_deg(angle):
		"""
		Convert theta angle in radians into latitude angle in degrees.

		float: angle
			Theta angle in radians.

		float:
			Latitude angle in degrees.
		"""
		return 90. - math.degrees(angle)


	@staticmethod
	def lat_deg_t_rad(angle):
		"""
		Convert latitude angle in degrees into theta angle in radians.

		float: angle
			Latitude angle in degrees.

		float:
			Theta angle in radians.
		"""
		return 0.5 * math.pi - math.radians(angle)


	@staticmethod
	def p_rad_lon_deg(angle):
		"""
		Convert phi angle in radians into longitude angle in degrees.

		float: angle
			Phi angle in radians.

		float:
			Longitude angle in degrees.
		"""
		return math.degrees(angle)


	@staticmethod
	def lon_deg_p_rad(angle):
		"""
		Convert longitude angle in degrees into phi angle in radians.

		float: angle
			Longitude angle in degrees.

		float:
			Phi angle in radians.
		"""
		return math.radians(angle)


	@staticmethod
	def spherical_cartesian(point):
		"""
		Convert a point in spherical coordinate to Cartesian coordinate.

		tuple of 3: point
			A point (r,t,p) in spherical coordinate.

		tuple of 3:
			A point (x,y,z) in Cartesian coordinate.
		"""
		sin_t = math.sin(point[1])

		return (point[0] * sin_t * math.cos(point[2]),
				point[0] * sin_t * math.sin(point[2]),
				point[0] * math.cos(point[1]))


	@staticmethod
	def cartesian_spherical(point):
		"""
		Convert a point in Cartesian coordinate to spherical coordinate.

		tuple of 3: point
			A point (x,y,z) in Cartesian coordinate.

		tuple of 3:
			A point (r,t,p) in spherical coordinate.
		"""
		x_2 = point[0] * point[0]
		y_2 = point[1] * point[1]
		z_2 = point[2] * point[2]

		r = math.sqrt(x_2 + y_2 + z_2)
		t = math.atan2(math.sqrt(x_2 + y_2), point[2])
		p = math.atan2(point[1], point[0])
		if p < 0.0:
			p += 2.0 * math.pi

		return (r, t, p)


	@staticmethod
	def distance(xy_point0, xy_point1):
		dx = xy_point1[0] - xy_point0[0]
		dy = xy_point1[1] - xy_point0[1]

		return math.sqrt(dx * dx + dy * dy)


	@staticmethod
	def get_nice_numbers(min_val, max_val, num_vals):
		gross_step = (max_val - min_val) / num_vals

		step = math.pow(10, math.floor(math.log(gross_step, 10)))
		if (5.0 * step) <= gross_step:
			step = 5.0 * step
		elif (2.5 * step) <= gross_step:
			step = 2.5 * step
		elif (2.0 * step) <= gross_step:
			step = 2.0 * step

		num_values = math.floor(
			math.floor(max_val / step) - math.ceil(min_val / step)) + 1
		min_value = math.ceil(min_val / step) * step
		max_value = math.floor(max_val / step) * step

		num_digits = round(abs(math.log(step, 10))) + 1
		values = []
		for i in range(num_values):
			values.append(round(min_value + i * step, num_digits))

		return values
