# Contains various membership function types
# ------------------------------------------

"""	 Membership function type
"""
abstract type MF end


"""	 Triangular membership function type

	 TriangularMF(l_vertex, center, r_vertex)

	 Properties
	 ----------
	 `l_vertex`, `center` and `r_vertex` are the vertices of the triangle, in order

	 `eval` function returns membership value at a point
	 `mean_at` function returns mean value at line clipped by given firing strength
"""
mutable struct TriangularMF<:MF

	l_vertex::Real
	center::Real
	r_vertex::Real

	eval::Function
	mean_at::Function
	get_n::Function

	function TriangularMF(l_vertex::Real, center::Real, r_vertex::Real)

		if l_vertex <= center <= r_vertex

			this = new()

			this.l_vertex = l_vertex
			this.center = center
			this.r_vertex = r_vertex

			this.eval = function eval(x)
				maximum([minimum([((x - this.l_vertex) / (this.center - this.l_vertex)), ((this.r_vertex - x) / (this.r_vertex - this.center))]), 0])
			end

			this.mean_at = function mean_at(firing_strength)

				if firing_strength != 1
					p1 = (this.center - this.l_vertex) * firing_strength + this.l_vertex
					p2 = (this.center - this.r_vertex) * firing_strength + this.r_vertex
					(p1 + p2) / 2
				elseif firing_strength == 1
					return this.center
				end

			end

			this.n_opt = function get_n(a, b, tol)
				return 0
			end

			this

		else

			error("invalid vertices")

		end

	end

end


"""
 Gaussian membership function type

	GaussianMF(center, sigma)

	 Properties
	 ----------
	 `center` is the center of the distribution
	 `sigma` determines width of the distribution

	 `eval` function returns membership value at a point
	 `mean_at` function returns mean value at line clipped by given firing strength

"""
mutable struct GaussianMF<:MF

	center::Real
	sigma::Real

	eval::Function
	mean_at::Function
	get_n::Function

	function GaussianMF(center::Real, sigma::Real)

		this = new()

		this.center = center
		this.sigma = sigma

		this.eval = function eval(x)
			exp( - 0.5 * ((x - this.center) / this.sigma) ^ 2)
		end

		this.mean_at = function mean_at(firing_strength)
			this.center
		end

		this.get_n = function get_n(a, b, tol)
			c = this.center
			s = this.sigma
			sdiff(x) = this.eval(x)*(c^2 - 2*c*x - s^2 + x^2)/s^4
			xs = collect(a:100:b)
			ys = sdiff.(xs)
			im = argmax(abs.(ys))
			ym = ys[im]
			xm = a + (b - a)/100*im
			sqrt((b - a)^3*abs(ym)/(12*tol))
		end

		this

	end

end

"""
	 Generalised Bell membership function type

	BellMF(a, b, c)

	 Properties
	 ----------
	 `a`, `b` and `c` the usual bell parameters with `c` being the center

	 `eval` function returns membership value at a point
	 `mean_at` function returns mean value at line clipped by given firing strength

"""
mutable struct BellMF<:MF

	a::Real
	b::Real
	c::Real

	eval::Function
	mean_at::Function

	function BellMF(a::Real, b::Real, c::Real)

		this = new()

		this.a = a
		this.b = b
		this.c = c

		this.eval = function eval(x)
			1 / (1 + abs((x - this.c) / this.a) ^ (2 * this.b))
		end

		this.mean_at = function mean_at(firing_strength)
			this.c
		end

		this

	end

end

"""
	 Trapezoidal membership function type

	TrapezoidalMF(l_bottom_vertex, l_top_vertex, r_top_vertex, r_bottom_vertex)

	 Properties
	 ----------
	 `l_bottom_vertex`, `l_top_vertex`, `r_top_vertex` and `r_bottom_vertex` are the vertices of the trapezoid, in order

	 `eval` function returns membership value at a point
	 `mean_at` function returns mean value at line clipped by given firing strength

"""
mutable struct TrapezoidalMF<:MF

	l_bottom_vertex::Real
	l_top_vertex::Real
	r_top_vertex::Real
	r_bottom_vertex::Real

	eval::Function
	mean_at::Function

	function TrapezoidalMF(l_bottom_vertex::Real, l_top_vertex::Real, r_top_vertex::Real, r_bottom_vertex::Real)

		if l_bottom_vertex <= l_top_vertex <= r_top_vertex <= r_bottom_vertex

			this = new()

			this.l_bottom_vertex = l_bottom_vertex
			this.l_top_vertex = l_top_vertex
			this.r_top_vertex = r_top_vertex
			this.r_bottom_vertex = r_bottom_vertex

			this.eval = function eval(x)
				maximum([minimum([((x - this.l_bottom_vertex) / (this.l_top_vertex - this.l_bottom_vertex)), 1, ((this.r_bottom_vertex - x) / (this.r_bottom_vertex - this.r_top_vertex))]), 0])
			end

			this.mean_at = function mean_at(firing_strength)
				p1 = (this.l_top_vertex - this.l_bottom_vertex) * firing_strength + this.l_bottom_vertex
				p2 = (this.r_top_vertex - this.r_bottom_vertex) * firing_strength + this.r_bottom_vertex
				(p1 + p2) / 2
			end

			this

		else

			error("invalid vertices")

		end

	end

end

"""
	 Sigmoid membership function type

	SigmoidMF(a, c, limit)

	 Properties
	 ----------
	 `a` controls slope
	 `c` is the crossover point
	 `limit` sets the extreme limit

	 `eval` function returns membership value at a point
	 `mean_at` function returns mean value at line clipped by given firing strength

"""
mutable struct SigmoidMF<:MF

	a::Real
	c::Real
	limit::Real

	eval::Function
	mean_at::Function

	function SigmoidMF(a::Real, c::Real, limit::Real)

		if (a > 0 && limit > c) || (a < 0 && limit < c)

			this = new()

			this.a = a
			this.c = c
			this.limit = limit

			this.eval = function eval(x)
				1 / (1 + exp(-this.a * (x - this.c)))
			end

			this.mean_at = function mean_at(firing_strength)

				if firing_strength == 1
					p_firing_strength = 0.999
				elseif firing_strength == 0
					p_firing_strength = 0.001
				else
					p_firing_strength = firing_strength
				end

				p1 = -log((1 / p_firing_strength) - 1) / this.a + this.c
				p2 = this.limit
				(p1 + p2) / 2

			end

			this

		else

			error("invalid parameters")

		end

	end
end

"""
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	Generalised Bell membership function type

	BellMF(a, b, c)

	 Properties
	 ----------
	 `a`, `b` and `c` the usual bell parameters with `c` being the center

	 `eval` function returns membership value at a point
	 `mean_at` function returns mean value at line clipped by given firing strength

"""
mutable struct CutMF<:MF

	x1::Real
	x2::Real
	α::Real
	toCutMF::MF

	eval::Function
	mean_at::Function

	function CutMF(x1::Real, x2::Real, α::Real, toCutMF::MF)

		this = new()

		this.x1 = x1
		this.x2 = x2
		this.α = α
		this.toCutMF = toCutMF

		this.eval = function eval(x)
			if x1 <= x <= x2
				return this.α
			else
				return this.toCutMF.eval(x)
			end
		end

		this.mean_at = function mean_at(firing_strength)
			if firing_strength >= α
				return this.toCutMF.mean_at(α)
			else
				return this.toCutMF.mean_at(firing_strength)
		end

		this

	end

end
