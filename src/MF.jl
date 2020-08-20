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
	cut::Function

	function TriangularMF(l_vertex::Real, center::Real, r_vertex::Real)

		if l_vertex <= center <= r_vertex

			this = new()

			this.l_vertex = l_vertex
			this.center = center
			this.r_vertex = r_vertex

			function inv(α)
				x1 = (this.center - this.l_vertex) * α + this.l_vertex
				x2 = (this.center - this.r_vertex) * α + this.r_vertex
				return x1, x2
			end
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

			this.n_opt = function get_n(tol)
				return 0
			end

			this.cut = function cut(α)
				x1, x2 = inv(α)
				return CutMF(x1, x2, α, this)
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
	cut::Function

	function GaussianMF(center::Real, sigma::Real)

		this = new()

		this.center = center
		this.sigma = sigma

		function inv(α)
			x1 = this.center - this.sigma*sqrt(-2*log(α))
			x2 = this.center + this.sigma*sqrt(-2*log(α))
			return x1, x2
		end

		this.eval = function eval(x)
			exp( - 0.5 * ((x - this.center) / this.sigma) ^ 2)
		end

		this.mean_at = function mean_at(firing_strength)
			this.center
		end

		this.get_n = function get_n(tol)
			l, u = inv(α)
			c = this.center
			s = this.sigma
			sdiff(x) = this.eval(x)*(c^2 - 2*c*x - s^2 + x^2)/s^4
			xs = collect(l:(u-l)/100:u)
			ys = sdiff.(xs)
			ym = maximum(abs.(ys))
			sqrt((b - a)^3*ym/(12*tol))
		end

		this.cut = function cut(α)
			x1, x2 = inv(α)
			return CutMF(x1, x2, α, this)
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
	get_n::Function
	cut::Function

	function BellMF(a::Real, b::Real, c::Real)

		this = new()

		this.a = a
		this.b = b
		this.c = c

		function inv(α)
			x1 = this.c - this.a*((1-α)/α)^(1/(2*b))
			x2 = this.c + this.a*((1-α)/α)^(1/(2*b))
			return x1, x2
		end
		this.eval = function eval(x)
			1 / (1 + abs((x - this.c) / this.a) ^ (2 * this.b))
		end

		this.mean_at = function mean_at(firing_strength)
			this.c
		end

		this.get_n = function get_n(tol)
			l, u = inv(α)
			return (u-l)^3/(12*tol)
		end
		this.cut = function cut(α)
			x1, x2 = inv(α)
			return CutMF(x1, x2, α, this)
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
	get_n::Function
	cut::Function

	function TrapezoidalMF(l_bottom_vertex::Real, l_top_vertex::Real, r_top_vertex::Real, r_bottom_vertex::Real)

		if l_bottom_vertex <= l_top_vertex <= r_top_vertex <= r_bottom_vertex

			this = new()

			this.l_bottom_vertex = l_bottom_vertex
			this.l_top_vertex = l_top_vertex
			this.r_top_vertex = r_top_vertex
			this.r_bottom_vertex = r_bottom_vertex

			function inv(α)
				x1 = (this.l_top_vertex - this.l_bottom_vertex) * α + this.l_bottom_vertex
				x2 = (this.r_top_vertex - this.r_bottom_vertex) * α + this.r_bottom_vertex
				return x1, x2
			end

			this.eval = function eval(x)
				maximum([minimum([((x - this.l_bottom_vertex) / (this.l_top_vertex - this.l_bottom_vertex)), 1, ((this.r_bottom_vertex - x) / (this.r_bottom_vertex - this.r_top_vertex))]), 0])
			end

			this.mean_at = function mean_at(firing_strength)
				p1 = (this.l_top_vertex - this.l_bottom_vertex) * firing_strength + this.l_bottom_vertex
				p2 = (this.r_top_vertex - this.r_bottom_vertex) * firing_strength + this.r_bottom_vertex
				(p1 + p2) / 2
			end

			this.get_n = function get_n(tol)
				return 0
			end

			this.cut = function cut(α)
				x1, x2 = inv(α)
				return CutMF(x1, x2, α, this)
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
	get_n::Function
	cut::Function

	function SigmoidMF(a::Real, c::Real, limit::Real)

		if (a > 0 && limit > c) || (a < 0 && limit < c)

			this = new()

			this.a = a
			this.c = c
			this.limit = limit

			function inv(α)
				return this.c - (1/this.a)*log((1-α)/α)
			end

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

			this.get_n = function get_n(tol)
				a = this.a
				c = this.c
				if a > 0
					u = this.limit
					l = inv(0.01)
				else
					l = this.limit
					u = inv(0.01)
				end
				aux1(x) = exp(-a*(x-c))
				aux2(x) = 2*a^2*(aux(x))^2/(aux(x) + 1)^3
				sdiff(x) = aux1(x) - a^2*aux(x)/(aux(x) + 1)^2
				xs = collect(l:(u-l)/100:u)
				ys = sdiff.(xs)
				ym = maximum(abs.(ys))
				sqrt((u - l)^3*ym/(12*tol))
			end

			this.cut = function cut(α)
				x = inv(α)
				if this.a > 0
					x1 = x
					x2 = this.limit
				else
					x1 = this.limit
					x2 = x
				end
				return CutMF(x1, x2, α, this)
			end

			this

		else

			error("invalid parameters")

		end

	end
end

"""
	Cut MF

	CutMF(a, b, c)

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
	cut::Function

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
		end

		this.get_n = function get_n(tol)
			return this.toCutMF.get_n(tol)
		end

		this.cut = function cut(α)
			if α >= this.α
				return this
			else
				return this.toCutMF.cut(α)
			end
		end

		this

	end

end
