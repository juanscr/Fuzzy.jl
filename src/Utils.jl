"""
	Estimates the integral of f between a and b by the trapezoidal rule

    trapz(f, a, b, n)

	 Parameters
	 ----------
     `f` is the function to integrate
     `a` lower bound for the interval
     `b` upper bound for the interval
     `n` number of partition of [a,b]
"""

function trapz(f, a, b, n)
    h = (b - a)/n
    integral = f(a) + f(b)
    i = 1
    while i < n
        x1 = a + i*h
        integral += 2*f(x1)
        i += 1
    end
    integral *= h/2
    return integral
end
