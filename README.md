# Fuzzy EAFIT Fork

Mamdani and Sugeno type Fuzzy Inference System in julia. This code is a fork of
[Phelipe](https://github.com/phelipe/Fuzzy.jl).

This fork adds some new features and corrects some deprecated functions.

Tested with Julia 1.5.0.

## QuickStart

```julia
julia> using Fuzzy
```

### Mamdani FIS

- Create input, output membership functions and rules

```julia
julia> input_a = Dict("small" => TriangularMF(1, 2, 3), "large" => TriangularMF(4, 5, 6))
julia> input_b = Dict("small" => TriangularMF(1, 2, 3))

julia> inputs = [input_a, input_b]
julia> output = Dict("small" => TriangularMF(1, 2, 3))

julia> rule = Rule(["large", "small"], "small")
julia> rules = [rule]
```

- Create FIS

```julia
julia> fis = FISMamdani(inputs, output, rules)
```

- Find output

```julia
julia> in_vals = [4.7, 2.3]
julia> eval_fis(fis, in_vals)
```

### Sugeno FIS

- Create input membership functions and rules with consequence coefficients

```julia
julia> input_a = Dict("small" => TriangularMF(1, 2, 3), "large" => TriangularMF(5, 6, 7))
julia> input_b = Dict("small" => TriangularMF(1, 2, 3))

julia> inputs = [input_a, input_b]

julia> rule1 = Rule(["large", "small"], [1.0, 1.0, 1.0])
julia> rule2 = Rule(["small", "small"], [0.0, 0.0, 5.0])
julia> rules = [rule]
```

- Create FIS

```julia
julia> fis = FISSugeno(inputs, rules)
```

- Find output

```julia
julia> in_vals = [2.3, 1.2]
julia> eval_fis(fis, in_vals)
```

## Features

- FIS

  - Mamdani
  - Sugeno

- Membership functions

  - Triangular
  - Gaussian
  - Bell
  - Trapezoidal
  - Sigmoid

- Defuzzification

  - Mean of Maximum
  - Weighted Average (default)
  - _Centroid (coming soon)_

- T-Norm

  - Minimum (MIN)
  - Algebraic product (A-PROD)
  - Bounded difference (B-DIF)
  - Drastic product (D-PROD)
  - Einstein product (E-PROD)
  - Hamacher product (H-PROD)

- S-Norm

  - Maximum (MAX)
  - Algebraic sum (A-SUM)
  - Bounded sum (B-SUM)
  - Drastic sum (D-SUM)
  - Einstein sum (E-SUM)
  - Hamacher sum (H-SUM)

## License

MIT
