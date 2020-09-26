# NumericalAnalysis

[![Build Status](https://travis-ci.com/ZhouZhuofei/NumericalAnalysis.jl.svg?branch=master)](https://travis-ci.com/ZhouZhuofei/NumericalAnalysis.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ZhouZhuofei/NumericalAnalysis.jl?svg=true)](https://ci.appveyor.com/project/ZhouZhuofei/NumericalAnalysis-jl)
[![Coverage](https://coveralls.io/repos/github/ZhouZhuofei/NumericalAnalysis.jl/badge.svg?branch=master)](https://coveralls.io/github/ZhouZhuofei/NumericalAnalysis.jl?branch=master)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://zhouzhuofei.github.io/NumericalAnalysis.jl/docs/build/index.html)
[![License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](https://github.com/ZhouZhuofei/NumericalAnalysis.jl/blob/master/LICENSE)


****

I am a beginner of Julia and numerical analysis, there may be many problems, 👏 to discuss with me.🤣

***

## Basic

here are some function:

```julia
julia> using NumericalAnalysis

```

Now, i just write two function:

- N the derivative, use ForwardDiff packageto calculate $\frac{dy}{dx}$, then recursive to get Nth derivative.(emmm, I feel a bit slow)
- Taylor Polynomial, get the value nth Taylor Ploynomial.



`NthDerivative(f, x, n)`

```julia
NthDerivative(x->x^3, 4, 4) # 0
NthDerivative(x->x^3, 4, 2) # 24
NthDerivative(cos, 1, 5) # NthDerivative(cos, 1, 5)
```

`TaylorPolynomials(f, x, x₀, n)`

```julia
TaylorPolynomials(cos, 0.1, 0, 6) # TaylorPolynomials(cos, 0.1, 0, 6)
TaylorPolynomials(x->x^3, 1.1, 1, 6) #1.331000..
```

.....
