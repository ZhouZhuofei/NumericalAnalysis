# Numerical Analysis
*****
[![Build Status](https://travis-ci.com/ZhouZhuofei/NumericalAnalysis.jl.svg?branch=master)](https://travis-ci.com/ZhouZhuofei/NumericalAnalysis.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ZhouZhuofei/NumericalAnalysis.jl?svg=true)](https://ci.appveyor.com/project/ZhouZhuofei/NumericalAnalysis-jl)
[![Coverage](https://coveralls.io/repos/github/ZhouZhuofei/NumericalAnalysis.jl/badge.svg?branch=master)](https://coveralls.io/github/ZhouZhuofei/NumericalAnalysis.jl?branch=master)
[![codecov](https://codecov.io/gh/ZhouZhuofei/NumericalAnalysis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ZhouZhuofei/NumericalAnalysis.jl)
[![License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](https://github.com/ZhouZhuofei/NumericalAnalysis.jl/blob/master/LICENSE)

****
It's a julia package for Numerical Analysis.

Now, i just write a little function in it.
I am a beginner of Julia and numerical analysis, there may be many problems, 👏 to discuss with me.🤣

***

## Basic

here are some function:



Now, i just write some functions:

- N the derivative, use ForwardDiff packageto calculate $\frac{dy}{dx}$, then recursive to get Nth derivative.(emmm, I feel a bit slow)
- Taylor Polynomial, get the value nth Taylor Ploynomial.
- Bisection function, find root
***
***




- NthDerivative(f, x, n)
```@example 1
using NumericalAnalysis
NthDerivative(x->x^3, 4, 4) # 0
```
```@example 1
NthDerivative(x->x^3, 4, 2) # 24
```
```@example 1
NthDerivative(cos, 1, 5) # NthDerivative(cos, 1, 5)
```

********************

- TaylorPolynomials(f, x, x₀, n)
```@example 1
TaylorPolynomials(cos, 0.1, 0, 6) # TaylorPolynomials(cos, 0.1, 0, 6)
```
```@example 1
TaylorPolynomials(x->x^3, 1.1, 1, 6) #1.331000..
```

********************

- Bisection(f, a, b)
```@example 1
using NumericalAnalysis: SEq1
SEq1.Bisection(sin, π/2, 3π/2)
```
```@example 1
SEq1.Bisection(x->x-1, 0, 3π/2, output_seq=true)
```
