# Interpolation and the Lagrange Polynomial


##  nth Lagrange interpolating polynomial

If $x_{0},x_{1},x_{2},...,x_{n}$ are n+1 distinct numbers and f is a function whose values are given at these numbers, then a unique polynomial $P(x)$ of degree at most $n$ exists with:

$f(x_{k})=P(x_{k}), for \quad each \quad k=0,1,2,...,n$

This polynomial is given by

$P(x)=f(x_{0})L_{n,0}(x)+...+f(x_{n})L_{n,n}(x)=\sum_{k=0}^{n}f(x_{k})L_{n,k}(x),$

Where, for each $k=0,1,2,...,n$

$L_{n,k}(x)=\frac{(x-x_{0})(x-x_{1})...(x-x_{k-1})(x-x_{x+1})...(x-x_{n})}{(x_{k}-x_{0})(x_{k}-x_{1})...(x_{k}-x_{k-1})(x_{k}-x_{k+1})...(x_{k}-x_{n})}=\prod_{i=0;i\ne k}^{n}\frac{(x-x_{i})}{(x_{k}-x_{i})}.$

```@example 3
using Plots
using NumericalAnalysis: Polynomial
x = 1:0.01:5
s = Polynomial.Lagrange([2, 11/4, 4], [1/2, 4/11, 1/4])
plot(x, [s.(x) 1 ./ x], title = "nth Lagrange polynomial", label = ["y = p(x)" "y = f(x)"], xlabel="x", ylabel="y")
```

## Method
```@autodocs
Modules = [NumericalAnalysis.Polynomial]
```
