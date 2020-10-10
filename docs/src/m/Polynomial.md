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
***
***
## Neville’s Iterated Interpolation
Let $f$ be defined at $x_{0},x_{1},...,x_{k}$ and let $x_{i}$ and $x_{j}$ be two distinct numbers in this set. Then

$P(x)=\frac{(x-x_{j})P_{0,1,...,j-1,j+1,...,k}(x)-(x-x_{i})P_{0,1,...,i-1,i+1,...,k}(x)}{(x_{i}-x_{j})}$

```@example 3
x = [1.0, 1.3, 1.6, 1.9, 2.2]
y = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
Polynomial.Neville(x, y, 1.5, tab=true)
```

## Newton’s Divided-Difference Formula
Suppose that $P_{n}(x)$ is the nth Lagrange polynomial that agrees with the function $f$ at the distinct numbers $x_0 , x_1 , . . . , x_n$ . Although this polynomial is unique, there are alternate algebraic representations that are useful in certain situations. The divided differences of f with respect to $x_0, x_1, . . . , x_n$ are used to express $P_n(x)$ in the form

$P_{n}(x)=a_{0}+a_{1}(x-x_{0})+a_{2}(x-x_{0})(x-x_{1})+...+a_{n}(x-x_{0})...(x-x_{n-1}),$

So

$a_{0}=P_{n}(x_{0})=f(x_{0})$

$P_{n}(x_{1})=a_{0}+a_{1}(x-x_{0})=f(x_{1})$

$a_{1} = \frac{f(x_{1})-f(x_{0})}{x_{1}-x_{0}}.$

```@example 3
x = [1.0, 1.3, 1.6, 1.9, 2.2]
y = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
print(Polynomial.NDDF(x, y, backward=true))
```


## Method
```@autodocs
Modules = [NumericalAnalysis.Polynomial]
```
