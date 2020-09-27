# Solutions of Equations in One variable



The most basic problems of Numerical approximation the **root-finding problem**.
$$
f(x) = 0
$$




## The Bisection Method



The Bisection technique based on the **Intermediate Value Theorem**, called the **Binary-search method**

First, suppose f is continuous function defined on the interval [a, b], with $f(a)\cdot f(b) <0$,The Intermediate Value Theorem tell us that  $\exist p \in (a, b)$ with $f(p) = 0$.

To begin:

- Set $a_1 = a$ and $b_1=b$, and $p_1 = \frac{a_1+b_1}{2}$
- if $f(p_1) = 0$, So done.
- If $f(p_1) \ne 0$, then $f(p_1)$ has the same sign as either $f(a)$ and $f(b)$
  - If $f(p_1)$ and $f(a)$ have the same sign, root in $(p_1, b_1)$, change: $a_2 = p_1$ and $b_2 = b_1$.
  - if $f(p_1)$ and $f(b)$ have the same sign, root in $(a_1, p_1)$, change: $a_2 = a_1$  and $b_2=p_1$.

```@example 1
using NumericalAnalysis: SEq1
SEq1.Bisection(x->cos(x)-x, 0.5, π/4)
```

you can see steps output:

```@example 1
SEq1.Bisection(x->cos(x)-x, 0.5, π/4, output_seq=true)
```
****
****

## fixed_point Iteration

Given a root-finding problem $f(p) = 0$, we can define functions $g$ with a fixed point at $p$ in a number of ways,  as:

$g(x) = x - f(x)$

Conversely, if the function $g$ has a fixed point at $p$, then the function defined by:

$f(x) = x - g(x)$

has a zero at $p$.

```@example 1
SEq1.fixed_point(x -> cos(x), 0.74)
```

```@example 1
SEq1.fixed_point(x -> cos(x), 0.74, output_seq=true)
```
***
***

## Newton's Method
Suppose that $f ∈ C^2[a, b]$. Let $p₀ ∈ [a, b]$ be an approximation to p such that $f'(p₀) ≠ 0$ and $|p − p₀|$ is “small.” Consider the first Taylor polynomial for $f(x)$ expanded about $p₀$ and evaluated at $x = p$.

$f(p) = f(p₀)+(p - p₀)f'(p₀)+\frac{(p-p₀)^2}{2}f''(ξ(p₀))$

where $ξ(p₀)$ lies between $p$ and $p₀$, since $f(p)=0$, Newton’s method is derived by assuming that since $|p − p₀|$ is small, the term involving $(p-p₀)^2$ is much smaller. So:

$0 \approx f(p₀)+(p-p₀)f'(p₀)$

So:

$p \approx p₀ - \frac{f(p₀)}{f'(p₀)}$

This sets the stage for Newton’s method, which starts with an initial approximation $p₀$, and generates the sequence $\{ p_n \}_{n=0}^{∞}$, by

$p_n = p_{n-1} - \frac{f(p_{n-1})}{f'(p_{n-1})}$

```@example 1
SEq1.Newton(x->cos(x) - x, 0.74)
```

```@example 1
SEq1.Newton(x->cos(x) - x, 0.74, output_seq=true)
```
***
***

## The Secant Method
To circumvent the problem of the derivative evaluation in Newton’s method, we introduce a slight variation. By definition,

$f'(p_{n-1}) = \lim_{x \to p_{n-1}} \frac{f(x) - f(p_{n-1})}{x-p_{n-1}}$

If $p_{n-2}$ is close to $p_{n-1}$, then

$f'(p_{n-1}) \approx \frac{f(p_{n-2}) - f(p_{n-1})}{p_{n-2} - p_{n-1}}$

So

$p_{n} = p_{n-1} - \frac{f(p_{n-1})(p_{n-1} - p_{n-2})}{f(p_{n-1})- f(p_{n-2})}$

```@example 1
SEq1.Secant(x -> cos(x) - x, 0.5, π/4)
```

```@example 1
SEq1.Secant(x -> cos(x) - x, 0.5, π/4, output_seq=true)
```
***
***

## False position
The method of False Position (also called Regula Falsi) generates approximations in the same manner as the Secant method, but it includes a test to ensure that the root is always bracketed between successive iterations

First choose initial approximations $p₀$ and $p₁$ with $f(p₀)·f(p₁) < 0$. The approximation $p₂$ is chosen in the same manner as in the Secant method, as the x-intercept of the line joining $(p₀, f(p₀))$ and $(p₁, f(p₁))$. To decide which secant line to use to compute $p₃$, consider $f(p₂) · f(p₁)$, or $sgn[f(p₂)]·sgn[f(p₁)]$.

- If $sgn[f(p₂)]·sgn[f(p₁)] < 0$, then $p₁$ and $p₂$ bracket a root. Choose $p₃$ as the x-intercept of the line joining $(p₁,f(p₁))$ and $(p₂,f(p₂))$.

- If not, choose $p₃$ as the x-intercept of the line joining $(p₀, f(p₀))$ and $(p₂, f(p₂))$, and then interchange the indices on $p₀$ and $p₁$.

```@example 1
SEq1.FalsePos(x->cos(x)-x, 0.5, π/4)
```

```@example 1
SEq1.FalsePos(x->cos(x)-x, 0.5, π/4, output_seq=true)
```
***
***





## Method
```@autodocs
Modules = [NumericalAnalysis.SEq1]
```
