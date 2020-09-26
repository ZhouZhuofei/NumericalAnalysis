# Solutions of Equations in One variable



The most basic problems of Numerical approximation the **root-finding problem**.
$$
f(x) = 0
$$



## Introduction
### The Bisection Method



The Bisection technique based on the **Intermediate Value Theorem**, called the **Binary-search method**



First, suppose f is continuous function defined on the interval [a, b], with $f(a)\cdot f(b) <0$,The Intermediate Value Theorem tell us that  $\exist p \in (a, b)$ with $f(p) = 0$.

To begin:

- Set $a_1 = a$ and $b_1=b$, and $p_1 = \frac{a_1+b_1}{2}$
- if $f(p_1) = 0$, So done.
- If $f(p_1) \ne 0$, then $f(p_1)$ has the same sign as either $f(a)$ and $f(b)$
  - If $f(p_1)$ and $f(a)$ have the same sign, root in $(p_1, b_1)$, change: $a_2 = p_1$ and $b_2 = b_1$.
  - if $f(p_1)$ and $f(b)$ have the same sign, root in $(a_1, p_1)$, change: $a_2 = a_1$  and $b_2=p_1$.

```@example
  using NumericalAnalysis: SEq1
  SEq1.Bisection(sin, π/2, 3π/2)
```
and the real answers
```@example
π
```

## Method
```@autodocs
Modules = [NumericalAnalysis.SEq1]
```
