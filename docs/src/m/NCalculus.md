# Numerical Differentiation and Integration
***
***


## Differentiation

The derivative of the function $f$ at $x_0$ is

$f'(x_0) = lim_{h \to 0}\frac{f(x_0+h) - f(x_0)}{h}$




### Three and Five Point formula.

when $x_1 = x_0+h$ and $x_2 = x_0 + 2h$ for some $h\ne0$ gives

**Three-Point Formulas**

$f'(x_0) = \frac{1}{h}[-\frac{3}{2}f(x_0) + 2f(x_1) - \frac{1}{2}f(x_2)] + \frac{h^2}{3}f^{(3)}(\xi_0)$



**Three-Point Endpoint Formula**

$f'(x_0) = \frac{1}{2h}[-3f(x_0)+4f(x_0+h) - f(x_0+2h)] + \frac{h^2}{3}f^{(3)}(\xi_0)$



**Three-Point Midpoint formula**

$f'(x_0) = \frac{1}{2h}[f(x_0+h) - f(x_0-h)] - \frac{h^2}{6}f^{(3)}(\xi_0)$

```@example 4
using NumericalAnalysis: NCalculus
NCalculus.ThreePoint(x->sin(x), 0, 0.1)
```

**Five-Point MidPoint formula**

$f'(x_0) = \frac{1}{12h}[f(x_0 - 2h) - 8f(x_0 -h)+8f(x_0+h)-f(x_0 + 2h)]+\frac{h^4}{30}f^{(5)}(\xi_0)$

**Five-Point EndPoint formula**

$f'(x_0) = \frac{1}{12h}[-25f(x_0)+48f(x_0 +h)-36f(x_0+2h)+16f(x_0+3h)-3f(x_0+4h)]+\frac{h^4}{5}f^{(5)}(\xi_0)$

```@example 4
NCalculus.FivePoint(x->sin(x), 0, 0.1)
```





## Integration

### Trapezoidal Rule

$\int_a^b f(x)dx = \frac{h}{2}[f(a)+f(b)] - \frac{h^3}{12}f''(\xi)$

```@example 4
NCalculus.Trapezoidal(x->x^2, 0, 2)
```

### Simpson's Rule

$\int_a^b f(x)dx = \frac{h}{3}[f(a) + f(a+h) + f(b)] - \frac{h^5}{90}f^{(4)}(\xi)$

```@example 4
NCalculus.Simpson(x->x^2, 0, 2)
```

### Newton-Cotes Formulas
will calculate n=1,2,3,4. formula rule. if you input n in $[1,4]$, will get a value.

```@example 4
NCalculus.Newton_cotes(x->sin(x),0,π/4)
```

```@example 4
NCalculus.Newton_cotes(x->sin(x), 0, π/4, n=4)
```

```@example 4
NCalculus.Newton_cotes(x->sin(x),0,π/4,n=3, method="Open")
```

### Romberg Integration

$R_{k,j} = R_{k,j-1} + \frac{1}{4^{j-1} - 1}(R_{k,j-1} - R_{k-1, j-1}), k = j,j+1,...$

```@example 4
NCalculus.Romberg(x->sin(x), 0, π/4, 5, seq_tab=true)
```
```@example 4
NCalculus.Romberg(x->sin(x), 0, π/4, 5)
```

### Gaussian Quadrature

Suppose we want to determine $c_1, c_2, x_1$ and $x_2$ so

$\int_{-1}^{1}f(x)dx=c_1f(x_1) + c_2f(x_2)$

we use a Ploynomial $f(x) = a_0 + a_1x+a_2x^2+a_3x^3$

we get:
a solution $c_1 = c_2=1, x_1 = -\frac{\sqrt{3}}{3}, x_2=\frac{\sqrt{3}}{3}$

So:

$\int_{-1}^{1}f(x)dx \approx f(\frac{-\sqrt{3}}{3}) + f(\frac{\sqrt{3}}{3})$

```@example 4
NCalculus.Gaussian_Quad(x->sin(x),n=5, a=0, b=π/4)
```

```@example 4
NCalculus.Gaussian_Quad(x->x^2)
```

## Multiple Integrals

In this section, Consider the double integral

$\iint_R f(x,y)dA$

where $R = {(x,y) | a \leq x \leq b, c \leq y \leq d}$. for some constants a, b, c, d, is a rectangular in the plane.

we apply the **Composite Trapezoidal rule** to calculate:

### Simpson Double Integral

$\int_a^b\int_{ydown(x)}^{yup(x)}f(x, y)dydx$

```@example 4
f(x, y) = log(x + 2y)
y_up(x) = 1.5
y_down(x) = 1.0
NCalculus.SimpsonDoubleIntegral(f, (1.4, 2.0), y_up, y_down)
```

### Gaussian Double Integral

```@example 4
NCalculus.GaussianDoubleIntegral(f, (1.4, 2.0), y_up, y_down)
```



## Method

```@autodocs
Modules = [NumericalAnalysis.NCalculus]
```
