# Numerical Differentiation and Integration
***
***


## Differentiation

The derivative of the function $f$ at $x_0$ is

$f'(x_0) = lim_{h \to 0}\frac{f(x_0+h) - f(x_0)}{h}$





when $x_1 = x_0+h$ and $x_2 = x_0 + 2h$ for some $h\ne0$ gives

**Three-Point Formulas**



$f'(x_0) = \frac{1}{h}[-\frac{3}{2}f(x_0) + 2f(x_1) - \frac{1}{2}f(x_2)] + \frac{h^2}{3}f^{(3)}(\xi_0)$



**Three-Point Endpoint Formula**

$f'(x_0) = \frac{1}{2h}[-3f(x_0)+4f(x_0+h) - f(x_0+2h)] + \frac{h^2}{3}f^{(3)}(\xi_0)$



**Three-Point Midpoint formula**





## Integration
