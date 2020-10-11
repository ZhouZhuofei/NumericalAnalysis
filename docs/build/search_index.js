var documenterSearchIndex = {"docs":
[{"location":"m/Polynomial/#Interpolation-and-the-Lagrange-Polynomial","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"","category":"section"},{"location":"m/Polynomial/#nth-Lagrange-interpolating-polynomial","page":"Interpolation and the Lagrange Polynomial","title":"nth Lagrange interpolating polynomial","text":"","category":"section"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"If x_0x_1x_2x_n are n+1 distinct numbers and f is a function whose values are given at these numbers, then a unique polynomial P(x) of degree at most n exists with:","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"f(x_k)=P(x_k) for quad each quad k=012n","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"This polynomial is given by","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"P(x)=f(x_0)L_n0(x)++f(x_n)L_nn(x)=sum_k=0^nf(x_k)L_nk(x)","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"Where, for each k=012n","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"L_nk(x)=frac(x-x_0)(x-x_1)(x-x_k-1)(x-x_x+1)(x-x_n)(x_k-x_0)(x_k-x_1)(x_k-x_k-1)(x_k-x_k+1)(x_k-x_n)=prod_i=0ine k^nfrac(x-x_i)(x_k-x_i)","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"using Plots\nusing NumericalAnalysis: Polynomial\nx = 1:0.01:5\ns = Polynomial.Lagrange([2, 11/4, 4], [1/2, 4/11, 1/4])\nplot(x, [s.(x) 1 ./ x], title = \"nth Lagrange polynomial\", label = [\"y = p(x)\" \"y = f(x)\"], xlabel=\"x\", ylabel=\"y\")","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"","category":"page"},{"location":"m/Polynomial/#Neville’s-Iterated-Interpolation","page":"Interpolation and the Lagrange Polynomial","title":"Neville’s Iterated Interpolation","text":"","category":"section"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"Let f be defined at x_0x_1x_k and let x_i and x_j be two distinct numbers in this set. Then","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"P(x)=frac(x-x_j)P_01j-1j+1k(x)-(x-x_i)P_01i-1i+1k(x)(x_i-x_j)","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"x = [1.0, 1.3, 1.6, 1.9, 2.2]\ny = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]\nPolynomial.Neville(x, y, 1.5, tab=true)","category":"page"},{"location":"m/Polynomial/#Newton’s-Divided-Difference-Formula","page":"Interpolation and the Lagrange Polynomial","title":"Newton’s Divided-Difference Formula","text":"","category":"section"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"Suppose that P_n(x) is the nth Lagrange polynomial that agrees with the function f at the distinct numbers x_0  x_1      x_n . Although this polynomial is unique, there are alternate algebraic representations that are useful in certain situations. The divided differences of f with respect to x_0 x_1     x_n are used to express P_n(x) in the form","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"P_n(x)=a_0+a_1(x-x_0)+a_2(x-x_0)(x-x_1)++a_n(x-x_0)(x-x_n-1)","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"So","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"a_0=P_n(x_0)=f(x_0)","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"P_n(x_1)=a_0+a_1(x-x_0)=f(x_1)","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"a_1 = fracf(x_1)-f(x_0)x_1-x_0","category":"page"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"x = [1.0, 1.3, 1.6, 1.9, 2.2]\ny = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]\nprint(Polynomial.NDDF(x, y, backward=true))","category":"page"},{"location":"m/Polynomial/#Method","page":"Interpolation and the Lagrange Polynomial","title":"Method","text":"","category":"section"},{"location":"m/Polynomial/","page":"Interpolation and the Lagrange Polynomial","title":"Interpolation and the Lagrange Polynomial","text":"Modules = [NumericalAnalysis.Polynomial]","category":"page"},{"location":"m/Polynomial/#NumericalAnalysis.Polynomial","page":"Interpolation and the Lagrange Polynomial","title":"NumericalAnalysis.Polynomial","text":"Polynomial\n\none of functions mapping the set of real numbers is the algebraic polynomials, the set of functions of the form:\n\np_n(x) = a_n x^n + a_n-1 x^n-1 +  +a_1x + a_0\n\n\n\n\n\n","category":"module"},{"location":"m/Polynomial/#NumericalAnalysis.Polynomial.Lagrange-Tuple{Array{T,1} where T,Array{T,1} where T}","page":"Interpolation and the Lagrange Polynomial","title":"NumericalAnalysis.Polynomial.Lagrange","text":"Lagrange(X::Vector, Y::Vector; a::String=\"x\")\n\nuse the (x_0 y_0) (x_1 y_1)  (x_n y_n), to caculate a polynomial. input X and Y, default use x to express functions. output: function.\n\n\n\n\n\n","category":"method"},{"location":"m/Polynomial/#NumericalAnalysis.Polynomial.NDDF-Tuple{Array{T,1} where T,Array{T,1} where T}","page":"Interpolation and the Lagrange Polynomial","title":"NumericalAnalysis.Polynomial.NDDF","text":"NDDF(x::Vector, y::Vector; a::String=\"x\", simple::Bool=true, tab::Bool=false, backward::Bool=false)\n\ninput x_1 x_2 x_3 x_n,valuesf(x_1), f(x_2), ...,f(x_n), set Argssimpel = trueoutput simplify ans, set tab to output a table. setbackward=true` use the Newton Backward–Difference Formula.\n\n\n\n\n\n","category":"method"},{"location":"m/Polynomial/#NumericalAnalysis.Polynomial.Neville-Tuple{Array{T,1} where T,Array{T,1} where T,Real}","page":"Interpolation and the Lagrange Polynomial","title":"NumericalAnalysis.Polynomial.Neville","text":"Neville(x::Vector, y::Vector, x₀::Real; tab::Bool=false)\n\ninput two vector, x and y. satisfy the function y=f(x), then input x₀, which you want to  approximate f(x₀). The default tab setting is false, just output the estimate. you can set true to output a result table.\n\n\n\n\n\n","category":"method"},{"location":"m/SEq1/#Solutions-of-Equations-in-One-variable","page":"SEq in one variable","title":"Solutions of Equations in One variable","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"The most basic problems of Numerical approximation the root-finding problem. $ f(x) = 0 $","category":"page"},{"location":"m/SEq1/#The-Bisection-Method","page":"SEq in one variable","title":"The Bisection Method","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"The Bisection technique based on the Intermediate Value Theorem, called the Binary-search method","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"First, suppose f is continuous function defined on the interval [a, b], with f(a)cdot f(b) 0,The Intermediate Value Theorem tell us that  exist p in (a b) with f(p) = 0.","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"To begin:","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"Set a_1 = a and b_1=b, and p_1 = fraca_1+b_12\nif f(p_1) = 0, So done.\nIf f(p_1) ne 0, then f(p_1) has the same sign as either f(a) and f(b)\nIf f(p_1) and f(a) have the same sign, root in (p_1 b_1), change: a_2 = p_1 and b_2 = b_1.\nif f(p_1) and f(b) have the same sign, root in (a_1 p_1), change: a_2 = a_1  and b_2=p_1.","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"using NumericalAnalysis: SEq1\nSEq1.Bisection(x->cos(x)-x, 0.5, π/4)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"you can see steps output:","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.Bisection(x->cos(x)-x, 0.5, π/4, output_seq=true)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/#fixed_point-Iteration","page":"SEq in one variable","title":"fixed_point Iteration","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"Given a root-finding problem f(p) = 0, we can define functions g with a fixed point at p in a number of ways,  as:","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"g(x) = x - f(x)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"Conversely, if the function g has a fixed point at p, then the function defined by:","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"f(x) = x - g(x)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"has a zero at p.","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.fixed_point(x -> cos(x), 0.74)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.fixed_point(x -> cos(x), 0.74, output_seq=true)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/#Newton's-Method","page":"SEq in one variable","title":"Newton's Method","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"Suppose that f  C^2a b. Let p₀  a b be an approximation to p such that f(p₀)  0 and p  p₀ is “small.” Consider the first Taylor polynomial for f(x) expanded about p₀ and evaluated at x = p.","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"f(p) = f(p₀)+(p - p₀)f(p₀)+frac(p-p₀)^22f(ξ(p₀))","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"where ξ(p₀) lies between p and p₀, since f(p)=0, Newton’s method is derived by assuming that since p  p₀ is small, the term involving (p-p₀)^2 is much smaller. So:","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"0 approx f(p₀)+(p-p₀)f(p₀)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"So:","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"p approx p₀ - fracf(p₀)f(p₀)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"This sets the stage for Newton’s method, which starts with an initial approximation p₀, and generates the sequence  p_n _n=0^, by","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"p_n = p_n-1 - fracf(p_n-1)f(p_n-1)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.Newton(x->cos(x) - x, 0.74)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.Newton(x->cos(x) - x, 0.74, output_seq=true)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/#The-Secant-Method","page":"SEq in one variable","title":"The Secant Method","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"To circumvent the problem of the derivative evaluation in Newton’s method, we introduce a slight variation. By definition,","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"f(p_n-1) = lim_x to p_n-1 fracf(x) - f(p_n-1)x-p_n-1","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"If p_n-2 is close to p_n-1, then","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"f(p_n-1) approx fracf(p_n-2) - f(p_n-1)p_n-2 - p_n-1","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"So","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"p_n = p_n-1 - fracf(p_n-1)(p_n-1 - p_n-2)f(p_n-1)- f(p_n-2)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.Secant(x -> cos(x) - x, 0.5, π/4)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.Secant(x -> cos(x) - x, 0.5, π/4, output_seq=true)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/#False-position","page":"SEq in one variable","title":"False position","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"The method of False Position (also called Regula Falsi) generates approximations in the same manner as the Secant method, but it includes a test to ensure that the root is always bracketed between successive iterations","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"First choose initial approximations p₀ and p₁ with f(p₀)f(p₁)  0. The approximation p₂ is chosen in the same manner as in the Secant method, as the x-intercept of the line joining (p₀ f(p₀)) and (p₁ f(p₁)). To decide which secant line to use to compute p₃, consider f(p₂)  f(p₁), or sgnf(p₂)sgnf(p₁).","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"If sgnf(p₂)sgnf(p₁)  0, then p₁ and p₂ bracket a root. Choose p₃ as the x-intercept of the line joining (p₁f(p₁)) and (p₂f(p₂)).\nIf not, choose p₃ as the x-intercept of the line joining (p₀ f(p₀)) and (p₂ f(p₂)), and then interchange the indices on p₀ and p₁.","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.FalsePos(x->cos(x)-x, 0.5, π/4)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.FalsePos(x->cos(x)-x, 0.5, π/4, output_seq=true)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/#Modified-Newton's-Way","page":"SEq in one variable","title":"Modified Newton's Way","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"Newton and ModifiedNewton, Both methods are rapidly convergent to the actual zero.","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"p_n = p_n-1 - fracf(p_n-1)f(p_n-1)f(p_n-1)^2 - f(p_n-1)f(p_n-1)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.ModifiedNewton(x -> cos(x) - x, 0.5)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.ModifiedNewton(x -> cos(x) - x, 0.5, output_seq=true)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"","category":"page"},{"location":"m/SEq1/#Muller-method","page":"SEq in one variable","title":"Muller method","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"find complex Zeros.","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"search in Google for more about Müller’s Method","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.Muller(x->x^4-3x^3+x^2+x+1, 0.5, -0.5, 0.0)","category":"page"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"SEq1.Muller(x->x^4-3x^3+x^2+x+1, 0.5, -0.5, 0.0, output_seq=true)","category":"page"},{"location":"m/SEq1/#Method","page":"SEq in one variable","title":"Method","text":"","category":"section"},{"location":"m/SEq1/","page":"SEq in one variable","title":"SEq in one variable","text":"Modules = [NumericalAnalysis.SEq1]","category":"page"},{"location":"m/SEq1/#NumericalAnalysis.SEq1","page":"SEq in one variable","title":"NumericalAnalysis.SEq1","text":"SEq1\n\nSolutions of Equations in One Variable here are functions to find root:\n\nBisection\nFixed-Point Iteration\nNewton Method\nThe Secant Method\nThe False Position Method\nMüller’s Method\n\n\n\n\n\n","category":"module"},{"location":"m/SEq1/#NumericalAnalysis.SEq1.Bisection-Tuple{Function,Real,Real}","page":"SEq in one variable","title":"NumericalAnalysis.SEq1.Bisection","text":"Bisection(f::Function, a::Real, b::Real; TOL::Float64=0.0000001, N::Int=40, output_seq::Bool=false)\n\nUse Bisection To find a solution to.\n\nf(x) = 0\n\nYou should use the function, and a interval a b, tolerance TOL, maximum number of iterations N.\n\nThen you will get approximate solution p or can't find in maximum N and return -1.\n\nyou can set output_seq=true to output sequence.\n\n\n\n\n\n","category":"method"},{"location":"m/SEq1/#NumericalAnalysis.SEq1.FalsePos-Tuple{Function,Real,Real}","page":"SEq in one variable","title":"NumericalAnalysis.SEq1.FalsePos","text":"FalsePos(f::Function, p₀::Real, p₁::Real; TOL::Float64=0.00000001, N::Int=40, output_seq::Bool=false)\n\nUse the false position method to find root. a function and interval a b, then you can set Args or use default.\n\n\n\n\n\n","category":"method"},{"location":"m/SEq1/#NumericalAnalysis.SEq1.ModifiedNewton-Tuple{Function,Real}","page":"SEq in one variable","title":"NumericalAnalysis.SEq1.ModifiedNewton","text":"ModifiedNewton(f::Function, p₀::Real; TOL::Float64=0.00000001, N::Int=20, output_seq::Bool=false)\n\nNewton and ModifiedNewton, Both methods are rapidly convergent to the actual zero.\n\n\n\n\n\n","category":"method"},{"location":"m/SEq1/#NumericalAnalysis.SEq1.Muller-Tuple{Function,Real,Real,Real}","page":"SEq in one variable","title":"NumericalAnalysis.SEq1.Muller","text":"Muller(f::Function, p₀::Real, p₁::Real, p₂::Real; TOL::Float64=0.00000001, N::Int=20, output_seq::Bool=false)\n\n\n\n\n\n","category":"method"},{"location":"m/SEq1/#NumericalAnalysis.SEq1.Newton-Tuple{Function,Real}","page":"SEq in one variable","title":"NumericalAnalysis.SEq1.Newton","text":"Newton(f::Function, p₀::Real; TOL::Float64=0.0000001, N::Int=20, output_seq=false)\n\nUse Newton method to find root.\n\nyou should use a function and initial p₀, then you can set Args or use default.\n\n\n\n\n\n","category":"method"},{"location":"m/SEq1/#NumericalAnalysis.SEq1.Secant-Tuple{Function,Real,Real}","page":"SEq in one variable","title":"NumericalAnalysis.SEq1.Secant","text":"Secant(f::Function, p₀::Real, )\n\nUse the Secant Method to find root A function and p₀, p₁(p₀<p₁), then you can set Args or use default.\n\n\n\n\n\n","category":"method"},{"location":"m/SEq1/#NumericalAnalysis.SEq1.fixed_point-Tuple{Function,Real}","page":"SEq in one variable","title":"NumericalAnalysis.SEq1.fixed_point","text":"fixed_point(f::Function, p₀::Real; TOL::Float64=0.0000001, N::Int=30, output_seq::Bool=false)\n\nUse fixed point Iteration to find root.\n\nThe number p is a fixed point for a given function g if g(p) = p.\n\ninput a function, initial p₀, optional Args: tolerancs TOL, maximum number of iterations N, output_seq.\n\nOutput a root p or sequence.\n\n\n\n\n\n","category":"method"},{"location":"m/basic/#Basic-Functions","page":"Basic function","title":"Basic Functions","text":"","category":"section"},{"location":"m/basic/","page":"Basic function","title":"Basic function","text":"some basic functions in numerical analysis, you could use.","category":"page"},{"location":"m/basic/#Usage","page":"Basic function","title":"Usage","text":"","category":"section"},{"location":"m/basic/#Functions","page":"Basic function","title":"Functions","text":"","category":"section"},{"location":"m/basic/","page":"Basic function","title":"Basic function","text":"Modules = [NumericalAnalysis]","category":"page"},{"location":"m/basic/#NumericalAnalysis.NthDerivative-Tuple{Function,Real,Int64}","page":"Basic function","title":"NumericalAnalysis.NthDerivative","text":"NthDrivative(f::Function, x::Real, n::Int)\n\nreturn the nth derivative value at x\n\n\n\n\n\n","category":"method"},{"location":"m/basic/#NumericalAnalysis.TaylorPolynomials-Tuple{Function,Real,Real,Int64}","page":"Basic function","title":"NumericalAnalysis.TaylorPolynomials","text":"TaylorPolynomials(f::Function, x::Real, x₀::Real, n::Int)\n\nreturn the value of nth Taylor Ploynomial at x.\n\n\n\n\n\n","category":"method"},{"location":"m/basic/#NumericalAnalysis.absolute_error-Tuple{Any,Any}","page":"Basic function","title":"NumericalAnalysis.absolute_error","text":"absolute_error(x, x̂)\n\nreturn the Absolute error of x and x̂\n\n\n\n\n\n","category":"method"},{"location":"m/basic/#NumericalAnalysis.relative_error-Tuple{Any,Any}","page":"Basic function","title":"NumericalAnalysis.relative_error","text":"relative_error(x, x̂)\n\nThe relative error between X1 and X2\n\n\n\n\n\n","category":"method"},{"location":"m/basic/#NumericalAnalysis.ξ-Tuple{Any,Any}","page":"Basic function","title":"NumericalAnalysis.ξ","text":"ξ(x₀, x)\n\nreturn the number ξ(x) bwteen x₀ and x.\n\n\n\n\n\n","category":"method"},{"location":"#Numerical-Analysis","page":"Home","title":"Numerical Analysis","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Coverage) (Image: codecov) (Image: License)","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"It's a julia package for Numerical Analysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Now, i just write a little function in it. I am a beginner of Julia and numerical analysis, there may be many problems, 👏 to discuss with me.🤣","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Basic","page":"Home","title":"Basic","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"here are some function:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Now, i just write some functions:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Basic\nN the derivative, use ForwardDiff packageto calculate fracdydx, then recursive to get Nth derivative.(emmm, I feel a bit slow)\nTaylor Polynomial, get the value nth Taylor Ploynomial.\nSolutions for equation in one Variable\nBisection function, find root\nfixed_point function.\nNewton's Method\nThe Secant Method\nThe False Position Method\nModified Newton's Method\nMüller’s Method\nInterpolation and the Lagrange Polynomial\nnth Larange interpolating polynomial\nNeville’s Iterated Interpolation\nNewton’s Divided-Difference Formula","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"NthDerivative(f, x, n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"using NumericalAnalysis\nNthDerivative(x->x^3, 4, 4) # 0","category":"page"},{"location":"","page":"Home","title":"Home","text":"NthDerivative(x->x^3, 4, 2) # 24","category":"page"},{"location":"","page":"Home","title":"Home","text":"NthDerivative(cos, 1, 5) # NthDerivative(cos, 1, 5)","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"TaylorPolynomials(f, x, x₀, n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"TaylorPolynomials(cos, 0.1, 0, 6) # TaylorPolynomials(cos, 0.1, 0, 6)","category":"page"},{"location":"","page":"Home","title":"Home","text":"TaylorPolynomials(x->x^3, 1.1, 1, 6) #1.331000..","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Bisection(f, a, b)","category":"page"},{"location":"","page":"Home","title":"Home","text":"using NumericalAnalysis: SEq1\nSEq1.Bisection(sin, π/2, 3π/2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"SEq1.Bisection(x->x-1, 0, 3π/2, output_seq=true)","category":"page"}]
}
