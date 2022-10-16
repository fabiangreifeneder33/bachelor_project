# Bachelor-Thesis Project

The goal of this project is to find interpolating Bézier-curves for some given
interpolation points.

To test our program we use the [main.py-script](main.py):
<br />First we have to set up some interpolation points, e.g. ten 
points on the Lemniscate of Girono equispaced in the parameter space. Then we call 
the init-method of the [class BZ](IP_Curve.py), which sets up a linear Matrix-Equation
Ax = b derived by the Bézier-representation formula $$\gamma (t) = ∑ Q~i~ * u~i~(t)$$, where 
the $$Q~i~$$'s are the yet unknown control-points of the interpolating curve $$\gamma (t)$$ and $$u~i~$$
is the i-th Bernstein-polynomial.

To solve this equation we use the Richardson-method. Depending on the parameter
method, which is set when calling the BZ-class function richardson_iter(), either 
simple, modified or cyclical Richardson-method is executed. We finish the iteration
when there is no more improvement in the max euclidean distance between the ip-points 
and their BZ-representation results (tolerance: 1e-4)

Finally, we plot our computed interpolating Bézier-curve by calling the BZ-class
function plot_results(), plot_ip_points() respectively. Further we also plot an
error-function which shows max_error by iteration