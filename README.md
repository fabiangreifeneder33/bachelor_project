# Bachelor-Thesis Project

The goal of this project is to find interpolating B-Spline/Bézier-curves for some given
interpolation points.

To test our program we use the [main.py-script](main.py):
<br />First we have to set up some interpolation points, e.g. ten 
points on the Lemniscate of Girono equispaced in the parameter space. Then we call 
the init-method of the [class IP_Curve](IP_Curve.py), which sets up a linear Matrix-Equation
Ax = b derived by the Bézier/B-Spline-representation formula ɣ(t) = ∑ Q_i * u_i(t), where 
the Q_i's are the yet unknown control-points of the interpolating curve ɣ and u_i 
is the i-th Bernstein-polynomial, B-Spline respectively.

To solve this equation we use the Richardson-method. Depending on the parameter
method, which is set when calling the IP-class function richardson_iter(), either 
simple, modified or cyclical Richardson-method is executed. We finish the iteration
when there is no more improvement in the max euclidean distance between the ip-points 
and their IP_Curve results (tolerance: 1e-4)

Finally, we plot our computed interpolating curve by calling the IP_Curve member
function plot_results(), plot_ip_points() respectively. Further we also plot an
error-function which shows max_error by iteration