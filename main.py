from matplotlib import pyplot as plt
from ip_functions import *
from math import pi
from IP_Curve import IP_Curve

# depending on the curve we use we have to choose different parameter-intervals and numbers of control-points
intervals = {helix: [0, 6 * pi], lemniscate: [-pi / 2, 3 * pi / 2], heart: [0, 2 * pi]}
cp_numbers = {helix: 19, lemniscate: 10, heart: 20}


if __name__ == "__main__":

    # select curve and set the number of control-points n_ and the parameter interval intvl_
    curve = lemniscate
    a, b = intervals[curve]

    # choose interpolation points where a, b determine the parameter-space
    interpolation_points = create_ip_points(curve, n=10, equispaced=True, sigma=0, a=a, b=b)

    # set-up curve-representation (with "B-Spline" or "Bezier" as method)
    ip_curve = IP_Curve(interpolation_points, method="B-Spline")

    # apply Richardson-method to our curve until there is no more improvement to the max-error
    while ip_curve.iter < 10:
        ip_curve.richardson_iter(method="simple")

    # plot computed interpolating curve
    ip_curve.plot_results(curve, a, b)
    ip_curve.plot_ip_points()
    plt.show()

    # plot error function
    ip_curve.plot_errors()
    plt.show()
