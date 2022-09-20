from matplotlib import pyplot as plt
from bz_functions import *
from math import pi
from BZ import BZ


if __name__ == "__main__":

    # depending on the curve we use we have to choose different parameter-intervals and numbers of control-points
    intervals = {helix: [0, 6 * pi], lemniscate: [-pi / 2, 3 * pi / 2], heart: [0, 2 * pi]}
    cp_numbers = {helix: 19, lemniscate: 10, heart: 20}

    # select curve and set the number of control-points n_ and the parameter interval intvl_
    curve = helix
    intvl_ = intervals[curve]
    n_ = cp_numbers[curve]

    # determine interpolation points equi-spaced in the parameter space
    interpolation_points = []
    S = np.linspace(intvl_[0], intvl_[1], n_)
    for s in S:
        interpolation_points += [curve(s)]

    # set-up Bézier-representation
    bz_repr = BZ(interpolation_points)

    # apply Richardson-method to Bézier curve until there is no more improvement to the max-error
    while bz_repr.is_improving():
        bz_repr.richardson_iter(method="cyclical")

    # plot computed Bézier curve
    bz_repr.plot_results(curve, *intvl_)
    bz_repr.plot_ip_points()
    plt.show()

    # plot error function
    bz_repr.plot_errors()
    plt.show()
