import scipy.special as special
import scipy.interpolate as sp
import numpy as np
import math


def bernstein_poly(k: int, n: int, t: float):
    return special.binom(n, k) * t ** k * (1 - t) ** (n - k)


def bezier_collocation_matrix(n, intvl_start=0., intvl_end=1.):
    A = np.zeros([n, n])
    decomp = np.linspace(intvl_start, intvl_end, num=n)

    for i, t in enumerate(decomp):
        for k in range(n):
            A[i, k] = bernstein_poly(k, n-1, t)
    return A


def b_spline_collocation_matrix(n, intvl_start=0., intvl_end=1., degree=3):
    A = np.zeros([n, n])
    knots = np.linspace(intvl_start, intvl_end, num=n)
    b_spline_knots = [intvl_start] * ((degree + 1) // 2) + list(knots) + [intvl_end] * ((degree + 1) // 2)

    for i, t in enumerate(knots):
        for k in range(n):
            temp = np.zeros(n)
            temp[k] = 1
            spl = sp.BSpline(np.array(b_spline_knots), temp, k=degree)
            A[i, k] = spl(t)
    return A


# heart-shaped curve as example for our target curve
def heart(t):
    return [16 * math.sin(t)**3, 13 * math.cos(t) - 5 * math.cos(2*t) - 2 * math.cos(3*t) - math.cos(4*t)]


# Lemniscate of Girono
def lemniscate(t):
    return [math.cos(t), math.cos(t) * math.sin(t)]


# Helix as a 3D-example
def helix(t):
    return [5 * math.cos(t), 5 * math.sin(t), t]
