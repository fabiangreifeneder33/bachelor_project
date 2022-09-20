from scipy import special
import numpy as np
import math


def bernstein_poly(k: int, n: int, t: float):
    return special.binom(n, k) * t ** k * (1 - t) ** (n - k)


def collocation_matrix(n: int, a=0, b=1):
    A = np.zeros([n + 1, n + 1])
    decomp = np.linspace(a, b, num=n + 1)

    for i, t in enumerate(decomp):
        for k in range(n + 1):
            A[i, k] = bernstein_poly(k, n, t)
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