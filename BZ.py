from matplotlib import pyplot as plt
from bz_functions import *
import math
from math import pi
import sys


# class BZ for Bézier representation of an interpolating curve
class BZ:

    def __init__(self, points: list):

        # degree of Bernstein-polynomials
        self.degree = len(points) - 1

        # dimension of the euclidean space which our curves lie in
        self.dim = len(points[0])

        # error: unsuitable dimension
        if self.dim not in [2, 3]:
            sys.exit("ERROR: UNSUITABLE DIMENSION")

        # interpolation-points (P_i's in Carnicer-Paper)
        self.ip_points = points

        # parameter interval (t_i's lie here)
        self.intvl = [0, 1]

        # ignore (for plotting)
        self.method_str = ""

        # derive matrix equation
        self.A = collocation_matrix(self.degree, *self.intvl)
        self.b = np.asarray(points)

        # solutions and errors of each iteration
        self.solutions = {0: np.asarray(points)}
        self.max_errors = {}
        self.iter = 0

        # parameters for modified, cyclical Richardson method respectively
        self.w, self.lambda_min, self.lambda_max = None, None, None
        self.m, self.w_iter = None, None

        self.l = 0.5

    # apply a single Richardson-Iteration where used method depends on argument method
    def richardson_iter(self, method="simple", m=10):

        # apply simple Richardson
        if method == "simple":
            self.simple_richardson()

        # apply modified (i.e. "gedämpftes") Richardson method
        elif method == "modified":
            self.modified_richardson()

        # apply cyclical Richardson method
        elif method == "cyclical":
            # length of a cycle
            self.m = m
            self.cyclical_richardson()

        else:
            sys.exit("ERROR: UNSUITABLE METHOD FOR RICHARDSON")

        self.iter += 1

        # store euclidean errors of current iteration
        self.update_errors()

    # Bézier curve of current iteration (input: float t, output: curve evaluated at t)
    def curve(self, t):
        control_points = self.solutions[self.iter]

        # coordinates of ɣ(t) get stored in list results
        result = []

        # compute ɣ(t) for each coordinate (formula: ɣ(t) = ∑ Q_i * u_i(t), where u_i is the i-th Bernstein-polynomial)
        for i in range(self.dim):
            result.append(sum(control_points[k, i] * bernstein_poly(k, self.degree, t) for k in range(self.degree+1)))

        return result

    def plot_results_2D(self, target_curve, t_0, t_n):

        # plot computed Bézier representation
        T = np.linspace(self.intvl[0], self.intvl[1], 100)
        plt_points = np.asarray([self.curve(t) for t in T])
        plt.plot(plt_points[:, 0], plt_points[:, 1], color='blue', linestyle="dashed", label="Bezier-representation")

        # plot target curve for comparison
        T = np.linspace(t_0, t_n, 100)
        plt_points = np.asarray([target_curve(t) for t in T])
        plt.plot(plt_points[:, 0], plt_points[:, 1], color='black', label="target-curve")

    def plot_results_3D(self, target_curve, t_0, t_n):

        ax = plt.axes(projection='3d')

        # plot computed Bézier representation
        T = np.linspace(self.intvl[0], self.intvl[1], 100)
        plt_points = np.asarray([self.curve(t) for t in T])
        ax.plot3D(plt_points[:, 0], plt_points[:, 1], plt_points[:, 2], color='blue', linestyle="dashed", label="Bezier-representation")

        # plot target curve for comparison
        T = np.linspace(t_0, t_n, 100)
        plt_points = np.asarray([target_curve(t) for t in T])
        ax.plot3D(plt_points[:, 0], plt_points[:, 1], plt_points[:, 2], color='black', label="target-curve")

    def plot_results(self, *args):

        if self.dim == 2:
            self.plot_results_2D(*args)

        elif self.dim == 3:
            self.plot_results_3D(*args)

        plt.legend()
        plt.title("Results after " + str(self.iter) + " iterations of " + self.method_str + " \nError: " + str(self.max_errors[self.iter]))

    def plot_control_points(self):
        # plot control points
        for p in self.solutions[self.iter]:
            plt.plot(*p, color="red", marker="o")

    def plot_ip_points(self):
        # plot interpolation points
        for p in self.ip_points:
            plt.plot(*p, color="red", marker="o")

    def plot_errors(self):
        # plot error function (error: max euclidean distance between the ip-points and their BZ-representation results)
        plt.plot(self.max_errors.keys(), self.max_errors.values(), color='blue')
        plt.xlabel("iterations")
        plt.ylabel("max euclidean error")
        plt.title("Max errors by iterations")

    def is_improving(self, tolerance):
        # returns false if there is no more improvement for the max error (tolerance: 5 iterations, eps: 1e-4)
        if self.iter < tolerance + 1:
            return True

        elif self.max_errors[self.iter - tolerance] - self.max_errors[self.iter] < 1e-4:
            return False

        return True

    def update_errors(self):
        # store max euclidean distance between the ip-points and their BZ-representation results
        euclidean_distances = [np.linalg.norm(e) for e in self.A.dot(self.solutions[self.iter]) - self.b]
        self.max_errors[self.iter] = max(euclidean_distances)

    def simple_richardson(self):
        # method-string for plotting results
        self.method_str = "Richardson"

        # compute iteration result
        Id = np.identity(self.A.shape[0])
        self.solutions[self.iter + 1] = self.b + (Id - self.A).dot(self.solutions[self.iter])

    def modified_richardson(self):
        # method-string for plotting results
        self.method_str = "modified Richardson"

        # if parameter w is not yet computed: w = 2 / (λ_min + λ_max)
        if self.w is None:
            eigen_values = np.linalg.eig(self.A)[0]
            self.w = 2 / (min(eigen_values) + max(eigen_values))

        # compute iteration result
        self.solutions[self.iter + 1] = self.solutions[self.iter] + self.w * (self.b - self.A.dot(self.solutions[self.iter]))

    def cyclical_richardson(self):
        # method-string for plotting results
        self.method_str = "cyclical Richardson"

        # if eigenvalues not yet computed, compute them
        if self.lambda_min is None:
            eigen_values = np.linalg.eig(self.A)[0]
            self.lambda_min, self.lambda_max = min(eigen_values), max(eigen_values)

        # compute w_k
        k = self.iter % self.m
        intermediate_res = self.lambda_min + self.lambda_max + (self.lambda_max - self.lambda_min) * math.cos((2 * k - 1) / (2 * self.m) * pi)
        self.w_iter = 2 / intermediate_res

        # compute iteration result
        self.solutions[self.iter + 1] = self.solutions[self.iter] + self.w_iter * (self.b - self.A.dot(self.solutions[self.iter]))

    # TODO: Herleitung
    def f(self, x):
        return 0.5 * np.dot(np.dot(self.A, x), x) - np.dot(self.b, x)

    def df(self, x):
        return np.dot(self.A, x) - self.b

    # gradient descent
    def gradient_descent_iter(self):
        # current iterate
        x_iter = self.solutions[self.iter]

        # apply gradient descent with parameter self.l
        self.solutions[self.iter + 1] = self.solutions[self.iter] - self.l * self.df(x_iter)

        self.iter += 1

        # store euclidean errors of current iteration
        self.update_errors()

        # if there is no improvement for 5 iterations we make lambda bigger
        if not self.is_improving(5):
            self.l *= 2
            print(self.iter)

