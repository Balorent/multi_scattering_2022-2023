#                                 MATHS.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    15 December 2022
# ------------------------------------------------------------------------


# Imports ----------------------------------------------------------------
import numpy as np
import math
import scipy.special
import random
import view


# Start ------------------------------------------------------------------
if __name__ == "__main__":
    print("main.py should be started instead")
    exit()


# Variables --------------------------------------------------------------

# Scatterers parameters #
N = 1                                       # Number of scatterers
radiusBall = 5                              # Radius of the ball in which the scatterers will be randomly placed
coordinates = [[0, 0] for i in range(N)]    # Coordinates of each scatterer [[x_1, y_1], [x_2, y_2], ..., [x_N, y_N]]

# Incident wave parameters #
k = 1                                      # Wave number [1/m]

# interaction potential parameters #
alpha = 1

# vector a #
a = np.zeros(N, dtype=complex)              # vector a from equation (3)


# Functions --------------------------------------------------------------
def phi_sph(x, y, k):
    """
    This function returns the value of a SPHERICAL wave given by the expression
        exp(ikr)/r
    at the point r = (x, y).
    """
    return G(k, np.sqrt(x*x + y*y))


def phi_pl(x, k):
    """
    This function returns the value of a PLANE wave given by the expression
        exp(ikx)
    at the point r = (x, ...).
    """
    return np.exp(1j*k*x)


def G(k, r): #Green function
    """
    This function return the value of the Green function G(k, r).
    """
    return -1/(2*math.pi) * scipy.special.kv(0, -1j*k*r)


def SurfaceUnitSphere(d):
    """
    This function return the value of the surface of a unitary sphere in d dimensions
    """
    return


def I(k, r):
    """
    This function return the value of the regularized Green function I(k, r) = -Im[G(k, r)].
    """
    if r != 0:
        res = -np.imag(G(k, r))
    else:
        res = 1/4
    return res


def invF(k, alpha):
    """
    This function return the value of the F used in the Foldy-Lax formalism
    """
    if view.model_type.get() == 1:  # Maximal model
        res = 1j/4
    elif view.model_type.get() == 0:  # Hard-sphere model  OR  delta-like potential if k*alpha << 1
        res = -I(k, 0) * G(k, alpha) / I(k, alpha)
    return res


def initialize_scatterers():
    """
    This function attributes to each of the N scatterers a random position inside
    a ball centered on the origin, and of radius radiusBall.
    """
    for i in range(N):
        xi = random.random() * 2*radiusBall - radiusBall
        y_max = radiusBall * np.sin(np.arccos(xi/radiusBall))
        yi = random.random() * 2*y_max - y_max
        coordinates[i] = [xi, yi]


def minimize_stddev():
    theta_contour = np.linspace(0, 2 * math.pi, 1000)
    x_contour = view.rc_value.get() * np.cos(theta_contour)
    y_contour = view.rc_value.get() * np.sin(theta_contour)
    x_grid, y_grid = np.linspace(-5, 5, 20), np.linspace(-5, 5, 20)
    min_stddev = None
    min_x, min_y = None, None
    for iy in range(len(y_grid)):
        for ix in range(len(x_grid)):
            x, y = x_grid[ix], y_grid[iy]
            if not(x == 0 and y == 0):
                coordinates[0] = [x, y]
                compute_a()
                psi = phi_sph(x_contour, y_contour, k)
                for i in range(N):
                    dx = x_contour - coordinates[i][0]
                    dy = y_contour - coordinates[i][1]
                    psi += a[i] * G(k, np.sqrt(dx * dx + dy * dy))
                psi /= phi_sph(x_contour, y_contour, k)
                psi = (np.abs(psi)) ** 2
                psi /= sum(psi)
                mean = sum(theta_contour * psi)
                stddev = np.sqrt(sum((theta_contour - mean) ** 2 * psi))
                print("(",  round(x, 2), ", ",   round(y, 2), ")  :  ", round(stddev, 4))
                if min_stddev is None or stddev < min_stddev:
                    min_stddev = stddev
                    min_x, min_y = x, y
    print("solution : ", min_x, min_y, "  :  ", min_stddev)




def compute_a():
    """
    This function computes the vector a containing the N values of the
    coefficients a[i] seen at equation (1). To do so, it uses the system
    of equations (3).
    """
    global a
    M = np.zeros((N, N), dtype=complex)  # Matrix M(k) from equation (3)
    vecPhi = np.zeros(N, dtype=complex)  # vector phi from equation (3)
    for i in range(N):
        if view.wave_type.get():
            vecPhi[i] = phi_sph(coordinates[i][0], coordinates[i][1], k)
        else:
            vecPhi[i] = phi_pl(coordinates[i][0], k)
        for j in range(N):
            if i != j:
                dx = coordinates[i][0] - coordinates[j][0]
                dy = coordinates[i][1] - coordinates[j][1]
                M[i, j] = -G(k, np.sqrt(dx**2 + dy**2))
            else:
                M[i, i] = invF(k, alpha)
    invM = np.linalg.inv(M)
    a = invM.dot(vecPhi)


def det_m(re_k, im_k):
    """
    This function computes the determinant of the matrix M
    for a certain k in the complex plane
    """
    M2 = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            if i != j:
                dx = coordinates[i][0] - coordinates[j][0]
                dy = coordinates[i][1] - coordinates[j][1]
                M2[i, j] = -G(re_k + 1j*im_k, np.sqrt(dx**2 + dy**2))
            else:
                M2[i, i] = invF(k, alpha)
    return np.linalg.det(M2)
