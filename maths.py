#                                 MATH.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    20 october 2022
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
N = 2                                       # Number of scatterers
radiusBall = 5                              # Radius of the ball in which the scatterers will be randomly placed
coordinates = [[0, 0] for i in range(N)]    # Coordinates of each scatterer [[x_1, y_1], [x_2, y_2], ..., [x_N, y_N]]

# Incident wave parameters #
k = 10                                      # Wave number [1/m]

# Scattering model parameter #
invFmax = 1j / 4                            # inverse of the scattering amplitude F(k)

# Matrices and vectors used in the computation #
M = np.zeros((N, N), dtype=complex)         # Matrix M(k) from equation (3)
invM = np.zeros((N, N), dtype=complex)      # inverse of the matrix M
vecPhi = np.zeros(N, dtype=complex)         # vector phi from equation (3)
a = np.zeros(N, dtype=complex)              # vector a from equation (3)


# Functions --------------------------------------------------------------
def phi_sph(x, y, k):
    """
    This function returns the value of a SPHERICAL wave given by the expression
        exp(ikr)/r
    at the point r = (x, y).
    """
    r = np.sqrt(x*x + y*y)
    return np.exp(1j*k*r)/r


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


def initialize_scatterers():
    """
    This function attributes to each of the N scatterers a random position inside
    a ball centered on the origin, and of radius radiusBall.
    """
    for i in range(N):
        xi = random.random() * 2*radiusBall - radiusBall
        yMax = radiusBall * np.sin(np.arccos(xi/radiusBall))
        yi = random.random() * 2*yMax - yMax
        coordinates[i] = [xi, yi]


def compute_a():
    """
    This function computes the vector a containing the N values of the
    coefficients a[i] seen at equation (1). To do so, it uses the system
    of equations (3).
    """
    global a
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
                M[i, i] = invFmax
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
                M2[i, i] = invFmax
    return np.linalg.det(M2)
