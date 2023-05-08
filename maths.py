#                                 MATHS.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    22 february 2023
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
A = 0


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


def incident_wave(x, y, k):
    if view.wave_type.get():
        res = phi_sph(x,y,k)
    else:
        res = phi_pl(x, k)
    return res


def G(k, r): #Green function
    """
    This function return the value of the Green function G(k, r).
    """
    return -1/(2*math.pi) * scipy.special.kv(0, -1j*k*r)


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

    # lim = 4
    # c = 1
    # y = np.linspace(-lim, lim, N)
    # for i, yi in enumerate(y):
    #     xi = yi**2 / (4*c) - c
    #     coordinates[i] = [xi, yi]


def compute_a():
    """
    This function computes the vector a containing the N values of the
    coefficients a[i] seen at equation (1). To do so, it uses the system
    of equations (3).
    """
    if N != 0:
        x, y = np.array(coordinates)[:,0], np.array(coordinates)[:,1]
        dx = x[:, None] - x
        dy = y[:, None] - y
        r = np.sqrt(dx ** 2 + dy ** 2) + np.identity(N)
        M = invF(k, alpha) * np.identity(N, dtype=complex) - G(k, r) * (1-np.identity(N, dtype=complex))
        vecPhi = incident_wave(x, y, k)
        res = np.linalg.solve(M, vecPhi)
    else:
        res = []
    return res


def compute_A():
    theta_res = 1000
    theta_contour = np.linspace(0, 2*np.pi, theta_res)
    r_contour = np.array([np.cos(theta_contour), np.sin(theta_contour)]).transpose()
    coord = np.array(coordinates).reshape((N, 2))
    s = np.dot(np.exp(-1j * k * np.dot(r_contour, coord.T)), np.array(a).reshape((N, 1)))
    integrand = np.abs(s)**2 + 2*np.real(s)
    A = 2*np.pi / (2*np.pi + sum(integrand * 2*np.pi/theta_res))
    return A


def compute_current_diff(theta_res):
    theta_contour = np.linspace(0, 2 * np.pi, theta_res)
    r_contour = np.array([np.cos(theta_contour), np.sin(theta_contour)]).transpose()
    coord = np.array(coordinates).reshape((N, 2))
    s = np.dot(np.exp(-1j * k * np.dot(r_contour, coord.T)), a.reshape((N, 1)))
    integrand = np.abs(s) ** 2 + 2 * np.real(s)
    current_diff = A + A * integrand - 1
    return np.array(current_diff).reshape(theta_res)


def compute_cross_section(theta_res):
    theta_contour = np.linspace(0, 2 * np.pi, theta_res)
    r_contour = np.array([np.cos(theta_contour), np.sin(theta_contour)]).transpose()
    coord = np.array(coordinates).reshape((N, 2))
    s = np.dot(np.exp(-1j * k * np.dot(r_contour, coord.T)), a.reshape((N, 1)))
    return np.array(1/k * np.abs(s)**2).reshape(theta_res)



def compute_psi(x_contour, y_contour, a, A):
    psi = incident_wave(x_contour, y_contour, k)
    for i in range(N):
        dx = x_contour - coordinates[i][0]
        dy = y_contour - coordinates[i][1]
        r = np.sqrt(dx**2 + dy**2)
        psi += a[i] * G(k, r)
    return psi


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


def directional_stat(theta, res):
    d_theta = theta[1] - theta[0]
    norm = sum(res * d_theta)
    mean = np.arctan2(sum(res * np.sin(theta) * d_theta), sum(res * np.cos(theta) * d_theta))
    mean_res_length = np.sqrt(sum(res * np.cos(theta) * d_theta)**2 + sum(res * np.sin(theta) * d_theta)**2) / norm
    variance = 1-mean_res_length
    ang_std_dev = np.sqrt(-2 * np.log(mean_res_length))
    return mean, mean_res_length, variance, ang_std_dev