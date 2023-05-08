import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds


def G(k, r): #Green function
    """
    This function return the value of the Green function G(k, r).
    """
    return -1/(2*np.pi) * scipy.special.kv(0, -1j*k*r)


def phi_sph(x, y, k):
    """
    This function returns the value of a SPHERICAL wave given by the expression
        exp(ikr)/r
    at the point r = (x, y).
    """
    return G(k, np.sqrt(x*x + y*y))


def invF(k):
    """
    This function return the value of the F used in the Foldy-Lax formalism
    """
    res = 1j/4
    return res


def compute_a(N, x):
    """
    This function computes the vector a containing the N values of the
    coefficients a[i] seen at equation (1). To do so, it uses the system
    of equations (3).
    """
    M = np.zeros((N, N), dtype=complex)  # Matrix M(k) from equation (3)
    vecPhi = [phi_sph(x[2*i], x[2*i+1], k) for i in range(N)]  # vector phi from equation (3)
    for i in range(N):
        for j in range(N):
            if i != j:
                dx = x[2*i] - x[2*j]
                dy = x[2*i+1] - x[2*j+1]
                M[i, j] = -G(k, np.sqrt(dx**2 + dy**2))
            else:
                M[i, i] = invF(k)
    return np.linalg.inv(M).dot(vecPhi)


def compute_psi(x_contour, y_contour, a, N, x):
    psi = phi_sph(x_contour, y_contour, k)
    for i in range(N):
        dx = x_contour - x[2*i]
        dy = y_contour - x[2*i+1]
        r = np.sqrt(dx**2 + dy**2)
        psi += a[i] * G(k, r)
    return psi


def directional_stat(theta, res):
    d_theta = theta[1] - theta[0]
    norm = sum(res * d_theta)
    mean = np.arctan2(sum(res * np.sin(theta) * d_theta), sum(res * np.cos(theta) * d_theta))
    mean_res_length = np.sqrt(sum(res * np.cos(theta) * d_theta)**2 + sum(res * np.sin(theta) * d_theta)**2) / norm
    # mean, mean resultant length, variance, angular standard deviation
    return mean, mean_res_length, 1-mean_res_length, np.sqrt(-2 * np.log(mean_res_length))


def opti_geometric(x, k, theta_contour, x_contour, y_contour, plot_bool):
    global plot_save
    if len(x)%2 == 1:
        print("ERROR : the array provided contains an odd number of entries, while it should contain an even number (2n coordinates)")
        return 0
    else:
        N = int(len(x)/2)
        a = compute_a(N, x)
        psi = compute_psi(x_contour, y_contour, a, N, x)
        f = (psi - phi_sph(x_contour, y_contour, k)) / phi_sph(x_contour, y_contour, k)
        current_diff = np.abs(f)**2 + 2*np.real(f)
        A = 1/4 / (1/4 + 1/(8*np.pi) * sum(current_diff*2*np.pi / theta_res))
        current_diff = A*current_diff + A - 1
        if plot_bool:
            plot_save = current_diff
        mean, mean_res_length, variance, ang_std_dev = directional_stat(theta_contour, np.maximum(current_diff, 0))
        print(ang_std_dev)
        print("\n")
        return ang_std_dev


k = 1
theta_res = 1000
r = 1000
theta_contour = np.linspace(0, 2*np.pi, theta_res)
x_contour = r * np.cos(theta_contour)
y_contour = r * np.sin(theta_contour)
plot_save = np.zeros(theta_res)

x = [-1, 0, -2, 0, -3, 0]
bounds = Bounds([-.5, -.5, -.5, -.5, -.5, -.5], [.5, .5, .5, .5, .5, .5])
minim = minimize(opti_geometric, x, args=(k, theta_contour, x_contour, y_contour, 0), method='trust-constr', bounds=bounds)
print("min = ", minim.x)
print("LOAD :")
print(k)
for n in range(int((len(x))/2)):
    print(minim.x[2*n], ", ", minim.x[2*n+1], sep='')

std_dev = opti_geometric(minim.x, k, theta_contour, x_contour, y_contour, 1)

plt.plot(np.linspace(0, 2*np.pi, theta_res), plot_save)
plt.xlim(0, 2*np.pi)
plt.ylim(min(plot_save) + min(plot_save) / 10, max(plot_save) + max(plot_save) / 10)
plt.show()


# Utiliser les containtes à bon escient : imposer que les atomes soient sur une ligne (toujours le meme y par exemple),
# et voir quelle réparition en x minimise la variance.