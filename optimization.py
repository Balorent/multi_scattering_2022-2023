import numpy as np
import scipy.special
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import Bounds
import os.path
import time


def directional_stat(theta, res):
    d_theta = theta[1] - theta[0]
    norm = sum(res * d_theta)
    mean = np.arctan2(sum(res * np.sin(theta) * d_theta), sum(res * np.cos(theta) * d_theta))
    mean_res_length = np.sqrt(sum(res * np.cos(theta) * d_theta)**2 + sum(res * np.sin(theta) * d_theta)**2) / norm
    # mean, mean resultant length, variance, angular standard deviation
    return mean, mean_res_length, 1-mean_res_length, np.sqrt(-2 * np.log(mean_res_length))


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


def compute_a_1(N, x):
    """
    This function computes the vector a containing the N values of the
    coefficients a[i] seen at equation (1). To do so, it uses the system
    of equations (3).
    """
    M = np.zeros((N, N), dtype=complex)  # Matrix M(k) from equation (3)
    vecPhi = [phi_sph(x[2*i], x[2*i+1], k_opti) for i in range(N)]  # vector phi from equation (3)
    for i in range(N):
        for j in range(N):
            if i != j:
                dx = x[2*i] - x[2*j]
                dy = x[2*i+1] - x[2*j+1]
                M[i, j] = -G(k_opti, np.sqrt(dx ** 2 + dy ** 2))
            else:
                M[i, i] = invF(k_opti)
    return np.linalg.inv(M).dot(vecPhi)


def compute_psi_1(x_contour, y_contour, a, N, x):
    psi = phi_sph(x_contour, y_contour, k_opti)
    for i in range(N):
        dx = x_contour - x[2*i]
        dy = y_contour - x[2*i+1]
        r = np.sqrt(dx**2 + dy**2)
        psi += a[i] * G(k_opti, r)
    return psi


def opti_geometric_1(x, k, theta_contour, x_contour, y_contour, plot_bool):
    # Optimize the position of the N scatterers
    global plot_save
    if len(x)%2 == 1:
        print("ERROR : the array provided contains an odd number of entries, while it should contain an even number (2n coordinates)")
        return 0
    else:
        N = int(len(x)/2) - 1
        a = compute_a_1(N, x)
        psi = compute_psi_1(x_contour, y_contour, a, N, x)
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


def compute_a_2(N, x):
    """
    This function computes the vector a containing the N values of the
    coefficients a[i] seen at equation (1). To do so, it uses the system
    of equations (3).
    """
    k = x[-1]
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


def compute_psi_2(x_contour, y_contour, a, N, x):
    k = x[-1]
    psi = phi_sph(x_contour, y_contour, k)
    for i in range(N):
        dx = x_contour - x[2*i]
        dy = y_contour - x[2*i+1]
        r = np.sqrt(dx**2 + dy**2)
        psi += a[i] * G(k, r)
    return psi


# def opti_geometric_2(x, theta_contour, x_contour, y_contour, plot_bool):
#     # Optimize the position of the N scatterers and of k
#     global plot_save1
#     if len(x)%2 != 1:
#         print("ERROR : the array provided contains an even number of entries, while it should contain an odd number (2n coordinates + k)")
#         return 0
#     else:
#         N = int(len(x)/2)
#         k = x[-1]
#         a = compute_a_2(N, x)
#         psi = compute_psi_2(x_contour, y_contour, a, N, x)
#         f = (psi - phi_sph(x_contour, y_contour, k)) / phi_sph(x_contour, y_contour, k)
#         current_diff = np.abs(f)**2 + 2*np.real(f)
#         A = 1/4 / (1/4 + 1/(8*np.pi) * sum(current_diff*2*np.pi / theta_res))
#         current_diff = A*current_diff + A - 1
#         if plot_bool:
#             plot_save1 = current_diff
#         mean, mean_res_length, variance, ang_std_dev = directional_stat(theta_contour, np.maximum(current_diff, 0))
#         print(N, count, ang_std_dev)
#         print("\n")
#         return ang_std_dev


def opti_geometric(x, plot_bool):
    # Optimize the position of the N scatterers and of k
    global plot_save
    if len(x)%2 != 1:
        print("ERROR : the array provided contains an even number of entries, while it should contain an odd number (2n coordinates + k)")
        return 0
    else:
        N = int(len(x)/2)
        k = x[-1]
        a = compute_a_2(N, x)
        s = np.dot(np.exp(-1j * k * np.dot(r_contour, coord.T)), a.reshape((N, 1)))
        integrand = np.abs(s) ** 2 + 2 * np.real(s)
        A = 2 * np.pi / (2 * np.pi + sum(integrand * 2 * np.pi / theta_res))
        current_diff = np.array(A + A * integrand - 1).reshape(theta_res)
        if plot_bool:
            plot_save = current_diff
        mean, mean_res_length, variance, ang_std_dev = directional_stat(theta_contour, np.maximum(current_diff, 0))
        print(N, count, ang_std_dev)
        print("\n")
        return ang_std_dev


for N in [5]:
    for count in [55]:
        theta_res = 1000
        r = 1000
        theta_contour = np.linspace(0, 2*np.pi, theta_res)
        x_contour = r * np.cos(theta_contour)
        y_contour = r * np.sin(theta_contour)
        plot_save = np.zeros(theta_res)

        x = np.random.rand(2*N)*10-5
        r_contour = np.array([np.cos(theta_contour), np.sin(theta_contour)]).transpose()
        coord = np.array(x).reshape((N, 2))
        x = np.append(x, 1)
        bounds_inf = [-5 for i in range(2*N)]
        bounds_inf = np.append(bounds_inf, 0)
        bounds_sup = [5 for i in range(2*N)]
        bounds_sup = np.append(bounds_sup, 10)
        bounds = Bounds(bounds_inf, bounds_sup)
        minim = minimize(opti_geometric, x, args=(False), method='trust-constr', bounds=bounds)
        k_opti = minim.x[-1]
        coord_opti = minim.x[:-1]

        print("min = ", minim.x)
        print("LOAD :")
        print(k_opti)

        for n in range(N):
            print(minim.x[2*n], ", ", minim.x[2*n+1], sep='')

        print("\n\n\n")

        std_dev = opti_geometric(minim.x, True)

        save_name = "opti/opti_positive_k_varies/N_" + str(N) + "/save_" + str(count)
        if not (os.path.exists(save_name)):
            f = open(save_name, "w")
        else:
            print("Error : File name already exists !")
            f = open("opti/opti_positive_k_varies/N_" + str(N) + "/temp_" + str(count), "w")

        f.write("N : " + str(N) + "\n")
        f.write("std_dev : " + str(std_dev) + "\n\n")
        f.write("\n")

        f.write("initial situation : " + "\n")
        f.write(str(x[-1]) + "\n")
        for i in range(N):
            f.write(str(x[2 * i]) + ", " + str(x[2 * i + 1]) + "\n")
        f.write("\n")

        f.write("optimized situation : " + "\n")
        f.write(str(k_opti) + "\n")
        for i in range(N):
            f.write(str(coord_opti[2 * i]) + ", " + str(coord_opti[2 * i + 1]) + "\n")
        f.write("\n")

        f.write("==THETA PLOT DATA==\n")
        for count in range(len(plot_save)):
            f.write("theta_dat : " + str("{:.5e}".format(plot_save[count])) + "\n")
        f.close()

        # plt.plot(np.linspace(0, 2 * np.pi, theta_res), plot_save)
        # section = np.linspace(0, 2 * np.pi, theta_res)
        # plt.fill_between(section, plot_save, alpha=.5)
        # plt.xlim(0, 2 * np.pi)
        # plt.ylim(min(plot_save) + min(plot_save) / 10, max(plot_save) + max(plot_save) / 10)
        # plt.show()
