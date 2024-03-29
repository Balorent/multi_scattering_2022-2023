## IMPORTS ##
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.special
import random
from mpl_toolkits.basemap import Basemap

def incident_wave(x, y, z, k):
    # return np.exp(1j * k * x)
    return G(k, np.sqrt(x*x + y*y + z*z))


def G(k, r): #Green function
    """
    This function return the value of the Green function G(k, r).
    """
    return (-1)/(4*np.pi)*np.exp(1j*k*r)/r


def I(k, r):
    """
    This function return the value of the regularized Green function I(k, r) = -Im[G(k, r)].
    """
    if r != 0:
        res = -np.imag(G(k, r))
    else:
        res = math.pi/2 * 4*math.pi / (2*math.pi)**3 * k
    return res


def invF(k, alpha):
    """
    This function return the value of the F used in the Foldy-Lax formalism
    """
    # Max
    res = 1/I(k, 0)
    # Hard-sphere
    # res = -I(k, 0) * G(k, alpha) / I(k, alpha)
    return res


def compute_a(k, alpha):
    """
    This function computes the vector a containing the N values of the
    coefficients a[i] seen at equation (1). To do so, it uses the system
    of equations (3).
    """
    M = np.zeros((N, N), dtype=complex)  # Matrix M(k) from equation (3)
    vecPhi = np.zeros(N, dtype=complex)  # vector phi from equation (3)
    for i in range(N):
        vecPhi[i] = incident_wave(coordinates[i][0], coordinates[i][1], coordinates[i][2], k)
        for j in range(N):
            if i != j:
                dx = coordinates[i][0] - coordinates[j][0]
                dy = coordinates[i][1] - coordinates[j][1]
                dz = coordinates[i][2] - coordinates[j][2]
                M[i, j] = -G(k, np.sqrt(dx**2 + dy**2 + dz**2))
            else:
                M[i, i] = invF(k, alpha)
    invM = np.linalg.inv(M)
    return invM.dot(vecPhi)


def directionnal_statistics(theta, phi, res):
    # RAPPEL : res[theta, phi]
    d_theta = theta[1] - theta[0]
    d_phi  = phi[1] - phi[0]
    norm, rho_x, rho_y, rho_z = 0, 0, 0, 0

    for t in range(len(theta)):
        for p in range(len(phi)):
            norm += d_theta * d_phi * res[t, p] * np.cos(theta[t])
            rho_x += d_theta * d_phi * res[t, p] * np.cos(theta[t])**2 * np.cos(phi[p])
            rho_y += d_theta * d_phi * res[t, p] * np.cos(theta[t])**2 * np.sin(phi[p])
            rho_z += d_theta * d_phi * res[t, p] * np.cos(theta[t]) * np.sin(theta[t])

    mean_res_length = np.sqrt(rho_x**2 + rho_y**2 + rho_z**2) / norm
    theta_mean = np.arcsin(rho_z / np.sqrt(rho_x**2 + rho_y**2 + rho_z**2))
    phi_mean = np.arctan2(rho_y, rho_x)

    # print("theta_mean = " + str(theta_mean))
    # print("phi_mean = " + str(phi_mean))
    return theta_mean, phi_mean, mean_res_length



## ATOMS POSITION ##
N = 1
k = 0.01
alpha = 1
radiusBall = 5
coordinates = []
for i in range(N):
    xi = random.random() * 2 * radiusBall - radiusBall
    yi = random.random() * 2 * radiusBall - radiusBall
    zi = random.random() * 2 * radiusBall - radiusBall
    coordinates.append([xi, yi, zi])
a = compute_a(k, alpha)

## SAMPLING ##
# ATTENTION : comme theta va de -pi à pi, il n'est pas défini comme en coordonnées sphériques usuelles,
# et les formules doivent donc etre adaptées (inverser sin et cos)
RAD = 180/np.pi
r = 10000

theta_res = 499
phi_res = 499
theta = np.linspace(-np.pi / 2, np.pi / 2, theta_res)
phi = np.linspace(-np.pi, np.pi, phi_res)
phi_mesh,theta_mesh = np.meshgrid(phi, theta)

d_theta = theta[1] - theta[0]
d_phi = phi[1] - phi[0]

x_contour = r * np.cos(theta_mesh) * np.cos(phi_mesh)
y_contour = r * np.cos(theta_mesh) * np.sin(phi_mesh)
z_contour = r * np.sin(theta_mesh)
psi = incident_wave(x_contour, y_contour, z_contour, k)
for i in range(N):
    dx = x_contour - coordinates[i][0]
    dy = y_contour - coordinates[i][1]
    dz = z_contour - coordinates[i][2]
    psi += a[i] * G(k, np.sqrt(dx*dx + dy*dy + dz*dz))
f = (psi - incident_wave(x_contour, y_contour, z_contour, k)) / G(k, r)
f_abs = np.abs(f)**2
psi_abs = np.abs(psi)**2

# # Onde sphérique (théorème optique)
dI = f_abs + 2*np.real(f)

A = k/(4*np.pi) / ( k/(4*np.pi) + k/(4*np.pi)**2 * sum(sum(dI * d_phi * d_theta * np.cos(theta_mesh))) )

psi_abs *= A
dI = A + A*dI - 1

to_plot = dI # k/(4*np.pi)**2 * (dI - 1)

print( sum(sum(dI * d_phi * d_theta * np.cos(theta_mesh))) )

# Onde plane (théorème optique)
# I_s = k*sum(sum(f_abs / (4*np.pi)**2 * d_phi * d_theta * np.cos(theta_mesh)))
# I_itf = np.imag(f[int(theta_res/2)][int(phi_res/2)])
# test = I_s + I_itf
# print("  ", I_s, "\n+ ", I_itf)
# print("= ", test)


## MOLLWEIDE 2d PLOT ##
fig = plt.figure()
ax = fig.add_subplot(1,2,2)
m = Basemap(projection='moll',lon_0=0,resolution='c') #projection='moll',lon_0=0,resolution='c'
min = to_plot.min()
max = to_plot.max() + to_plot.max()/1000
cont = m.contourf(phi_mesh * RAD, theta_mesh * RAD, to_plot, np.arange(min, max, (max - min) / 200), cmap=plt.cm.jet, latlon=True)
cbar = m.colorbar(cont)
parallels = np.arange(-90.,90.,30.)
m.drawparallels(parallels,labels=[False,False,False,False])
meridians = np.arange(-180.,180.,30.)
m.drawmeridians(meridians,labels=[False,False,False,False])


# DIRECTIONAL STATISTICS
# theta_mean, phi_mean = mean(theta, phi, res)
theta_mean2, phi_mean2, mean_res_length = directionnal_statistics(theta, phi, to_plot)
m.scatter(phi_mean2 * RAD, theta_mean2 * RAD, marker='o', color='w', latlon=True, s=50)
print("mean_res_length = " + str(mean_res_length))



## SCATTERS 3d PLOT ##
ax2 = fig.add_subplot(1,2,1, projection='3d')
for i in range(N):
    size = 5 + 100*(np.abs(a[i]) - np.abs(a).min()) / np.abs(a).max()
    ax2.scatter(coordinates[i][0], coordinates[i][1], coordinates[i][2], marker='o', color='k', s=10)
ax2.scatter(0, 0, 0, marker='x', color='r')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.axis('equal')



## PLOT SHOW ##
plt.show()