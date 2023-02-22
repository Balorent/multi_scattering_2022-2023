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
    return -1/(2*math.pi) * (-1j*k / (2*math.pi*r))**(1/2) * scipy.special.kv(1/2, -1j*k*r)


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
    # res = 1/I(k, 0)
    # Hard-sphere
    res = -I(k, 0) * G(k, alpha) / I(k, alpha)
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


## SETUP ##
N = 30
k = 1
alpha = 1
radiusBall = 5
coordinates = [[0, 0, 0] for i in range(N)]
for i in range(N):
    xi = random.random() * 2 * radiusBall - radiusBall
    yi = random.random() * 2 * radiusBall - radiusBall
    zi = random.random() * 2 * radiusBall - radiusBall
    coordinates[i] = [xi, yi, zi]
a = compute_a(k, alpha)


## PLOT 1D ##
# theta_mesh = np.linspace(0, 2*math.pi, 1000)
# r = 1000000
# x_contour = r*np.cos(theta_mesh)
# y_contour = r*np.sin(theta_mesh)
# res = incident_wave(x_contour, k)
# for i in range(N):
#     dx = x_contour - coordinates[i][0]
#     dy = y_contour - coordinates[i][1]
#     res += a[i] * G(k, np.sqrt(dx * dx + dy * dy))
# res -= incident_wave(x_contour, k)
# res /= np.exp(1j*k*r)/r
# res = np.abs(res)**2
# # PLOTS
# fig = plt.figure()
# ax = fig.add_subplot(111)
# plt.plot(theta_mesh, res, 'k')
# plt.ylim([0, max(res)*1.1])
# plt.show()



## PLOT 2D ##
ra = np.linspace(-np.pi, np.pi, 400)
dec= np.linspace(-np.pi/2, np.pi/2, 400)
phi,theta = np.meshgrid(ra,dec)
RAD = 180/np.pi
r = 10000
# ATTENTION : comme theta va de -pi à pi, il n'est pas défini comme en coordonnées sphériques usuelles,
# et les formules doivent donc etre adaptées (inverser sin et cos)
x_contour = r*np.cos(theta)*np.cos(phi)
y_contour = r*np.cos(theta)*np.sin(phi)
z_contour = r*np.sin(theta)
res = incident_wave(x_contour, y_contour, z_contour, k)
for i in range(N):
    dx = x_contour - coordinates[i][0]
    dy = y_contour - coordinates[i][1]
    dz = z_contour - coordinates[i][2]
    res += a[i] * G(k, np.sqrt(dx*dx + dy*dy + dz*dz))
res -= incident_wave(x_contour, y_contour, z_contour, k)
res /= np.exp(1j*k*r)/r
res = np.abs(res)**2

# MOLLWEIDE (basemap)
fig = plt.figure()
ax = fig.add_subplot(1,2,2)
m = Basemap(projection='moll',lon_0=0,resolution='c') #projection='moll',lon_0=0,resolution='c'
min = (res.min()>10**(-10)) * (res.min()-10**(-10))
max = res.max()+res.max()/100
cont = m.contourf(phi*RAD, theta*RAD, res, np.arange(min, max, (max - min)/200), cmap=plt.cm.jet, latlon=True)
cbar = m.colorbar(cont)

parallels = np.arange(-90.,90.,30.)
m.drawparallels(parallels,labels=[False,False,False,False])
meridians = np.arange(-180.,180.,30.)
m.drawmeridians(meridians,labels=[False,False,False,False])

# SCATTERS (basemap)
ax2 = fig.add_subplot(1,2,1, projection='3d')
for i in range(N):
    size = 5 + 100*(np.abs(a[i]) - np.abs(a).min()) / np.abs(a).max()
    ax2.scatter(coordinates[i][0], coordinates[i][1], coordinates[i][2], marker='o', color='k', s=10)
ax2.scatter(0, 0, 0, marker='x', color='r')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.axis('equal')

plt.show()