#                                 MAIN.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    22 february 2023
# Description  :    This program aims at simulating the wave function of a
#                   spherical wave in a cloud chamber, and the effect of the
#                   scattering of this wave by a relatively small number of
#                   scatterers.
#
#                   In particular, emphasis was placed on providing the visual
#                   tools to determine how such a situation can lead to the
#                   well-known linear tracks observed in the experimental
#                   cloud chambers.
#
#                   The method used to compute the multi-scattered wave function
#                   is the Foldy-Lax method, in which the wave function is given
#                   by
#
#                                          N
#                     psi(r) = phi(r) + sum    [ a[i] * G(k, r | x_i) ],                (1)
#                                          i=1
#
#                   where :
#                       - psi is the wave function
#                       - phi is the incident wave
#                       - G is the green function
#                       - k is the wave number of the incident wave
#                       - x_i is the position vector of the ith scatterer
#                       - N is the number of scatterers
#
#                   The different coefficients a[i] can be found to be
#                   solutions of the equation :
#
#                                                   N
#                     a[i] = F(k) * [ phi(x_i) + sum    [ a[j] * G(k, x_i | x_j) ]].    (2)
#                                                   j!=i
#
#                   This is equivalent to say that the vector a containing
#                   the coefficients a[i] is solution of the system :
#
#                     M(k) * a = phi,                                                   (3)
#
#                   where :
#                       - M(k) is a matrix N*N
#                       - phi is the N*1 vector [phi(x_1), ..., phi(n_N)]^T
# ------------------------------------------------------------------------
import view

if __name__ == "__main__":
    view.initialise()
