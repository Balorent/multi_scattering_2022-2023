#                                 PLOTS.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    22 february 2023
# ------------------------------------------------------------------------

# Imports ----------------------------------------------------------------
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import math
import matplotlib.colors as colors
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb

import view
import maths
import controller
import complexTools

# Start ------------------------------------------------------------------
if __name__ == "__main__":
    print("main.py should be started instead")
    exit()


class PlotXY:
    def __init__(self, root, x_res, x_min, x_max, y_res, y_min, y_max):
        # Set the attributes
        self.root = root
        self.x_res = x_res
        self.x_min = x_min
        self.x_max = x_max
        self.y_res = y_res
        self.y_min = y_min
        self.y_max = y_max

        # Create the mesh
        self.x_mesh, self.y_mesh = np.meshgrid(np.linspace(x_min, x_max, x_res), np.linspace(y_min, y_max, y_res))
        self.dx = abs(self.x_mesh[1][1] - self.x_mesh[0][0])
        self.dy = abs(self.y_mesh[1][1] - self.y_mesh[0][0])

        # Create the axis
        self.ax = root.add_subplot()
        self.ax.set(adjustable='box', aspect='equal')
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_xlabel('x [d]')
        self.ax.set_ylabel('y [d]')
        view.plot_canvas_XY.draw()
        self.psi = np.zeros((x_res, y_res), dtype=complex)
        self.psi_abs = np.zeros((x_res, y_res))
        self.pcm = self.ax.pcolormesh(self.x_mesh, self.y_mesh, abs(self.psi_abs), alpha=0, shading='auto')
        # self.root.colorbar(self.pcm, ax=self.ax, fraction=0.046, pad=0.04)
        self.background = view.plot_canvas_XY.copy_from_bbox(self.ax.bbox)

        # Create the min and max bounds for the color legend
        self.scale_min = None
        self.scale_max = None
        self.cmap_type = 'YlOrRd'

        # Create the circles representing the scatterers
        self.coordinates = None
        self.scatterer_list = []

        # Create the circle representing the contour
        self.contour = None

    def first_plot(self):
        self.coordinates = maths.coordinates

        # Compute a and A
        maths.a = maths.compute_a()
        maths.A = maths.compute_A()

        # compute the square modulus of the wave function
        self.psi = maths.compute_psi(self.x_mesh, self.y_mesh, maths.a, maths.A)
        self.psi_abs = (np.abs(self.psi)) ** 2
        # self.psi_real /= (self.dx * self.dy * sum(sum(self.psi_real)))  # Normalization such that the integral(|psi|²) = 1

        # Define the bounds of the min and max bounds for the color legend
        self.scale_min = self.psi_abs.min()
        self.scale_max = self.psi_abs.max()
        controller.update_textbox(view.scale_min_textbox, round(self.scale_min, 5))
        controller.update_textbox(view.scale_max_textbox, round(self.scale_max, 5))

        # Display the pcolormesh
        self.pcm = self.ax.imshow(self.psi_abs, norm=colors.LogNorm(vmin=self.scale_min, vmax=self.scale_max),
                                  cmap=self.cmap_type,
                                  origin="lower",
                                  extent=[self.x_min, self.x_max, self.y_min, self.y_max])

        # Display the scatterers
        for i in range(maths.N):
            circle = patches.Circle((self.coordinates[i][0], self.coordinates[i][1]), radius=0.1, color='black')
            self.scatterer_list.append(circle)
            self.ax.add_patch(circle)

        # Add the circular contour
        self.contour = patches.Circle((0, 0), radius=view.rc_value.get(), color='black', linewidth=0.3, fill=None)
        self.ax.add_patch(self.contour)

    def update_plot(self, update_scatterers, update_res_bound, new_x=0, new_y=0):
        self.psi_abs = np.zeros((self.x_res, self.y_res))

        scattering_amp_bool = controller.scattering_amp_bool
        view_type = view.view_type_2

        # Update scatterers position if necessary
        if update_scatterers >= 0:
            self.scatterer_list[update_scatterers].center = new_x, new_y
            maths.coordinates[update_scatterers][0], maths.coordinates[update_scatterers][1] = new_x, new_y
            self.coordinates[update_scatterers][0], maths.coordinates[update_scatterers][1] = new_x, new_y

        # Compute a and A
        maths.a = maths.compute_a()
        maths.A = maths.compute_A()

        # Update scatterers size if necessary
        if scattering_amp_bool:
            for i in range(maths.N):
                self.scatterer_list[i].radius = abs(maths.a[i]) / max(abs(maths.a)) * (self.x_max - self.x_min) / 20

        # compute the square modulus of the wave function
        if view.view_type_1.get() == 0:
            self.psi = maths.compute_psi(self.x_mesh, self.y_mesh, maths.a, maths.A)
        elif view.view_type_1.get() == 1:
            self.psi = maths.incident_wave(self.x_mesh, self.y_mesh, maths.k)
        else:
            self.psi = maths.compute_psi(self.x_mesh, self.y_mesh, maths.a, maths.A) - maths.incident_wave(self.x_mesh, self.y_mesh, maths.k)

        if view_type.get() == 0:
            self.psi_abs = (np.abs(self.psi)) ** 2
        elif view_type.get() == 1:
            self.psi_abs = np.angle(self.psi) + 2*np.pi*(np.angle(self.psi) < 0)
        # else:
        #     self.psi_abs /= (self.dx * self.dy * sum(sum(np.abs(self.psi)**2)))  # Normalization such that the integral(|psi|²) = 1

        if not math.isnan(np.abs(self.psi_abs[0][0])):
            if not update_res_bound and view_type.get() != 2:
                self.pcm.set_array(self.psi_abs)
            elif view_type.get() == 0:
                self.pcm = self.ax.imshow(self.psi_abs,
                                          cmap='YlOrRd',
                                          origin="lower",
                                          extent=[self.x_min, self.x_max, self.y_min, self.y_max])
                if view.scale_type.get() == 0:
                    self.pcm.set_norm(colors.LogNorm(vmin=self.scale_min, vmax=self.scale_max))
                elif view.scale_type.get() == 1:
                    self.pcm.set_norm(
                        colors.PowerNorm(gamma=view.pow_scale_value, vmin=self.scale_min, vmax=self.scale_max))
                elif view.scale_type.get() == 2:
                    self.pcm.set_norm(colors.BoundaryNorm(boundaries=[self.scale_min,
                                                                      float(view.step_scale_textbox.get()),
                                                                      self.scale_max],
                                                          ncolors=256))

            elif view_type.get() == 1:
                self.pcm = self.ax.imshow(self.psi_abs, norm=colors.Normalize(vmin=0, vmax=2*math.pi),
                                          cmap='hsv',
                                          origin="lower",
                                          extent=[self.x_min, self.x_max, self.y_min, self.y_max])
            elif view_type.get() == 2:
                self.pcm = self.ax.imshow(complexTools.domaincol_c(self.psi, 0.9),
                                          origin="lower",
                                          extent=[self.x_min, self.x_max, self.y_min, self.y_max])
            view.plot_canvas_XY.restore_region(self.background)
            self.ax.draw_artist(self.pcm)
            for i in range(maths.N):
                self.ax.draw_artist(self.scatterer_list[i])
                self.scatterer_list[i].remove()
                self.ax.add_patch(self.scatterer_list[i])
            self.ax.draw_artist(self.contour)
            view.plot_canvas_XY.blit(self.ax.bbox)

    def update_mesh(self, update_x_res, update_x_min, update_x_max, update_y_res, update_y_min, update_y_max):
        self.x_mesh, self.y_mesh = np.meshgrid(np.linspace(update_x_min, update_x_max, update_x_res),
                                               np.linspace(update_y_min, update_y_max, update_y_res))
        self.dx = abs(self.x_mesh[1][1] - self.x_mesh[0][0])
        self.dy = abs(self.y_mesh[1][1] - self.y_mesh[0][0])
        for i in range(maths.N):
            self.scatterer_list[i].radius = 0.1 * (update_x_max - update_x_min) / 20
        self.ax.set_xlim(update_x_min, update_x_max)
        self.ax.set_ylim(update_y_min, update_y_max)
        self.x_res = update_x_res
        self.x_min = update_x_min
        self.x_max = update_x_max
        self.y_res = update_y_res
        self.y_min = update_y_min
        self.y_max = update_y_max


class PlotK:
    def __init__(self, root, re_k_res, re_k_min, re_k_max, im_k_res, im_k_min, im_k_max):
        # Set the attributes
        self.root = root
        self.re_k_res = re_k_res
        self.re_k_min = re_k_min
        self.re_k_max = re_k_max
        self.im_k_res = im_k_res
        self.im_k_min = im_k_min
        self.im_k_max = im_k_max

        # Create the axis
        self.ax = root.add_subplot()
        self.ax.set(adjustable='box')
        self.ax.set_xlim(re_k_min, re_k_max)
        self.ax.set_ylim(im_k_min, im_k_max)
        self.ax.set_xlabel('Re(k) [1/nm]')
        self.ax.set_ylabel('Im(k) [1/nm]')

        # Create the mesh
        self.re_k_mesh, self.im_k_mesh = np.meshgrid(np.linspace(re_k_min, re_k_max, re_k_res), np.linspace(im_k_min, im_k_max, im_k_res))

    def first_plot(self):
        print("plotK")


class PlotDetM:
    def __init__(self, root, im_k_res, im_k_min, im_k_max):
        # Set the attributes
        self.root = root
        self.im_k_res = im_k_res
        self.im_k_min = im_k_min
        self.im_k_max = im_k_max

        # Create the mesh
        self.im_k_mesh = np.linspace(im_k_min, im_k_max, im_k_res)
        self.det_m = np.zeros(im_k_res)

        # Create the axis
        self.ax = root.add_subplot()
        self.ax.set(adjustable='box')
        self.ax.set_xlim(im_k_min, im_k_max)
        self.ax.set_xlabel('$\Im(k) [1/nm]$')
        self.ax.set_ylabel('det(M(k))')
        view.plot_canvas_DetM.draw()
        self.background = view.plot_canvas_DetM.copy_from_bbox(self.ax.bbox)
        self.line = None

    def first_plot(self):
        re_k = maths.k
        for i in range(self.im_k_res):
            self.det_m[i] = np.abs(maths.det_m(re_k, self.im_k_mesh[i]))
        self.ax.set_ylim(self.det_m.min()-self.det_m.min()/10, self.det_m.max()+self.det_m.max()/10)
        self.line = self.ax.semilogy(self.im_k_mesh, self.det_m, color='red', linewidth=1)
        self.ax.add_line(self.line[0])

    def update_plot(self):
        re_k = maths.k
        for i in range(self.im_k_res):
            self.det_m[i] = np.abs(maths.det_m(re_k, self.im_k_mesh[i]))
        if not math.isnan(self.det_m[0]):
            self.ax.set_ylim(self.det_m.min()-self.det_m.min()/10, self.det_m.max()+self.det_m.max()/10)
            self.line[0].set_data(self.im_k_mesh, self.det_m)
            view.plot_canvas_DetM.restore_region(self.background)
            self.ax.draw_artist(self.line[0])
            view.plot_canvas_DetM.blit(self.ax.bbox)
            view.plot_canvas_DetM.draw()

    def update_mesh(self, update_im_k_res, update_im_k_min, update_im_k_max):
        self.im_k_mesh = np.linspace(update_im_k_min, update_im_k_max, update_im_k_res)
        self.ax.set_xlim(update_im_k_min, update_im_k_max)
        self.im_k_res = update_im_k_res
        self.im_k_min = update_im_k_min
        self.im_k_max = update_im_k_max
        self.det_m = np.zeros(update_im_k_res)


class PlotTheta:
    def __init__(self, root, theta_res, theta_min, theta_max):
        # Set the attributes
        self.root = root
        self.theta_res = theta_res
        self.theta_min = theta_min
        self.theta_max = theta_max

        # Create the meshes
        self.theta_contour = np.linspace(0, 2 * math.pi, theta_res)
        self.x_contour = view.rc_value.get() * np.cos(self.theta_contour)
        self.y_contour = view.rc_value.get() * np.sin(self.theta_contour)
        self.dtheta = abs(self.theta_contour[1] - self.theta_contour[0])
        self.f, self.psi = np.zeros(theta_res, dtype=complex), np.zeros(theta_res, dtype=complex)
        self.current_diff, self.cross_section = np.zeros(theta_res), np.zeros(theta_res)
        self.to_plot = np.zeros(theta_res)

        # Create the axis
        self.ax = root.add_subplot()
        self.ax.set(adjustable='box')
        self.ax.set_xlim(theta_min, theta_max)
        self.ax.set_ylim(0, 2)
        self.ax.set_xlabel('\u03B8 [rad]')
        self.ax.set_ylabel('$\\frac{d\u03C3}{d\u03B8}$ or $\u0394 \\, J$')
        view.plot_canvas_theta.draw()
        self.background = view.plot_canvas_theta.copy_from_bbox(self.ax.bbox)
        self.line = None

    def first_plot(self):
        self.psi = maths.compute_psi(self.x_contour, self.y_contour, maths.a, maths.A)
        self.f = ( self.psi - maths.incident_wave(self.x_contour, self.y_contour, maths.k) ) / maths.phi_sph(self.x_contour, self.y_contour, maths.k)

        if view.wave_type.get(): # Spherical wave
            self.current_diff = (np.abs(self.f))**2 + 2*np.real(self.f)
            A = 1/4 / ( 1/4 + 1/(8*np.pi) * sum(self.current_diff * 2*np.pi/self.theta_res) )
            self.current_diff = maths.A*self.current_diff + maths.A - 1
            self.ax.set_ylim(min(self.current_diff), max(self.current_diff))
            self.to_plot = self.current_diff
            # optical theorem (spherical)
            # print("Spherical optical theorem :", 1/(8*np.pi) * sum(self.current_diff * 2*np.pi/self.theta_res), "=", 0)
        else: # Plane wave
            self.cross_section = np.abs(self.f)**2
            self.ax.set_ylim(min(self.cross_section), max(self.cross_section))
            self.to_plot = self.cross_section
            # optical theorem (plane)
            # I_s = maths.k * sum(np.abs(self.f) ** 2 * 2 * np.pi / self.theta_res / (8 * np.pi * maths.k))
            # I_itf = np.imag(self.f[0])
            # print("    Plane optical theorem :", I_s, "+", I_itf, "=", 0)

        if not math.isnan(self.to_plot[0]):
            self.line = self.ax.plot(self.theta_contour, self.to_plot, color='blue', linewidth=1)
            self.ax.add_line(self.line[0])
        controller.update_textbox(view.variance_textbox, 1)
        controller.update_textbox(view.stddev_textbox, 1)

        (mean, mean_res_length, variance, ang_std_dev) = maths.directional_stat(self.theta_contour, self.to_plot)
        controller.update_textbox(view.variance_textbox, round(variance, 5))
        controller.update_textbox(view.stddev_textbox, round(ang_std_dev, 5))

    def update_plot(self):
        if not controller.asymptotic_bool:
            self.psi = maths.compute_psi(self.x_contour, self.y_contour, maths.a, maths.A)
            self.f = (self.psi - maths.incident_wave(self.x_contour, self.y_contour, maths.k)) / maths.phi_sph(self.x_contour, self.y_contour, maths.k)
            if view.wave_type.get():  # Spherical wave
                self.current_diff = np.abs(self.f) ** 2 + 2 * np.real(self.f)
                A = 1/4 / (1/4 + 1/(8*np.pi) * sum(self.current_diff * 2*np.pi/self.theta_res))
                self.current_diff = A*self.current_diff + A - 1
                self.ax.set_ylim(min(self.current_diff), max(self.current_diff))
                self.to_plot = self.current_diff
            else:  # Plane wave
                self.cross_section = np.abs(self.f) ** 2 / maths.k
                self.ax.set_ylim(min(self.cross_section), max(self.cross_section))
                self.to_plot = self.cross_section

        else:
            if view.wave_type.get(): # Spherical wave
                self.current_diff = maths.compute_current_diff(self.theta_res)
                self.ax.set_ylim(min(self.current_diff), max(self.current_diff))
                self.to_plot = self.current_diff
                # optical theorem (spherical)
                # print("Spherical optical theorem :", 1/(8*np.pi) * sum(self.current_diff * 2*np.pi/self.theta_res), "=", 0)
            else: # Plane wave
                self.cross_section = maths.compute_cross_section(self.theta_res)
                self.ax.set_ylim(min(self.cross_section), max(self.cross_section))
                self.to_plot = self.cross_section
                # optical theorem (plane)
                # I_s = maths.k * sum(np.abs(self.f)**2 * 2*np.pi/self.theta_res / (8*np.pi*maths.k))
                # I_itf = np.imag(self.f[0])
                # print("    Plane optical theorem :", I_s, "+", I_itf, "=", 0)

        view.plot_canvas_theta.draw()
        if not math.isnan(self.to_plot[0]):
            view.plot_canvas_theta.restore_region(self.background)
            self.line[0].set_data(self.theta_contour, self.to_plot)
            view.plot_canvas_theta.restore_region(self.background)
            self.ax.draw_artist(self.line[0])
            view.plot_canvas_theta.blit(self.ax.bbox)

        (mean, mean_res_length, variance, ang_std_dev) = maths.directional_stat(self.theta_contour, self.to_plot)
        controller.update_textbox(view.variance_textbox, round(variance, 5))
        controller.update_textbox(view.stddev_textbox, round(ang_std_dev, 5))

    def update_mesh(self, update_res):
        if update_res:
            self.theta_contour = np.linspace(0, 2 * math.pi, view.theta_res_value)
        self.x_contour = float(view.rc_textbox.get()) * np.cos(self.theta_contour)
        self.y_contour = float(view.rc_textbox.get()) * np.sin(self.theta_contour)
        self.dtheta = abs(self.theta_contour[1] - self.theta_contour[0])
        self.theta_res = view.theta_res_value
        self.f = np.zeros(self.theta_res, dtype=complex)
