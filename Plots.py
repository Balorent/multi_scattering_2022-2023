#                                 PLOTS.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    4 November 2022
# ------------------------------------------------------------------------

# Imports ----------------------------------------------------------------
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import math
import matplotlib.colors as colors
import matplotlib.patches as patches

import view
import maths
import controller

# Start ------------------------------------------------------------------
if __name__ == "__main__":
    print("main.py should be started instead")
    exit()


class PlotXY:
    def __init__(self, root, x_res, x_min, x_max, y_res, y_min, y_max, root_grid):
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
        self.ax = root.add_subplot(root_grid[:, 0])
        self.ax.set(adjustable='box', aspect='equal')
        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_xlabel('x [nm]')
        self.ax.set_ylabel('y [nm]')
        view.plot_canvas.draw()
        self.psi = np.zeros((x_res, y_res))
        self.pcm = self.ax.pcolormesh(self.x_mesh, self.y_mesh, self.psi, alpha=0, shading='auto')
        # self.root.colorbar(self.pcm, ax=self.ax, fraction=0.046, pad=0.04)
        self.background = view.plot_canvas.copy_from_bbox(self.ax.bbox)

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

        # Compute a
        maths.compute_a()
        k = maths.k
        a = maths.a

        # compute the square modulus of the wave function
        if view.wave_type.get():
            self.psi = maths.phi_sph(self.x_mesh, self.y_mesh, k)
        else:
            self.psi = maths.phi_pl(self.x_mesh, k)
        for i in range(maths.N):
            self.psi += a[i] * \
                 maths.G(k, np.sqrt((self.x_mesh - self.coordinates[i][0]) * (self.x_mesh - self.coordinates[i][0]) +
                                    (self.y_mesh - self.coordinates[i][1]) * (self.y_mesh - self.coordinates[i][1])))
        self.psi = (np.abs(self.psi)) ** 2
        self.psi /= (self.dx * self.dy * sum(sum(self.psi)))  # Normalization such that the integral(|psi|²) = 1

        # Define the bounds of the min and max bounds for the color legend
        self.scale_min = self.psi.min()
        self.scale_max = self.psi.max()
        controller.update_textbox(view.scale_min_textbox, round(self.scale_min, 5))
        controller.update_textbox(view.scale_max_textbox, round(self.scale_max, 5))

        # Display the pcolormesh
        self.pcm = self.ax.pcolormesh(self.x_mesh, self.y_mesh, self.psi,
                                      norm=colors.LogNorm(vmin=self.scale_min, vmax=self.scale_max),
                                      cmap=self.cmap_type,
                                      shading='auto')

        # Display the scatterers
        for i in range(maths.N):
            circle = patches.Circle((self.coordinates[i][0], self.coordinates[i][1]), radius=0.1, color='black')
            self.scatterer_list.append(circle)
            self.ax.add_patch(circle)

        # Add the circular contour
        self.contour = patches.Circle((0, 0), radius=view.rc_value.get(), color='black', linewidth=0.3, fill=None)
        self.ax.add_patch(self.contour)

    def update_plot(self, update_scatterers, update_res_bound, new_x=0, new_y=0):
        scattering_amp_bool = controller.scattering_amp_bool
        phase_view_bool = controller.phase_view_bool

        # Update scatterers position if necessary
        if update_scatterers >= 0:
            self.scatterer_list[update_scatterers].center = new_x, new_y
            maths.coordinates[update_scatterers][0], maths.coordinates[update_scatterers][1] = new_x, new_y
            self.coordinates[update_scatterers][0], maths.coordinates[update_scatterers][1] = new_x, new_y

        # Compute a
        maths.compute_a()
        k = maths.k
        a = maths.a

        # Update scatterers size if necessary
        if scattering_amp_bool:
            for i in range(maths.N):
                self.scatterer_list[i].radius = abs(a[i]) / max(abs(a)) * (self.x_max - self.x_min) / 20

        # compute the square modulus of the wave function
        if view.wave_type.get():
            self.psi = maths.phi_sph(self.x_mesh, self.y_mesh, k)
        else:
            self.psi = maths.phi_pl(self.x_mesh, k)
        for i in range(maths.N):
            self.psi += a[i] * \
                 maths.G(k, np.sqrt((self.x_mesh - self.coordinates[i][0]) * (self.x_mesh - self.coordinates[i][0]) +
                                    (self.y_mesh - self.coordinates[i][1]) * (self.y_mesh - self.coordinates[i][1])))

        if not phase_view_bool:
            self.psi = (np.abs(self.psi)) ** 2
            self.psi /= (self.dx * self.dy * sum(sum(self.psi)))  # Normalization such that the integral(|psi|²) = 1
            if self.cmap_type != 'YlOrRd':
                self.pcm.set(norm=colors.LogNorm(vmin=self.scale_min, vmax=self.scale_max), cmap='YlOrRd')
                self.cmap_type = 'YlOrRd'
        else:
            self.psi = np.angle(self.psi)
            if self.cmap_type != 'twilight':
                self.pcm.set(norm=colors.Normalize(vmin=-math.pi, vmax=math.pi), cmap='twilight')
                self.cmap_type = 'twilight'

        if not math.isnan(self.psi[0][0]):
            if not update_res_bound:
                self.pcm.set_array(self.psi)
            elif not phase_view_bool:
                self.pcm = self.ax.pcolormesh(self.x_mesh, self.y_mesh, self.psi,
                                              norm=colors.LogNorm(vmin=self.scale_min, vmax=self.scale_max),
                                              cmap='YlOrRd',
                                              shading='auto')
            else:
                self.pcm = self.ax.pcolormesh(self.x_mesh, self.y_mesh, self.psi,
                                              norm=colors.Normalize(vmin=-math.pi, vmax=math.pi),
                                              cmap='twilight',
                                              shading='auto')
            view.plot_canvas.restore_region(self.background)
            self.ax.draw_artist(self.pcm)
            for i in range(maths.N):
                self.ax.draw_artist(self.scatterer_list[i])
                self.scatterer_list[i].remove()
                self.ax.add_patch(self.scatterer_list[i])
            self.ax.draw_artist(self.contour)
            view.plot_canvas.blit(self.ax.bbox)

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
    def __init__(self, root, re_k_res, re_k_min, re_k_max, im_k_res, im_k_min, im_k_max, root_grid):
        # Set the attributes
        self.root = root
        self.re_k_res = re_k_res
        self.re_k_min = re_k_min
        self.re_k_max = re_k_max
        self.im_k_res = im_k_res
        self.im_k_min = im_k_min
        self.im_k_max = im_k_max

        # Create the axis
        self.ax = root.add_subplot(root_grid[1, 1])
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
    def __init__(self, root, im_k_res, im_k_min, im_k_max, root_grid):
        # Set the attributes
        self.root = root
        self.im_k_res = im_k_res
        self.im_k_min = im_k_min
        self.im_k_max = im_k_max

        # Create the mesh
        self.im_k_mesh = np.linspace(im_k_min, im_k_max, im_k_res)
        self.det_m = np.zeros(im_k_res)

        # Create the axis
        self.ax = root.add_subplot(root_grid[1, 1])
        self.ax.set(adjustable='box')
        self.ax.set_xlim(im_k_min, im_k_max)
        self.ax.set_xlabel('Im(k) [1/nm]')
        self.ax.set_ylabel('det(M(k))')
        view.plot_canvas.draw()
        self.background = view.plot_canvas.copy_from_bbox(self.ax.bbox)
        self.line = None

    def first_plot(self):
        re_k = maths.k
        for i in range(self.im_k_res):
            self.det_m[i] = np.abs(maths.det_m(re_k, self.im_k_mesh[i]))
        self.ax.set_ylim(0, 2*self.det_m[-1])
        self.line = self.ax.plot(self.im_k_mesh, self.det_m, color='red', linewidth=1)
        self.ax.add_line(self.line[0])

    def update_plot(self):
        re_k = maths.k
        for i in range(self.im_k_res):
            self.det_m[i] = np.abs(maths.det_m(re_k, self.im_k_mesh[i]))
        if not math.isnan(self.det_m[0]):
            self.ax.set_ylim(0, 2*self.det_m[-1])
            self.line[0].set_data(self.im_k_mesh, self.det_m)
            view.plot_canvas.restore_region(self.background)
            self.ax.draw_artist(self.line[0])
            view.plot_canvas.blit(self.ax.bbox)

    def update_mesh(self, update_im_k_res, update_im_k_min, update_im_k_max):
        self.im_k_mesh = np.linspace(update_im_k_min, update_im_k_max, update_im_k_res)
        self.ax.set_xlim(update_im_k_min, update_im_k_max)
        self.im_k_res = update_im_k_res
        self.im_k_min = update_im_k_min
        self.im_k_max = update_im_k_max
        self.det_m = np.zeros(update_im_k_res)


class PlotTheta:
    def __init__(self, root, theta_res, theta_min, theta_max, root_grid):
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
        self.psi = np.zeros(theta_res)

        # Create the axis
        self.ax = root.add_subplot(root_grid[0, 1])
        self.ax.set(adjustable='box')
        self.ax.set_xlim(theta_min, theta_max)
        self.ax.set_ylim(0, 1)
        self.ax.set_xlabel('\u03B8 [rad]')
        self.ax.set_ylabel('|\u03A8|²')
        view.plot_canvas.draw()
        self.background = view.plot_canvas.copy_from_bbox(self.ax.bbox)
        self.line = None

    def first_plot(self):
        k = maths.k
        a = maths.a
        coordinates = maths.coordinates
        if view.wave_type:
            self.psi = maths.phi_sph(self.x_contour, self.y_contour, k)
        else:
            self.psi = maths.phi_pl(self.x_contour, k)
        for i in range(maths.N):
            dx = self.x_contour - coordinates[i][0]
            dy = self.y_contour - coordinates[i][1]
            self.psi += a[i] * maths.G(k, np.sqrt(dx * dx + dy * dy))
        self.psi = (np.abs(self.psi)) ** 2
        self.psi /= max(self.psi)
        if not math.isnan(self.psi[0]):
            self.line = self.ax.plot(self.theta_contour, self.psi, color='blue', linewidth=1)
            self.ax.add_line(self.line[0])
        controller.update_textbox(view.entropy_textbox, 1)
        controller.update_textbox(view.stddev_textbox, 1)

        self.psi /= sum(self.psi)
        entropy = -sum(self.psi * np.log(self.psi))
        mean = sum(self.theta_contour * self.psi)
        stddev = np.sqrt(sum((self.theta_contour - mean) ** 2 * self.psi))
        controller.update_textbox(view.entropy_textbox, round(entropy, 5))
        controller.update_textbox(view.stddev_textbox, round(stddev, 5))

    def update_plot(self):
        phase_view_bool = controller.phase_view_bool
        k = maths.k
        a = maths.a
        coordinates = maths.coordinates
        if view.wave_type.get():
            self.psi = maths.phi_sph(self.x_contour, self.y_contour, k)
        else:
            self.psi = maths.phi_pl(self.x_contour, k)
        for i in range(maths.N):
            dx = self.x_contour - coordinates[i][0]
            dy = self.y_contour - coordinates[i][1]
            self.psi += a[i] * maths.G(k, np.sqrt(dx * dx + dy * dy))
        if not phase_view_bool:
            self.psi = (np.abs(self.psi)) ** 2
            self.psi /= max(self.psi)
            if self.ax.get_ylim() != (0, 1):
                self.ax.set_ylim(0, 1)
                view.plot_canvas.draw()
        else:
            self.psi = np.angle(self.psi)
            if self.ax.get_ylim() == (0, 1):
                self.ax.set_ylim(-math.pi, math.pi)
                view.plot_canvas.draw()
        if not math.isnan(self.psi[0]):
            view.plot_canvas.restore_region(self.background)
            self.line[0].set_data(self.theta_contour, self.psi)
            view.plot_canvas.restore_region(self.background)
            self.ax.draw_artist(self.line[0])
            view.plot_canvas.blit(self.ax.bbox)

        self.psi /= sum(self.psi)
        entropy = -sum(self.psi * np.log(self.psi))
        polar_dots = self.psi * np.exp(1j * self.theta_contour)
        mean = np.angle(sum(polar_dots))
        theta_contour_2 = self.theta_contour - (mean + 2*math.pi*(mean<0)) \
                          - [2*math.pi*(0 < mean < math.pi and 0 < self.theta_contour[i] > mean + math.pi) for i in range(len(self.theta_contour))] \
                          + [2*math.pi*(math.pi < mean+2*math.pi < 2*math.pi and 0 <= self.theta_contour[i] < mean + math.pi) for i in range(len(self.theta_contour))]
        stddev = np.sqrt(sum(theta_contour_2**2 * self.psi))
        controller.update_textbox(view.entropy_textbox, round(entropy, 5))
        controller.update_textbox(view.stddev_textbox, round(stddev, 5))

    def update_mesh(self, update_res):
        if update_res:
            self.theta_contour = np.linspace(0, 2 * math.pi, view.theta_res_value)
        self.x_contour = float(view.rc_textbox.get()) * np.cos(self.theta_contour)
        self.y_contour = float(view.rc_textbox.get()) * np.sin(self.theta_contour)
        self.dtheta = abs(self.theta_contour[1] - self.theta_contour[0])
        self.theta_res = view.theta_res_value
        self.psi = np.zeros(update_res)

# mean = sum(self.theta_contour*self.psi)
# stddev = np.sqrt(sum((self.theta_contour - mean)**2 * self.psi))
