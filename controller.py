#                             CONTROLLER.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    22 february 2023
# ------------------------------------------------------------------------

# Imports ----------------------------------------------------------------
import matplotlib
matplotlib.use('TkAgg')
import tkinter as tk
import matplotlib.colors as colors
import matplotlib.patches as patches
import numpy as np
import math
import random
from tkinter import NSEW
import os.path
import json

import view
import maths

# Start ------------------------------------------------------------------
if __name__ == "__main__":
    print("main.py should be started instead")
    exit()

# Variables --------------------------------------------------------------
coordinates = maths.coordinates
in_range = -1
selected_scatterer = -1

asymptotic_bool = 0
scattering_amp_bool = 0
hide_scatterer_bool = 0


# Events handlers --------------------------------------------------------
def button_press_callback(event):
    """ Detects clicks and checks if in range of an editable plot """
    global in_range, selected_scatterer
    if event.inaxes == view.xy_plot.ax:
        test = 0
        for i in range(maths.N):
            x_i, y_i = event.xdata, event.ydata
            if (x_i - coordinates[i][0]) ** 2 + (y_i - coordinates[i][1]) ** 2 \
                <= view.xy_plot.scatterer_list[i].radius ** 2:
                r_i, theta_i = np.sqrt(coordinates[i][0] ** 2 + coordinates[i][1] ** 2), \
                               np.arctan2(coordinates[i][1], coordinates[i][0]) + \
                               (np.arctan2(coordinates[i][1], coordinates[i][0]) < 0) * 2 * math.pi
                in_range = i
                if selected_scatterer < maths.N:
                    view.xy_plot.scatterer_list[selected_scatterer].set(color='black')
                selected_scatterer = in_range
                view.xy_plot.scatterer_list[i].set(color='dodgerblue')
                set_rightpanel(coordinates[i][0], coordinates[i][1], r_i, theta_i)
                test = 1
        if not test and selected_scatterer >= 0:
            if selected_scatterer < maths.N:
                view.xy_plot.scatterer_list[selected_scatterer].set(color='black')
            selected_scatterer = -1
            set_rightpanel(0, 0, 0, 0)
    view.plot_canvas_XY.draw()
    view.plot_canvas_theta.draw()
    view.plot_canvas_DetM.draw()


def button_release_callback(event):
    """ Detects press release and removes any interaction with the plot """
    global in_range
    in_range = -1


def motion_notify_callback(event):
    """ Detects movement while clicking and updates the corresponding values """
    global in_range

    if event.button and event.inaxes == view.xy_plot.ax and in_range >= 0:
        x_i, y_i = event.xdata, event.ydata
        r_i, theta_i = np.sqrt(x_i ** 2 + y_i ** 2), np.arctan2(y_i, x_i) + (np.arctan2(y_i, x_i) < 0) * 2 * math.pi
        view.xy_plot.update_plot(in_range, 0, x_i, y_i)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()
        set_rightpanel(x_i, y_i, r_i, theta_i)


def set_rightpanel(x, y, r, theta):
    view.x_value.set(x)
    update_textbox(view.x_textbox, round(x, 3))
    view.y_value.set(y)
    update_textbox(view.y_textbox, round(y, 3))
    view.r_value.set(r)
    update_textbox(view.r_textbox, round(r, 3))
    view.theta_value.set(theta)
    update_textbox(view.theta_textbox, round(theta, 3))


# Text-boxes -------------------------------------------------------------
def update_textbox(textbox, value):
    textbox.delete(0, tk.END)
    textbox.insert(0, value)


def update_x_from_tb(event):
    global selected_scatterer
    if selected_scatterer >= 0:
        # new position variables
        new_x, y = float(view.x_textbox.get()), maths.coordinates[selected_scatterer][1]
        new_r, new_theta = np.sqrt(y ** 2 + new_x ** 2), np.arctan2(y, new_x) + (np.arctan2(y, new_x) < 0) * 2 * math.pi

        # update the x slider
        view.x_value.set(new_x)
        print(view.x_value.get())
        # update the r slider and textbox
        view.r_value.set(new_r)
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the theta slider and textbox
        view.theta_value.set(new_theta)
        update_textbox(view.theta_textbox, round(new_theta, 3))

        # update the plots
        view.xy_plot.update_plot(selected_scatterer, 0, new_x, y)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_y_from_tb(event):
    global selected_scatterer
    if selected_scatterer >= 0:
        # new position variables
        x, new_y = maths.coordinates[selected_scatterer][0], float(view.y_textbox.get())
        new_r, new_theta = np.sqrt(new_y ** 2 + x ** 2), np.arctan2(new_y, x) + (np.arctan2(new_y, x) < 0) * 2 * math.pi

        # update the y slider
        view.y_value.set(new_y)
        # update the r slider and textbox
        view.r_value.set(new_r)
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the theta slider and textbox
        view.theta_value.set(new_theta)
        update_textbox(view.theta_textbox, round(new_theta, 3))

        # update the plots
        view.xy_plot.update_plot(selected_scatterer, 0, x, new_y)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_r_from_tb(event):
    global selected_scatterer
    if selected_scatterer >= 0:
        # new position variables
        new_r, theta = float(view.r_textbox.get()), float(view.theta_textbox.get())
        new_x, new_y = new_r * np.cos(theta), new_r * np.sin(theta)

        # update the r slider
        view.r_value.set(new_r)
        # update the x slider and textbox
        view.x_value.set(new_x)
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the y slider and textbox
        view.y_value.set(new_y)
        update_textbox(view.y_textbox, round(new_y, 3))

        # update the plots
        view.xy_plot.update_plot(selected_scatterer, 0, new_x, new_y)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_theta_from_tb(event):
    global selected_scatterer
    if selected_scatterer >= 0:
        # new position variables
        r, new_theta = float(view.r_textbox.get()), float(view.theta_textbox.get())
        new_x, new_y = r * np.cos(new_theta), r * np.sin(new_theta)

        # update the theta slider
        view.theta_value.set(new_theta)
        # update the x slider and textbox
        view.x_value.set(new_x)
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the y slider and textbox
        view.y_value.set(new_y)
        update_textbox(view.y_textbox, round(new_y, 3))

        # update the plots
        view.xy_plot.update_plot(selected_scatterer, 0, new_x, new_y)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_k_from_tb(event):
    if float(view.k_textbox.get()) != 0:
        new_k, new_lambda = float(view.k_textbox.get()), 2*math.pi/float(view.k_textbox.get())
        view.k_value.set(new_k)
        view.lambda_value.set(new_lambda)
        maths.k = new_k
        update_textbox(view.lambda_textbox, round(new_lambda, 7))
        view.xy_plot.update_plot(-1, 0)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_lambda_from_tb(event):
    if float(view.lambda_textbox.get()) != 0:
        new_k, new_lambda = 2 * math.pi / float(view.lambda_textbox.get()), float(view.lambda_textbox.get())
        view.k_value.set(new_k)
        view.lambda_value.set(new_lambda)
        maths.k = new_k
        update_textbox(view.k_textbox, round(new_k, 3))
        view.xy_plot.update_plot(-1, 0)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_alpha_from_tb(event):
    if float(view.alpha_textbox.get()) != 0:
        new_alpha = float(view.alpha_textbox.get())
        view.alpha_value.set(new_alpha)
        maths.alpha = new_alpha
        view.xy_plot.update_plot(-1, 0)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_rc_from_tb(event):
    new_rc = float(view.rc_textbox.get())
    view.rc_value.set(new_rc)
    view.xy_plot.contour.radius = new_rc
    view.theta_plot.update_mesh(0)
    view.theta_plot.update_plot()
    view.xy_plot.update_plot(-1, 0)


def update_pow_arg_from_tb(event):
    new_arg = float(view.pow_scale_textbox.get())
    view.pow_scale_value = new_arg
    if view.scale_type.get() == 1:
        refresh_scale()


def update_step_arg_from_tb(event):
    new_arg = float(view.step_scale_textbox.get())
    view.step_scale_value.set(new_arg)
    if view.scale_type.get() == 2:
        refresh_scale()


def update_scale_min_from_tb(event):
    new_scale_min = float(view.scale_min_textbox.get())
    view.xy_plot.scale_min = new_scale_min
    refresh_scale()


def update_scale_max_from_tb(event):
    new_scale_max = float(view.scale_max_textbox.get())
    view.xy_plot.scale_max = new_scale_max
    refresh_scale()


# Sliders ----------------------------------------------------------------
def update_x_from_slider(value):
    global selected_scatterer
    if selected_scatterer >= 0:
        # new position variables
        new_x, y = float(value), maths.coordinates[selected_scatterer][1]
        new_r, new_theta = np.sqrt(y ** 2 + new_x ** 2), np.arctan2(y, new_x) + (np.arctan2(y, new_x) < 0) * 2 * math.pi

        # update the x textbox
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the r slider and textbox
        view.r_value.set(new_r)
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the theta slider and textbox
        view.theta_value.set(new_theta)
        update_textbox(view.theta_textbox, round(new_theta, 3))

        # update the plot
        view.xy_plot.update_plot(selected_scatterer, 0, new_x, y)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_y_from_slider(value):
    global selected_scatterer
    if selected_scatterer >= 0:
        # new position variables
        x, new_y = maths.coordinates[selected_scatterer][0], float(value)
        new_r, new_theta = np.sqrt(new_y ** 2 + x ** 2), np.arctan2(new_y, x) + (np.arctan2(new_y, x) < 0) * 2 * math.pi

        # update the y textbox
        update_textbox(view.y_textbox, round(new_y, 3))
        # update the r slider and textbox
        view.r_value.set(new_r)
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the theta slider and textbox
        view.theta_value.set(new_theta)
        update_textbox(view.theta_textbox, round(new_theta, 3))

        # update the plot
        view.xy_plot.update_plot(selected_scatterer, 0, x, new_y)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_r_from_slider(value):
    global selected_scatterer
    if selected_scatterer >= 0:
        # new position variables
        new_r, theta = float(value), float(view.theta_textbox.get())
        new_x, new_y = new_r * np.cos(theta), new_r * np.sin(theta)

        # update the r textbox
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the x slider and textbox
        view.x_value.set(new_x)
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the y slider and textbox
        view.y_value.set(new_y)
        update_textbox(view.y_textbox, round(new_y, 3))

        # update the plot
        view.xy_plot.update_plot(selected_scatterer, 0, new_x, new_y)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_theta_from_slider(value):
    global selected_scatterer
    if selected_scatterer >= 0:
        # new position variables
        r, new_theta = float(view.r_textbox.get()), float(value)
        new_x, new_y = r * np.cos(new_theta), r * np.sin(new_theta)

        # update the theta textbox
        update_textbox(view.theta_textbox, round(new_theta, 3))
        # update the x slider and textbox
        view.x_value.set(new_x)
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the y slider and textbox
        view.y_value.set(new_y)
        update_textbox(view.y_textbox, round(new_y, 3))

        # update the plot
        view.xy_plot.update_plot(selected_scatterer, 0, new_x, new_y)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_k_from_slider(value):
    if float(value) != 0:
        new_k, new_lambda = float(value), 2*math.pi/float(value)
        maths.k = new_k
        update_textbox(view.k_textbox, new_k)
        update_textbox(view.lambda_textbox, round(new_lambda, 7))
        view.lambda_value.set(new_lambda)
        view.xy_plot.update_plot(-1, 0, )
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_lambda_from_slider(value):
    if float(value) != 0:
        new_k, new_lambda = 2 * math.pi / float(value), float(value)
        maths.k = new_k
        update_textbox(view.k_textbox, round(new_k, 3))
        update_textbox(view.lambda_textbox, new_lambda)
        view.k_value.set(new_k)
        view.xy_plot.update_plot(-1, 0, )
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_alpha_from_slider(value):
    if float(value) != 0:
        new_alpha = float(value)
        maths.alpha = new_alpha
        update_textbox(view.alpha_textbox, new_alpha)
        view.xy_plot.update_plot(-1, 0, )
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def update_step_arg_from_slider(value):
    new_arg = float(value)
    update_textbox(view.step_scale_textbox, new_arg)
    if view.scale_type.get() == 2:
        refresh_scale()


def update_rc_from_slider(value):
    new_rc = float(value)
    update_textbox(view.rc_textbox, new_rc)
    view.xy_plot.contour.radius = new_rc
    view.theta_plot.update_mesh(0)
    view.theta_plot.update_plot()
    view.xy_plot.update_plot(-1, 0)


# Buttons ----------------------------------------------------------------
def add_scatterer():
    maths.N += 1
    maths.a = np.zeros(maths.N, dtype=complex)

    # Add the coordinates to maths.coordinates
    xi = random.random() * 2 * maths.radiusBall - maths.radiusBall
    y_max = maths.radiusBall * np.sin(np.arccos(xi / maths.radiusBall))
    yi = random.random() * 2 * y_max - y_max
    maths.coordinates.append([xi, yi])

    # Add the circle to the ax of xy_plot
    circle = patches.Circle((xi, yi), radius=0.1, color='black')
    view.xy_plot.scatterer_list.append(circle)
    view.xy_plot.ax.add_patch(circle)
    view.xy_plot.update_plot(-1, 0)
    view.theta_plot.update_plot()
    view.resonances_plot.update_plot()


def remove_scatterer():
    global selected_scatterer
    if selected_scatterer >= 0:
        maths.N -= 1
        maths.a = np.zeros(maths.N, dtype=complex)
        maths.coordinates.remove(maths.coordinates[selected_scatterer])
        view.xy_plot.scatterer_list[selected_scatterer].remove()
        view.xy_plot.scatterer_list.remove(view.xy_plot.scatterer_list[selected_scatterer])
        view.xy_plot.update_plot(-1, 0)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()


def remove_all():
    for i in reversed(range(maths.N)):
        maths.N = 0
        maths.a = np.zeros(maths.N, dtype=complex)
        maths.coordinates.remove(maths.coordinates[i])
        view.xy_plot.scatterer_list[i].remove()
        view.xy_plot.scatterer_list.remove(view.xy_plot.scatterer_list[i])
    view.xy_plot.update_plot(-1, 0)
    view.theta_plot.update_plot()
    view.resonances_plot.update_plot()


def plane_wave():
    view.wave_type.set(0)
    view.xy_plot.update_plot(-1, 0)
    view.theta_plot.update_plot()


def spherical_wave():
    view.wave_type.set(1)
    view.xy_plot.update_plot(-1, 0)
    view.theta_plot.update_plot()


def max_model():
    view.model_type.set(1)
    view.xy_plot.update_plot(-1, 0)
    view.theta_plot.update_plot()
    view.resonances_plot.update_plot()


def hs_model():
    view.model_type.set(0)
    view.xy_plot.update_plot(-1, 0)
    view.theta_plot.update_plot()
    view.resonances_plot.update_plot()


def scale_auto_refresh():
    if view.xy_plot.pcm is not None:
        view.xy_plot.scale_min, view.xy_plot.scale_max = view.xy_plot.pcm.get_array().min(), view.xy_plot.pcm.get_array().max()
        refresh_scale()


def refresh_scale():
    scale_min, scale_max = view.xy_plot.scale_min, view.xy_plot.scale_max
    view.xy_plot.update_plot(-1, 1)
    update_textbox(view.scale_min_textbox, round(scale_min, 5))
    update_textbox(view.scale_max_textbox, round(scale_max, 5))
    view.plot_canvas_XY.draw()


def log_scale():
    refresh_scale()


def pow_scale():
    refresh_scale()


def step_scale():
    refresh_scale()


def asymptotic_value():
    global asymptotic_bool
    asymptotic_bool = not asymptotic_bool
    view.theta_plot.update_plot()


def change_scattering_amp_visu():
    global scattering_amp_bool
    a = maths.a
    scattering_amp_bool = not scattering_amp_bool
    for i in range(maths.N):
        if scattering_amp_bool:
            view.xy_plot.scatterer_list[i].radius = abs(a[i])/max(abs(a)) * (view.xy_plot.x_max - view.xy_plot.x_min)/20
        else:
            view.xy_plot.scatterer_list[i].radius = 0.1*(view.xy_plot.x_max - view.xy_plot.x_min)/20
    view.plot_canvas_XY.draw()


def hide_scatterers():
    global hide_scatterer_bool
    hide_scatterer_bool = not hide_scatterer_bool
    for i in range(maths.N):
        if hide_scatterer_bool:
            view.xy_plot.scatterer_list[i].remove()
        else:
            view.xy_plot.ax.add_patch(view.xy_plot.scatterer_list[i])
    view.plot_canvas_XY.draw()


def total_view():
    view.xy_plot.update_plot(-1, 1)


def incident_view():
    view.xy_plot.update_plot(-1, 1)


def scattered_view():
    view.xy_plot.update_plot(-1, 1)


def modulus_view():
    view.xy_plot.update_plot(-1, 1)
    view.theta_plot.update_plot()


def phase_view():
    view.xy_plot.update_plot(-1, 1)
    view.theta_plot.update_plot()


def complex_view():
    view.xy_plot.update_plot(-1, 1)
    view.theta_plot.update_plot()


def apply_res_bounds():
    view.x_res_value = int(view.x_res_textbox.get())
    view.x_min_value = float(view.x_min_textbox.get())
    view.x_max_value = float(view.x_max_textbox.get())
    view.y_res_value = int(view.y_res_textbox.get())
    view.y_min_value = float(view.y_min_textbox.get())
    view.y_max_value = float(view.y_max_textbox.get())
    view.im_k_res_value = int(view.im_k_res_textbox.get())
    view.im_k_min_value = float(view.im_k_min_textbox.get())
    view.im_k_max_value = float(view.im_k_max_textbox.get())
    view.theta_res_value = int(view.theta_res_textbox.get())

    view.xy_plot.update_mesh(view.x_res_value, view.x_min_value, view.x_max_value,
                             view.y_res_value, view.y_min_value, view.y_max_value)
    view.theta_plot.update_mesh(view.theta_res_value)
    view.resonances_plot.update_mesh(view.im_k_res_value, view.im_k_min_value, view.im_k_max_value)

    view.xy_plot.update_plot(-1, 1)
    if view.scale_type.get() == 0:
        view.xy_plot.pcm.set_norm(colors.LogNorm(vmin=view.xy_plot.scale_min, vmax=view.xy_plot.scale_max))
    elif view.scale_type.get() == 1:
        view.xy_plot.pcm.set_norm(colors.PowerNorm(gamma=view.pow_scale_value, vmin=view.xy_plot.scale_min, vmax=view.xy_plot.scale_max))
    elif view.scale_type.get() == 2:
        view.xy_plot.pcm.set_norm(colors.BoundaryNorm(boundaries=[view.xy_plot.scale_min,
                                                                  float(view.step_scale_textbox.get()),
                                                                  view.xy_plot.scale_max], ncolors=256))
    view.theta_plot.update_plot()
    view.resonances_plot.update_plot()
    view.plot_canvas_XY.draw()
    view.plot_canvas_theta.draw()
    view.plot_canvas_DetM.draw()


def lattice_1d():
    N = int(view.n_1d_textbox.get())
    d = float(view.d_1d_textbox.get())
    angle = float(view.angle_1d_textbox.get())
    x0 = view.x0_1d_textbox.get()
    y0 = view.y0_1d_textbox.get()
    r0 = view.r0_1d_textbox.get()
    theta0 = view.theta0_1d_textbox.get()
    if ((x0 == "") and (y0 == "")) or ((r0 == "") and (theta0 == "")):
        if (x0 == "") and (y0 == ""):
            x0 = float(r0) * np.cos(float(theta0))
            y0 = float(r0) * np.sin(float(theta0))
        elif (r0 == "") and (theta0 == ""):
            x0 = float(x0)
            y0 = float(y0)
        maths.N += N
        maths.a = np.zeros(maths.N, dtype=complex)
        for i in range(N):
            xi = x0 + i*d*np.cos(angle)
            yi = y0 + i*d*np.sin(angle)
            maths.coordinates.append([xi, yi])
            circle = patches.Circle((xi, yi), radius=0.1, color='black')
            view.xy_plot.scatterer_list.append(circle)
            view.xy_plot.ax.add_patch(circle)
        view.xy_plot.update_plot(-1, 0)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()
    else:
        print("Wrong position for the 1D lattice (either provide valid x0 and y0, or valid r0 and theta0, not more "
              "neither less)")


def lattice_2d():
    Nx = int(view.nx_2d_textbox.get())
    Ny = int(view.ny_2d_textbox.get())
    dx = float(view.dx_2d_textbox.get())
    dy = float(view.dy_2d_textbox.get())
    x0 = float(view.x0_2d_textbox.get())
    y0 = float(view.y0_2d_textbox.get())
    angle = float(view.angle_2d_textbox.get())
    maths.N += Nx*Ny
    maths.a = np.zeros(maths.N, dtype=complex)
    for i in range(Ny):
        for j in range(Nx):
            xi = x0 + j * dx * np.cos(angle) - i * dy * np.sin(angle)
            yi = y0 + j * dx * np.sin(angle) + i * dy * np.cos(angle)
            maths.coordinates.append([xi, yi])
            circle = patches.Circle((xi, yi), radius=0.1, color='black')
            view.xy_plot.scatterer_list.append(circle)
            view.xy_plot.ax.add_patch(circle)
    view.xy_plot.update_plot(-1, 0)
    view.theta_plot.update_plot()
    view.resonances_plot.update_plot()


def parabola():
    # set up the parameters
    N = int(view.n_parabola_textbox.get())
    angle = float(view.angle_parabola_textbox.get())
    d = float(view.d_from_origin_parabola_textbox.get())
    step = float(view.step_parabola_textbox.get())
    maths.N += N
    maths.a = np.zeros(maths.N, dtype=complex)

    # fill the y-axis
    y = np.zeros((N))
    ind = int((N - 1) / 2)
    if N%2 != 0:
        y_0 = 0
        y[ind] = y_0
        ind_center_right = ind
        ind_center_left = ind
    else:
        y_0 = step/2
        y[ind] = -y_0
        y[ind+1] = y_0
        ind_center_left = ind
        ind_center_right = ind+1

    for i in range(int((N-1)/2)):
        roots = np.roots([1, 0, 16*d**2 - 2*y_0**2, -32*d**2 * y_0, 16*d**2 * y_0**2 + y_0**4 - 16*d**2 * step**2])
        y_1 = np.max(roots)
        y[ind_center_right+i+1] = y_1
        y[ind_center_left - i - 1] = -y_1
        y_0 = y_1

    # fill the x-axis
    list = []
    for i, yi in enumerate(y):
        xi = yi ** 2 / (4 * d) - d
        list.append([xi, yi])

    # tilt the parabola + plot it
    for i in range(N):
        xi = list[i][0]*np.cos(angle) - list[i][1]*np.sin(angle)
        yi = list[i][0]*np.sin(angle) + list[i][1]*np.cos(angle)
        maths.coordinates.append([xi, yi])
        circle = patches.Circle((xi, yi), radius=0.1, color='black')
        view.xy_plot.scatterer_list.append(circle)
        view.xy_plot.ax.add_patch(circle)
    view.xy_plot.update_plot(-1, 0)
    view.theta_plot.update_plot()
    view.resonances_plot.update_plot()



def save():
    save_name = "saves/" + view.save_title_textbox.get()
    if not (os.path.exists(save_name)):
        f = open(save_name, "w")
        f.write("==GENERAL INFORMATIONS==\n")
        f.write("Number of atoms : " + str(maths.N) + "\n")
        f.write("Positions of the atoms : " + str(maths.coordinates) + "\n")
        if view.wave_type:
            wave_str = "spherical"
        else:
            wave_str = "plane"
        f.write("Wave type : " + wave_str + "\n")
        f.write("Wave number k : " + str(maths.k) + "\n")
        if view.model_type:
            model_str = "maximal"
        else:
            model_str = "hard-sphere"
        f.write("Potential model : " + model_str + "\n")
        f.write("Scattering length alpha : " + str(maths.alpha) + "\n")
        f.write("Circle radius : " + str(view.xy_plot.contour.radius) + "\n")
        f.write("x discretization : " + str(view.xy_plot.x_min) +
                ":" + str(view.xy_plot.x_res) +
                ":" + str(view.xy_plot.x_max) + "\n")
        f.write("y discretization : " + str(view.xy_plot.y_min) +
                ":" + str(view.xy_plot.y_res) +
                ":" + str(view.xy_plot.y_max) + "\n")
        f.write("theta discretization : " + str(view.theta_plot.theta_min) +
                ":" + str(view.theta_plot.theta_res) +
                ":" + str(view.theta_plot.theta_max) + "\n")
        f.write("im(k) discretization : " + str(view.resonances_plot.im_k_min) +
                ":" + str(view.resonances_plot.im_k_res) +
                ":" + str(view.resonances_plot.im_k_max) + "\n\n")

        f.write("==XY PLOT DATA== [ format : x/y/psi(x, y) ]\n")
        for i in range(len(view.xy_plot.x_mesh)):
            for j in range(len(view.xy_plot.y_mesh[0])):
                f.write("xy_dat : " + str("{:.2e}".format(view.xy_plot.psi_abs[i, j])) + "\n")
        f.write("\n")

        f.write("==THETA PLOT DATA== [ format : theta/psi(R, theta ]\n")
        for i in range(len(view.theta_plot.theta_contour)):
            f.write("theta_dat : " + str("{:.2e}".format(view.theta_plot.to_plot[i])) + "\n")
        f.write("\n")

        f.write("==DET(M) PLOT DATA== [ format : im(k)/det(M(k)) ]\n")
        for i in range(len(view.resonances_plot.im_k_mesh)):
            f.write("detm_dat : " + str("{:.2e}".format(view.resonances_plot.det_m[i])) + "\n")

        f.close()
    else:
        print("Error : File name already exists !")


def load():
    remove_all()
    new_str = view.load_textbox.get("1.0", "end-1c")
    list = new_str.split("\n")
    new_k = float(list[0])
    new_coord = []
    for i in range(len(list)-1):
        if list[i+1] != "":
            new_coord.append(list[i+1].split(", "))
            new_coord[i] = [float(elem) for elem in new_coord[i]]

    # update k
    if new_k != 0:
        new_lambda = 2*math.pi/new_k
        maths.k = new_k
        update_textbox(view.k_textbox, new_k)
        update_textbox(view.lambda_textbox, round(new_lambda, 7))
        view.lambda_value.set(new_lambda)

    # update scatterers
    for elem in new_coord:
        maths.N += 1
        maths.a = np.zeros(maths.N, dtype=complex)

        # Add the coordinates to maths.coordinates
        xi = elem[0]
        yi = elem[1]
        maths.coordinates.append([xi, yi])

        # Add the circle to the ax of xy_plot
        circle = patches.Circle((xi, yi), radius=0.1, color='black')
        view.xy_plot.scatterer_list.append(circle)
        view.xy_plot.ax.add_patch(circle)
        view.xy_plot.update_plot(-1, 0)
        view.theta_plot.update_plot()
        view.resonances_plot.update_plot()

    print("Loaded !")