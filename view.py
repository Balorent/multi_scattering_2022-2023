#                                 VIEW.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    22 february 2023
# ------------------------------------------------------------------------

# Imports ----------------------------------------------------------------
import matplotlib
matplotlib.use('TkAgg')
import tkinter as tk
from tkinter import ttk, HORIZONTAL, LEFT, BOTH, RIGHT, Y, VERTICAL
import matplotlib.gridspec as gridspec
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import platform

from Plots import PlotXY
from Plots import PlotK
from Plots import PlotDetM
from Plots import PlotTheta
import maths
import controller

# Start ------------------------------------------------------------------
if __name__ == "__main__":
    print("main.py should be started instead")
    exit()

# Variables --------------------------------------------------------------

# Main panels #
window = tk.Tk()
plot_panel_XY = tk.Frame(window, background="white")
plot_canvas_XY = None
plot_panel_theta = tk.Frame(window, background="white")
plot_canvas_theta = None
plot_panel_DetM = tk.Frame(window, background="white")
plot_canvas_DetM = None
control_panel = tk.Frame(window, width=2000)

# Control panel widgets #
x_slider = None
x_textbox = None
x_value = tk.DoubleVar()
y_slider = None
y_textbox = None
y_value = tk.DoubleVar()
r_slider = None
r_textbox = None
r_value = tk.DoubleVar()
theta_slider = None
theta_textbox = None
theta_value = tk.DoubleVar()
add_scatterer_button = None
remove_scatterer_button = None
remove_all_button = None
plane_wave_button = None
spherical_wave_button = None
wave_type = tk.IntVar()  # 1 if spherical, 0 if plane
wave_type.set(1)
k_slider = None
k_textbox = None
k_value = tk.DoubleVar()
k_value.set(maths.k)
lambda_slider = None
lambda_textbox = None
lambda_value = tk.DoubleVar()
lambda_value.set(2*math.pi/maths.k)
max_model_button = None
hs_model_button = None
model_type = tk.IntVar()  # 0 if hard-sphere, 1 if maximal
model_type.set(1)
alpha_slider = None
alpha_textbox = None
alpha_value = tk.DoubleVar()
alpha_value.set(1)
refresh_scale_button = None
scale_max_textbox = None
scale_min_textbox = None
scale_min_value = None
scale_max_value = None
log_scale_button = None
pow_scale_button = None
step_scale_button = None
scale_type = tk.IntVar()  # 0 if log, 1 if power, 2 if step
scale_type.set(0)
pow_scale_textbox = None
pow_scale_value = 1
step_scale_textbox = None
step_scale_slider = None
step_scale_value = tk.DoubleVar()
step_scale_value.set(0.001)
rc_slider = None
rc_textbox = None
rc_value = tk.DoubleVar()
rc_value.set(8)
asymptotic_button = None
variance_textbox = None
stddev_textbox = None
total_view_button = None
incident_view_button = None
scattered_view_button = None
scattering_amp_button = None
view_type_1 = tk.IntVar()  # 0 if total, 1 if incident, 2 if scattered
view_type_1.set(0)
hide_scatterers_button = None
modulus_view_button = None
phase_view_button = None
complex_view_button = None
view_type_2 = tk.IntVar()  # 0 if modulus, 1 if phase, 2 if complex
view_type_2.set(0)
width_cv_slider = None
center_cv_slider = None
x_res_textbox = None
x_res_value = 100
x_max_textbox = None
x_max_value = 10
x_min_textbox = None
x_min_value = -10
y_res_textbox = None
y_res_value = 100
y_max_textbox = None
y_max_value = 10
y_min_textbox = None
y_min_value = -10
im_k_res_textbox = None
im_k_res_value = 100
im_k_max_textbox = None
im_k_max_value = 0
im_k_min_textbox = None
im_k_min_value = -2
theta_res_textbox = None
theta_res_value = 1000
apply_res_bounds_button = None
lattice_1d_button = None
n_1d_textbox = None
d_1d_textbox = None
x0_1d_textbox = None
y0_1d_textbox = None
r0_1d_textbox = None
theta0_1d_textbox = None
angle_1d_textbox = None
lattice_2d_button = None
nx_2d_textbox = None
ny_2d_textbox = None
dx_2d_textbox = None
dy_2d_textbox = None
x0_2d_textbox = None
y0_2d_textbox = None
parabola_button = None
n_parabola_textbox = None
angle_parabola_textbox = None
d_from_origin_parabola_textbox = None
step_parabola_textbox = None
angle_2d_textbox = None
save_button = None
save_title_textbox = None
load_button = None
load_textbox = None

# Plots #
xy_plot = None
theta_plot = None
resonances_plot = None


def initialise():
    global plot_canvas_XY, plot_canvas_theta, plot_canvas_DetM, xy_plot, theta_plot, resonances_plot
    initialize_window()
    initialize_plot_panel()
    initialize_control_panel()
    # maths.minimize_stddev()
    maths.initialize_scatterers()
    xy_plot.first_plot()
    theta_plot.first_plot()
    resonances_plot.first_plot()

    # Loop
    plot_canvas_XY.draw()
    plot_canvas_XY.draw()
    plot_canvas_XY.draw()
    window.mainloop()


def initialize_window():
    if platform.system() == "Linux":
        window.attributes('-zoomed', True)
    elif platform.system() == "Windows":
        window.state('zoomed')
    window['bg'] = 'white'
    window.title("Alpha wave function")
    window.iconbitmap("")
    window.protocol("WM_DELETE_WINDOW", exit)


def initialize_plot_panel():
    global xy_plot, theta_plot, resonances_plot, plot_canvas_XY, plot_canvas_theta, plot_canvas_DetM

    # plot XY
    plot_figure_XY = plt.Figure(figsize=(window.winfo_screenheight() * 8 / 1080, window.winfo_screenheight() * 8 / 1080))
    plot_panel_XY = tk.Frame(window, background="white")
    plot_panel_XY.grid(row=1, column=0, rowspan=2)
    plot_canvas_XY = FigureCanvasTkAgg(plot_figure_XY, plot_panel_XY)
    plot_canvas_XY.get_tk_widget().pack(side=tk.LEFT)
    plot_canvas_XY.draw()
    xy_plot = PlotXY(plot_figure_XY, 100, -10, 10, 100, -10, 10)
    plot_figure_XY.canvas.mpl_connect('button_release_event', controller.button_release_callback)
    plot_figure_XY.canvas.mpl_connect('motion_notify_event', controller.motion_notify_callback)
    plot_figure_XY.canvas.mpl_connect('button_press_event', controller.button_press_callback)

    # plot theta
    plot_figure_theta = plt.Figure(figsize=(window.winfo_screenheight() * 7 / 1080, window.winfo_screenheight() * 2.5 / 1080))
    plot_panel_theta = tk.Frame(window, background="white")
    plot_panel_theta.grid(row=1, column=1)
    plot_canvas_theta = FigureCanvasTkAgg(plot_figure_theta, plot_panel_theta)
    plot_canvas_theta.get_tk_widget().pack(side=tk.LEFT)
    plot_canvas_theta.draw()
    theta_plot = PlotTheta(plot_figure_theta, 1000, 0, 2 * math.pi)

    # plot detM
    plot_figure_DetM = plt.Figure(figsize=(window.winfo_screenheight() * 7 / 1080, window.winfo_screenheight() * 5 / 1080))
    plot_panel_DetM = tk.Frame(window, background="white")
    plot_panel_DetM.grid(row=2, column=1)
    plot_canvas_DetM = FigureCanvasTkAgg(plot_figure_DetM, plot_panel_DetM)
    plot_canvas_DetM.get_tk_widget().pack(side=tk.LEFT)
    plot_canvas_DetM.draw()
    resonances_plot = PlotDetM(plot_figure_DetM, im_k_res_value, im_k_min_value, im_k_max_value)


def initialize_control_panel():
    global x_slider, x_textbox, x_value, y_slider, y_textbox, y_value, r_slider, r_textbox, r_value, theta_slider, \
        theta_textbox, theta_value, add_scatterer_button, remove_scatterer_button, remove_all_button, \
        plane_wave_button, spherical_wave_button, k_slider, k_textbox, lambda_slider, lambda_textbox, \
        max_model_button, hs_model_button, alpha_slider, alpha_textbox, refresh_scale_button, scale_max_textbox, \
        scale_min_textbox, log_scale_button, pow_scale_button, step_scale_button, scale_type, pow_scale_textbox, \
        step_scale_textbox, step_scale_slider, rc_slider, rc_textbox, rc_value, asymptotic_button, variance_textbox, stddev_textbox, \
        total_view_button, incident_view_button, scattered_view_button, scattering_amp_button, hide_scatterers_button, \
        modulus_view_button, phase_view_button, complex_view_button, x_res_textbox, x_res_value, x_max_textbox, \
        x_max_value, x_min_textbox, x_min_value, y_res_textbox, y_res_value, y_max_textbox, y_max_value, \
        y_min_textbox, y_min_value, theta_res_textbox, theta_res_value, apply_res_bounds_button, im_k_res_textbox, \
        im_k_max_textbox, im_k_min_textbox, lattice_1d_button, n_1d_textbox, d_1d_textbox, x0_1d_textbox, \
        y0_1d_textbox, r0_1d_textbox, theta0_1d_textbox, angle_1d_textbox, lattice_2d_button, nx_2d_textbox, \
        ny_2d_textbox, dx_2d_textbox, dy_2d_textbox, x0_2d_textbox, y0_2d_textbox, angle_2d_textbox, parabola_button, \
        save_button, n_parabola_textbox, angle_parabola_textbox, d_from_origin_parabola_textbox, \
        step_parabola_textbox, xy_plot, theta_plot, resonances_plot, save_title_textbox, load_button, load_textbox

    control_panel.grid(row=0, column=2, rowspan=4, columnspan=2, sticky=tk.NSEW)
    control_canvas = tk.Canvas(control_panel, width=375, height=1000)
    control_canvas.pack(side=LEFT, fill=BOTH)
    scroll_bar = ttk.Scrollbar(control_panel, orient=VERTICAL, command=control_canvas.yview)
    scroll_bar.pack(side=RIGHT, fill=Y)
    control_canvas.configure(yscrollcommand=scroll_bar.set)
    control_canvas.bind('<Configure>', lambda e: control_canvas.configure(scrollregion=control_canvas.bbox("all")))
    control_panel_utilities = tk.Frame(control_canvas)
    control_canvas.create_window((0, 0), window=control_panel_utilities, anchor="nw")

    row = 0

    # Coordinates
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Selected atom").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    x_slider = tk.Scale(control_panel_utilities, from_=-10, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="x [nm]", showvalue=False, variable=x_value, command=controller.update_x_from_slider)
    x_textbox = tk.Entry(control_panel_utilities, width=10)
    x_slider.grid(row=row, column=0, columnspan=2)
    x_textbox.grid(row=row, column=2)
    row += 1

    y_slider = tk.Scale(control_panel_utilities, from_=-10, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="y [nm]", showvalue=False, variable=y_value, command=controller.update_y_from_slider)
    y_textbox = tk.Entry(control_panel_utilities, width=10)
    y_slider.grid(row=row, column=0, columnspan=2)
    y_textbox.grid(row=row, column=2)
    row += 1

    r_slider = tk.Scale(control_panel_utilities, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="r [m]", showvalue=False, variable=r_value, command=controller.update_r_from_slider)
    r_textbox = tk.Entry(control_panel_utilities, width=10)
    r_slider.grid(row=row, column=0, columnspan=2)
    r_textbox.grid(row=row, column=2)
    row += 1

    theta_slider = tk.Scale(control_panel_utilities, from_=0, to=2 * math.pi, resolution=0.01, orient=tk.HORIZONTAL,
                            length=200, label="\u03B8 [deg]", showvalue=False, variable=theta_value,
                            command=controller.update_theta_from_slider)
    theta_textbox = tk.Entry(control_panel_utilities, width=10)
    theta_textbox = tk.Entry(control_panel_utilities, width=10)
    theta_slider.grid(row=row, column=0, columnspan=2)
    theta_textbox.grid(row=row, column=2)
    row += 1

    # Add/remove atoms
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Add/remove atom").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    add_scatterer_button = tk.Button(control_panel_utilities, text="Add a scatterer", command=controller.add_scatterer)
    add_scatterer_button.grid(row=row, column=0, columnspan=3)
    row += 1

    remove_scatterer_button = tk.Button(control_panel_utilities, text="Remove a scatterer",
                                        command=controller.remove_scatterer)
    remove_scatterer_button.grid(row=row, column=0, columnspan=3)
    row += 1

    remove_all_button = tk.Button(control_panel_utilities, text="Remove all", command=controller.remove_all)
    remove_all_button.grid(row=row, column=0, columnspan=3)
    row += 1

    # Incident wave
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Incident wave parameters").grid(row=row, column=0, sticky="nsew",
                                                                            columnspan=3)
    row += 1

    plane_wave_button = tk.Radiobutton(control_panel_utilities, text="Plane wave", variable=wave_type, value=0,
                                       command=controller.plane_wave)
    plane_wave_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    spherical_wave_button = tk.Radiobutton(control_panel_utilities, text="Spherical wave", variable=wave_type, value=1,
                                           command=controller.spherical_wave)
    spherical_wave_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    k_slider = tk.Scale(control_panel_utilities, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="k [1/nm]", showvalue=False, variable=k_value, command=controller.update_k_from_slider)
    k_textbox = tk.Entry(control_panel_utilities, width=10)
    k_slider.grid(row=row, column=0, columnspan=2)
    k_slider.set(10)
    k_textbox.grid(row=row, column=2)
    row += 1
    controller.update_textbox(k_textbox, round(maths.k, 5))

    lambda_slider = tk.Scale(control_panel_utilities, from_=0.5, to=10, resolution=0.1, orient=tk.HORIZONTAL, length=200,
                             label="\u03BB [nm]", showvalue=False, variable=lambda_value,
                             command=controller.update_lambda_from_slider)
    lambda_textbox = tk.Entry(control_panel_utilities, width=10)
    lambda_slider.grid(row=row, column=0, columnspan=2)
    lambda_slider.set(2*math.pi/maths.k)
    lambda_textbox.grid(row=row, column=2)
    row += 1
    controller.update_textbox(lambda_textbox, round(2*math.pi/maths.k, 5))

    # Potential parameters
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Potential model parameters").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    max_model_button = tk.Radiobutton(control_panel_utilities, text="Maximal model", variable=model_type, value=1,
                                       command=controller.max_model)
    max_model_button.grid(row=row, column=0, sticky=tk.W)
    row += 1
    hs_model_button = tk.Radiobutton(control_panel_utilities, text="Hard-sphere model", variable=model_type, value=0,
                                       command=controller.hs_model)
    hs_model_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    alpha_slider = tk.Scale(control_panel_utilities, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="\u03B1 [nm]", showvalue=False, variable=alpha_value, command=controller.update_alpha_from_slider)
    alpha_textbox = tk.Entry(control_panel_utilities, width=10)
    alpha_slider.grid(row=row, column=0, columnspan=2)
    alpha_slider.set(1)
    alpha_textbox.grid(row=row, column=2)
    row += 1
    controller.update_textbox(alpha_textbox, maths.alpha)

    # Colorbar
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Colorbar parameters").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    refresh_scale_button = tk.Button(control_panel_utilities, text="Refresh legend scale",
                                     command=controller.scale_auto_refresh)
    refresh_scale_button.grid(row=row, column=0, columnspan=3)
    row += 1

    tk.Label(control_panel_utilities, text="scale max = ").grid(row=row, column=0, sticky=tk.E)
    scale_max_textbox = tk.Entry(control_panel_utilities, width=10)
    scale_max_textbox.grid(row=row, column=2)
    row += 1

    tk.Label(control_panel_utilities, text="scale min = ").grid(row=row, column=0, sticky=tk.E)
    scale_min_textbox = tk.Entry(control_panel_utilities, width=10)
    scale_min_textbox.grid(row=row, column=2)
    row += 1

    log_scale_button = tk.Radiobutton(control_panel_utilities, text="logarithmic scale", variable=scale_type, value=0,
                                      command=controller.log_scale)
    log_scale_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    pow_scale_button = tk.Radiobutton(control_panel_utilities, text="power scale", variable=scale_type, value=1,
                                      command=controller.pow_scale)
    pow_scale_textbox = tk.Entry(control_panel_utilities, width=10)
    pow_scale_button.grid(row=row, column=0, sticky=tk.W)
    pow_scale_textbox.grid(row=row, column=2)
    row += 1
    controller.update_textbox(pow_scale_textbox, round(pow_scale_value, 5))

    step_scale_button = tk.Radiobutton(control_panel_utilities, text="step scale", variable=scale_type, value=2,
                                       command=controller.step_scale)
    step_scale_textbox = tk.Entry(control_panel_utilities, width=10)
    step_scale_slider = tk.Scale(control_panel_utilities, from_=0.0001, to=0.005, resolution=0.0001,
                                 orient=tk.HORIZONTAL, length=200, label="step", showvalue=False,
                                 variable=step_scale_value, command=controller.update_step_arg_from_slider)
    step_scale_button.grid(row=row, column=0, sticky=tk.W)
    step_scale_textbox.grid(row=row, column=2)
    row += 1
    step_scale_slider.grid(row=row, column=0, columnspan=2)
    row += 1
    controller.update_textbox(step_scale_textbox, round(step_scale_value.get(), 5))

    # Contour
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Contour parameters").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    rc_slider = tk.Scale(control_panel_utilities, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL,
                         length=200, label="Contour Radius [m]", showvalue=False, variable=rc_value,
                         command=controller.update_rc_from_slider)
    rc_textbox = tk.Entry(control_panel_utilities, width=10)
    rc_slider.grid(row=row, column=0, columnspan=2)
    rc_textbox.grid(row=row, column=2)
    row += 1
    controller.update_textbox(rc_textbox, round(rc_value.get(), 5))
    asymptotic_button = tk.Checkbutton(control_panel_utilities, text="Asymptotic value",
                                           command=controller.asymptotic_value)
    asymptotic_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    # Functions to optimize
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Directional statistics").grid(row=row, column=0, sticky="nsew",
                                                                       columnspan=3)
    row += 1

    tk.Label(control_panel_utilities, text="variance V = ").grid(row=row, column=0, pady=2, columnspan=2, sticky=tk.E)
    variance_textbox = tk.Entry(control_panel_utilities, width=10)
    variance_textbox.grid(row=row, column=2)
    row += 1

    tk.Label(control_panel_utilities, text="std dev S = ").grid(row=row, column=0, pady=2, columnspan=2, sticky=tk.E)
    stddev_textbox = tk.Entry(control_panel_utilities, width=10)
    stddev_textbox.grid(row=row, column=2)
    row += 1


    # Other visual parameters
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Other visual parameters").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    # total_view_button = None
    # incident_view_button = None
    # scattered_view_button = None

    # total_view_button, incident_view_button, scattered_view_button
    total_view_button = tk.Radiobutton(control_panel_utilities, text="Total wave", variable=view_type_1,
                                       value=0, command=controller.total_view)
    incident_view_button = tk.Radiobutton(control_panel_utilities, text="Incident wave", variable=view_type_1,
                                          value=1, command=controller.incident_view)
    scattered_view_button = tk.Radiobutton(control_panel_utilities, text="Scattered wave", variable=view_type_1,
                                           value=2, command=controller.scattered_view)
    total_view_button.grid(row=row, column=0, sticky=tk.W)
    row += 1
    incident_view_button.grid(row=row, column=0, sticky=tk.W)
    row += 1
    scattered_view_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    scattering_amp_button = tk.Checkbutton(control_panel_utilities, text="Scattering amplitude visualization",
                                           command=controller.change_scattering_amp_visu)
    scattering_amp_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    hide_scatterers_button = tk.Checkbutton(control_panel_utilities, text="Hide scatterers", command=controller.hide_scatterers)
    hide_scatterers_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    # modulus_view_button, phase_view_button, complex_view_button
    modulus_view_button = tk.Radiobutton(control_panel_utilities, text="Square modulus view", variable=view_type_2,
                                         value=0, command=controller.modulus_view)
    phase_view_button = tk.Radiobutton(control_panel_utilities, text="Phase view", variable=view_type_2,
                                       value=1, command=controller.phase_view)
    complex_view_button = tk.Radiobutton(control_panel_utilities, text="Complex view", variable=view_type_2,
                                         value=2, command=controller.complex_view)
    modulus_view_button.grid(row=row, column=0, sticky=tk.W)
    row += 1
    phase_view_button.grid(row=row, column=0, sticky=tk.W)
    row += 1
    complex_view_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    res_bound_panel = tk.Frame(control_panel_utilities)
    res_bound_panel.grid(row=row, column=0, columnspan=3)
    tk.Label(res_bound_panel, text="x res").grid(row=0, column=0, pady=2)
    x_res_textbox = tk.Entry(res_bound_panel, width=5)
    x_res_textbox.grid(row=0, column=1)
    controller.update_textbox(x_res_textbox, round(xy_plot.x_res, 5))
    tk.Label(res_bound_panel, text="x min").grid(row=0, column=2, ipadx=10)
    x_min_textbox = tk.Entry(res_bound_panel, width=5)
    x_min_textbox.grid(row=0, column=3)
    controller.update_textbox(x_min_textbox, round(xy_plot.x_min, 5))
    tk.Label(res_bound_panel, text="x max").grid(row=0, column=4, ipadx=10)
    x_max_textbox = tk.Entry(res_bound_panel, width=5)
    x_max_textbox.grid(row=0, column=5)
    controller.update_textbox(x_max_textbox, round(xy_plot.x_max, 5))

    tk.Label(res_bound_panel, text="y res").grid(row=1, column=0, pady=2)
    y_res_textbox = tk.Entry(res_bound_panel, width=5)
    y_res_textbox.grid(row=1, column=1)
    controller.update_textbox(y_res_textbox, round(xy_plot.y_res, 5))
    tk.Label(res_bound_panel, text="y min").grid(row=1, column=2, ipadx=10)
    y_min_textbox = tk.Entry(res_bound_panel, width=5)
    y_min_textbox.grid(row=1, column=3)
    controller.update_textbox(y_min_textbox, round(xy_plot.y_min, 5))
    tk.Label(res_bound_panel, text="y max").grid(row=1, column=4, ipadx=10)
    y_max_textbox = tk.Entry(res_bound_panel, width=5)
    y_max_textbox.grid(row=1, column=5)
    controller.update_textbox(y_max_textbox, round(xy_plot.y_max, 5))

    tk.Label(res_bound_panel, text="imk res").grid(row=2, column=0, pady=2)
    im_k_res_textbox = tk.Entry(res_bound_panel, width=5)
    im_k_res_textbox.grid(row=2, column=1)
    controller.update_textbox(im_k_res_textbox, round(im_k_res_value, 3))
    tk.Label(res_bound_panel, text="imk min").grid(row=2, column=2, ipadx=10)
    im_k_min_textbox = tk.Entry(res_bound_panel, width=5)
    im_k_min_textbox.grid(row=2, column=3)
    controller.update_textbox(im_k_min_textbox, round(im_k_min_value, 3))
    tk.Label(res_bound_panel, text="imk max").grid(row=2, column=4, ipadx=10)
    im_k_max_textbox = tk.Entry(res_bound_panel, width=5)
    im_k_max_textbox.grid(row=2, column=5)
    controller.update_textbox(im_k_max_textbox, round(im_k_max_value, 3))

    tk.Label(res_bound_panel, text="\u03B8 res").grid(row=3, column=0, pady=2)
    theta_res_textbox = tk.Entry(res_bound_panel, width=5)
    theta_res_textbox.grid(row=3, column=1)
    row += 1
    controller.update_textbox(theta_res_textbox, round(theta_plot.theta_res, 5))

    apply_res_bounds_button = tk.Button(control_panel_utilities, text="Apply", command=controller.apply_res_bounds)
    apply_res_bounds_button.grid(row=row, column=0, columnspan=3)
    row += 1

    # Presets
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Presets").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    presets_panel = tk.Frame(control_panel_utilities)
    presets_panel.grid(row=row, column=0, columnspan=3)
    row_presets_panel = 0

    lattice_1d_button = tk.Button(presets_panel, text="Lattice 1D", command=controller.lattice_1d)
    lattice_1d_button.grid(row=row_presets_panel, column=0, columnspan=2)
    row_presets_panel += 1
    tk.Label(presets_panel, text="N").grid(row=row_presets_panel, column=0, pady=2)
    n_1d_textbox = tk.Entry(presets_panel, width=5)
    n_1d_textbox.grid(row=row_presets_panel, column=1)
    row_presets_panel += 1
    tk.Label(presets_panel, text="d").grid(row=row_presets_panel, column=0, pady=2)
    d_1d_textbox = tk.Entry(presets_panel, width=5)
    d_1d_textbox.grid(row=row_presets_panel, column=1)
    row_presets_panel += 1
    tk.Label(presets_panel, text="x0").grid(row=row_presets_panel, column=0, pady=2)
    x0_1d_textbox = tk.Entry(presets_panel, width=5)
    x0_1d_textbox.grid(row=row_presets_panel, column=1)
    row_presets_panel += 1
    tk.Label(presets_panel, text="y0").grid(row=row_presets_panel, column=0, pady=2)
    y0_1d_textbox = tk.Entry(presets_panel, width=5)
    y0_1d_textbox.grid(row=row_presets_panel, column=1)
    row_presets_panel += 1
    tk.Label(presets_panel, text="r0").grid(row=row_presets_panel, column=0, pady=2)
    r0_1d_textbox = tk.Entry(presets_panel, width=5)
    r0_1d_textbox.grid(row=row_presets_panel, column=1)
    row_presets_panel += 1
    tk.Label(presets_panel, text="\u03B80").grid(row=row_presets_panel, column=0, pady=2)
    theta0_1d_textbox = tk.Entry(presets_panel, width=5)
    theta0_1d_textbox.grid(row=row_presets_panel, column=1)
    row_presets_panel += 1
    tk.Label(presets_panel, text="\u03B1").grid(row=row_presets_panel, column=0, pady=2)
    angle_1d_textbox = tk.Entry(presets_panel, width=5)
    angle_1d_textbox.grid(row=row_presets_panel, column=1)
    row_presets_panel += 1

    row_presets_panel = 0
    lattice_2d_button = tk.Button(presets_panel, text="Lattice 2D", command=controller.lattice_2d)
    lattice_2d_button.grid(row=row_presets_panel, column=2, columnspan=2)
    row_presets_panel += 1
    tk.Label(presets_panel, text="Nx").grid(row=row_presets_panel, column=2, pady=2)
    nx_2d_textbox = tk.Entry(presets_panel, width=5)
    nx_2d_textbox.grid(row=row_presets_panel, column=3)
    row_presets_panel += 1
    tk.Label(presets_panel, text="Ny").grid(row=row_presets_panel, column=2, pady=2)
    ny_2d_textbox = tk.Entry(presets_panel, width=5)
    ny_2d_textbox.grid(row=row_presets_panel, column=3)
    row_presets_panel += 1
    tk.Label(presets_panel, text="dx").grid(row=row_presets_panel, column=2, pady=2)
    dx_2d_textbox = tk.Entry(presets_panel, width=5)
    dx_2d_textbox.grid(row=row_presets_panel, column=3)
    row_presets_panel += 1
    tk.Label(presets_panel, text="dy").grid(row=row_presets_panel, column=2, pady=2)
    dy_2d_textbox = tk.Entry(presets_panel, width=5)
    dy_2d_textbox.grid(row=row_presets_panel, column=3)
    row_presets_panel += 1
    tk.Label(presets_panel, text="x0").grid(row=row_presets_panel, column=2, pady=2)
    x0_2d_textbox = tk.Entry(presets_panel, width=5)
    x0_2d_textbox.grid(row=row_presets_panel, column=3)
    row_presets_panel += 1
    tk.Label(presets_panel, text="y0").grid(row=row_presets_panel, column=2, pady=2)
    y0_2d_textbox = tk.Entry(presets_panel, width=5)
    y0_2d_textbox.grid(row=row_presets_panel, column=3)
    row_presets_panel += 1
    tk.Label(presets_panel, text="\u03B1").grid(row=row_presets_panel, column=2, pady=2)
    angle_2d_textbox = tk.Entry(presets_panel, width=5)
    angle_2d_textbox.grid(row=row_presets_panel, column=3)
    row_presets_panel += 1

    row_presets_panel = 0
    parabola_button = tk.Button(presets_panel, text="Parabola", command=controller.parabola)
    parabola_button.grid(row=row_presets_panel, column=4, columnspan=2)
    row_presets_panel += 1
    tk.Label(presets_panel, text="N").grid(row=row_presets_panel, column=4, pady=2)
    n_parabola_textbox = tk.Entry(presets_panel, width=5)
    n_parabola_textbox.grid(row=row_presets_panel, column=5)
    row_presets_panel += 1
    tk.Label(presets_panel, text="\u03B1").grid(row=row_presets_panel, column=4, pady=2)
    angle_parabola_textbox = tk.Entry(presets_panel, width=5)
    angle_parabola_textbox.grid(row=row_presets_panel, column=5)
    row_presets_panel += 1
    tk.Label(presets_panel, text="d").grid(row=row_presets_panel, column=4, pady=2)
    d_from_origin_parabola_textbox = tk.Entry(presets_panel, width=5)
    d_from_origin_parabola_textbox.grid(row=row_presets_panel, column=5)
    row_presets_panel += 1
    tk.Label(presets_panel, text="step").grid(row=row_presets_panel, column=4, pady=2)
    step_parabola_textbox = tk.Entry(presets_panel, width=5)
    step_parabola_textbox.grid(row=row_presets_panel, column=5)

    # Saving options
    row+=1
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Saving options").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    save_button = tk.Button(control_panel_utilities, text="Save", command=controller.save)
    save_button.grid(row=row, column=0, columnspan=3)
    row += 1
    save_title_textbox = tk.Entry(control_panel_utilities, width=20)
    save_title_textbox.grid(row=row, column=0, columnspan=3)
    row += 1

    # Loading options
    ttk.Separator(control_panel_utilities, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel_utilities, text="Loading options").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    load_textbox = tk.Text(control_panel_utilities, width=40, height=10)
    load_textbox.grid(row=row, column=0, columnspan=3)
    row += 1
    load_button = tk.Button(control_panel_utilities, text="Load", command=controller.load)
    load_button.grid(row=row, column=0, columnspan=3)
    row += 1

    # Bind textbox actions
    x_textbox.bind("<Return>", controller.update_x_from_tb)
    y_textbox.bind("<Return>", controller.update_y_from_tb)
    r_textbox.bind("<Return>", controller.update_r_from_tb)
    theta_textbox.bind("<Return>", controller.update_theta_from_tb)
    k_textbox.bind("<Return>", controller.update_k_from_tb)
    lambda_textbox.bind("<Return>", controller.update_lambda_from_tb)
    alpha_textbox.bind("<Return>", controller.update_alpha_from_tb)
    rc_textbox.bind("<Return>", controller.update_rc_from_tb)
    pow_scale_textbox.bind("<Return>", controller.update_pow_arg_from_tb)
    step_scale_textbox.bind("<Return>", controller.update_step_arg_from_tb)
    scale_min_textbox.bind("<Return>", controller.update_scale_min_from_tb)
    scale_max_textbox.bind("<Return>", controller.update_scale_max_from_tb)