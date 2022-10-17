#                                 VIEW.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    16 october 2022
# ------------------------------------------------------------------------

# Imports ----------------------------------------------------------------
import tkinter as tk
from tkinter import ttk, HORIZONTAL, LEFT, BOTH, RIGHT, Y, VERTICAL
import matplotlib.gridspec as gridspec
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from Plots import PlotXY
from Plots import PlotK
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
plot_panel = tk.Frame(window, background="white")
plot_canvas = None
control_panel = tk.Frame(window)

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
plane_wave_button = None
spherical_wave_button = None
wave_type = tk.IntVar()  # 1 if spherical, 0 if plane
wave_type.set(1)
k_slider = None
k_textbox = None
k_value = tk.DoubleVar()
k_value.set(maths.k)
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
scattering_amp_button = None
hide_scatterers_button = None
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
theta_res_textbox = None
theta_res_value = 1500
apply_res_bounds_button = None
save_button = None

# Plots #
xy_plot = None
theta_plot = None
resonances_plot = None


def initialise():
    initialize_window()
    initialize_plot_panel()
    initialize_control_panel()
    maths.initialize_scatterers()
    xy_plot.first_plot()
    theta_plot.first_plot()
    resonances_plot.first_plot()

    # Loop
    window.mainloop()


def initialize_window():
    window.state('zoomed')
    window['bg'] = 'white'
    window.title("Alpha wave function")
    window.iconbitmap("")
    window.protocol("WM_DELETE_WINDOW", exit)


def initialize_plot_panel():
    global xy_plot, theta_plot, resonances_plot, plot_canvas

    plot_panel.pack(side=tk.LEFT, fill=BOTH, expand=1)
    plot_figure = plt.Figure(figsize=(16, 8))
    plot_canvas = FigureCanvasTkAgg(plot_figure, plot_panel)
    plot_canvas.get_tk_widget().pack(side=tk.LEFT)
    plot_grid = gridspec.GridSpec(2, 2, height_ratios=[1, 2], width_ratios=[3, 2])

    xy_plot = PlotXY(plot_figure, 100, -10, 10, 100, -10, 10, plot_grid)
    theta_plot = PlotTheta(plot_figure, 1000, 0, 2 * math.pi, plot_grid)
    resonances_plot = PlotK(plot_figure, 100, 0, 5, 100, -5, 0, plot_grid)

    plot_figure.canvas.mpl_connect('button_press_event', controller.button_press_callback)
    plot_figure.canvas.mpl_connect('button_release_event', controller.button_release_callback)
    plot_figure.canvas.mpl_connect('motion_notify_event', controller.motion_notify_callback)


def initialize_control_panel():
    global x_slider, x_textbox, x_value, y_slider, y_textbox, y_value, r_slider, r_textbox, r_value, theta_slider, \
        theta_textbox, theta_value, plane_wave_button, spherical_wave_button, k_slider, k_textbox, \
        refresh_scale_button, scale_max_textbox, scale_min_textbox, log_scale_button, pow_scale_button, \
        step_scale_button, scale_type, pow_scale_textbox, step_scale_textbox, step_scale_slider, rc_slider, \
        rc_textbox, rc_value, scattering_amp_button, hide_scatterers_button, x_res_textbox, x_res_value, \
        x_max_textbox, x_max_value, x_min_textbox, x_min_value, y_res_textbox, y_res_value, y_max_textbox, \
        y_max_value, y_min_textbox, y_min_value, theta_res_textbox, theta_res_value, apply_res_bounds_button, \
        save_button

    control_panel.pack(side=tk.RIGHT, fill=BOTH, expand=1)

    row = 0

    # Coordinates
    ttk.Separator(control_panel, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel, text="Coordinates of the selected atom").grid(row=row, column=0, sticky="nsew",
                                                                          columnspan=3)
    row += 1

    x_slider = tk.Scale(control_panel, from_=-10, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="x [nm]", showvalue=0, variable=x_value, command=controller.update_x_from_slider)
    x_textbox = tk.Entry(control_panel, width=10)
    x_slider.grid(row=row, column=0, columnspan=2)
    x_textbox.grid(row=row, column=2)
    row += 1

    y_slider = tk.Scale(control_panel, from_=-10, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="y [nm]", showvalue=0, variable=y_value, command=controller.update_y_from_slider)
    y_textbox = tk.Entry(control_panel, width=10)
    y_slider.grid(row=row, column=0, columnspan=2)
    y_textbox.grid(row=row, column=2)
    row += 1

    r_slider = tk.Scale(control_panel, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="r [m]", showvalue=0, variable=r_value, command=controller.update_r_from_slider)
    r_textbox = tk.Entry(control_panel, width=10)
    r_slider.grid(row=row, column=0, columnspan=2)
    r_textbox.grid(row=row, column=2)
    row += 1

    theta_slider = tk.Scale(control_panel, from_=0, to=2 * math.pi, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                            label="\u03B8 [rad]", showvalue=0, variable=theta_value,
                            command=controller.update_theta_from_slider)
    theta_textbox = tk.Entry(control_panel, width=10)
    theta_textbox = tk.Entry(control_panel, width=10)
    theta_slider.grid(row=row, column=0, columnspan=2)
    theta_textbox.grid(row=row, column=2)
    row += 1

    # Incident wave
    ttk.Separator(control_panel, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel, text="Incident wave parameters").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    plane_wave_button = tk.Radiobutton(control_panel, text="Plane wave", variable=wave_type, val=0,
                                       command=controller.plane_wave)
    plane_wave_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    spherical_wave_button = tk.Radiobutton(control_panel, text="Spherical wave", variable=wave_type, val=1,
                                           command=controller.spherical_wave)
    spherical_wave_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    k_slider = tk.Scale(control_panel, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL, length=200,
                        label="k [1/m]", showvalue=0, variable=k_value, command=controller.update_k_from_slider)
    k_textbox = tk.Entry(control_panel, width=10)
    k_slider.grid(row=row, column=0, columnspan=2)
    k_slider.set(10)
    k_textbox.grid(row=row, column=2)
    row += 1
    controller.update_textbox(k_textbox, round(maths.k, 5))

    # Colorbar
    ttk.Separator(control_panel, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel, text="Colorbar parameters").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    refresh_scale_button = tk.Button(control_panel, text="Refresh legend scale", command=controller.scale_auto_refresh)
    refresh_scale_button.grid(row=row, column=0, columnspan=3)
    row += 1

    tk.Label(control_panel, text="scale max").grid(row=row, column=0, sticky=tk.W)
    scale_max_textbox = tk.Entry(control_panel, width=10)
    scale_max_textbox.grid(row=row, column=2)
    row += 1

    tk.Label(control_panel, text="scale min").grid(row=row, column=0, sticky=tk.W)
    scale_min_textbox = tk.Entry(control_panel, width=10)
    scale_min_textbox.grid(row=row, column=2)
    row += 1

    log_scale_button = tk.Radiobutton(control_panel, text="logarithmic scale", variable=scale_type, value=0,
                                      command=controller.log_scale)
    log_scale_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    pow_scale_button = tk.Radiobutton(control_panel, text="power scale", variable=scale_type, value=1,
                                      command=controller.pow_scale)
    pow_scale_textbox = tk.Entry(control_panel, width=10)
    pow_scale_button.grid(row=row, column=0, sticky=tk.W)
    pow_scale_textbox.grid(row=row, column=2)
    row += 1
    controller.update_textbox(pow_scale_textbox, round(pow_scale_value, 5))

    step_scale_button = tk.Radiobutton(control_panel, text="step scale", variable=scale_type, value=2,
                                       command=controller.step_scale)
    step_scale_textbox = tk.Entry(control_panel, width=10)
    step_scale_slider = tk.Scale(control_panel, from_=0.0001, to=0.005, resolution=0.0001, orient=tk.HORIZONTAL,
                                 length=200, label="step", showvalue=0, variable=step_scale_value,
                                 command=controller.update_step_arg_from_slider)
    step_scale_button.grid(row=row, column=0, sticky=tk.W)
    step_scale_textbox.grid(row=row, column=2)
    row += 1
    step_scale_slider.grid(row=row, column=0, columnspan=2)
    row += 1
    controller.update_textbox(step_scale_textbox, round(step_scale_value.get(), 5))

    # Contour
    ttk.Separator(control_panel, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel, text="Contour parameters").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    rc_slider = tk.Scale(control_panel, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL,
                         length=200, label="Contour Radius [m]", showvalue=0, variable=rc_value,
                         command=controller.update_rc_from_slider)
    rc_textbox = tk.Entry(control_panel, width=10)
    rc_slider.grid(row=row, column=0, columnspan=2)
    rc_textbox.grid(row=row, column=2)
    row += 1
    controller.update_textbox(rc_textbox, round(rc_value.get(), 5))

    # Other visual parameters
    ttk.Separator(control_panel, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel, text="Other visual parameters").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    scattering_amp_button = tk.Checkbutton(control_panel, text="Scattering amplitude visualization",
                                           command=controller.change_scattering_amp_visu)
    scattering_amp_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    hide_scatterers_button = tk.Checkbutton(control_panel, text="Hide scatterers", command=controller.hide_scatterers)
    hide_scatterers_button.grid(row=row, column=0, sticky=tk.W)
    row += 1

    res_bound_panel = tk.Frame(control_panel)
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

    tk.Label(res_bound_panel, text="\u03B8 res").grid(row=2, column=0, pady=2)
    theta_res_textbox = tk.Entry(res_bound_panel, width=5)
    theta_res_textbox.grid(row=2, column=1)
    row += 1
    controller.update_textbox(theta_res_textbox, round(theta_plot.theta_res, 5))

    apply_res_bounds_button = tk.Button(control_panel, text="Apply", command=controller.apply_res_bounds)
    apply_res_bounds_button.grid(row=row, column=0, columnspan=3)
    row += 1

    # Saving options
    ttk.Separator(control_panel, orient=HORIZONTAL).grid(row=row, column=0, ipadx=150, pady=10, columnspan=3)
    row += 1
    tk.Label(control_panel, text="Saving options").grid(row=row, column=0, sticky="nsew", columnspan=3)
    row += 1

    save_button = tk.Button(control_panel, text="Save", command=controller.save)
    save_button.grid(row=row, column=0, columnspan=3)
    row += 1

    # Bind textbox actions
    x_textbox.bind("<Return>", controller.update_x_from_tb)
    y_textbox.bind("<Return>", controller.update_y_from_tb)
    r_textbox.bind("<Return>", controller.update_r_from_tb)
    theta_textbox.bind("<Return>", controller.update_theta_from_tb)
    k_textbox.bind("<Return>", controller.update_k_from_tb)
    rc_textbox.bind("<Return>", controller.update_rc_from_tb)
    pow_scale_textbox.bind("<Return>", controller.update_pow_arg_from_tb)
    step_scale_textbox.bind("<Return>", controller.update_step_arg_from_tb)
    scale_min_textbox.bind("<Return>", controller.update_scale_min_from_tb)
    scale_max_textbox.bind("<Return>", controller.update_scale_max_from_tb)
