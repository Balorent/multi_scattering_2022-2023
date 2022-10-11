#                                 VIEW.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    08 october 2022
# ------------------------------------------------------------------------

if __name__ == "__main__":
    print("main.py should be started instead")
    exit()


import tkinter as tk
from tkinter import ttk, HORIZONTAL, LEFT, BOTH, RIGHT, Y, VERTICAL
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import math

import maths
import controller
import view

x_res = 100
y_res = 100
x_min = -10
x_max = 10
y_min = -10
y_max = 10
xMesh = None
yMesh = None

t_res = 1500 #radian
rContour = 8
thetaContour = None
xContour = None
yContour = None

window = tk.Tk()
radiusScatterers = 0.1 * (x_max - x_min)/20

left_frame = tk.Frame()
figure = plt.Figure(figsize=(5, 5.9))
canvas = FigureCanvasTkAgg(figure, left_frame)
gsLP = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
axMesh = figure.add_subplot(gsLP[0])
axContour = figure.add_subplot(gsLP[1])
contour = None
contour_plot = None

right_frame = tk.Frame(window)
right_frame.pack(side=tk.RIGHT, fill=BOTH, expand=1)
right_canvas = tk.Canvas(right_frame, width=300)
right_canvas.pack(side=LEFT, fill=BOTH, expand=1)
scroll_bar = ttk.Scrollbar(right_frame, orient=VERTICAL, command=right_canvas.yview)
scroll_bar.pack(side=RIGHT, fill=Y)
right_canvas.configure(yscrollcommand=scroll_bar.set)
right_canvas.bind('<Configure>', lambda e: right_canvas.configure(scrollregion = right_canvas.bbox("all")))
right_frame_utilities = tk.Frame(right_canvas)
right_canvas.create_window((0, 0), window=right_frame_utilities, anchor="nw")
right_frame_up = tk.Frame(right_frame_utilities)
right_frame_up.grid(row=0, column=0)
right_frame_res_bound = tk.Frame(right_frame_utilities)
right_frame_res_bound.grid(row=1, column=0)
right_frame_save = tk.Frame(right_frame_utilities)
right_frame_save.grid(row=2, column=0)

scatterList = []
patchList = []

init = 0

WF = None

scale_min = None
scale_max = None

#rightPanel
x_value = tk.DoubleVar()
x_slider = None
x_textbox = None
y_value = tk.DoubleVar()
y_slider = None
y_textbox = None
r_value = tk.DoubleVar()
r_slider = None
r_textbox = None
t_value = tk.DoubleVar()
t_slider = None
t_textbox = None

plane_wave = None
spherical_wave = None
k_value = tk.DoubleVar()
k_value.set(maths.k)
k_slider = None
k_textbox = None

rc_value = tk.DoubleVar()
rc_value.set(rContour)
rc_slider = None
rc_textbox = None

scale_min_textbox = None
scale_max_textbox = None
refresh_scale = None
scale_type = "log"
logScale = None
powScale = None
powScale_textbox = None
powArg = 1
stepScale = None
stepScale_textbox = None
stepScale_value = tk.DoubleVar()
stepScale_slider = None
stepArg = 0.002
stepScale_value.set(stepArg)

scattering_amp = None
hide_scatterer = None

x_res_textbox = None
x_min_textbox = None
x_max_textbox = None

y_res_textbox = None
y_min_textbox = None
y_max_textbox = None

t_res_textbox = None

apply_res_bounds = None

save_button = None


def initialise():
    # Initialize the view
    view.initializeWindow()
    view.initializeLeftPanel()
    view.initializeRightPanel()
    maths.initialize_scatterers()
    initializeMesh()
    initializeContour()
    maths.compute_a()
    view.plotWF()
    view.plotScatterers()
    view.init = 1
    view.axMesh.add_patch(view.contour)

    # Bind plot actions with corresponding functions
    figure.canvas.mpl_connect('button_press_event', controller.button_press_callback)
    figure.canvas.mpl_connect('button_release_event', controller.button_release_callback)
    figure.canvas.mpl_connect('motion_notify_event', controller.motion_notify_callback)

    # Bind command to sliders
    x_slider.configure(command=controller.update_X_from_slider)
    y_slider.configure(command=controller.update_Y_from_slider)
    r_slider.configure(command=controller.update_R_from_slider)
    t_slider.configure(command=controller.update_T_from_slider)
    k_slider.configure(command=controller.update_K_from_slider)
    rc_slider.configure(command=controller.update_RC_from_slider)
    stepScale_slider.configure(command=controller.update_stepArg_from_slider)

    # # # Bind radio button actions
    plane_wave.configure(command=controller.wavePl)
    spherical_wave.configure(command=controller.waveSph)

    # Bind textbox actions
    x_textbox.bind("<Return>", controller.update_X_from_tb)
    y_textbox.bind("<Return>", controller.update_Y_from_tb)
    r_textbox.bind("<Return>", controller.update_R_from_tb)
    t_textbox.bind("<Return>", controller.update_T_from_tb)
    k_textbox.bind("<Return>", controller.update_K_from_tb)
    rc_textbox.bind("<Return>", controller.update_RC_from_tb)
    powScale_textbox.bind("<Return>", controller.update_powArg)
    stepScale_textbox.bind("<Return>", controller.update_stepArg_from_tb)
    view.scale_min_textbox.bind("<Return>", controller.update_scale_min_from_tb)
    view.scale_max_textbox.bind("<Return>", controller.update_scale_max_from_tb)

    # Loop
    window.mainloop()


def initializeWindow():
    window.title("Alpha wave function")
    window.iconbitmap("")
    window.resizable(False, False)
    window.protocol("WM_DELETE_WINDOW", exit)

    left_frame.pack(side=tk.LEFT)


def initializeLeftPanel():
    # Add figure to window
    canvas.get_tk_widget().pack(side=tk.TOP)

    # Create and initialise plots
    axMesh.set_xlim(x_min, x_max)
    axMesh.set_ylim(y_min, y_max)
    axMesh.set_xlabel('x [m]')
    axMesh.set_ylabel('y [m]')

    # Create and initialise contour
    axContour.set_xlim(0, 2 * math.pi)
    axContour.set_xlabel('\u03B8 [rad]')
    axContour.set_ylabel('|\u03A8|²')

    # Add contour to the plot
    view.contour = patches.Circle((0, 0), radius=rContour, color='black', linewidth=0.25, fill=None)
    axMesh.add_patch(view.contour)



def initializeRightPanel():
    row = 0
    #################
    ## coordinates ##
    #################
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=0, ipadx=110, pady=10)
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=1, ipadx=35, pady=10)
    row += 1
    title = tk.Label(right_frame_up, text="Coordinates of the selected atom")
    title.grid(row=row, column=0, sticky="nsew")
    row += 1

    #x-y sliders
    view.x_slider = tk.Scale(right_frame_up, from_=-10, to=10, resolution=0.01, orient=tk.HORIZONTAL,
                             length=200, label="x [m]", showvalue=0, variable=view.x_value)
    view.x_textbox = tk.Entry(right_frame_up, width=10)
    view.y_slider = tk.Scale(right_frame_up, from_=-10, to=10, resolution=0.01, orient=tk.HORIZONTAL,
                             length=200, label="y [m]", showvalue=0, variable=view.y_value)
    view.y_textbox = tk.Entry(right_frame_up, width=10)
    view.x_slider.grid(row=row, column=0)
    view.x_textbox.grid(row=row, column=1)
    row += 1
    view.y_slider.grid(row=row, column=0)
    view.y_textbox.grid(row=row, column=1)
    row += 1

    #r-theta sliders
    view.r_slider = tk.Scale(right_frame_up, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL,
                             length=200, label="r [m]", showvalue=0, variable=view.r_value)
    view.r_textbox = tk.Entry(right_frame_up, width=10)
    view.t_slider = tk.Scale(right_frame_up, from_=0, to=2 * math.pi, resolution=0.01, orient=tk.HORIZONTAL,
                             length=200, label="\u03B8 [rad]", showvalue=0, variable=view.t_value)
    view.t_textbox = tk.Entry(right_frame_up, width=10)
    view.r_slider.grid(row=row, column=0)
    view.r_textbox.grid(row=row, column=1)
    row += 1
    view.t_slider.grid(row=row, column=0)
    view.t_textbox.grid(row=row, column=1)
    row += 1

    ###################
    ## Incident wave ##
    ###################
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=0, ipadx=110, pady=10)
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=1, ipadx=35, pady=10)
    row += 1
    title = tk.Label(right_frame_up, text="Incident wave parameters")
    title.grid(row=row, column=0, sticky="nsew")
    row += 1

    # Type of wave
    view.plane_wave = tk.Radiobutton(right_frame_up, text="Plane wave", variable=maths.waveType, val=False)
    view.plane_wave.grid(row=row, column=0, sticky=tk.W)
    row += 1
    view.spherical_wave = tk.Radiobutton(right_frame_up, text="Spherical wave", variable=maths.waveType, val=True)
    view.spherical_wave.grid(row=row, column=0, sticky=tk.W)
    row += 1

    #k value
    view.k_slider = tk.Scale(right_frame_up, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL,
                             length=200, label="k [1/m]", showvalue=0, variable=view.k_value)
    view.k_textbox = tk.Entry(right_frame_up, width=10)
    view.k_slider.grid(row=row, column=0)
    view.k_slider.set(maths.k)
    view.k_textbox.grid(row=row, column=1)
    controller.update_textbox(view.k_textbox, maths.k)
    row += 1

    ##############
    ## Colorbar ##
    ##############
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=0, ipadx=110, pady=10)
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=1, ipadx=35, pady=10)
    row += 1
    title = tk.Label(right_frame_up, text="Colorbar parameters")
    title.grid(row=row, column=0, sticky="nsew")
    row += 1

    # refresh_scale button
    view.refresh_scale = tk.Button(right_frame_up, text="Refresh legend scale", command=view.auto_refresh)
    view.refresh_scale.grid(row=row, column=0, sticky=tk.W)
    row += 1

    # min and max Scale value
    x_max = tk.Label(right_frame_up, text="scale max")
    x_max.grid(row=row, column=0, sticky=tk.W)
    view.scale_max_textbox = tk.Entry(right_frame_up, width=10)
    view.scale_max_textbox.grid(row=row, column=1)
    row += 1

    x_min = tk.Label(right_frame_up, text="scale min")
    x_min.grid(row=row, column=0, sticky=tk.W)
    view.scale_min_textbox = tk.Entry(right_frame_up, width=10)
    view.scale_min_textbox.grid(row=row, column=1)
    row += 1

    # Scale function radio button
    view.logScale = tk.Radiobutton(right_frame_up, text="logarithmic scale", variable=view.scale_type, value="log",
                                   command=controller.logScale)
    view.powScale = tk.Radiobutton(right_frame_up, text="power scale", variable=view.scale_type, value="pow",
                                   command=controller.powScale)
    view.stepScale = tk.Radiobutton(right_frame_up, text="step scale", variable=view.scale_type, value="step",
                                    command=controller.stepScale)
    view.powScale_textbox = tk.Entry(right_frame_up, width=10)
    view.stepScale_textbox = tk.Entry(right_frame_up, width=10)
    view.stepScale_slider = tk.Scale(right_frame_up, from_=0.00009, to=0.005, resolution=0.0001, orient=tk.HORIZONTAL,
                                     length=200, label="step", showvalue=0, variable=stepScale_value)

    view.logScale.grid(row=row, column=0, sticky=tk.W)
    row += 1
    view.powScale.grid(row=row, column=0, sticky=tk.W)
    view.powScale_textbox.grid(row=row, column=1)
    controller.update_textbox(view.powScale_textbox, powArg)
    row += 1
    view.stepScale.grid(row=row, column=0, sticky=tk.W)
    view.stepScale_textbox.grid(row=row, column=1)
    controller.update_textbox(view.stepScale_textbox, stepArg)
    row += 1

    view.stepScale_slider.grid(row=row, column=0, sticky=tk.W)
    row += 1


    #############
    ## Contour ##
    #############
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=0, ipadx=110, pady=10)
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=1, ipadx=35, pady=10)
    row += 1
    title = tk.Label(right_frame_up, text="Contour parameters")
    title.grid(row=row, column=0, sticky="nsew")
    row += 1

    # Radius of the countour
    view.rc_slider = tk.Scale(right_frame_up, from_=0, to=10, resolution=0.01, orient=tk.HORIZONTAL,
                              length=200, label="Contour Radius [m]", showvalue=0, variable=view.rc_value)
    view.rc_textbox = tk.Entry(right_frame_up, width=10)
    view.rc_slider.grid(row=row, column=0)
    view.rc_slider.set(view.rContour)
    view.rc_textbox.grid(row=row, column=1)
    row += 1
    controller.update_textbox(view.rc_textbox, rContour)


    #############################
    ## Other visual parameters ##
    #############################
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=0, ipadx=110, pady=10)
    line = ttk.Separator(right_frame_up, orient=HORIZONTAL)
    line.grid(row=row, column=1, ipadx=35, pady=10)
    row += 1
    title = tk.Label(right_frame_up, text="Other visual parameters")
    title.grid(row=row, column=0, sticky="nsew")
    row += 1

    # button
    view.scattering_amp = tk.Checkbutton(right_frame_up, text="Scattering amplitude visualization", command=controller.change_scattering_amp_visu)
    view.scattering_amp.grid(row=row, column=0, sticky=tk.W)
    row += 1

    view.hide_scatterer = tk.Checkbutton(right_frame_up, text="Hide scatterers", command=controller.hide_scatterer)
    view.hide_scatterer.grid(row=row, column=0, sticky=tk.W)
    row += 1

    # resolution and limits
    title_x_res = tk.Label(right_frame_res_bound, text="x res")
    title_x_min = tk.Label(right_frame_res_bound, text="x min")
    title_x_max = tk.Label(right_frame_res_bound, text="x max")
    view.x_res_textbox = tk.Entry(right_frame_res_bound, width=5)
    view.x_min_textbox = tk.Entry(right_frame_res_bound, width=5)
    view.x_max_textbox = tk.Entry(right_frame_res_bound, width=5)
    title_x_res.grid(row=0, column=0, pady=2)
    view.x_res_textbox.grid(row=0, column=1)
    controller.update_textbox(view.x_res_textbox, view.x_res)
    title_x_min.grid(row=0, column=2, ipadx=10)
    view.x_min_textbox.grid(row=0, column=3)
    controller.update_textbox(view.x_min_textbox, view.x_min)
    title_x_max.grid(row=0, column=4, ipadx=10)
    view.x_max_textbox.grid(row=0, column=5)
    controller.update_textbox(view.x_max_textbox, view.x_max)

    title_y_res = tk.Label(right_frame_res_bound, text="y res")
    title_y_min = tk.Label(right_frame_res_bound, text="y min")
    title_y_max = tk.Label(right_frame_res_bound, text="y max")
    view.y_res_textbox = tk.Entry(right_frame_res_bound, width=5)
    view.y_min_textbox = tk.Entry(right_frame_res_bound, width=5)
    view.y_max_textbox = tk.Entry(right_frame_res_bound, width=5)
    title_y_res.grid(row=1, column=0, pady=2)
    view.y_res_textbox.grid(row=1, column=1)
    controller.update_textbox(view.y_res_textbox, view.y_res)
    title_y_min.grid(row=1, column=2, ipadx=10)
    view.y_min_textbox.grid(row=1, column=3)
    controller.update_textbox(view.y_min_textbox, view.y_min)
    title_y_max.grid(row=1, column=4, ipadx=10)
    view.y_max_textbox.grid(row=1, column=5)
    controller.update_textbox(view.y_max_textbox, view.y_max)

    title_t_res = tk.Label(right_frame_res_bound, text="\u03B8 res")
    view.t_res_textbox = tk.Entry(right_frame_res_bound, width=5)
    title_t_res.grid(row=2, column=0, pady=2)
    view.t_res_textbox.grid(row=2, column=1)
    controller.update_textbox(view.t_res_textbox, view.t_res)

    view.apply_res_bounds = tk.Button(right_frame_res_bound, text="Apply", command=view.apply_res_bounds_method)
    apply_res_bounds.grid(row=3, column=0)


    ####################
    ## Saving options ##
    ####################
    line = ttk.Separator(right_frame_save, orient=HORIZONTAL)
    line.grid(row=row, column=0, ipadx=110, pady=10)
    line = ttk.Separator(right_frame_save, orient=HORIZONTAL)
    line.grid(row=row, column=1, ipadx=35, pady=10)
    row += 1
    title = tk.Label(right_frame_save, text="Saving options")
    title.grid(row=row, column=0, sticky="nsew")
    row += 1

    view.save_button = tk.Button(right_frame_save, text="Save", command=controller.save)
    view.save_button.grid(row=row, column=0, sticky="nsew")
    row += 1


    ###################
    ## End separator ##
    ###################
    line = ttk.Separator(right_frame_save, orient=HORIZONTAL)
    line.grid(row=row, column=0, ipadx=110, pady=10)
    line = ttk.Separator(right_frame_save, orient=HORIZONTAL)
    line.grid(row=row, column=1, ipadx=35, pady=10)
    row += 1


def initializeMesh():
    view.yMesh, view.xMesh = np.meshgrid(np.linspace(y_min, y_max, y_res), np.linspace(x_min, x_max, x_res))


def initializeContour():
    view.thetaContour = np.linspace(0, 2*math.pi, view.t_res)
    view.xContour = rContour * np.cos(view.thetaContour)
    view.yContour = rContour * np.sin(view.thetaContour)


def plotWF():
    k = maths.k
    a = maths.a
    coordinates = maths.coordinates
    if maths.waveType == "spherical":
        z = maths.phi_sph(xMesh, yMesh, k)
    else:
        z = maths.phi_pl(xMesh, k)
    for i in range(maths.N):
        z += a[i] * maths.G(k, np.sqrt((xMesh - coordinates[i][0]) * (xMesh - coordinates[i][0]) + (yMesh - coordinates[i][1]) * (yMesh - coordinates[i][1])))
    z = (np.abs(z)) ** 2
    z /= normalize_z2D(z)
    view.scale_min = z.min()
    view.scale_max = z.max()
    controller.update_textbox(view.scale_min_textbox, round(view.scale_min, 5))
    controller.update_textbox(view.scale_max_textbox, round(view.scale_max, 5))
    view.WF = axMesh.pcolormesh(xMesh, yMesh, z, norm=colors.LogNorm(vmin=view.scale_min, vmax=view.scale_max), cmap='YlOrRd', shading='auto')
    figure.colorbar(WF, ax=axMesh)
    compute_contour()


def updateWF():
    k = maths.k
    a = maths.a
    coordinates = maths.coordinates
    if maths.waveType == "spherical":
        z = maths.phi_sph(xMesh, yMesh, maths.k)
    else:
        z = maths.phi_pl(xMesh, maths.k)
    for i in range(maths.N):
        z += a[i] * maths.G(k, np.sqrt((xMesh - coordinates[i][0]) * (xMesh - coordinates[i][0]) + (yMesh - coordinates[i][1]) * (yMesh - coordinates[i][1])))
    z = (np.abs(z)) ** 2
    z /= normalize_z2D(z)
    if not math.isnan(z[0][0]):
        view.WF.set_array(z)
    else:
        view.WF.set_array(np.zeros((maths.N, maths.N)))
    compute_contour()


def updateWF_change_res_bound():
    k = maths.k
    a = maths.a
    coordinates = maths.coordinates
    if maths.waveType == "spherical":
        z = maths.phi_sph(xMesh, yMesh, maths.k)
    else:
        z = maths.phi_pl(xMesh, maths.k)
    for i in range(maths.N):
        z += a[i] * maths.G(k, np.sqrt((xMesh - coordinates[i][0]) * (xMesh - coordinates[i][0]) + (yMesh - coordinates[i][1]) * (yMesh - coordinates[i][1])))
    z = (np.abs(z)) ** 2
    z /= normalize_z2D(z)
    view.WF = axMesh.pcolormesh(xMesh, yMesh, z, cmap='YlOrRd', shading='auto')
    if scale_type == "log":
        view.WF.set_norm(colors.LogNorm(vmin=scale_min, vmax=scale_max))
    elif scale_type == "pow":
        view.WF.set_norm(colors.PowerNorm(gamma=powArg, vmin=scale_min, vmax=scale_max))
    elif scale_type == "step":
        view.WF.set_norm(colors.BoundaryNorm(boundaries=[scale_min, stepArg, scale_max], ncolors=256))
    compute_contour()


def normalize_z2D(z):
    dS = abs(xMesh[1][1]-xMesh[0][0]) * abs(yMesh[1][1]-yMesh[0][0])
    return sum(sum(z))*dS


def normalize_z1D(z):
    dT = abs(thetaContour[1] - thetaContour[0])
    return sum(z)*dT


def compute_contour():
    global contour_plot
    if contour_plot is not None:
        contour_plot.remove()
        contour_plot = None
    k = maths.k
    a = maths.a
    coordinates = maths.coordinates
    if maths.waveType == "spherical":
        z = maths.phi_sph(xContour, yContour, k)
    else:
        z = maths.phi_pl(xContour, k)
    for i in range(maths.N):
        dx = xContour - coordinates[i][0]
        dy = yContour - coordinates[i][1]
        z += a[i] * maths.G(k, np.sqrt(dx*dx + dy*dy))
    z = (np.abs(z)) ** 2
    z /= normalize_z1D(z)
    if not math.isnan(z[0]):
        axContour.set_ylim(0, z.max()+0.01)
        contour_plot, = axContour.plot(thetaContour, z, color='blue', linewidth=1)


def updateContour():
    view.xContour = rContour * np.cos(view.thetaContour)
    view.yContour = rContour * np.sin(view.thetaContour)


def plotScatterers():
    coordinates = maths.coordinates
    for i in range(maths.N):
        circleI = patches.Circle((coordinates[i][0], coordinates[i][1]), radius=radiusScatterers, color='black')
        patchI = axMesh.add_patch(circleI)
        scatterList.append(circleI)
        patchList.append(patchI)


def updateScatterers(i):
    patchList[i].remove()
    axMesh.add_patch(scatterList[i])


def updateContourCircle():
    view.contour.remove()
    axMesh.add_patch(view.contour)


def current_click():
    for i in range(maths.N):
        if controller.selected_scatter == i:
            scatterList[i].set(color='dodgerblue')
            view.updateScatterers(i)
        else:
            scatterList[i].set(color='black')
            view.updateScatterers(i)


def refresh_scale_method():
    global scale_type, powArg, stepArg, scale_min, scale_max
    if scale_type == "log":
        view.WF.set_norm(colors.LogNorm(vmin=scale_min, vmax=scale_max))
    elif scale_type == "pow":
        view.WF.set_norm(colors.PowerNorm(gamma=powArg, vmin=scale_min, vmax=scale_max))
    elif scale_type == "step":
        view.WF.set_norm(colors.BoundaryNorm(boundaries=[scale_min, stepArg, scale_max], ncolors=256))
    controller.update_textbox(view.scale_min_textbox, round(scale_min, 5))
    controller.update_textbox(view.scale_max_textbox, round(scale_max, 5))
    view.plt.draw()
    view.canvas.draw()


def auto_refresh():
    global scale_min, scale_max
    if WF is not None:
        scale_min, scale_max = WF.get_array().min(), WF.get_array().max()
        refresh_scale_method()


def apply_res_bounds_method():
    view.x_res = int(view.x_res_textbox.get())
    view.x_min = float(view.x_min_textbox.get())
    view.x_max = float(view.x_max_textbox.get())
    view.y_res = int(view.y_res_textbox.get())
    view.y_min = float(view.y_min_textbox.get())
    view.y_max = float(view.y_max_textbox.get())
    view.t_res = int(view.t_res_textbox.get())
    initializeMesh()
    initializeContour()
    view.updateWF_change_res_bound()
    view.update_bounds()
    view.plt.draw()
    view.canvas.draw()


def update_bounds():
    axMesh.set_xlim(x_min, x_max)
    axMesh.set_ylim(y_min, y_max)
    view.radiusScatterers = 0.1 * (x_max - x_min) / 20
    for i in range(maths.N):
        scatterList[i].radius = view.radiusScatterers
        updateScatterers(i)
        updateContourCircle()
