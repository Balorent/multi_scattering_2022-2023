#                              CONTROLLER.PY
# ------------------------------------------------------------------------
# Author       :    Baptiste Lorent
# Last edition :    08 october 2022
# ------------------------------------------------------------------------


# Imports ----------------------------------------------------------------

# External libraries
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Other files
import view
import maths


# Start ------------------------------------------------------------------
if __name__ == "__main__":
    print("main.py should be started instead")
    exit()


# Variables --------------------------------------------------------------
N = maths.N
coordinates = maths.coordinates
in_range = -1
selected_scatter = -1

scattering_amp_visu_state = 0
hide_scatterer_bool = 0


# Functions --------------------------------------------------------------
def initialise():
    view.initialise()


def button_press_callback(event):
    """ Detects clicks and checks if in range of an editable plot """
    global in_range, selected_scatter
    test = 0
    if event.inaxes:
        for i in range(N):
            x_i, y_i = event.xdata, event.ydata
            if (event.xdata-coordinates[i][0])**2 + (event.ydata-coordinates[i][1])**2 <= view.scatterList[i].radius**2:
                r_i, theta_i = np.sqrt(x_i ** 2 + y_i ** 2), np.arctan2(y_i, x_i) + (np.arctan2(y_i, x_i) < 0)*2*math.pi
                test = 1
                in_range = i
                selected_scatter = i
                view.current_click()
                set_rightpanel(x_i, y_i, r_i, theta_i)
    if not test:
        selected_scatter = -1
        set_rightpanel(0, 0, 0, 0)
        view.current_click()
    view.plt.draw()
    view.canvas.draw()


def button_release_callback(event):
    """ Detects press release and removes any interaction with the plot """
    global in_range
    if selected_scatter == -1:
        set_rightpanel(0, 0, 0, 0)
        view.current_click()
    in_range = -1
    view.plt.draw()
    view.canvas.draw()


def motion_notify_callback(event):
    """ Detects movement while clicking and updates the corresponding values """
    global in_range

    if event.button and event.inaxes and in_range >= 0:
        x_i, y_i = event.xdata, event.ydata
        r_i, theta_i = np.sqrt(x_i ** 2 + y_i ** 2), np.arctan2(y_i, x_i) + (np.arctan2(y_i, x_i) < 0) * 2 * math.pi
        update_scatterer_move(x_i, y_i, in_range)
        set_rightpanel(x_i, y_i, r_i, theta_i)
        view.plt.draw()
        view.canvas.draw()

def update_scatterer_move(x, y, i):
    view.scatterList[i].center = x, y
    coordinates[i][0], coordinates[i][1] = x, y
    maths.compute_a()
    view.updateWF()
    view.updateScatterers(i)
    update_scattering_amp_visu()


def update_k_change(newK):
    maths.k = newK
    maths.compute_a()
    view.updateWF()
    update_scattering_amp_visu()


def update_phi_change():
    maths.compute_a()
    view.updateWF()
    update_scattering_amp_visu()


def set_rightpanel(x, y, r, theta):
    view.x_value.set(x)
    update_textbox(view.x_textbox, round(x, 3))
    view.y_value.set(y)
    update_textbox(view.y_textbox, round(y, 3))
    view.r_value.set(r)
    update_textbox(view.r_textbox, round(r, 3))
    view.t_value.set(theta)
    update_textbox(view.t_textbox, round(theta, 3))


def update_textbox(textbox, value):
    textbox.delete(0, view.tk.END)
    textbox.insert(0, value)


def update_X_from_tb(event):
    if selected_scatter >= 0:
        # new position variables
        new_x, y = float(view.x_textbox.get()), maths.coordinates[selected_scatter][1]
        new_r, new_theta = np.sqrt(y ** 2 + new_x ** 2), np.arctan2(y, new_x) + (np.arctan2(y, new_x) < 0) * 2 * math.pi

        # update the x slider
        view.x_value.set(new_x)
        # update the r slider and textbox
        view.r_value.set(new_r)
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the theta slider and textbox
        view.t_value.set(new_theta)
        update_textbox(view.t_textbox, round(new_theta, 3))

        # update the plot
        update_scatterer_move(new_x, y, selected_scatter)

        # refresh the view
        view.plt.draw()
        view.canvas.draw()


def update_Y_from_tb(event):
    if selected_scatter >= 0:
        # new position variables
        x, new_y = maths.coordinates[selected_scatter][0], float(view.y_textbox.get())
        new_r, new_theta = np.sqrt(new_y ** 2 + x ** 2), np.arctan2(new_y, x) + (np.arctan2(new_y, x) < 0) * 2 * math.pi

        # update the y slider
        view.y_value.set(new_y)
        # update the r slider and textbox
        view.r_value.set(new_r)
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the theta slider and textbox
        view.t_value.set(new_theta)
        update_textbox(view.t_textbox, round(new_theta, 3))

        # update the plot
        update_scatterer_move(x, new_y, selected_scatter)

        # refresh the view
        view.plt.draw()
        view.canvas.draw()


def update_R_from_tb(event):
    if selected_scatter >= 0:
        # new position variables
        new_r, theta = float(view.r_textbox.get()), float(view.t_textbox.get())
        new_x, new_y = new_r * np.cos(theta), new_r * np.sin(theta)

        # update the r slider
        view.r_value.set(new_r)
        # update the x slider and textbox
        view.x_value.set(new_x)
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the y slider and textbox
        view.y_value.set(new_y)
        update_textbox(view.y_textbox, round(new_y, 3))

        # update the plot
        update_scatterer_move(new_x, new_y, selected_scatter)

        # refresh the view
        view.plt.draw()
        view.canvas.draw()


def update_T_from_tb(event):
    if selected_scatter >= 0:
        # new position variables
        r, new_theta = float(view.r_textbox.get()), float(view.t_textbox.get())
        new_x, new_y = r * np.cos(new_theta), r * np.sin(new_theta)

        # update the theta slider
        view.t_value.set(new_theta)
        # update the x slider and textbox
        view.x_value.set(new_x)
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the y slider and textbox
        view.y_value.set(new_y)
        update_textbox(view.y_textbox, round(new_y, 3))

        # update the plot
        update_scatterer_move(new_x, new_y, selected_scatter)

        # refresh the view
        view.plt.draw()
        view.canvas.draw()


def update_X_from_slider(value):
    if selected_scatter >= 0:
        # new position variables
        new_x, y = float(value), maths.coordinates[selected_scatter][1]
        new_r, new_theta = np.sqrt(y ** 2 + new_x ** 2), np.arctan2(y, new_x) + (np.arctan2(y, new_x) < 0) * 2 * math.pi

        # update the x textbox
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the r slider and textbox
        view.r_value.set(new_r)
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the theta slider and textbox
        view.t_value.set(new_theta)
        update_textbox(view.t_textbox, round(new_theta, 3))

        # update the plot
        update_scatterer_move(new_x, y, selected_scatter)

        # refresh the view
        view.plt.draw()
        view.canvas.draw()


def update_Y_from_slider(value):
    if selected_scatter >= 0:
        # new position variables
        x, new_y = maths.coordinates[selected_scatter][0], float(value)
        new_r, new_theta = np.sqrt(new_y ** 2 + x ** 2), np.arctan2(new_y, x) + (np.arctan2(new_y, x) < 0) * 2 * math.pi

        # update the y textbox
        update_textbox(view.y_textbox, round(new_y, 3))
        # update the r slider and textbox
        view.r_value.set(new_r)
        update_textbox(view.r_textbox, round(new_r, 3))
        # update the theta slider and textbox
        view.t_value.set(new_theta)
        update_textbox(view.t_textbox, round(new_theta, 3))

        # update the plot
        update_scatterer_move(x, new_y, selected_scatter)

        # refresh the view
        view.plt.draw()
        view.canvas.draw()


def update_R_from_slider(value):
    if selected_scatter >= 0:
        # new position variables
        new_r, theta = float(value), float(view.t_textbox.get())
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
        update_scatterer_move(new_x, new_y, selected_scatter)

        # refresh the view
        view.plt.draw()
        view.canvas.draw()


def update_T_from_slider(value):
    if selected_scatter >= 0:
        # new position variables
        r, new_theta = float(view.r_textbox.get()), float(value)
        new_x, new_y = r * np.cos(new_theta), r * np.sin(new_theta)

        # update the theta textbox
        update_textbox(view.t_textbox, round(new_theta, 3))
        # update the x slider and textbox
        view.x_value.set(new_x)
        update_textbox(view.x_textbox, round(new_x, 3))
        # update the y slider and textbox
        view.y_value.set(new_y)
        update_textbox(view.y_textbox, round(new_y, 3))

        # update the plot
        update_scatterer_move(new_x, new_y, selected_scatter)

        # refresh the view
        view.plt.draw()
        view.canvas.draw()


def update_K_from_tb(event):
    new_k = float(view.k_textbox.get())
    view.k_value.set(new_k)
    update_k_change(new_k)
    view.plt.draw()
    view.canvas.draw()


def update_K_from_slider(value):
    new_k = float(value)
    update_k_change(new_k)
    update_textbox(view.k_textbox, new_k)
    view.plt.draw()
    view.canvas.draw()


def update_RC_from_tb(event):
    new_rc = float(view.rc_textbox.get())
    view.rContour = new_rc
    view.contour.radius = new_rc
    view.rc_value.set(new_rc)
    view.updateContour()
    view.compute_contour()
    view.plt.draw()
    view.canvas.draw()


def update_RC_from_slider(value):
    new_rc = float(value)
    view.rContour = new_rc
    view.contour.radius = new_rc
    update_textbox(view.rc_textbox, new_rc)
    view.updateContour()
    view.compute_contour()
    view.plt.draw()
    view.canvas.draw()


def waveSph():
    maths.waveType = "spherical"
    update_phi_change()
    view.plt.draw()
    view.canvas.draw()


def wavePl():
    maths.waveType = "plane"
    update_phi_change()
    view.plt.draw()
    view.canvas.draw()


def logScale():
    view.scale_type = "log"
    view.refresh_scale_method()


def powScale():
    view.scale_type = "pow"
    view.refresh_scale_method()


def stepScale():
    view.scale_type = "step"
    view.refresh_scale_method()


def update_powArg(event):
    newArg = float(view.powScale_textbox.get())
    view.powArg = newArg
    if view.scale_type == "pow":
        view.refresh_scale_method()


def update_stepArg_from_tb(event):
    newArg = float(view.stepScale_textbox.get())
    view.stepArg = newArg
    view.stepScale_value.set(newArg)
    if view.scale_type == "step":
        view.refresh_scale_method()
    view.plt.draw()
    view.canvas.draw()


def update_stepArg_from_slider(value):
    newArg = float(value)
    view.stepArg = newArg
    update_textbox(view.stepScale_textbox, newArg)
    if view.scale_type == "step":
        view.refresh_scale_method()
    view.plt.draw()
    view.canvas.draw()


def update_scale_min_from_tb(event):
    new_scale_min = float(view.scale_min_textbox.get())
    view.scale_min = new_scale_min
    view.refresh_scale_method()

def update_scale_max_from_tb(event):
    new_scale_max = float(view.scale_max_textbox.get())
    view.scale_max = new_scale_max
    view.refresh_scale_method()


def change_scattering_amp_visu():
    global scattering_amp_visu_state
    a = maths.a
    scattering_amp_visu_state = not scattering_amp_visu_state
    for i in range(maths.N):
        if scattering_amp_visu_state:
            view.scatterList[i].radius = abs(a[i])/max(abs(a)) * (view.x_max - view.x_min)/20
        else:
            view.scatterList[i].radius = view.radiusScatterers
    view.plt.draw()
    view.canvas.draw()


def hide_scatterer():
    global hide_scatterer_bool
    hide_scatterer_bool = not hide_scatterer_bool
    for i in range(maths.N):
        if hide_scatterer_bool:
            view.scatterList[i].remove()
        else:
            view.axMesh.add_patch(view.scatterList[i])
        view.plt.draw()
        view.canvas.draw()


def update_scattering_amp_visu():
    global scattering_amp_visu_state
    a = maths.a
    if scattering_amp_visu_state:
        for i in range(maths.N):
            view.scatterList[i].radius = abs(a[i]) / max(abs(a)) * (view.x_max - view.x_min)/20
            view.plt.draw()
            view.canvas.draw()


def save():
    yMesh, xMesh = np.meshgrid(np.linspace(view.y_min, view.y_max, 500), np.linspace(view.x_min, view.x_max, 500))
    fig, ax = plt.subplots()

    k = maths.k
    a = maths.a
    coordinates = maths.coordinates
    if maths.waveType == "spherical":
        z = maths.phi_sph(xMesh, yMesh, k)
    else:
        z = maths.phi_pl(xMesh, k)
    for i in range(maths.N):
        z += a[i] * maths.G(k, np.sqrt(
            (xMesh - coordinates[i][0]) * (xMesh - coordinates[i][0]) + (yMesh - coordinates[i][1]) * (
                        yMesh - coordinates[i][1])))
    z = (np.abs(z)) ** 2
    z /= view.normalize_z2D(z)
    view.scale_min = z.min()
    view.scale_max = z.max()
    WF = ax.pcolormesh(xMesh, yMesh, z, norm=colors.LogNorm(vmin=view.scale_min, vmax=view.scale_max),
                                cmap='YlOrRd', shading='auto')
    fig.colorbar(WF, ax=ax)
    plt.axis('scaled')
    plt.savefig('xy_plot.png')
    plt.close(fig)

    fig2, ax2 = plt.subplots()
    thetaContour = np.linspace(0, 2*math.pi, 10000)
    xContour = view.rContour * np.cos(thetaContour)
    yContour = view.rContour * np.sin(thetaContour)
    coordinates = maths.coordinates
    if maths.waveType == "spherical":
        z = maths.phi_sph(xContour, yContour, k)
    else:
        z = maths.phi_pl(xContour, k)
    for i in range(maths.N):
        dx = xContour - coordinates[i][0]
        dy = yContour - coordinates[i][1]
        z += a[i] * maths.G(k, np.sqrt(dx * dx + dy * dy))
    z = (np.abs(z)) ** 2
    z /= view.normalize_z1D(z)
    if not math.isnan(z[0]):
        ax2.set_ylim(0, z.max() + 0.01)
        ax2.plot(thetaContour, z, color='blue', linewidth=1)
    ax2.set_aspect(30)
    plt.savefig('theta_plot.png')

    print("Image Saved !")