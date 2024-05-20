import numpy as np
import matplotlib.pyplot as plt
import sympy
import re
import string


def f_x(x, y, t=0):
    return x*y - 4

def f_y(x, y, t=0):
    return (x-4)*(y-x)

def Runge_kutt4(x_0, y_0, a, b, h):
    points_x, points_y = [], []
    step = a
    while step < b:
        k1 = h * f_x(x_0, y_0)
        q1 = h * f_y(x_0, y_0)
        k2 = h * f_x(x_0 + k1 / 3, y_0 + q1 / 3)
        q2 = h * f_y(x_0 + k1 / 3, y_0 + q1 / 2)
        k3 = h * f_x(x_0 - k1 / 2 + k2, y_0 - q1 / 2 + q2)
        q3 = h * f_y(x_0 - k1 / 2 + k2, y_0 - q1 / 2 + q2)
        k4 = h * f_x(x_0 + k1 - k2 + k3, y_0 + q1 - q2 + q3)
        q4 = h * f_y(x_0 + k1 - k2 + k3, y_0 + q1 - q2 + q3)
        x_0 += (k1 + 3 * k2 + 3 * k3 + k4) / 8
        y_0 += (q1 + 3 * q2 + 3 * q3 + q4) / 8
        points_x.append(x_0)
        points_y.append(y_0)
        step += h
    plt.plot(points_x, points_y, color='blue')

def solveWequilibrium():
    xy_range = [-5, 7]
    h = 0.05
    a, b = 0, 200
    shift = 1
    equilibrium_points = np.array([[-2, -2], [2, 2], [4, 1]])
    fig = plt.figure(1, (10, 10))
    plt.xlim(*xy_range)
    plt.ylim(*xy_range)
    X, Y = np.meshgrid(
        np.linspace(*xy_range, 100),
        np.linspace(*xy_range, 100)
    )
    plt.streamplot(
        X, Y,
        f_x(X, Y), f_y(X, Y),
        color='cyan',
        density=1,
        arrowstyle='->',
        arrowsize=1.5
    )

    for i in range(len(equilibrium_points)):
        plt.pause(1)
        Runge_kutt4(
            equilibrium_points[i][0], 
            equilibrium_points[i][1], 
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0],
            equilibrium_points[i][1] - shift, 
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0],
            equilibrium_points[i][1] + shift, 
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] - shift, 
            equilibrium_points[i][1], a, b, h)
        Runge_kutt4(
            equilibrium_points[i][0] + shift, 
            equilibrium_points[i][1], 
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] - shift, 
            equilibrium_points[i][1] - shift, 
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] - shift, 
            equilibrium_points[i][1] + shift, 
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] + shift, 
            equilibrium_points[i][1] - shift, 
            a, b, h
            )
        Runge_kutt4(
            equilibrium_points[i][0] + shift, 
            equilibrium_points[i][1] + shift, 
            a, b, h
        )
        shift -= 0.4

    plt.plot(
        *np.hsplit(equilibrium_points, 2),
        c='red',
        linestyle='None',
        marker='o'
    )
    plt.show()

def solveWOUTequilibrium(dotx, doty):
    xy_range = [-5, 7]
    h = 0.05
    a, b = 0, 200
    shift = 1
    points = np.array(sympy.solve([dotx, doty], ['x', 'y'], check=False))
    equilibrium_points = np.array([[*p] for p in points])
    fig = plt.figure(1, (10, 10))
    plt.xlim(*xy_range)
    plt.ylim(*xy_range)
    X, Y = np.meshgrid(
        np.linspace(*xy_range, 100),
        np.linspace(*xy_range, 100)
    )
    plt.streamplot(
        X, Y,
        f_x(X, Y), f_y(X, Y),
        color='cyan',
        density=1,
        arrowstyle='->',
        arrowsize=1.5
    )

    for i in range(len(equilibrium_points)):
        Runge_kutt4(
            equilibrium_points[i][0],
            equilibrium_points[i][1],
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0],
            equilibrium_points[i][1] - shift,
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0],
            equilibrium_points[i][1] + shift,
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] - shift,
            equilibrium_points[i][1], a, b, h)
        Runge_kutt4(
            equilibrium_points[i][0] + shift,
            equilibrium_points[i][1],
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] - shift,
            equilibrium_points[i][1] - shift,
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] - shift,
            equilibrium_points[i][1] + shift,
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] + shift,
            equilibrium_points[i][1] - shift,
            a, b, h
        )
        Runge_kutt4(
            equilibrium_points[i][0] + shift,
            equilibrium_points[i][1] + shift,
            a, b, h
        )
        shift -= 0.4

    plt.plot(
        *np.hsplit(equilibrium_points, 2),
        c='red',
        linestyle='None',
        marker='o'
    )
    plt.show()


solveWequilibrium()
# solveWOUTequilibrium('x*y - 4', '(x-4)*(y-x)')
