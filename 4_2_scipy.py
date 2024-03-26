from scipy.optimize import root
import numpy as np


def equations(xy):
    x, y = xy
    f = y - np.sin(x) - 1
    g = y - x**2 + 1
    return [f, g]


def jacobian(xy):
    x, y = xy
    df_dx = -np.cos(x)
    df_dy = 1
    dg_dx = -2 * x
    dg_dy = 1
    return [[df_dx, df_dy], [dg_dx, dg_dy]]


initial_guess = [1, 1]

solution = root(equations, initial_guess, jac=jacobian, method='hybr')

if solution.success:
    x, y = solution.x
    print(f"Решение: x = {x}, y = {y}")
else:
    print("Решение не найдено. Метод Ньютона не сошелся.")


initial_guess = [-1, -1]

solution = root(equations, initial_guess, jac=jacobian, method='hybr')

if solution.success:
    x, y = solution.x
    print(f"Решение: x = {x}, y = {y}")
else:
    print("Решение не найдено. Метод Ньютона не сошелся.")
