

# Определение функций для системы уравнений


from matplotlib import pyplot as plt
import numpy as np


def dx_dt(x, y):
    return x*y-4


def dy_dt(x, y):
    return (x-4)*(y-x)


# Создание сетки значений x и y
x = np.linspace(-7, 7, 100)
y = np.linspace(-7, 7, 100)
X, Y = np.meshgrid(x, y)

# Вычисление значений производных dx/dt и dy/dt на сетке
dx = dx_dt(X, Y)
dy = dy_dt(X, Y)

# Построение фазового портрета
plt.figure(figsize=(8, 8))
plt.streamplot(X, Y, dx, dy, color='b', linewidth=1, density=1.5, arrowstyle='->',
               arrowsize=1.5)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Линии фазового портрета')
plt.grid(True)
plt.axvline(0, color='black')
plt.axhline(0, color='black')
plt.scatter(-2, -2, s=50, c='red')
plt.scatter(2, 2, s=50, c='red')
plt.scatter(4, 1, s=50, c='red')
plt.show()
