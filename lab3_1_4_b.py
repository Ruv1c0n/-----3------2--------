from matplotlib import pyplot as plt
import numpy as np
from scipy.special import roots_legendre
from prettytable.colortable import ColorTable, Themes

'''Входные данные'''
A = 1
B = 2
C = 1.5
D = 2.5
M = 25
TAU = (D - C) / M
EPS = 0.001


'''Входная функция'''


def f(x, t):
    return np.sin(t / (1 + x ** 2) + 0.001 * x)


'''Квадратурная формула трапеций'''


def Trapezoid(t, N):
    h = (B - A) / N
    result = (A + B) / 2
    for i in range(1, N):
        x_i = A + i * h
        result += f(x_i, t)
    return h * result


'''Метод удвоения'''


def integrateDoubleSteps():
    double_steps_array = np.ones((0, 2), dtype='float64')
    N_arr = np.empty(0, dtype='int')

    # Перебор точек t
    for i in range(M):
        N = 1
        t = C + i * TAU
        integral_prev = Trapezoid(t, N)
        integral_curr = Trapezoid(t, 2 * N)

        # Поиск подходящего числа N
        while abs(integral_curr - integral_prev) > EPS:
            integral_prev = integral_curr
            N *= 2
            integral_curr = Trapezoid(t, 2 * N)

        double_steps_array = np.append(
            double_steps_array,
            [[integral_prev, t]],
            axis=0
        )
        N_arr = np.append(
            N_arr,
            [N],
            axis=0
        )

    return double_steps_array, N_arr


'''Метод с использованием квадратурной формулы Гаусса'''


def quadratureGaussForm():
    gauss_w_3_nodes = np.empty((0, 2), dtype='float64')
    gauss_w_4_nodes = np.empty((0, 2), dtype='float64')

    # Создание узлов
    nodes_3, weights_3 = roots_legendre(3)
    nodes_4, weights_4 = roots_legendre(4)

    # Перевод узлов к нужной форме
    nodes_3_transformed = ((B - A) / 2) * nodes_3 + ((A + B) / 2)
    weights_3_transformed = ((B - A) / 2) * weights_3

    nodes_4_transformed = ((B - A) / 2) * nodes_4 + ((A + B) / 2)
    weights_4_transformed = ((B - A) / 2) * weights_4

    # Подсчет интеграла по узлам
    for i in range(M):
        t = C + i * TAU
        result_3_nodes = sum(weights_3_transformed * f(nodes_3_transformed, t))
        result_4_nodes = sum(weights_4_transformed * f(nodes_4_transformed, t))
        gauss_w_3_nodes = np.append(
            gauss_w_3_nodes,
            [[result_3_nodes, t]],
            axis=0
        )
        gauss_w_4_nodes = np.append(
            gauss_w_4_nodes,
            [[result_4_nodes, t]],
            axis=0
        )

    return gauss_w_3_nodes, gauss_w_4_nodes


# Запуск алгоритмов
double_steps_array, N_iteration = integrateDoubleSteps()
gauss_w_3_nodes, gauss_w_4_nodes = quadratureGaussForm()

# Создание таблицы для вывода
aproximation = ColorTable(theme=Themes.OCEAN)
aproximation.field_names = [
    "Gauss with 3 nodes",
    "Gauss with 4 nodes",
    "Doubling method",
    "N"
]
for i in range(len(double_steps_array)):
    aproximation.add_row([
        gauss_w_3_nodes[i][0],
        gauss_w_4_nodes[i][0],
        double_steps_array[i][0],
        N_iteration[i]
    ])
aproximation.align = 'l'
aproximation.align["N"] = 'c'
print(aproximation)

# Вывод графика
fig, axs = plt.subplots(1, 4)
fig.canvas.manager.set_window_title('Aproximation')
fig.set_figwidth(22)
fig.set_figheight(7)

axs[0].set_xlabel("Точка t")
axs[0].set_ylabel("Значение функции F")
axs[0].errorbar(
    gauss_w_3_nodes[:, 1],
    gauss_w_3_nodes[:, 0],
    label='Gauss with 3 nodes'
)
axs[0].errorbar(
    gauss_w_4_nodes[:, 1],
    gauss_w_4_nodes[:, 0],
    label='Gauss with 4 nodes'
)
axs[0].errorbar(
    double_steps_array[:, 1],
    double_steps_array[:, 0],
    label='Doubling method'
)
axs[0].legend(loc='lower right')

axs[1].axis('off')
axs[3].axis('off')

aproximation_table_data = []
for i in range(len(gauss_w_3_nodes)):
    aproximation_table_data.append(
        [gauss_w_3_nodes[i][0],
         gauss_w_4_nodes[i][0],
         double_steps_array[i][0],
         N_iteration[i]]
    )
aproximation_table = axs[2].table(
    cellText=aproximation_table_data,
    colLabels=[
        "Gauss with 3 nodes",
        "Gauss with 4 nodes",
        "Doubling method",
        "N"
    ],
    rowLoc='center',
    colLoc='center',
    colColours=["palegreen"] * 4,
    loc='center'
)
aproximation_table.set_fontsize(30)
aproximation_table.scale(3.5, 1.5)
axs[2].axis('off')

plt.show()


# fig, ax = plt.subplots()
# fig.canvas.manager.set_window_title('Aproximation')
# ax.set_xlabel("Точка t")
# ax.set_ylabel("Значение функции F")
# ax.errorbar(
#     gauss_w_3_nodes[:, 1],
#     gauss_w_3_nodes[:, 0],
#     label='Gauss with 3 nodes'
# )
# ax.errorbar(
#     gauss_w_4_nodes[:, 1],
#     gauss_w_4_nodes[:, 0],
#     label='Gauss with 4 nodes'
# )
# ax.errorbar(
#     double_steps_array[:, 1],
#     double_steps_array[:, 0],
#     label='Doubling method'
# )
# ax.legend(loc='lower right')

# fig3 = plt.figure(figsize=(10, 10))
# ax3 = fig3.add_subplot()
# fig3.canvas.manager.set_window_title('Aproximation')
# aproximation_table_data = []
# for i in range(len(gauss_w_3_nodes)):
#     aproximation_table_data.append(
#         [gauss_w_3_nodes[i][0],
#          gauss_w_4_nodes[i][0],
#          double_steps_array[i][0],
#          N_iteration[i]]
#     )
# aproximation_table = ax3.table(
#     cellText=aproximation_table_data,
#     colLabels=[
#         "Gauss with 3 nodes",
#         "Gauss with 4 nodes",
#         "Doubling method",
#         "N"
#     ],
#     rowLoc='center',
#     colLoc='center',
#     colColours=["palegreen"] * 4,
#     loc='center'
# )
# aproximation_table.set_fontsize(50)
# aproximation_table.set_fontsize(30)
# aproximation_table.scale(1, 2)
# ax3.axis('off')

# plt.show()
