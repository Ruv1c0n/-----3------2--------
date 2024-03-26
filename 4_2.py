from scipy.optimize import root
import numpy as np


def Function(x, y):
    return np.array([
        y - np.sin(x) - 1,
        y - x**2 + 1
    ])


def Phi(F):
    return np.sum(np.square(F))


def Jacobian(x, y):
    return np.array([[-np.cos(x), -2 * x], [1, 1]])


def gradient_descent(X, eps=1e-16, max_iter=1000):
    delta_k = X
    for _ in range(max_iter):

        gradient = np.dot(Jacobian(*delta_k), Function(*delta_k))
        '''Вариант А'''
        # n = Phi(Function(*delta_k))
        # d = np.dot(gradient, gradient)
        '''Вариант В'''
        n = np.dot(Function(*delta_k),
                   np.dot(np.transpose(Jacobian(*delta_k)), gradient))
        d = np.dot(np.dot(np.transpose(Jacobian(*delta_k)), gradient),
                   np.dot(np.transpose(Jacobian(*delta_k)), gradient))
        lambda_k = n / d

        delta_k1 = delta_k - lambda_k * gradient

        if np.linalg.norm(delta_k1 - delta_k) / np.linalg.norm(delta_k) < eps:
            return delta_k1

        delta_k = delta_k1
    return delta_k


X = np.array([-1, -1], dtype=float)
print(gradient_descent(X))
X = np.array([1, 1], dtype=float)
print(gradient_descent(X))
