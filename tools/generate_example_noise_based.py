import matplotlib.pyplot as plt
import numpy as np


def generate(sample_size):
    x = np.random.uniform(0, 1, sample_size)
    e = np.random.uniform(-0.1, 0.1, sample_size)

    y = 2*x + e
    return x, y, e


if __name__ == "__main__":
    n = 1000

    x, y, e = generate(n)

    a = (np.sum((x- np.mean(x))) * np.sum((y- np.mean(y)))) / np.sum((x- np.mean(x)) ** 2)
    b = np.mean(y) - a * np.mean(x)
    yhat = a*x + b

    plt.plot(x, y)
    plt.plot(x, yhat)
    plt.show()
