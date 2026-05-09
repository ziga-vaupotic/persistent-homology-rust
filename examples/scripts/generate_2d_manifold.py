import numpy as np
import csv

def circle(n=200, noise=0.05, radius=1.0):
    t = np.linspace(0, 2*np.pi, n)
    x = radius * np.cos(t) + noise * np.random.randn(n)
    y = radius * np.sin(t) + noise * np.random.randn(n)
    return np.column_stack((x, y))

def spiral(n=300, noise=0.05):
    t = np.linspace(0, 4*np.pi, n)
    r = t
    x = r * np.cos(t) + noise * np.random.randn(n)
    y = r * np.sin(t) + noise * np.random.randn(n)
    return np.column_stack((x, y))


def save_csv(points, filename="data.csv"):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(points)

data = spiral(n=100, noise=0.3)
save_csv(data, "../../data.csv")