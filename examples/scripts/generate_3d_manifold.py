import numpy as np
import csv
import matplotlib.pyplot as plt

def torus(n=1000, R=3.0, r=1.0, noise=0.02):
    u = np.random.uniform(0, 2*np.pi, n)
    v = np.random.uniform(0, 2*np.pi, n)

    x = (R + r*np.cos(v)) * np.cos(u)
    y = (R + r*np.cos(v)) * np.sin(u)
    z = r*np.sin(v)

    pts = np.column_stack((x, y, z))

    if noise:
        pts += noise * np.random.randn(n, 3)

    return pts


def double_torus(n=2000,
                 R=2.0,
                 r=0.7,
                 separation=3.0,
                 noise=0.02):

    n1 = n // 2
    n2 = n - n1

    T1 = torus(n1, R, r, noise=0)
    T2 = torus(n2, R, r, noise=0)

    T1[:, 0] -= separation
    T2[:, 0] += separation

    bridge_n = max(100, n // 10)

    t = np.random.uniform(-1, 1, bridge_n)

    bridge = np.column_stack([
        separation * t,
        0.3 * np.random.randn(bridge_n),
        0.3 * np.random.randn(bridge_n)
    ])

    pts = np.vstack([T1, T2, bridge])

    if noise:
        pts += noise * np.random.randn(*pts.shape)

    return pts


def mobius(n=2000, width=0.4, noise=0.01):
    u = np.random.uniform(0, 2*np.pi, n)
    v = np.random.uniform(-width, width, n)

    x = (1 + v*np.cos(u/2)) * np.cos(u)
    y = (1 + v*np.cos(u/2)) * np.sin(u)
    z = v*np.sin(u/2)

    pts = np.column_stack((x, y, z))

    if noise:
        pts += noise * np.random.randn(n, 3)

    return pts


def save_csv(points, filename="data.csv"):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(points)

save_csv(mobius(n=500, noise=0), "../data/mobius.csv")
save_csv(double_torus(n=500, noise=0.2), "../data/double_torus.csv")
save_csv(torus(n=500, noise=0.2), "../data/torus.csv")
