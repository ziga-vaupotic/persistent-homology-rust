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



def ring(n=300, r_inner=2.0, r_outer=4.0, noise=0.05, uniform_area=True):
    theta = np.linspace(0, 2*np.pi, n, endpoint=False)

    if uniform_area:
        # uniform in area: r^2 is uniform
        u = np.random.rand(n)
        r = np.sqrt(u * (r_outer**2 - r_inner**2) + r_inner**2)
    else:
        # uniform in radius
        r = np.random.uniform(r_inner, r_outer, size=n)

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    # add noise
    x += noise * np.random.randn(n)
    y += noise * np.random.randn(n)

    return np.column_stack((x, y))

def save_csv(points, filename="data.csv"):
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(points)

data = spiral(n=100, noise=0.3)

#data = ring(n=500, r_inner=2.0, r_outer=4.0, noise=0.3)
save_csv(data, "../../data.csv")