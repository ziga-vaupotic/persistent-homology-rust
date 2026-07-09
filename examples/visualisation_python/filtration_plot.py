import csv
import matplotlib.pyplot as plt


def load_points(path):
    points = []

    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            x = float(row[0])
            y = float(row[1])
            points.append((x, y))

    return points


def load_filtration(path):
    simplices = []

    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            eps = float(row[0])
            vertices = list(map(int, row[1:]))
            simplices.append((eps, vertices))

    return simplices


def at_scale(simplices, eps):
    return [s for s in simplices if s[0] <= eps]


def split_simplices(simplices):
    points = set()
    edges = []
    triangles = []

    for _, v in simplices:
        if len(v) == 1:
            points.add(v[0])
        elif len(v) == 2:
            edges.append(v)
        elif len(v) == 3:
            triangles.append(v)

    return sorted(points), edges, triangles


def plot(points, simplices_pts, edges, triangles):

    fig, ax = plt.subplots()

    for i in simplices_pts:
        x, y = points[i]
        ax.scatter(x, y, color='tab:green')


    for i, j in edges:
        x = [points[i][0], points[j][0]]
        y = [points[i][1], points[j][1]]
        ax.plot(x, y, color='tab:blue')

    for i, j, k in triangles:
        xs = [points[i][0], points[j][0], points[k][0]]
        ys = [points[i][1], points[j][1], points[k][1]]
        ax.fill(xs, ys, alpha=0.2, color='grey')

    ax.axis('off')
    ax.set_aspect('equal')
    #plt.title("Vietoris–Rips filtration at fixed ε")
    #plt.show()

    plt.show()


if __name__ == "__main__":
    points = load_points("../../data.csv")
    simplices = load_filtration("../../filtration.csv")

    eps = 0.3

    filtered = at_scale(simplices, eps)
    p_ids, edges, triangles = split_simplices(filtered)

    plot(points, p_ids, edges, triangles)
