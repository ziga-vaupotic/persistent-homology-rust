import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def load_points(path):
    points = []

    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            x = float(row[0])
            y = float(row[1])
            z = float(row[2])
            points.append((x, y, z))

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
    tetrahedra = []

    for _, v in simplices:
        if len(v) == 1:
            points.add(v[0])
        elif len(v) == 2:
            edges.append(v)
        elif len(v) == 3:
            triangles.append(v)
        elif len(v) == 4:
            tetrahedra.append(v)



    return sorted(points), edges, triangles, tetrahedra


def plot(points, simplices_pts, edges, triangles, tetrahedra):

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    print(len(simplices_pts), len(points))

    for i in simplices_pts:
        x, y, z = points[i]
        ax.scatter(x, y, z, color='tab:green', s=1)

    for i, j in edges:
        x = [points[i][0], points[j][0]]
        y = [points[i][1], points[j][1]]
        z = [points[i][2], points[j][2]]

        ax.plot(x, y, z, color='tab:blue')

    for i, j, k in triangles:
        verts = [[
            points[i],
            points[j],
            points[k]
        ]]

        poly = Poly3DCollection(
            verts,
            alpha=0.2,
            facecolor='red'
        )
        ax.add_collection3d(poly)


    for i, j, k, l in tetrahedra:
        faces = [
            [points[i], points[j], points[k]],
            [points[i], points[j], points[l]],
            [points[i], points[k], points[l]],
            [points[j], points[k], points[l]],
        ]

        poly = Poly3DCollection(
            faces,
            alpha=0.1,
            facecolor='pink',
        )
        ax.add_collection3d(poly)

    ax.axis('off')
    ax.set_aspect('equal')
    #plt.title("Vietoris–Rips filtration at fixed ε")
    #plt.show()

    plt.show()
    plt.savefig("filtration_plot_3d.png", dpi=500)


if __name__ == "__main__":
    points = load_points("../../data.csv")
    simplices = load_filtration("../../filtration.csv")

    eps = 0.5

    filtered = at_scale(simplices, eps)
    p_ids, edges, triangles, tetrahedra = split_simplices(filtered)

    plot(points, p_ids, edges, triangles, tetrahedra)