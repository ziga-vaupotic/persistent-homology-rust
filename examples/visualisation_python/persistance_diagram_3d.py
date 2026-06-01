import csv
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def load_filtration(path):
    simplices = []
    with open(path, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            eps = float(row[0])
            vertices = list(map(int, row[1:]))
            simplices.append((eps, vertices))
    return simplices


def load_persistence(path, filtration):
    births = []
    deaths = []
    dimensions = []
    with open(path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row) >= 2:
                dim = int(row[0]) if len(row) > 2 else 0
                birth_idx = int(row[1])
                birth = filtration[birth_idx][0] if birth_idx < len(filtration) else None

                if len(row) > 2 and row[2] != '':
                    death_idx = int(row[2])
                    death = filtration[death_idx][0] if death_idx < len(filtration) else None
                else:
                    death = None

                if death is not None and abs(birth - death) < 1e-10:
                    continue

                births.append(birth)
                deaths.append(death)
                dimensions.append(dim)
    return births, deaths, dimensions


def load_points(path):
    points = []
    with open(path, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 3:
                points.append((float(row[0]), float(row[1]), float(row[2])))
    return points


def get_filtration_at_eps(filtration, eps):
    filtered = [s for s in filtration if s[0] <= eps]

    points = set()
    edges = []
    triangles = []
    tetrahedra = []

    for _, v in filtered:
        if len(v) == 1:
            points.add(v[0])
        elif len(v) == 2:
            edges.append(v)
        elif len(v) == 3:
            triangles.append(v)
        elif len(v) == 4:
            tetrahedra.append(v)

    return sorted(points), edges, triangles, tetrahedra


def plot_complex(ax, points, point_ids, edges, triangles, tetrahedra, max_eps):
    ax.clear()
    ax.set_title(f"3D Rips Complex at ε = {max_eps:.3f}")
    ax.axis('off')

    if point_ids:
        xs, ys, zs = zip(*[points[i] for i in point_ids])
        ax.scatter(xs, ys, zs, color='tab:orange', edgecolors='black', s=30)

    for i, j in edges:
        x = [points[i][0], points[j][0]]
        y = [points[i][1], points[j][1]]
        z = [points[i][2], points[j][2]]
        ax.plot(x, y, z, color='tab:blue', alpha=0.7)

    for i, j, k in triangles:
        verts = [[points[i], points[j], points[k]]]
        poly = Poly3DCollection(verts, alpha=0.16, facecolor='gray', edgecolor='gray')
        ax.add_collection3d(poly)

    for i, j, k, l in tetrahedra:
        faces = [
            [points[i], points[j], points[k]],
            [points[i], points[j], points[l]],
            [points[i], points[k], points[l]],
            [points[j], points[k], points[l]],
        ]
        poly = Poly3DCollection(faces, alpha=0.08, facecolor='pink', edgecolor='gray')
        ax.add_collection3d(poly)

    ax.set_box_aspect([1, 1, 1])


if __name__ == '__main__':

    type = "torus"

    root_dir = Path(__file__).resolve().parent.parent.parent
    filtration = load_filtration(root_dir / f'examples/data/{type}_filtration.csv')
    births, deaths, dimensions = load_persistence(root_dir / f'examples/data/{type}_persistence.csv', filtration)
    points = load_points(root_dir / f'examples/data/{type}.csv')

    finite_deaths = [d for d in deaths if d is not None]
    max_val = max([b for b in births if b is not None] + finite_deaths)
    slider_eps = 0.0

    fig1, ax1 = plt.subplots(figsize=(6, 6))
    fig1.canvas.manager.set_window_title('Persistence Diagram')
    plt.subplots_adjust(left=0.18, bottom=0.2)

    colors = ['red', 'blue', 'green', 'purple', 'brown', 'orange']
    dims = sorted(set(dimensions))
    dim_labels = [f"H{dim}" for dim in dims]
    for idx, dim in enumerate(dims):
        xs = [b for b, d in zip(births, dimensions) if d == dim]
        ys = [deaths[i] if deaths[i] is not None else max_val * 1.05 for i, d in enumerate(dimensions) if d == dim]
        ax1.scatter(xs, ys, label=dim_labels[idx], color=colors[idx % len(colors)])

    ax1.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    threshold_line = ax1.axvline(slider_eps, color='gray', linestyle='--', alpha=0.5)

    ax1.set_xlabel('Birth ε')
    ax1.set_ylabel('Death ε')
    ax1.set_title('Persistence Diagram')
    ax1.legend(loc='upper left')
    ax1.set_xlim(0, max_val * 1)
    ax1.set_ylim(0, max_val * 1)

    slider_ax = fig1.add_axes([0.18, 0.08, 0.7, 0.04])
    eps_slider = Slider(slider_ax, 'Max ε', 0.0, max_val, valinit=slider_eps, valstep=max_val / 200)

    fig2 = plt.figure(figsize=(8, 8))
    ax2 = fig2.add_subplot(projection='3d')
    fig2.canvas.manager.set_window_title('3D Rips Complex')

    point_ids, edges, triangles, tetrahedra = get_filtration_at_eps(filtration, slider_eps)
    plot_complex(ax2, points, point_ids, edges, triangles, tetrahedra, slider_eps)

    def update(val):
        current_eps = eps_slider.val
        threshold_line.set_xdata([current_eps, current_eps])

        point_ids, edges, triangles, tetrahedra = get_filtration_at_eps(filtration, current_eps)
        plot_complex(ax2, points, point_ids, edges, triangles, tetrahedra, current_eps)
        fig2.canvas.draw_idle()
        fig1.canvas.draw_idle()

    eps_slider.on_changed(update)
    plt.show()
