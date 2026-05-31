
import csv
import math
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons
from matplotlib.patches import Circle

def load_filtration(path):
    simplices = []

    with open(path) as f:
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

                if death != None and (abs(birth - death) < 1e-10):
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
            if len(row) >= 2:
                points.append((float(row[0]), float(row[1])))
    return points


def get_filtration_at_eps(filtration, eps):
    simp = [x for x in filtration if x[0] < eps]

    points = set()
    edges = []
    triangles = []

    for _, v in simp:
        if len(v) == 1:
            points.add(v[0])
        elif len(v) == 2:
            edges.append(v)
        elif len(v) == 3:
            triangles.append(v)

    return sorted(points), edges, triangles


def plot_complex(ax, points, edges, triangles, max_eps):
    ax.clear()
    ax.set_title(f"Rips Complex at ε = {max_eps:.3f}")
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    ax.scatter(xs, ys, color='tab:orange', edgecolors='black', zorder=0)
    # Draw balls first
    for x, y in zip(xs, ys):
        circle = Circle(
            (x, y),
            radius=max_eps/2,
            color='tab:orange',
            alpha=0.2,
            zorder=-1
        )
        ax.add_patch(circle)

    for i, j in edges:
        x = [points[i][0], points[j][0]]
        y = [points[i][1], points[j][1]]
        ax.plot(x, y, color='tab:blue', alpha=0.7, zorder=1)

    for i, j, k in triangles:
        xs = [points[i][0], points[j][0], points[k][0]]
        ys = [points[i][1], points[j][1], points[k][1]]
        ax.fill(xs, ys, alpha=0.2, color='grey')

    #ax.relim()
    ax.autoscale_view()


if __name__ == '__main__':
    root_dir = Path(__file__).resolve().parent.parent.parent
    filtration = load_filtration(root_dir / 'filtration.csv')
    births, deaths, dimensions = load_persistence(root_dir / 'persistence.csv', filtration)
    points = load_points(root_dir / 'data.csv')

    finite_deaths = [d for d in deaths if d is not None]
    max_val = max([b for b in births if b is not None] + finite_deaths)
    slider_eps = 9

    fig1, ax1 = plt.subplots(figsize=(6, 6))
    fig1.canvas.manager.set_window_title('Persistence Diagram')
    plt.subplots_adjust(left=0.18, bottom=0.2)

    colors = ['red', 'blue', 'green', 'black', 'purple', 'brown']
    dims = sorted(set(dimensions))
    dim_labels = [f"H{dim}" for dim in dims]
    scatter_plots = []
    for idx, dim in enumerate(dims):
        xs = [b for b, d in zip(births, dimensions) if d == dim]
        ys = [deaths[i] if deaths[i] is not None else max_val * 1.05 for i in range(len(deaths)) if dimensions[i] == dim]
        scatter_plots.append(ax1.scatter(xs, ys, label=dim_labels[idx], color=colors[idx % len(colors)]))

    diag_line, = ax1.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    threshold_line = ax1.axvline(slider_eps, color='gray', linestyle='--', alpha=0.5)

    ax1.set_xlabel('Birth ε')
    ax1.set_ylabel('Death ε')
    ax1.set_title('Persistence Diagram')
    legend = [ax1.legend(loc='upper left')]
    ax1.set_xlim(0, max_val * 1)
    ax1.set_ylim(0, max_val * 1)

    checkbox_ax = fig1.add_axes([0.02, 0.4, 0.12, 0.15])
    check = CheckButtons(checkbox_ax, dim_labels, [True] * len(dim_labels))

    def toggle_visibility(label):
        index = dim_labels.index(label)
        scatter_plots[index].set_visible(not scatter_plots[index].get_visible())
        legend[0].remove()
        legend[0] = ax1.legend(loc='upper left')
        fig1.canvas.draw_idle()

    check.on_clicked(toggle_visibility)

    slider_ax = fig1.add_axes([0.18, 0.08, 0.7, 0.04])
    eps_slider = Slider(slider_ax, 'Max ε', 0.0, max_val / 4, valinit=slider_eps, valstep=max_val / 200)

    fig2, ax2 = plt.subplots(figsize=(6, 6))
    fig2.canvas.manager.set_window_title('Rips Complex')

    p_ids, edges, triangles = get_filtration_at_eps(filtration,0)
    plot_complex(ax2, points, edges, triangles,  slider_eps)

    def update(val):
        current_eps = eps_slider.val
        threshold_line.set_xdata([current_eps, current_eps])

        p_ids, edges, triangles = get_filtration_at_eps(filtration, current_eps)
        plot_complex(ax2, points, edges, triangles,  val)
        fig2.canvas.draw_idle()
        fig1.canvas.draw_idle()

    eps_slider.on_changed(update)

    plt.show()
