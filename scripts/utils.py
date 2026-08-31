import numpy as np
from vedo import Points, Sphere, Plotter, Axes


def generate_92_points(radius = 1, rings = 18, points_per_ring = 5):
    """
    Génère 92 points sur une sphère unité :
        - 2 pôles
        - 18 anneaux
        - 5 points par anneau

    Total = 2 + 18*5 = 92
    """

    points = [[0.0, 0.0, radius], [0.0, 0.0, -radius]]

    for k in range(1, rings + 1):

        z = 1.0 - 2.0 * k / (rings + 1)

        # Rayon du cercle correspondant sur la sphère
        rho = np.sqrt(1.0 - z**2)

        # 5-fold symmetry
        for j in range(points_per_ring):

            theta = 2.0 * np.pi * j / points_per_ring

            x = rho * np.cos(theta)
            y = rho * np.sin(theta)

            points.append([x, y, z])

    return np.array(points)


points = generate_92_points()

print("Nombre de points :", len(points))

# Affichage
for i, (x, y, z) in enumerate(points, start=1):
    print(f"{i:3d}  {x: .6f}  {y: .6f}  {z: .6f}")


# ============================================
# Points
# ============================================

pts = Points(
    points,
    r=10,
    c="red"
)


# ============================================
# Sphère unité transparente
# ============================================

sphere = Sphere(
    pos=(0, 0, 0),
    r=1
).alpha(0.15).c("lightblue")


# ============================================
# Axes
# ============================================

axes = Axes(
    xtitle="X",
    ytitle="Y",
    ztitle="Z"
)


# ============================================
# Affichage interactif
# ============================================

plotter = Plotter(
    axes=axes,
    bg="white"
)

plotter.show(
    sphere,
    pts,
    axes,
    viewup="z",
    interactive=True
)
