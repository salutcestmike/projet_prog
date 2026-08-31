import numpy as np
from vedo import Points, Sphere, Plotter, Axes


RADIUS = 5

def fibonacci_sphere(n, radius = RADIUS):
    points = []

    golden_angle = np.pi * (3 - np.sqrt(5))

    for i in range(n):

        z = (1 - 2 * (i + 0.5) / n)
        r = np.sqrt(1 - z*z) * RADIUS 

        theta = golden_angle * i

        x = r * np.cos(theta)
        y = r * np.sin(theta)

        points.append([x, y, z * RADIUS])

    return np.array(points)



points = fibonacci_sphere(92)
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
    r=RADIUS
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
