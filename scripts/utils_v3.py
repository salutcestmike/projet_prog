import numpy as np
from vedo import Points, Sphere, Plotter, Axes


RADIUS = 1

def sphere_points(n, radius = RADIUS):
    thetas = []
    phis = []
    points = []


    for k in range(1, n + 1):
        hk = -1 + 2 * ((k - 1) / (n - 1))
        theta = np.arccos(hk)

        if ((k == 1) or (k == n)):
            phi = 0
        else:
            phi = (phis[k - 2] + (3.6 / np.sqrt(n)) * (1 / np.sqrt(1 - (hk ** 2)))) % (2*np.pi)

        thetas.append(theta)
        phis.append(phi)
        points.append([np.sin(theta) * np.cos(phi) * radius,
                       np.sin(theta) * np.sin(phi) * radius, 
                       np.cos(theta) * radius])
    print(thetas)
    print(len(phis))
    return np.array(points)



points = sphere_points(92)
print(points)
print("Nombre de points :", len(points))

# # Affichage
# for i, (x, y, z) in enumerate(points, start=1):
#     print(f"{i:3d}  {x: .6f}  {y: .6f}  {z: .6f}")


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
