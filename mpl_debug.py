import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def set_axs(axs, elev, azim,focal_length):
    axs.set_proj_type('persp', focal_length=focal_length)
    axs.set_title(f" testtitle", fontsize=10)
    axs.view_init(elev=elev, azim=azim)
    print(axs.get_position())
    print(axs.get_proj())


def plot_planes(axs, alpha=0.3, min_=-200, max_=200):
    # show the x=0, y=0 and z=0 plane as transparent polygons
    minx, miny, minz = min_, min_, min_
    maxx, maxy, maxz = max_, max_, max_
    colorx, colory, colorz = 'red', 'blue', 'green'
    x0_plane = [(0, miny, minz), (0, miny, maxz), (0, maxy, maxz), (0, maxy, minz)]
    y0_plane = [(minx, 0, minz), (minx, 0, maxz), (maxx, 0, maxz), (maxx, 0, minz)]
    z0_plane = [(minx, miny, 0), (minx, maxy, 0), (maxx, maxy, 0), (maxx, miny, 0)]
    poly_x = Poly3DCollection([x0_plane], color=colorx, alpha=alpha)
    poly_y = Poly3DCollection([y0_plane], color=colory, alpha=alpha)
    poly_z = Poly3DCollection([z0_plane], color=colorz, alpha=alpha)
    axs.add_collection3d(poly_x)
    axs.add_collection3d(poly_y)
    axs.add_collection3d(poly_z)


def main():
    fig, axs = plt.subplots(1, 1, subplot_kw={'projection': '3d'})
    set_axs(axs, 30, 120, 10)  # define projection and view

    points = [
        (131.00, 0.00, 0.00),
        (29.10, 6.76, 6.76),
        (23.11, 118.42, 118.42),
        (- 3.53, 37.73, 37.73),
        (- 106.35, 47.04, 47.04),
        (- 149.96, -130.32, -30.32),
        (- 142.39, -34.84, -134.84),
        (1.57, -35.27, -35.27),
        (9.43, -30.68, -30.68),
        (25.51, -14.96, -14.96),
        (28.33, -9.20, -9.20),
    ]
    x, y, z = zip(*points)

    axs.scatter(x, y, z)  # Plot the points
    axs.plot(x, y, z, color='black')  # connect the points
    plot_planes(axs)  # show the x=0, y=0 and z=0 plane as transparent polygons


    plt.show()


if __name__ == "__main__":
    main()
