import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def set_axs(axs, ax_index, focal_length):
    axs[ax_index].set_proj_type('persp', focal_length=focal_length)
    axs[ax_index].set_title(f" testtitle {ax_index}", fontsize=10)
    axs[ax_index].view_init(elev=10*ax_index, azim=90+30*ax_index)
    print(ax_index)
    print(axs[ax_index].get_position())
    print(axs[ax_index].get_proj())


def plot_planes(axs):
    minx, miny, minz = -200, -200, -200
    maxx, maxy, maxz = 200, 200, 200
    colorx, colory, colorz = 'red', 'blue', 'green'
    alpha = 0.3
    # Plot a blue transparent (alpha=0.3) rectangle which represents the y=0 plane. The 3D coordinate
    x0_plane = [(0, miny, minz), (0, miny, maxz), (0, maxy, maxz), (0, maxy, minz)]
    y0_plane = [(minx, 0, minz), (minx, 0, maxz), (maxx, 0, maxz), (maxx, 0, minz)]
    z0_plane = [(minx, miny, 0), (minx, maxy, 0), (maxx, maxy, 0), (maxx, miny, 0)]
    for ax in axs:
        # show the y=0 plane as a blue, transparent polygon in the x-z plane
        poly_x = Poly3DCollection([x0_plane], color=colorx, alpha=alpha)
        poly_y = Poly3DCollection([y0_plane], color=colory, alpha=alpha)
        poly_z = Poly3DCollection([z0_plane], color=colorz, alpha=alpha)
        ax.add_collection3d(poly_x)
        ax.add_collection3d(poly_y)
        ax.add_collection3d(poly_z)


def main():
    fig, axs = plt.subplots(1, 3, subplot_kw={'projection': '3d'})

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

    # Plot the data
    for ax in axs:
        ax.plot(x, y, z, color='black')  # Add a black line connecting the points
        ax.scatter(x, y, z)

    plot_planes(axs)

    for ax_index in range(3):
        set_axs(axs, ax_index, 10)

    plt.show()


if __name__ == "__main__":
    main()
