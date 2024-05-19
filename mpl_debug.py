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

def read_points_from_file(filename):
    with open(filename, encoding='utf-8') as file:
        lines = file.readlines()
    x_lists = []
    y_lists = []
    z_lists = []
    params_list = []
    for line in lines:
        if line[0] == "a":  # header line
            params = line.strip()
            params_list.append(params)
            x_list = []
            y_list = []
            z_list = []
            x_lists.append(x_list)
            y_lists.append(y_list)
            z_lists.append(z_list)
        elif line[0] == "L":  # point line
            x, y, z = [float(i) for i in line.strip().split(",")[1:]]  # [1:] because the first value is the Label
            x_list.append(x)
            y_list.append(y)
            z_list.append(z)
            # print(f"{x=} {y=} {z=} ")
    return x_lists, y_lists, z_lists, params_list

def main():
    fig, axs = plt.subplots(1, 1, subplot_kw={'projection': '3d'})
    set_axs(axs, 30, 200, 10)  # define projection and view
    x_lists, y_lists, z_lists, params_list = read_points_from_file("debug_file.txt")
    axs.scatter(x_lists[0], y_lists[0], z_lists[0])  # Plot the points
    axs.plot(x_lists[0], y_lists[0], z_lists[0], color='black')  # connect the points
    plot_planes(axs)  # show the x=0, y=0 and z=0 plane as transparent polygons
    plt.tight_layout()
    fig.savefig(params_list[0] + ".png")
    plt.show()


if __name__ == "__main__":
    main()
