import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def set_axs(axs, ax_index, focal_length):
    axs[ax_index].set_proj_type('persp', focal_length=focal_length)
    axs[ax_index].set_title(f" testtitle {ax_index}", fontsize=10)
    # print camera position
    print(ax_index)
    print(axs[ax_index].get_position())
    print(axs[ax_index].get_proj())

def main():
    fig, axs = plt.subplots(1, 3, subplot_kw={'projection': '3d'})

    points = [
        (130.00, 0.00, 0.00),
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

    for ax_index in range(3):
        # axs[ax_index].set_proj_type('persp', focal_length=10)
        set_axs(axs, ax_index, 10)

    plt.show()


if __name__ == "__main__":
    main()

# --------------------
# In the view_init method, elev controls the elevation angle in the z plane,
# and azim controls the azimuth angle in the x, y plane. Adjust these values
# accordingly to set your desired viewpoint and viewing direction.
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import numpy as np
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# # Set viewpoint to (-1000, -500, 2000) and direct view towards (0, 0, 0)
# ax.view_init(elev=10, azim=120)
#
# # Plot your 3D data here
# # For example, plotting a point at (0, 0, 0)
# ax.scatter(0, 0, 0, color='r')
#
# plt.show()
