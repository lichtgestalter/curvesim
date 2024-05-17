import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def set_axs(axs, ax_index, focal_length):
    axs[ax_index].set_proj_type('persp', focal_length=focal_length)
    axs[ax_index].set_title(f" testtitle {ax_index}", fontsize=10)
    axs[ax_index].view_init(elev=10*ax_index, azim=90+30*ax_index)
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
# To show the y=0 plane and z=0 plane transparently in each plot, you can
# use the axvline and axhline functions in Matplotlib. You can set the
# transparency level using the alpha parameter. Here is an example of how
# you can add transparent planes to your plots:
# In this code snippet, the axhline function is used to plot the y=0 plane,
# and the axvline function is used to plot the z=0 plane. The alpha parameter
# in these functions controls the transparency level of the planes. You can
# adjust the alpha value to make the planes more or less transparent.
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
#
# def set_axs(axs, ax_index, focal_length):
#     axs[ax_index].set_proj_type('persp', focal_length=focal_length)
#     axs[ax_index].set_title(f" testtitle {ax_index}", fontsize=10)
#
#     # Plot transparent y=0 plane
#     axs[ax_index].axhline(0, color='blue', alpha=0.3)
#
#     # Plot transparent z=0 plane
#     axs[ax_index].axvline(0, color='green', alpha=0.3)
#
#     print(ax_index)
#     print(axs[ax_index].get_position())
#     print(axs[ax_index].get_proj())
#
# def main():
#     fig, axs = plt.subplots(1, 3, subplot_kw={'projection': '3d'})
#
#     for i in range(3):
#         set_axs(axs, i, focal_length=1.0)
#
#     plt.show()