import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def set_axs(axs, ax_index, focal_length, testy):
    axs[ax_index].set_proj_type('persp', focal_length=focal_length)
    axs[ax_index].set_title(f" testtitle {ax_index}", fontsize=10)
    # print camera position
    print(axs[ax_index].get_proj())
    print(testy)

def main():
    fig, axs = plt.subplots(1, 3, subplot_kw={'projection': '3d'})

    # Get the test data
    x, y, z = axes3d.get_test_data(0.09)

    points = [
        (30.00, 0.00, 0.00),
        (29.10, 6.76, 6.76),
        (23.11, 18.42, 18.42),
        (- 3.53, 37.73, 37.73),
        (- 106.35, 47.04, 47.04),
        (- 149.96, -30.32, -30.32),
        (- 142.39, -34.84, -34.84),
        (1.57, -35.27, -35.27),
        (9.43, -30.68, -30.68),
        (25.51, -14.96, -14.96),
        (28.33, -9.20, -9.20),
    ]

    # Plot the data
    for ax in axs:
        ax.plot_wireframe(x, y, z, rstride=10, cstride=10)

    for ax_index in range(3):
        set_axs(axs, ax_index, 10, 42 + ax_index * 5)

    plt.show()


if __name__ == "__main__":
    main()
