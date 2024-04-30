import matplotlib.pyplot as plt

myfile = "debug_file.txt"
with open(myfile) as file:
    lines = file.readlines()
for line in lines:
    if "a=" in line and "kleinomegaquer" in line:
        print(f"{line.strip()=}")
        continue
    while " = (" in line:
        print(f"L  {line.strip()=}")
print()
exit(1)


# Define the points
L0 = (96.98, 17.10, -17.36)
L22 = (96.15, 16.95, 21.64)
L45 = (80.67, 14.22, 57.36)
L90 = (17.10, 3.02, 98.48)
L135 = (-56.49, -9.96, 81.91)
L180 = (-96.98, -17.10, 17.36)
L225 = (-80.67, -14.22, -57.36)
L270 = (-17.10, -3.02, -98.48)
L315 = (56.49, 9.96, -81.91)

# Create the plot
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Left plot
axs[0].scatter([L0[0], L22[0], L45[0], L90[0], L135[0], L180[0], L225[0], L270[0], L315[0]],
               [L0[2], L22[2], L45[2], L90[2], L135[2], L180[2], L225[2], L270[2], L315[2]], color='green')

axs[0].plot([L0[0], L22[0], L45[0], L90[0], L135[0], L180[0], L225[0], L270[0], L315[0], L0[0]],
            [L0[2], L22[2], L45[2], L90[2], L135[2], L180[2], L225[2], L270[2], L315[2], L0[2]], color='grey')

axs[0].set_title('left/top a=100 e=0.00 i=90 O=10 koq=0')
axs[0].set_xlabel('x')
axs[0].set_ylabel('z')
axs[0].set_facecolor('black')
axs[0].grid(False)
axs[0].spines['bottom'].set_color('grey')
axs[0].spines['top'].set_color('grey')
axs[0].spines['left'].set_color('grey')
axs[0].spines['right'].set_color('grey')

# Right plot
axs[1].scatter([L0[0], L22[0], L45[0], L90[0], L135[0], L180[0], L225[0], L270[0], L315[0]],
               [L0[1], L22[1], L45[1], L90[1], L135[1], L180[1], L225[1], L270[1], L315[1]], color='green')

axs[1].plot([L0[0], L22[0], L45[0], L90[0], L135[0], L180[0], L225[0], L270[0], L315[0], L0[0]],
            [L0[1], L22[1], L45[1], L90[1], L135[1], L180[1], L225[1], L270[1], L315[1], L0[1]], color='grey')

axs[1].set_title('right/edge a=100 e=0.00 i=90 O=10 koq=0')
axs[1].set_xlabel('x')
axs[1].set_ylabel('y')
axs[1].set_facecolor('black')
axs[1].grid(False)
axs[1].spines['bottom'].set_color('grey')
axs[1].spines['top'].set_color('grey')
axs[1].spines['left'].set_color('grey')
axs[1].spines['right'].set_color('grey')

axs[0].set_xlim(-200, 200)
axs[0].set_ylim(-200, 200)
axs[1].set_xlim(-200, 200)
axs[1].set_ylim(-200, 200)

plt.tight_layout()
plt.show()

# Save the plot
fig.savefig('ellipse.png')