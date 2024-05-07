import tkinter as tk
import math

class PointViewer(tk.Tk):
    def __init__(self, points):
        super().__init__()
        self.title("3D Point Viewer")
        self.geometry("800x600")

        self.points = points
        self.viewer_position = [100, 100, 1000]

        self.canvas = tk.Canvas(self, bg="white", width=600, height=600)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.btn_frame = tk.Frame(self)
        self.btn_frame.pack(side=tk.RIGHT, fill=tk.Y)

        self.create_buttons()
        self.draw_axes()
        self.draw_points()

    def create_buttons(self):
        btn_texts = ["x-100", "x+100", "y-100", "y+100", "z-100", "z+100"]
        btn_commands = [self.move_x_minus, self.move_x_plus, self.move_y_minus,
                        self.move_y_plus, self.move_z_minus, self.move_z_plus]

        for text, command in zip(btn_texts, btn_commands):
            btn = tk.Button(self.btn_frame, text=text, command=command)
            btn.pack(pady=10)

    def draw_axes(self):
        # Clear canvas
        self.canvas.delete("axes")

        # Draw x-axis
        self.canvas.create_line(300, 300, 600, 300, arrow=tk.LAST, tags="axes")
        self.canvas.create_text(600, 300, text="x", anchor=tk.CENTER, tags="axes")

        # Draw y-axis
        self.canvas.create_line(300, 300, 300, 0, arrow=tk.LAST, tags="axes")
        self.canvas.create_text(300, 0, text="y", anchor=tk.CENTER, tags="axes")

        # Draw z-axis
        self.canvas.create_line(300, 300, 300, 600, arrow=tk.LAST, tags="axes")
        self.canvas.create_text(300, 600, text="z", anchor=tk.CENTER, tags="axes")

    def draw_points(self):
        # Clear canvas
        self.canvas.delete("points")

        # Draw points
        for point in self.points:
            _, x, y, z = point.split(",")
            x = float(x) - self.viewer_position[0]
            y = float(y) - self.viewer_position[1]
            z = float(z) - self.viewer_position[2]

            # Skip point if z is zero
            if z == 0:
                continue

            # Calculate 2D coordinates using perspective projection
            scaling_factor = 1000 / (1000 + z + 0.01)
            x_2d = 300 + scaling_factor * x
            y_2d = 300 + scaling_factor * y

            # Draw point as green circle
            self.canvas.create_oval(x_2d - 5, y_2d - 5, x_2d + 5, y_2d + 5, fill="green", tags="points")

        # Draw lines connecting points
        for i in range(len(self.points) - 1):
            point1 = self.points[i]
            point2 = self.points[i + 1]
            _, x1, y1, z1 = point1.split(",")
            _, x2, y2, z2 = point2.split(",")
            x1 = float(x1) - self.viewer_position[0]
            y1 = float(y1) - self.viewer_position[1]
            z1 = float(z1) - self.viewer_position[2]
            x2 = float(x2) - self.viewer_position[0]
            y2 = float(y2) - self.viewer_position[1]
            z2 = float(z2) - self.viewer_position[2]

            # Skip line if z of any point is zero
            if z1 == 0 or z2 == 0:
                continue

            # Calculate 2D coordinates using perspective projection
            scaling_factor1 = 1000 / (1000 + z1 + 0.01)
            scaling_factor2 = 1000 / (1000 + z2 + 0.01)
            x1_2d = 300 + scaling_factor1 * x1
            y1_2d = 300 + scaling_factor1 * y1
            x2_2d = 300 + scaling_factor2 * x2
            y2_2d = 300 + scaling_factor2 * y2

            # Draw line
            self.canvas.create_line(x1_2d, y1_2d, x2_2d, y2_2d, fill="black", tags="points")

    def move_x_minus(self):
        self.viewer_position[0] -= 100
        self.draw_axes()
        self.draw_points()

    def move_x_plus(self):
        self.viewer_position[0] += 100
        self.draw_axes()
        self.draw_points()

    def move_y_minus(self):
        self.viewer_position[1] -= 100
        self.draw_axes()
        self.draw_points()

    def move_y_plus(self):
        self.viewer_position[1] += 100
        self.draw_axes()
        self.draw_points()

    def move_z_minus(self):
        self.viewer_position[2] -= 100
        self.draw_axes()
        self.draw_points()

    def move_z_plus(self):
        self.viewer_position[2] += 100
        self.draw_axes()
        self.draw_points()

if __name__ == "__main__":
    points = [
        "L0,30.00,0.00,0.00",
        "L7,29.10,6.76,6.76",
        "L21,23.11,18.42,18.42",
        "L49,-3.53,37.73,37.73",
        "L112,-106.35,47.04,47.04",
        "L217,-149.96,-30.32,-30.32",
        "L224,-142.39,-34.84,-34.84",
        "L315,1.57,-35.27,-35.27",
        "L322,9.43,-30.68,-30.68",
        "L343,25.51,-14.96,-14.96",
        "L350,28.33,-9.20,-9.20"
    ]

    app = PointViewer(points)
    app.mainloop()
