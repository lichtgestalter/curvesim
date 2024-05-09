import tkinter as tk

class PointViewer(tk.Frame):
    def __init__(self, master, points):
        super().__init__(master)
        self.master = master
        self.pack(fill=tk.BOTH, expand=True)
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

    def draw_axes(self):
        # Clear canvas
        self.canvas.delete("axes")

        # Draw x-axis
        self.canvas.create_line(300, 300, 600, 300, arrow=tk.LAST, tags="axes")
        self.canvas.create_text(600, 300, text="x", anchor=tk.CENTER, tags="axes")

        # Draw y-axis
        self.canvas.create_line(300, 300, 300, 0, arrow=tk.LAST, tags="axes")
        self.canvas.create_text(300, 20, text="y", anchor=tk.CENTER, tags="axes")

        # Draw z-axis
        self.canvas.create_line(300, 300, 300, 600, arrow=tk.LAST, tags="axes")
        self.canvas.create_text(20, 300, text="z", anchor=tk.CENTER, tags="axes")

    def draw_points(self):
        # Clear canvas
        self.canvas.delete("points")

        # Draw points
        for point in self.points:
            label, x, y, z = point
            projected_x = (x - self.viewer_position[0]) * 300 / (z - self.viewer_position[2]) + 300
            projected_y = (y - self.viewer_position[1]) * 300 / (z - self.viewer_position[2]) + 300
            self.canvas.create_oval(projected_x - 5, projected_y - 5, projected_x + 5, projected_y + 5, fill="green", tags="points")

        # Connect points with lines
        for i in range(len(self.points) - 1):
            x1, y1, z1 = self.points[i][1:]
            x2, y2, z2 = self.points[i+1][1:]
            projected_x1 = (x1 - self.viewer_position[0]) * 300 / (z1 - self.viewer_position[2]) + 300
            projected_y1 = (y1 - self.viewer_position[1]) * 300 / (z1 - self.viewer_position[2]) + 300
            projected_x2 = (x2 - self.viewer_position[0]) * 300 / (z2 - self.viewer_position[2]) + 300
            projected_y2 = (y2 - self.viewer_position[1]) * 300 / (z2 - self.viewer_position[2]) + 300
            self.canvas.create_line(projected_x1, projected_y1, projected_x2, projected_y2, fill="black", tags="points")


if __name__ == "__main__":
    root = tk.Tk()
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
    viewer_start_position = (100, 100, 1000)
    canvas_size = (600, 600)
    axes_lenght = 200
    position_delta = 100
    app = PointViewer(root, points, viewer_start_position, canvas_size, axes_lenght, position_delta)
    root.mainloop()