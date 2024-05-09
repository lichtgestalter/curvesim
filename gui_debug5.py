import tkinter as tk
import math


class PointViewer(tk.Canvas):
    def __init__(self, master, points, viewer_start_position, canvas_size, axes_length, position_delta):
        super().__init__(master, width=canvas_size[0], height=canvas_size[1])
        self.points = points
        self.viewer_position = list(viewer_start_position)
        self.canvas_size = canvas_size
        self.axes_length = axes_length
        self.position_delta = position_delta

        self.draw_axes()
        self.draw_points()
        self.display_viewer_position()

        self.pack()
        self.create_buttons()

    def draw_axes(self):
        # x-axis
        self.create_line(0, self.canvas_size[1] // 2, self.canvas_size[0], self.canvas_size[1] // 2, fill="black")
        self.create_text(self.canvas_size[0] - 10, self.canvas_size[1] // 2 - 10, text="x", fill="black")

        # y-axis
        self.create_line(self.canvas_size[0] // 2, 0, self.canvas_size[0] // 2, self.canvas_size[1], fill="black")
        self.create_text(self.canvas_size[0] // 2 + 10, 10, text="y", fill="black")

        # z-axis
        self.create_line(self.canvas_size[0] // 2, self.canvas_size[1], self.canvas_size[0] // 2, 0, fill="black")
        self.create_text(self.canvas_size[0] // 2 - 10, self.canvas_size[1] - 10, text="z", fill="black")

    def draw_points(self):
        for point in self.points:
            self.draw_point(point)

    def draw_point(self, point):
        perspective_point = self.project_point(point)
        x, y = PointViewer.round_screen_coords(perspective_point)
        self.create_oval(x - 5, y - 5, x + 5, y + 5, fill="green")

    def project_point(self, point):
        x_diff = point[0] - self.viewer_position[0]
        y_diff = point[1] - self.viewer_position[1]
        z_diff = point[2] - self.viewer_position[2]

        if z_diff == 0:
            z_diff = 0.0001  # Prevent division by zero error

        scale_factor = self.axes_length / z_diff
        x_perspective = self.canvas_size[0] // 2 + x_diff * scale_factor
        y_perspective = self.canvas_size[1] // 2 - y_diff * scale_factor  # Invert y-axis for tkinter
        return x_perspective, y_perspective

    @staticmethod
    def round_screen_coords(point):
        return round(point[0]), round(point[1])

    def display_viewer_position(self):
        position_text = f"Pos: {self.viewer_position[0]}, {self.viewer_position[1]}, {self.viewer_position[2]}"
        self.create_text(10, 10, text=position_text, anchor="nw")

    def create_buttons(self):
        button_texts = ["x-", "x+", "y-", "y+", "z-", "z+"]
        button_commands = [self.move_x_neg, self.move_x_pos, self.move_y_neg, self.move_y_pos, self.move_z_neg, self.move_z_pos]
        for command, text in zip(button_commands, button_texts):
            button = tk.Button(self.master, text=text, command=command)
            button.pack(side="left")

    def update_canvas(self):
        self.delete("all")
        self.draw_axes()
        self.draw_points()
        self.display_viewer_position()

    def move_x_neg(self):
        self.viewer_position[0] -= self.position_delta
        self.update_canvas()

    def move_x_pos(self):
        self.viewer_position[0] += self.position_delta
        self.update_canvas()

    def move_y_neg(self):
        self.viewer_position[1] -= self.position_delta
        self.update_canvas()

    def move_y_pos(self):
        self.viewer_position[1] += self.position_delta
        self.update_canvas()

    def move_z_neg(self):
        self.viewer_position[2] -= self.position_delta
        self.update_canvas()

    def move_z_pos(self):
        self.viewer_position[2] += self.position_delta
        self.update_canvas()


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
    axes_length = 200
    position_delta = 100
    app = PointViewer(root, points, viewer_start_position, canvas_size, axes_length, position_delta)
    root.mainloop()
