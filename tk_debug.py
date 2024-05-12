import tkinter as tk

class PointViewer(tk.Canvas):
    def __init__(self, master, points, viewer_start_position, canvas_size, axes_length, position_delta, point_radius):
        super().__init__(master, width=canvas_size[0], height=canvas_size[1])
        self.points = points
        self.viewer_position = viewer_start_position
        self.canvas_size = canvas_size
        self.axes_length = axes_length
        self.position_delta = position_delta
        self.point_radius = point_radius
        self.canvas_center = [self.canvas_size[0] // 2, self.canvas_size[1] // 2]
        # self.offset = self.canvas_center.copy()   # Used to keep the origin (0, 0, 0) in the middle of the canvas
        self.offset = [0, 0]
        print(f"Canvas center: {self.canvas_center[0]:.0f}, {self.canvas_center[1]:.0f}")
        self.pack()
        self.create_buttons()
        self.update_canvas()

    def update_canvas(self):
        self.delete("all")
        self.update_offset()
        self.draw_axes()
        self.draw_points()
        self.display_viewer_position()

    def update_offset(self):
        self.origin = self.project_point([0, 0, 0])
        print(f"Origin before update: {self.origin[0]:.0f}, {self.origin[1]:.0f}")
        print(f"Offset before update: {self.offset[0]:.0f}, {self.offset[1]:.0f}")
        # self.offset[0] -= self.canvas_center[0]
        # self.offset[1] -= self.canvas_center[1]
        self.offset[0] += self.origin[0] - self.canvas_center[0]
        self.offset[1] += self.origin[1] - self.canvas_center[1]
        self.origin = self.project_point([0, 0, 0])
        print(f"Origin after update: {self.origin[0]:.0f}, {self.origin[1]:.0f}")
        print(f"Offset after update: {self.offset[0]:.0f}, {self.offset[1]:.0f}")

    def draw_axes(self):
        pad = 10
        # origin = self.project_point([0, 0, 0])
        x_arrow = self.project_point([self.axes_length, 0, 0])
        y_arrow = self.project_point([0, self.axes_length, 0])
        z_arrow = self.project_point([0, 0, self.axes_length])
        self.create_line(self.origin[0], self.origin[1], x_arrow[0], x_arrow[1], fill="red")
        self.create_line(self.origin[0], self.origin[1], y_arrow[0], y_arrow[1], fill="green")
        self.create_line(self.origin[0], self.origin[1], z_arrow[0], z_arrow[1], fill="blue")
        self.create_text(x_arrow[0], x_arrow[1] + pad, text="x", fill="red")
        self.create_text(y_arrow[0], y_arrow[1] + pad, text="y", fill="green")
        self.create_text(z_arrow[0], z_arrow[1] + pad, text="z", fill="blue")

    def draw_points(self):
        for point in self.points:
            self.draw_point(point, self.point_radius)
        coords = [self.project_point(point) for point in self.points]
        self.create_polygon(coords, fill="", outline="black", tags="points")

    def draw_point(self, point, radius):
        perspective_point = self.project_point(point)
        x, y = PointViewer.round_screen_coords(perspective_point)
        self.create_oval(x - radius, y - radius, x + radius, y + radius, fill="green")

    def project_point(self, point):
        x_diff = point[0] - self.viewer_position[0]
        y_diff = point[1] - self.viewer_position[1]
        z_diff = point[2] - self.viewer_position[2]

        if z_diff == 0:
            z_diff = 0.01  # Prevent division by zero error

        scale_factor = self.axes_length / z_diff
        x_perspective = self.canvas_size[0] // 2 + x_diff * scale_factor
        x_perspective -= self.offset[0]
        y_perspective = self.canvas_size[1] // 2 - y_diff * scale_factor  # Invert y-axis for tkinter
        y_perspective -= self.offset[1]
        return [x_perspective, y_perspective]

    @staticmethod
    def round_screen_coords(point):
        return round(point[0]), round(point[1])

    def display_viewer_position(self):
        position_text = f"Pos: {self.viewer_position[0]}, {self.viewer_position[1]}, {self.viewer_position[2]} Offset: {self.offset[0]:.0f}, {self.offset[1]:.0f}"
        self.create_text(10, 10, text=position_text, anchor="nw")

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

    def create_buttons(self):
        button_texts = ["x-", "x+", "y-", "y+", "z-", "z+"]
        button_commands = [self.move_x_neg, self.move_x_pos, self.move_y_neg, self.move_y_pos, self.move_z_neg, self.move_z_pos]
        for command, text in zip(button_commands, button_texts):
            button = tk.Button(self.master, text=text, command=command)
            button.pack(side="left")


def main():
    main_window = tk.Tk()
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
    viewer_start_position = [100, 100, 1000]
    canvas_size = (1300, 900)
    axes_length = 200
    position_delta = 100
    point_radius = 2
    PointViewer(main_window, points, viewer_start_position, canvas_size, axes_length, position_delta, point_radius)
    main_window.mainloop()


if __name__ == "__main__":
    main()
