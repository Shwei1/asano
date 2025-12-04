import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
import random
import tkinter as tk
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from partition import *
from graph_algorithms import *
import ast



class ImaiAsanoGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Imai–Asano Rectilinear Polygon Partition")

        self.outer = []         
        self.holes = []         
        self.chosen_H = []      
        self.chosen_V = []      

        control_frame = tk.Frame(master)
        control_frame.pack(side=tk.TOP, fill=tk.X)

        canvas_frame = tk.Frame(master)
        canvas_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        btn_load = tk.Button(
            control_frame,
            text="Load polygon from file",
            command=self.load_polygon_from_file,
        )
        btn_load.pack(side=tk.LEFT, padx=5, pady=5)

        btn_random = tk.Button(
            control_frame,
            text="Random rectilinear polygon",
            command=self.generate_random_polygon,
        )
        btn_random.pack(side=tk.LEFT, padx=5, pady=5)

        btn_partition = tk.Button(
            control_frame,
            text="Partition",
            command=self.partition_current_polygon,
        )
        btn_partition.pack(side=tk.LEFT, padx=5, pady=5)

        btn_clear = tk.Button(
            control_frame,
            text="Clear",
            command=self.clear_all,
        )
        btn_clear.pack(side=tk.LEFT, padx=5, pady=5)

        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect("equal", "box")
        self.ax.grid(False)

        self.canvas = FigureCanvasTkAgg(self.fig, master=canvas_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)

        self.redraw()

    def load_polygon_from_file(self):
        filename = filedialog.askopenfilename(
            title="Select polygon file",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if not filename:
            return

        try:
            with open(filename, "r") as f:
                contents = f.read()
        except OSError as e:
            messagebox.showerror("Error", f"Cannot open file:\n{e}")
            return

        if "outer: " not in contents:
            messagebox.showerror(
                "Error",
                'File must contain outer polygon as\nouter: [[x1, y1], ..., [xn, yn]]'
            )
            return

        outer = []
        holes = []
        scanning_holes = False
        holes_description = ""

        for line in contents.splitlines():
            if "outer: " in line:
                content = line.replace("outer: ", "")
                outer = ast.literal_eval(content)
            elif "holes: " in line:
                scanning_holes = True
                content = line.replace("holes: ", "")
                holes_description += content
            elif scanning_holes:
                holes_description += line

        if holes_description.strip():
            holes = ast.literal_eval(holes_description)
        else:
            holes = []

        outer = [list(v) for v in outer]
        holes = [[list(v) for v in hole] for hole in holes]

        self.outer = outer
        self.holes = holes
        self.chosen_H = []
        self.chosen_V = []

        self.redraw()



    def generate_random_polygon(self,
                                target_vertices: int = 100,
                                width: int = 60,
                                height: int = 30,
                                max_holes: int = 5):
        """
        Генерує більш складний ортогональний багатокутник із ~target_vertices вершинами
        (x-монотонний "гістограмний" багатокутник над віссю y=0) та кількома прямокутними отворами.

        - outer: приблизно 2 * n_cols + 3 вершин ≈ target_vertices
        - holes: до max_holes прямокутних отворів, що:
            * повністю лежать усередині outer,
            * не перетинаються і не торкаються одне одного.
        """

        import random

        n_cols = max(2, (target_vertices - 3) // 2)

        width = max(width, n_cols + 2)  # 0, ..., width; внутрішні x = 1..width-1

        xs_inner = sorted(random.sample(range(1, width), n_cols - 1))
        xs = [0] + xs_inner + [width]  
        m = len(xs)                    

        hs = [random.randint(height // 3, height) for _ in range(m)]

        verts = []
        verts.append([xs[0], 0])

        for i in range(m):
            verts.append([xs[i], hs[i]])
            if i < m - 1:
                verts.append([xs[i + 1], hs[i]])

        verts.append([xs[-1], 0])

        outer = []
        for v in verts:
            if not outer or outer[-1] != v:
                outer.append(v)

        def roof_height_at(x):
            if x <= xs[0]:
                return hs[0]
            if x >= xs[-1]:
                return hs[-1]
            for i in range(m - 1):
                if xs[i] <= x < xs[i + 1]:
                    return hs[i]
            return hs[-1]

        holes = []
        hole_rects = [] 

        def rect_overlap_or_touch(r1, r2):
            x1a, x2a, y1a, y2a = r1
            x1b, x2b, y1b, y2b = r2
            return not (x2a <= x1b or x1a >= x2b or y2a <= y1b or y1a >= y2b)

        max_trials = 200 

        for _ in range(max_holes):
            placed = False
            for _t in range(max_trials):
                max_hole_w = min(6, width - 2)  
                if max_hole_w < 2:
                    break
                w_hole = random.randint(2, max_hole_w)

                x1 = random.randint(1, width - 1 - w_hole)
                x2 = x1 + w_hole

                h_min = min(roof_height_at(x) for x in range(x1, x2))

                if h_min <= 3:
                    continue

                y2 = random.randint(2, h_min - 1)
                if y2 <= 2:
                    continue
                y1 = random.randint(1, y2 - 1)

                rect = (x1, x2, y1, y2)

                bad = False
                for r in hole_rects:
                    if rect_overlap_or_touch(rect, r):
                        bad = True
                        break
                if bad:
                    continue

                hole_rects.append(rect)
                hole = [
                    [x1, y1],
                    [x1, y2],
                    [x2, y2],
                    [x2, y1],
                ]
                holes.append(hole)
                placed = True
                break

            if not placed:
                break

        self.outer = outer
        self.holes = holes
        self.chosen_H = []
        self.chosen_V = []
        self.redraw()


    def partition_current_polygon(self):
        if not self.outer:
            messagebox.showinfo("Info", "No polygon loaded or generated.")
            return

        try:
            outer = [v[:] for v in self.outer]      
            holes = [[p[:] for p in hole] for hole in self.holes]

            outer, holes, chosen_H, chosen_V = compute_partition(outer, holes)

        except Exception as e:
            messagebox.showerror("Error during partition", str(e))
            return

        self.outer = outer
        self.holes = holes
        self.chosen_H = chosen_H
        self.chosen_V = chosen_V

        self.redraw()

    def clear_all(self):
        self.outer = []
        self.holes = []
        self.chosen_H = []
        self.chosen_V = []
        self.redraw()

    def redraw(self):
        self.ax.clear()
        self.ax.set_aspect("equal", "box")
        self.ax.grid(False)

        if self.outer:
            xs = [v[0] for v in self.outer] + [self.outer[0][0]]
            ys = [v[1] for v in self.outer] + [self.outer[0][1]]
            self.ax.plot(xs, ys, "k-")

        for hole in self.holes:
            xs = [v[0] for v in hole] + [hole[0][0]]
            ys = [v[1] for v in hole] + [hole[0][1]]
            self.ax.plot(xs, ys, "k-")

        for seg in self.chosen_H:
            x1, x2 = seg["x_from"], seg["x_to"]
            y = seg["y"]
            self.ax.plot([x1, x2], [y, y], linewidth=1)

        for seg in self.chosen_V:
            x = seg["x"]
            y1, y2 = seg["y_from"], seg["y_to"]
            self.ax.plot([x, x], [y1, y2], linewidth=1)

        self.canvas.draw()


if __name__ == "__main__":
    root = tk.Tk()
    app = ImaiAsanoGUI(root)
    root.mainloop()



