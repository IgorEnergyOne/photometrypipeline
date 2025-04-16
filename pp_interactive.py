import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from tkinter import filedialog, messagebox, simpledialog
from matplotlib.patches import Patch
import os
import numpy as np
from astropy.time import Time

DEFAULT_COLORS = ["blue", "red", "green", "orange", "purple", "brown", "pink", "gray", "olive", "cyan"]

class LightCurveData:
    def __init__(self):
        self.df = None
        self.filename = None

    def load(self, filepath):
        self.filename = filepath
        self.df = pd.read_csv(filepath)
        if 'rejected' not in self.df.columns:
            self.df['rejected'] = False

    def save(self):
        if self.df is not None and self.filename:
            self.df.to_csv(self.filename, index=False)

    def toggle_rejection(self, index):
        self.df.loc[index, 'rejected'] = not self.df.loc[index, 'rejected']

class LightCurvePlot:
    def __init__(self, figure, ax, master_frame):
        self.fig = figure
        self.ax = ax
        self.master_frame = master_frame
        self.selected_index = None
        self.xlabel = "Julian Date"
        self.ylabel = "Magnitude"
        self.title = "Lightcurve"
        self.x_data = None
        self.valid_color = "blue"
        self.rejected_color = "red"
        self.legend_frame = None

    def clear(self):
        self.ax.clear()

    def set_labels(self):
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(self.title)
        self.ax.title.set_picker(True)
        self.ax.xaxis.label.set_picker(True)
        self.ax.yaxis.label.set_picker(True)
        self.ax.invert_yaxis()

    def draw(self):
        self.fig.canvas.draw()

    def update(self, df, mode, time_mode, show_rejected, update_legend=True):
        self.clear()
        jd = df['julian_date']
        if time_mode == 'julian_date':
            x = jd
        elif time_mode == 'mjd':
            x = jd - 2400000.5
        elif time_mode == 'minutes':
            x = (jd - jd.min()) * 24 * 60
        self.x_data = x.to_numpy()

        self.xlabel = {
            'minutes': f"Minutes from {Time(jd.min(), format='jd').to_value('iso', subfmt='date_hm')} UT",
            'julian_date': "Julian Date",
            'mjd': "Modified Julian Date (MJD)"
        }[time_mode]

        y = df[f'{"mag" if mode == "target" else "inst_mag" if mode == "instrumental" else "mag_control"}'].to_numpy()
        yerr = df[f'{"sig" if mode == "target" else "inst_sig" if mode == "instrumental" else "sig_control"}'].to_numpy()

        mask_valid = ~df['rejected']
        mask_rejected = df['rejected']

        self.ax.errorbar(self.x_data[mask_valid], y[mask_valid], yerr=yerr[mask_valid], fmt='o', color=self.valid_color, picker=5)
        if show_rejected:
            self.ax.errorbar(self.x_data[mask_rejected], y[mask_rejected], yerr=yerr[mask_rejected], fmt='o', color=self.rejected_color, picker=5)

        if self.selected_index is not None:
            sx, sy = self.x_data[self.selected_index], y[self.selected_index]
            self.ax.plot(sx, sy, 'o', color='orange', markersize=10)
            label = os.path.splitext(os.path.basename(df.loc[self.selected_index, 'filename']))[0] if 'filename' in df.columns else str(self.selected_index)
            self.ax.text(sx, sy + 0.1, label, color='orange')

        if update_legend:
            self.update_legend()

        self.set_labels()
        self.draw()

    def update_legend(self):
        if self.legend_frame and self.legend_frame.winfo_exists():
            self.legend_frame.destroy()
        self.legend_frame = ttk.Frame(self.master_frame)
        self.legend_frame.pack(before=self.fig.canvas.get_tk_widget(), side="top", pady=5)
        legend_fig = plt.figure(figsize=(4, 0.5))
        legend_ax = legend_fig.add_subplot(111)
        legend_ax.axis('off')
        legend_ax.legend(handles=[
            Patch(color=self.valid_color, label="Valid"),
            Patch(color=self.rejected_color, label="Rejected"),
            Patch(color="orange", label="Selected")
        ], loc='center', ncol=3, frameon=False)
        legend_canvas = FigureCanvasTkAgg(legend_fig, master=self.legend_frame)
        legend_canvas.draw()
        legend_canvas.get_tk_widget().pack()

    def select_point(self, x_click, y_click, df, mode):
        y = df[f'{"mag" if mode == "target" else "inst_mag" if mode == "instrumental" else "mag_control"}'].to_numpy()
        distances = np.sqrt((self.x_data - x_click)**2 + (y - y_click)**2)
        self.selected_index = int(np.argmin(distances)) if distances.min() < 0.1 else None

    def move_selection(self, direction):
        if self.x_data is None:
            return
        if self.selected_index is None:
            self.selected_index = 0 if direction > 0 else len(self.x_data) - 1
        else:
            self.selected_index = max(0, min(len(self.x_data) - 1, self.selected_index + direction))

class LightCurveGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Lightcurve Viewer and Editor")

        self.data = LightCurveData()
        self.mode = 'target'
        self.time_mode = 'minutes'
        self.show_rejected = True

        self.create_widgets()
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.root.bind("q", self.confirm_exit)
        self.root.bind("r", self.toggle_rejection)
        self.root.bind("<Left>", self.move_left)
        self.root.bind("<Right>", self.move_right)

    def create_widgets(self):
        self.master_frame = ttk.Frame(self.root)
        self.master_frame.pack(fill=BOTH, expand=True)

        control_frame = ttk.Frame(self.master_frame, padding=10)
        control_frame.pack(side=TOP, fill=X)

        ttk.Button(control_frame, text="Open CSV", command=self.open_csv).pack(side=LEFT, padx=5)
        ttk.Button(control_frame, text="Save CSV", command=self.save_csv).pack(side=LEFT, padx=5)
        ttk.Button(control_frame, text="Save Plot", command=self.save_plot).pack(side=LEFT, padx=5)

        ttk.Label(control_frame, text="Mode:").pack(side=LEFT)
        self.mode_var = ttk.StringVar(value='target')
        ttk.Combobox(control_frame, textvariable=self.mode_var, values=['target', 'instrumental', 'control'], state='readonly').pack(side=LEFT)
        self.mode_var.trace_add('write', lambda *_: self.set_mode())

        ttk.Label(control_frame, text="Time Axis:").pack(side=LEFT)
        self.time_var = ttk.StringVar(value='minutes')
        ttk.Combobox(control_frame, textvariable=self.time_var, values=['minutes', 'julian_date', 'mjd'], state='readonly').pack(side=LEFT)
        self.time_var.trace_add('write', lambda *_: self.set_time_mode())

        ttk.Label(control_frame, text="Valid Color:").pack(side=LEFT)
        self.valid_color_var = ttk.StringVar(value="blue")
        ttk.Combobox(control_frame, textvariable=self.valid_color_var, values=DEFAULT_COLORS, state='readonly').pack(side=LEFT)
        self.valid_color_var.trace_add('write', lambda *_: self.update_colors())

        ttk.Label(control_frame, text="Rejected Color:").pack(side=LEFT)
        self.rejected_color_var = ttk.StringVar(value="red")
        ttk.Combobox(control_frame, textvariable=self.rejected_color_var, values=DEFAULT_COLORS, state='readonly').pack(side=LEFT)
        self.rejected_color_var.trace_add('write', lambda *_: self.update_colors())

        self.toggle_rejected_var = ttk.BooleanVar(value=True)
        ttk.Checkbutton(control_frame, text="Show Rejected", variable=self.toggle_rejected_var, command=self.set_show_rejected).pack(side=LEFT)

        plot_frame = ttk.Frame(self.master_frame)
        plot_frame.pack(fill=BOTH, expand=True, padx=10, pady=5)
        self.fig, self.ax = plt.subplots(figsize=(8, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, plot_frame).update()

        self.plot = LightCurvePlot(self.fig, self.ax, self.master_frame)

        self.canvas.mpl_connect("button_press_event", self.on_click)
        self.canvas.mpl_connect("pick_event", self.on_pick)

        help_frame = ttk.Frame(self.master_frame, padding=10)
        help_frame.pack(side=BOTTOM, fill=X)
        ttk.Label(help_frame, text="Hotkeys: [r] toggle rejection | [q] quit | [arrow left]/[arrow right] move selection | click = select/unselect or edit title/label").pack()



    def update_colors(self):
        self.plot.valid_color = self.valid_color_var.get()
        self.plot.rejected_color = self.rejected_color_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, update_legend=True)

    def open_csv(self):
        file = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
        if file:
            self.data.load(file)
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected)

    def save_csv(self):
        self.data.save()

    def save_plot(self):
        if self.data.df is not None:
            filename = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG Files", "*.png")])
            if filename:
                self.plot.fig.savefig(filename, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Saved", f"Plot saved as {filename}")

    def set_mode(self, _=None):
        self.mode = self.mode_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected)

    def set_time_mode(self, _=None):
        self.time_mode = self.time_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected)

    def set_show_rejected(self):
        self.show_rejected = self.toggle_rejected_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected)

    def on_click(self, event):
        if event.inaxes != self.ax or self.data.df is None:
            return
        if event.xdata and event.ydata:
            self.plot.select_point(event.xdata, event.ydata, self.data.df, self.mode)
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, update_legend=False)

    def on_pick(self, event):
        artist = event.artist
        if artist == self.ax.title:
            new = simpledialog.askstring("Edit Title", "New Title:", initialvalue=self.plot.title)
            if new: self.plot.title = new
        elif artist == self.ax.xaxis.label:
            new = simpledialog.askstring("Edit X Label", "New X Label:", initialvalue=self.plot.xlabel)
            if new: self.plot.xlabel = new
        elif artist == self.ax.yaxis.label:
            new = simpledialog.askstring("Edit Y Label", "New Y Label:", initialvalue=self.plot.ylabel)
            if new: self.plot.ylabel = new
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, update_legend=False)

    def toggle_rejection(self, _=None):
        if self.plot.selected_index is not None:
            self.data.toggle_rejection(self.plot.selected_index)
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, update_legend=False)

    def move_left(self, _=None):
        self.plot.move_selection(-1)
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, update_legend=False)

    def move_right(self, _=None):
        self.plot.move_selection(1)
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, update_legend=False)

    def confirm_exit(self, _=None):
        self.on_close()

    def on_close(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.root.quit()
            self.root.destroy()

if __name__ == '__main__':
    root = ttk.Window(themename="flatly")
    app = LightCurveGUI(root)
    root.mainloop()
