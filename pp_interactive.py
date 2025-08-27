#!/usr/bin/env python3

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
MARKERS = [ ('o', 'circle'), ('s', 'square'), ('p', 'pentagon'), ('x', 'x'), ('D', 'diamond'),
            ('*', 'star'), ('v', 'triangle_down'), ('^', 'triangle_up'), ('<', 'triangle_left'),
            ('>', 'triangle_right'), ('+', 'plus'), ('d', 'thin_diamond'),
            ]

def next_version(path: str) -> str:
    """
    Return path like originalname_1.ext, originalname_2.ext …,
    choosing the first number that is not on disk.
    """
    base, ext = os.path.splitext(path)
    i = 1
    while os.path.exists(f"{base}_{i}{ext}"):
        i += 1
    return f"{base}_{i}{ext}"


class LightCurveData:
    def __init__(self):
        self.df = None
        self.filename = None

    def load(self, filepath):
        self.filename = filepath
        self.df = pd.read_csv(filepath)
        if 'rejected' not in self.df.columns:
            self.df['rejected'] = False
        if 'sextractor_flags' not in self.df.columns:
            self.df['sextractor_flags'] = 0

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
        self.auto_title = True  # Start with auto title enabled
        self.x_data = None
        self.valid_color = "blue"
        self.rejected_color = "red"
        self.flagged_color = "orange"
        self.legend_frame = None

        # Marker and error bar properties
        self.marker_size = 5.0
        self.marker_style = 'o'  # Default circle marker
        self.errorbar_capsize = 3.0
        self.errorbar_capthick = 1.0
        self.errorbar_linewidth = 1.0

        # Available marker styles
        self.available_markers = MARKERS
        self.marker_dict = {name: marker for marker, name in self.available_markers}

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

    def target_name_parser(self, df):
        """parse target name from csv file"""
        try:
            target_name = df['target'].iloc[0]

        except KeyError:
            target_name = 'Lightcurve'
        return target_name

    def draw(self):
        self.fig.canvas.draw()

    def update(self, df, mode, time_mode, show_rejected, errorbar_type='calibrated', update_legend=True):
        self.clear()
        try:
            jd = df['julian_date']
        except KeyError:
            # show massage if no data for control star is available
            messagebox.showwarning("Warning", "No date is available (wrong file?).")
            return None
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

        if mode == 'target':
            y = df['mag'].to_numpy()
        elif mode == 'instrumental':
            y = df['inst_mag'].to_numpy()
        elif mode == 'control':
            try:
                y = df['mag_control'].to_numpy()
            except KeyError:
                # show massage if no data for control star is available
                messagebox.showwarning("Warning", "No control star data available.")
                return None
        else:
            raise ValueError(f"Invalid mode: {mode}")


        # Determine which error bars to use based on the selected type
        if errorbar_type == 'instrumental':
            yerr = df['inst_sig'].to_numpy()
        elif errorbar_type == 'calibrated':
            yerr = df['sig' if mode != 'control' else 'sig_control'].to_numpy()
        else:  # 'none'
            yerr = None

        mask_valid = ~df['rejected'] & (df['sextractor_flags'] == 0)
        mask_rejected = df['rejected']
        mask_flagged = (df['sextractor_flags'] > 0) & (~df['rejected'])

        # Plot with or without error bars based on selection
        errorbar_kwargs = {
            'fmt': self.marker_style,
            'markersize': self.marker_size,
            'capsize': self.errorbar_capsize,
            'capthick': self.errorbar_capthick,
            'elinewidth': self.errorbar_linewidth,
            'picker': 5
        }

        plot_kwargs = {
            'marker': self.marker_style,
            'markersize': self.marker_size,
            'linestyle': 'None',
            'picker': 5
        }

        if yerr is not None:
            self.ax.errorbar(self.x_data[mask_valid], y[mask_valid], yerr=yerr[mask_valid],
                           color=self.valid_color, **errorbar_kwargs)
            self.ax.errorbar(self.x_data[mask_flagged], y[mask_flagged], yerr=yerr[mask_flagged],
                           color=self.flagged_color, **errorbar_kwargs)
            if show_rejected:
                self.ax.errorbar(self.x_data[mask_rejected], y[mask_rejected],
                               yerr=yerr[mask_rejected], color=self.rejected_color, **errorbar_kwargs)
        else:
            self.ax.plot(self.x_data[mask_valid], y[mask_valid], color=self.valid_color, **plot_kwargs)
            self.ax.plot(self.x_data[mask_flagged], y[mask_flagged], color=self.flagged_color, **plot_kwargs)
            if show_rejected:
                self.ax.plot(self.x_data[mask_rejected], y[mask_rejected],
                           color=self.rejected_color, **plot_kwargs)

        if self.selected_index is not None:
            sx, sy = self.x_data[self.selected_index], y[self.selected_index]
            self.ax.plot(sx, sy, 'o', color='orange', markersize=10)
            label = os.path.splitext(os.path.basename(df.loc[self.selected_index, 'filename']))[0] if 'filename' in df.columns else str(self.selected_index)
            # plot at least 5 yticks on the plot
            self.ax.text(sx, sy + 0.1, label, color='orange')

        if update_legend:
            self.update_legend()

        # Update title only if in auto mode
        if self.auto_title:
            self.title = self.target_name_parser(df)

        # TODO: make at least 5 yticks on the plot
        self.set_labels()
        self.draw()

    def update_legend(self):
        if self.legend_frame and self.legend_frame.winfo_exists():
            self.legend_frame.destroy()
        self.legend_frame = ttk.Frame(self.master_frame)
        self.legend_frame.pack(before=self.fig.canvas.get_tk_widget(), side="top", pady=5)

        # Create a frame for the legend items
        legend_items_frame = ttk.Frame(self.legend_frame)
        legend_items_frame.pack()

        # Create legend items with clickable patches
        self.legend_handles = {
            'valid': self.valid_color,
            'rejected': self.rejected_color,
            'flagged': self.flagged_color
        }

        # Create a list of (label, color_var) pairs
        legend_items = [
            ("Valid", 'valid'),
            ("Rejected", 'rejected'),
            ("Flagged", 'flagged')
        ]

        for i, (label, color_key) in enumerate(legend_items):
            frame = ttk.Frame(legend_items_frame)
            frame.pack(side=LEFT, padx=10)

            # Create a colored patch using a simple colored frame
            color = self.legend_handles[color_key]
            patch_frame = ttk.Frame(frame, width=20, height=20)
            patch_frame.pack(side=LEFT, padx=2)
            patch_frame.pack_propagate(False)  # Prevent frame from resizing

            # Create a colored label inside the frame
            patch = ttk.Label(patch_frame, background=color, borderwidth=1, relief='solid')
            patch.pack(fill='both', expand=True)

            # Store the patch for later updates
            if not hasattr(self, 'legend_patches'):
                self.legend_patches = {}
            self.legend_patches[color_key] = patch

            # Add label
            ttk.Label(frame, text=label).pack(side=LEFT, padx=2)

            # Make the patch and label clickable
            def on_patch_click(event, key=color_key):
                self.on_legend_click(key)

            patch.bind('<Button-1>', on_patch_click)
            patch_frame.bind('<Button-1>', on_patch_click)

            # Also make the label clickable
            for widget in frame.winfo_children():
                if isinstance(widget, ttk.Label):
                    widget.bind('<Button-1>', on_patch_click)

    def on_legend_click(self, color_key):
        """Handle click on legend patch with a dropdown color picker"""
        # Create a new top-level window
        color_dialog = ttk.Toplevel()
        color_dialog.title(f"Select {color_key.capitalize()} Color")
        color_dialog.transient(self.master_frame)  # Set to be on top of the main window
        color_dialog.grab_set()  # Make the dialog modal

        # Position the dialog near the mouse click
        x = self.master_frame.winfo_pointerx()
        y = self.master_frame.winfo_pointery()
        color_dialog.geometry(f"+{x}+{y}")

        # Create a frame for the color selection
        frame = ttk.Frame(color_dialog, padding="10")
        frame.pack(fill=BOTH, expand=True)

        # Label
        ttk.Label(frame, text=f"Select color for {color_key} points:").pack(pady=5)

        # Create a combobox with color options
        color_var = ttk.StringVar(value=self.legend_handles[color_key])
        color_combo = ttk.Combobox(
            frame,
            textvariable=color_var,
            values=DEFAULT_COLORS,
            state='readonly',
            width=15
        )
        color_combo.pack(pady=5)

        # Preview the selected color
        preview_frame = ttk.Frame(frame, height=30, width=100)
        preview_frame.pack_propagate(False)
        preview_frame.pack(pady=5)
        preview = ttk.Frame(preview_frame, style=f"{color_var.get().title()}.TFrame")
        preview.pack(fill=BOTH, expand=True)

        def update_preview(event=None):
            color = color_var.get()
            preview.configure(style=f"{color.title()}.TFrame")

        color_var.trace_add('write', lambda *_: update_preview())

        # Create styles for color preview
        for color in DEFAULT_COLORS:
            ttk.Style().configure(f"{color.title()}.TFrame", background=color)

        # OK and Cancel buttons
        button_frame = ttk.Frame(frame)
        button_frame.pack(pady=5)

        def apply_color():
            color = color_var.get()
            if color and color in DEFAULT_COLORS:
                # Update the corresponding color variable
                if color_key == 'valid':
                    self.valid_color = color
                elif color_key == 'rejected':
                    self.rejected_color = color
                elif color_key == 'flagged':
                    self.flagged_color = color

                # Update the legend patch color
                if hasattr(self, 'legend_patches') and color_key in self.legend_patches:
                    self.legend_patches[color_key].configure(background=color)

                # Get the current plot data from the parent GUI
                if hasattr(self, 'parent_gui') and hasattr(self.parent_gui, 'data') and self.parent_gui.data.df is not None:
                    # Update the plot with new colors
                    self.update(
                        self.parent_gui.data.df,
                        self.parent_gui.mode,
                        self.parent_gui.time_mode,
                        self.parent_gui.show_rejected,
                        self.parent_gui.errorbar_type,
                        update_legend=False
                    )
                    self.fig.canvas.draw_idle()
            color_dialog.destroy()

        ttk.Button(button_frame, text="OK", command=apply_color).pack(side=LEFT, padx=5)
        ttk.Button(button_frame, text="Cancel", command=color_dialog.destroy).pack(side=LEFT, padx=5)

        # Bind Enter key to apply color
        color_dialog.bind('<Return>', lambda e: apply_color())
        color_dialog.bind('<Escape>', lambda e: color_dialog.destroy())

        # Initial preview update
        update_preview()

    def select_point(self, x_click, y_click, df, mode):
        y = df[f'{"mag" if mode == "target" else "inst_mag" if mode == "instrumental" else "mag_control"}'].to_numpy()
        # check the extent of the data
        x_ext = max(self.x_data) - min(self.x_data)
        y_ext = max(y) - min(y)
        extent = max(x_ext, y_ext)
        distances = np.sqrt((self.x_data - x_click)**2 + (y - y_click)**2)
        self.selected_index = int(np.argmin(distances)) if distances.min() < max(extent*0.01, 0.1) else None

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
        self.errorbar_type = 'calibrated'  # Can be 'instrumental', 'calibrated', or 'none'

        self.create_widgets()
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.root.bind("q", self.confirm_exit)
        self.root.bind("r", self.toggle_rejection)
        self.root.bind("a", self.cancel_selection)
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

        # Add error bar type selection
        ttk.Label(control_frame, text="Error Bars:").pack(side=LEFT, padx=(10, 0))
        self.errorbar_var = ttk.StringVar(value='calibrated')
        ttk.Combobox(control_frame, textvariable=self.errorbar_var,
                    values=['calibrated', 'instrumental', 'none'],
                    state='readonly', width=12).pack(side=LEFT)
        self.errorbar_var.trace_add('write', lambda *_: self.set_errorbar_type())

        ttk.Label(control_frame, text="Time Axis:").pack(side=LEFT)
        self.time_var = ttk.StringVar(value='minutes')
        ttk.Combobox(control_frame, textvariable=self.time_var, values=['minutes', 'julian_date', 'mjd'], state='readonly').pack(side=LEFT)
        self.time_var.trace_add('write', lambda *_: self.set_time_mode())

        self.toggle_rejected_var = ttk.BooleanVar(value=True)
        ttk.Checkbutton(control_frame, text="Show Rejected", variable=self.toggle_rejected_var, command=self.set_show_rejected).pack(side=LEFT)

        # Add marker settings button
        ttk.Button(control_frame, text="Marker Settings", command=self.show_marker_settings).pack(side=LEFT, padx=5)

        plot_frame = ttk.Frame(self.master_frame)
        plot_frame.pack(fill=BOTH, expand=True, padx=10, pady=5)
        self.fig, self.ax = plt.subplots(figsize=(8, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, plot_frame).update()

        self.plot = LightCurvePlot(self.fig, self.ax, self.master_frame)
        self.plot.parent_gui = self  # Set reference to parent GUI

        self.canvas.mpl_connect("button_press_event", self.on_click)
        self.canvas.mpl_connect("pick_event", self.on_pick)

        help_frame = ttk.Frame(self.master_frame, padding=10)
        help_frame.pack(side=BOTTOM, fill=X)
        ttk.Label(help_frame, text="Hotkeys: [r] toggle rejection "
                                   "| [a] cancel selection "
                                   "| [q] quit | "
                                   "[arrow left]/[arrow right] move selection "
                                   "| click = select/unselect or edit title/label").pack()


    def update_colors(self):
        self.plot.valid_color = self.valid_color_var.get()
        self.plot.rejected_color = self.rejected_color_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, update_legend=True)

    def open_csv(self):
        file = filedialog.askopenfilename(initialdir=os.getcwd(),
                                          filetypes=[("CSV Files", "*.csv")])
        if file:
            self.data.load(file)
            #  ⬇️  put filename in the title bar
            self.root.title(f"Lightcurve Viewer – {os.path.basename(file)}")
            self.plot.update(self.data.df, self.mode, self.time_mode,
                             self.show_rejected)

    def save_csv(self):
        if self.data.df is None or self.data.filename is None:
            messagebox.showerror("Nothing to save", "Load a CSV first.")
            return

        # ask the user; empty string means they pressed “Cancel”
        file = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(self.data.filename),
            initialfile=os.path.basename(self.data.filename),
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv")]
        )

        if not file:  # <- user hit Cancel
            file = next_version(self.data.filename)

        self.data.filename = file
        self.data.save()
        messagebox.showinfo("Saved", f"CSV saved as {file}")

    def save_plot(self):
        if self.data.df is None:
            return

        default_png = os.path.splitext(self.data.filename or "plot")[0] + ".png"

        file = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(default_png),
            initialfile=os.path.basename(default_png),
            defaultextension=".png",
            filetypes=[("PNG Files", "*.png")]
        )

        if not file:
            file = next_version(default_png)

        self.plot.fig.savefig(file, dpi=300, bbox_inches="tight")
        messagebox.showinfo("Saved", f"Plot saved as {file}")

    def set_mode(self, _=None):
        self.mode = self.mode_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected)

    def set_time_mode(self, _=None):
        self.time_mode = self.time_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected)

    def set_show_rejected(self):
        self.show_rejected = self.toggle_rejected_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type)

    def set_errorbar_type(self, _=None):
        self.errorbar_type = self.errorbar_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type)

    def on_click(self, event):
        if event.inaxes != self.ax or self.data.df is None:
            return
        if event.xdata and event.ydata:
            self.plot.select_point(event.xdata, event.ydata, self.data.df, self.mode)
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, update_legend=False)

    def on_pick(self, event):
        artist = event.artist
        if artist == self.ax.title:
            new = simpledialog.askstring("Edit Title", "New Title:", initialvalue=self.plot.title)
            if new: self.plot.title = new
            self.plot.auto_title = False
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
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, update_legend=False)

    def move_left(self, _=None):
        self.plot.move_selection(-1)
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, update_legend=False)

    def move_right(self, _=None):
        self.plot.move_selection(1)
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, update_legend=False)
    
    def cancel_selection(self, _=None):
        """Cancel the current point selection."""
        if self.plot.selected_index is not None:
            self.plot.selected_index = None
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, update_legend=False)

    def confirm_exit(self, _=None):
        self.on_close()

    def show_marker_settings(self):
        """Show dialog for marker and error bar settings"""
        settings_dialog = ttk.Toplevel()
        settings_dialog.title("Marker and Error Bar Settings")
        settings_dialog.transient(self.root)
        settings_dialog.grab_set()

        # Position the dialog near the main window
        x = self.root.winfo_x() + 50
        y = self.root.winfo_y() + 50
        settings_dialog.geometry(f"+{x}+{y}")

        frame = ttk.Frame(settings_dialog, padding=10)
        frame.pack(fill=BOTH, expand=True)

        # Marker style
        ttk.Label(frame, text="Marker Style:").grid(row=0, column=0, sticky=W, pady=2)
        marker_names = [name for marker, name in self.plot.available_markers]
        current_marker_name = next((name for marker, name in self.plot.available_markers
                                 if marker == self.plot.marker_style), 'circle')
        marker_var = ttk.StringVar(value=current_marker_name)
        marker_combo = ttk.Combobox(frame, textvariable=marker_var, values=marker_names, state='readonly')
        marker_combo.grid(row=0, column=1, sticky=EW, pady=2, padx=5)

        # Marker size
        ttk.Label(frame, text="Marker Size:").grid(row=1, column=0, sticky=W, pady=2)
        size_var = ttk.DoubleVar(value=self.plot.marker_size)
        size_scale = ttk.Scale(frame, from_=1, to=30, variable=size_var, orient=HORIZONTAL)
        size_scale.grid(row=1, column=1, sticky=EW, pady=2, padx=5)
        size_entry = ttk.Entry(frame, textvariable=size_var, width=5)
        size_entry.grid(row=1, column=2, sticky=W, pady=2, padx=5)

        # Error bar cap size
        ttk.Label(frame, text="Error Cap Size:").grid(row=2, column=0, sticky=W, pady=2)
        capsize_var = ttk.DoubleVar(value=self.plot.errorbar_capsize)
        capsize_scale = ttk.Scale(frame, from_=0, to=20, variable=capsize_var, orient=HORIZONTAL)
        capsize_scale.grid(row=2, column=1, sticky=EW, pady=2, padx=5)
        capsize_entry = ttk.Entry(frame, textvariable=capsize_var, width=5)
        capsize_entry.grid(row=2, column=2, sticky=W, pady=2, padx=5)

        # Error bar cap thickness
        ttk.Label(frame, text="Cap Thickness:").grid(row=3, column=0, sticky=W, pady=2)
        capthick_var = ttk.DoubleVar(value=self.plot.errorbar_capthick)
        capthick_scale = ttk.Scale(frame, from_=0.0, to=10, variable=capthick_var, orient=HORIZONTAL)
        capthick_scale.grid(row=3, column=1, sticky=EW, pady=2, padx=5)
        capthick_entry = ttk.Entry(frame, textvariable=capthick_var, width=5)
        capthick_entry.grid(row=3, column=2, sticky=W, pady=2, padx=5)

        # Error bar line width
        ttk.Label(frame, text="Error Bar Width:").grid(row=4, column=0, sticky=W, pady=2)
        linewidth_var = ttk.DoubleVar(value=self.plot.errorbar_linewidth)
        linewidth_scale = ttk.Scale(frame, from_=0.0, to=10, variable=linewidth_var, orient=HORIZONTAL)
        linewidth_scale.grid(row=4, column=1, sticky=EW, pady=2, padx=5)
        linewidth_entry = ttk.Entry(frame, textvariable=linewidth_var, width=5)
        linewidth_entry.grid(row=4, column=2, sticky=W, pady=2, padx=5)

        # Preview frame
        preview_frame = ttk.LabelFrame(frame, text="Preview", padding=5)
        preview_frame.grid(row=0, column=3, rowspan=5, padx=10, sticky=N+S)

        fig, ax = plt.subplots(figsize=(3, 2), dpi=80)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

        # Add sample points with error bars
        x = [0.2, 0.5, 0.8]
        y = [0.5, 0.5, 0.5]
        yerr = [0.2, 0.2, 0.2]

        preview_line = ax.errorbar(x, y, yerr=yerr, fmt='o', color='blue',
                                 markersize=size_var.get(),
                                 capsize=capsize_var.get(),
                                 capthick=capthick_var.get(),
                                 elinewidth=linewidth_var.get())

        canvas = FigureCanvasTkAgg(fig, master=preview_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=BOTH, expand=True)

        def update_preview(*args):
            try:
                # Get the current marker style
                marker_style = next((marker for marker, name in self.plot.available_markers
                                   if name == marker_var.get()), 'o')

                # Clear the current plot
                ax.clear()
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)

                # Add sample points with error bars using current settings
                x = [0.2, 0.5, 0.8]
                y = [0.5, 0.5, 0.5]
                yerr = [0.2, 0.2, 0.2]

                # Redraw the errorbar with current settings
                global preview_line
                preview_line = ax.errorbar(x, y, yerr=yerr,
                                         fmt=marker_style,
                                         color='blue',
                                         markersize=size_var.get(),
                                         capsize=capsize_var.get(),
                                         capthick=capthick_var.get(),
                                         elinewidth=linewidth_var.get())

                canvas.draw_idle()
            except Exception as e:
                print(f"Error updating preview: {e}")

        # Bind variables to update preview
        marker_var.trace_add('write', update_preview)
        size_var.trace_add('write', update_preview)
        capsize_var.trace_add('write', update_preview)
        capthick_var.trace_add('write', update_preview)
        linewidth_var.trace_add('write', update_preview)

        # Add Apply and Close buttons
        button_frame = ttk.Frame(frame)
        button_frame.grid(row=5, column=0, columnspan=4, pady=10)

        def apply_settings():
            # Update plot settings
            self.plot.marker_style = next((marker for marker, name in self.plot.available_markers
                                         if name == marker_var.get()), 'o')
            self.plot.marker_size = size_var.get()
            self.plot.errorbar_capsize = capsize_var.get()
            self.plot.errorbar_capthick = capthick_var.get()
            self.plot.errorbar_linewidth = linewidth_var.get()

            # Update the plot
            if self.data.df is not None:
                self.plot.update(
                    self.data.df,
                    self.mode,
                    self.time_mode,
                    self.show_rejected,
                    self.errorbar_type
                )

        ttk.Button(button_frame, text="Apply", command=apply_settings).pack(side=LEFT, padx=5)
        ttk.Button(button_frame, text="Close", command=settings_dialog.destroy).pack(side=LEFT, padx=5)

        # Make the window resizable
        settings_dialog.resizable(True, False)

        # Set focus to the dialog
        settings_dialog.focus_set()

        # Make the dialog modal
        settings_dialog.wait_window()

    def on_close(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            # Close any matplotlib figures
            plt.close('all')
            # Destroy the root window
            self.root.destroy()
            # Exit the application
            self.root.quit()

if __name__ == '__main__':
    root = ttk.Window(themename="flatly")
    app = LightCurveGUI(root)
    root.mainloop()
