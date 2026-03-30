# -*- coding: utf-8 -*-
"""
AsteroidImageViewer – popup window that displays the image for the currently
selected photometry point.
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import TYPE_CHECKING, Optional
import pandas as pd
import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from .utils import _safe
try:
    from PIL import Image, ImageTk
except ImportError:
    Image = None
    ImageTk = None
if TYPE_CHECKING:
    from .gui import LightCurveGUI
class AsteroidImageViewer:
    """Manages the asteroid image window: loading, zooming, overlay toggle,
    and persisting window state/zoom level."""
    def __init__(self, root: tk.Tk, data_provider: "LightCurveGUI"):
        self.root = root
        self.data_provider = data_provider
        self.zoom_level = 1.0
        self.window_geometry: Optional[str] = None
        self.image_window: Optional[tk.Toplevel] = None
        self._current_image_original = None
        self._current_overlay_original = None
        self.show_overlay_var = tk.BooleanVar(value=True)
    # ——— public API ———
    def toggle_visibility(self):
        """Toggle the visibility of the asteroid image window."""
        if self.image_window is not None and self.image_window.winfo_exists():
            if self.image_window.state() == 'withdrawn':
                self._restore_window()
            else:
                self._hide_window()
        else:
            self._create_window()
    # ——— window lifecycle ———
    def _restore_window(self):
        if self.window_geometry:
            _safe(self.image_window.geometry, self.window_geometry)
        self.image_window.deiconify()
        self.image_window.lift()
        self.root.after(100, lambda: self.root.focus_force())
    def _hide_window(self):
        self.window_geometry = _safe(self.image_window.geometry)
        self.image_window.withdraw()
    def _create_window(self):
        top = tk.Toplevel(self.root)
        self.image_window = top
        top.title("Asteroid image")
        top.transient(self.root)
        if self.window_geometry:
            _safe(top.geometry, self.window_geometry)
        self._bind_events(top)
        self.update_image()
        chk = ttk.Checkbutton(top, text="Toggle Overlay (O)",
                               variable=self.show_overlay_var,
                               command=self.refresh_display)
        chk.pack(side=BOTTOM, pady=5)
        self.root.focus_set()
        def return_focus(event=None):
            if self.root.winfo_exists():
                self.root.focus_set()
        top.bind("<Map>", return_focus)
        self.root.after(50, lambda: self.root.focus_force())
        self.root.after(200, lambda: self.root.focus_force())
    def _bind_events(self, win: tk.Toplevel):
        def _on_configure(evt=None):
            if evt and evt.widget == win:
                self.window_geometry = win.geometry()
        win.bind("<Configure>", _on_configure, add="+")
        win.bind("<MouseWheel>", self._on_zoom)
        win.bind("<Button-4>", self._on_zoom)
        win.bind("<Button-5>", self._on_zoom)
        win.bind("<o>", self._toggle_overlay_hotkey)
        win.bind("<O>", self._toggle_overlay_hotkey)
    # ——— image loading & display ———
    def update_image(self):
        """Load and display the image for the current selection."""
        if not self.image_window or not self.image_window.winfo_exists():
            return
        for widget in self.image_window.winfo_children():
            if widget.winfo_class() != 'TCheckbutton':
                widget.destroy()
        alias = self.data_provider.current_lc_alias
        idx = self.data_provider.current_point_index
        if alias is None or idx is None:
            ttk.Label(self.image_window, text="No point selected.").pack(padx=20, pady=20)
            return
        lc = self.data_provider.lightcurves.get(alias)
        if not lc:
            return
        try:
            row_data = lc.df.loc[idx]
        except (KeyError, IndexError):
            return
        raw_filename = None
        if 'filename' in row_data and pd.notna(row_data['filename']):
            raw_filename = str(row_data['filename']).strip()
        if not raw_filename:
            ttk.Label(self.image_window, text="No filename in data.").pack(padx=20, pady=20)
            return
        base_dir = os.getcwd()
        if lc.filename:
            base_dir = os.path.dirname(os.path.abspath(lc.filename))
        candidate_paths = []
        if os.path.isabs(raw_filename):
            candidate_paths.append(Path(raw_filename))
        else:
            bp = Path(base_dir)
            candidate_paths += [bp / raw_filename, bp.parent / raw_filename, bp.parent.parent / raw_filename]
        found_path = next((p for p in candidate_paths if p.exists()), None)
        if not found_path:
            msg = ("Image not found.\nSearched for: " + raw_filename + "\nIn:\n"
                   + "\n".join(str(p.parent) for p in candidate_paths))
            ttk.Label(self.image_window, text=msg).pack(padx=20, pady=20)
            return
        try:
            if Image is None:
                ttk.Label(self.image_window, text="Pillow not installed.").pack(padx=20, pady=20)
                return
            self._current_image_original = Image.open(str(found_path)).convert("RGBA")
            self._current_overlay_original = None
            potential_overlays = []
            if found_path.suffix.lower() != '.png':
                potential_overlays.append(found_path.with_suffix('.png'))
            potential_overlays.append(found_path.with_name(f"{found_path.stem}_overlay.png"))
            for ov_p in potential_overlays:
                if ov_p.exists():
                    self._current_overlay_original = Image.open(str(ov_p)).convert("RGBA")
                    break
            self.refresh_display()
            self.image_window.title(f"Image: {found_path.name}")
        except Exception as e:
            ttk.Label(self.image_window, text=f"Error loading image:\n{e}").pack(padx=20, pady=20)
    def refresh_display(self):
        """Resize and display the stored image at the current zoom level."""
        if not self.image_window or not self.image_window.winfo_exists():
            return
        if self._current_image_original is None:
            return
        img_label = None
        for widget in self.image_window.winfo_children():
            if isinstance(widget, ttk.Label) and hasattr(widget, 'image'):
                img_label = widget
                break
        if img_label is None:
            for widget in self.image_window.winfo_children():
                if widget.winfo_class() != 'TCheckbutton':
                    widget.destroy()
            img_label = ttk.Label(self.image_window)
            img_label.pack(padx=10, pady=10)
        orig_w, orig_h = self._current_image_original.size
        new_w, new_h = int(orig_w * self.zoom_level), int(orig_h * self.zoom_level)
        resized_img = self._current_image_original.resize((new_w, new_h), Image.LANCZOS)
        if self._current_overlay_original and self.show_overlay_var.get():
            ov = self._current_overlay_original.resize((new_w, new_h), Image.LANCZOS)
            resized_img.paste(ov, (0, 0), ov)
        photo = ImageTk.PhotoImage(resized_img)
        img_label.configure(image=photo)
        img_label.image = photo
    def _on_zoom(self, event):
        if self._current_image_original is None:
            return
        scale_factor = 1.1
        if event.num == 4 or event.delta > 0:
            self.zoom_level *= scale_factor
        elif event.num == 5 or event.delta < 0:
            self.zoom_level /= scale_factor
        self.zoom_level = max(0.1, min(20.0, self.zoom_level))
        self.refresh_display()
    def _toggle_overlay_hotkey(self, event=None):
        self.show_overlay_var.set(not self.show_overlay_var.get())
        self.refresh_display()
