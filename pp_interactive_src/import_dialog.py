# -*- coding: utf-8 -*-
"""
Interactive column mapping and preview dialog for custom delimited text files.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd
import tkinter as tk
import ttkbootstrap as ttk
from tkinter import messagebox
from ttkbootstrap.constants import *

from .config import AppConfig, get_config, save_user_config
from .loaders import CustomTextLoader, TextImportConfig


class TextImportDialog:
    """Dialog for previewing and configuring custom delimited text files.

    Allows the user to interactively:
    - Inspect first ~20 lines of the text file in a tabular preview
    - Configure delimiter, header presence, comment prefix, and rows to skip
    - Map columns to canonical fields (time/JD, mag, sig, band, rejected, etc.)
    - Select time representation (Julian Date, MJD, ISO-8601 UTC, or relative time with zero epoch)
    - Save, load, and manage named presets or use the last-used settings
    """

    def __init__(self, parent: tk.Widget, filepath: str, cfg: Optional[AppConfig] = None):
        self.parent = parent
        self.filepath = filepath
        self.cfg = cfg or get_config()
        self.result_config: Optional[TextImportConfig] = None

        self.top = tk.Toplevel(parent)
        self.top.title(f"Import Custom Text File – {os.path.basename(filepath)}")
        self.top.geometry("920x650")
        self.top.minsize(800, 550)
        self.top.transient(parent)
        self.top.grab_set()

        # UI State Variables
        self.delimiter_var = tk.StringVar(value="whitespace")
        self.custom_delim_var = tk.StringVar(value="")
        self.has_header_var = tk.BooleanVar(value=True)
        self.skip_rows_var = tk.IntVar(value=0)
        self.comment_char_var = tk.StringVar(value="#")

        self.time_format_var = tk.StringVar(value="julian_date")
        self.zero_time_var = tk.StringVar(value="")
        self.default_band_var = tk.StringVar(value="V")
        self.rejected_val_var = tk.StringVar(value="!")

        # Column mappings (field -> combobox choice)
        self.col_vars: Dict[str, tk.StringVar] = {
            "julian_date": tk.StringVar(),
            "mag": tk.StringVar(),
            "sig": tk.StringVar(),
            "inst_mag": tk.StringVar(),
            "inst_sig": tk.StringVar(),
            "band": tk.StringVar(),
            "rejected": tk.StringVar(),
            "batch_number": tk.StringVar(),
            "filename": tk.StringVar(),
        }

        self.preview_cols: List[str] = []
        self.preset_var = tk.StringVar(value="[None / Last Used]")

        self._init_initial_values()
        self._build_ui()
        self._update_preview()

        self.top.wait_window()

    def _init_initial_values(self):
        """Seed UI variables from last_used config or automatic detection."""
        ti = self.cfg.text_import
        last = ti.last_used or {}

        # Auto-detect as baseline
        auto_cfg = CustomTextLoader.auto_detect_config(self.filepath)

        self.delimiter_var.set(last.get("delimiter", auto_cfg.delimiter))
        self.skip_rows_var.set(last.get("skip_rows", auto_cfg.skip_rows))
        self.comment_char_var.set(last.get("comment_char", auto_cfg.comment_char))
        self.has_header_var.set(last.get("has_header", auto_cfg.has_header))
        self.time_format_var.set(last.get("time_format", auto_cfg.time_format))
        self.default_band_var.set(last.get("default_band", auto_cfg.default_band))

        # Check for zero time in headers if raw text
        try:
            with open(self.filepath, "r", encoding="utf-8", errors="ignore") as f:
                head = [f.readline() for _ in range(15)]
            import re
            for line in head:
                if "ZERO TIME" in line:
                    m = re.search(r"(\d+\.\d+)", line)
                    if m:
                        self.zero_time_var.set(m.group(1))
                        self.time_format_var.set("relative_days")
                        break
        except Exception:
            pass

    def _build_ui(self):
        """Construct the dialog layout."""
        main_frame = ttk.Frame(self.top, padding=10)
        main_frame.pack(fill=BOTH, expand=True)

        # Top row: Presets
        preset_frame = ttk.LabelFrame(main_frame, text="Preset Configurations", padding=6)
        preset_frame.pack(fill=X, pady=(0, 6))

        ttk.Label(preset_frame, text="Preset:").pack(side=LEFT, padx=(0, 4))
        self.preset_combo = ttk.Combobox(
            preset_frame,
            textvariable=self.preset_var,
            state="readonly",
            width=26,
        )
        self.preset_combo.pack(side=LEFT, padx=4)
        self._refresh_presets_combo()
        self.preset_combo.bind("<<ComboboxSelected>>", self._on_preset_selected)

        ttk.Button(preset_frame, text="Save Preset...", command=self._save_preset, bootstyle="outline-primary").pack(side=LEFT, padx=4)
        ttk.Button(preset_frame, text="Delete Preset", command=self._delete_preset, bootstyle="outline-danger").pack(side=LEFT, padx=4)

        # Settings row (Delimiters & parsing)
        parse_frame = ttk.LabelFrame(main_frame, text="File Parsing Options", padding=6)
        parse_frame.pack(fill=X, pady=(0, 6))

        # Delimiter selector
        ttk.Label(parse_frame, text="Delimiter:").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        delim_combo = ttk.Combobox(
            parse_frame,
            textvariable=self.delimiter_var,
            values=["whitespace", ",", "\\t (tab)", ";", "|", "custom"],
            state="readonly",
            width=12,
        )
        delim_combo.grid(row=0, column=1, sticky="w", padx=4, pady=2)
        delim_combo.bind("<<ComboboxSelected>>", lambda e: self._on_delim_changed())

        self.custom_delim_entry = ttk.Entry(parse_frame, textvariable=self.custom_delim_var, width=6)
        self.custom_delim_entry.grid(row=0, column=2, sticky="w", padx=4, pady=2)
        if self.delimiter_var.get() != "custom":
            self.custom_delim_entry.configure(state="disabled")

        # Skip rows
        ttk.Label(parse_frame, text="Skip Rows:").grid(row=0, column=3, sticky="w", padx=(12, 4), pady=2)
        skip_spin = ttk.Spinbox(parse_frame, from_=0, to=100, textvariable=self.skip_rows_var, width=5)
        skip_spin.grid(row=0, column=4, sticky="w", padx=4, pady=2)
        skip_spin.bind("<FocusOut>", lambda e: self._update_preview())
        skip_spin.bind("<Return>", lambda e: self._update_preview())

        # Comment char
        ttk.Label(parse_frame, text="Comment Char:").grid(row=0, column=5, sticky="w", padx=(12, 4), pady=2)
        comment_entry = ttk.Entry(parse_frame, textvariable=self.comment_char_var, width=5)
        comment_entry.grid(row=0, column=6, sticky="w", padx=4, pady=2)
        comment_entry.bind("<FocusOut>", lambda e: self._update_preview())
        comment_entry.bind("<Return>", lambda e: self._update_preview())

        # Has Header checkbox
        header_check = ttk.Checkbutton(
            parse_frame,
            text="First line is header",
            variable=self.has_header_var,
            command=self._update_preview,
        )
        header_check.grid(row=0, column=7, sticky="w", padx=(16, 4), pady=2)

        ttk.Button(parse_frame, text="Refresh Preview", command=self._update_preview, bootstyle="secondary").grid(
            row=0, column=8, sticky="e", padx=(12, 4), pady=2
        )

        # Middle: Preview Table (Treeview)
        preview_frame = ttk.LabelFrame(main_frame, text="Data Preview (first 20 rows)", padding=6)
        preview_frame.pack(fill=BOTH, expand=True, pady=(0, 6))

        tree_scroll_y = ttk.Scrollbar(preview_frame, orient="vertical")
        tree_scroll_y.pack(side=RIGHT, fill=Y)
        tree_scroll_x = ttk.Scrollbar(preview_frame, orient="horizontal")
        tree_scroll_x.pack(side=BOTTOM, fill=X)

        self.tree = ttk.Treeview(
            preview_frame,
            selectmode="none",
            yscrollcommand=tree_scroll_y.set,
            xscrollcommand=tree_scroll_x.set,
            height=8,
        )
        self.tree.pack(fill=BOTH, expand=True)
        tree_scroll_y.config(command=self.tree.yview)
        tree_scroll_x.config(command=self.tree.xview)

        # Bottom Frame: Column Mappings & Time Format
        bottom_box = ttk.Frame(main_frame)
        bottom_box.pack(fill=X, pady=(0, 6))

        # Column Mapping Frame
        col_frame = ttk.LabelFrame(bottom_box, text="Column Mappings", padding=6)
        col_frame.pack(side=LEFT, fill=BOTH, expand=True, padx=(0, 4))

        fields = [
            ("Time / JD *:", "julian_date", 0, 0),
            ("Magnitude *:", "mag", 0, 2),
            ("Uncertainty (sig):", "sig", 1, 0),
            ("Filter / Band:", "band", 1, 2),
            ("Rejected Flag:", "rejected", 2, 0),
            ("Batch ID:", "batch_number", 2, 2),
            ("Instrumental Mag:", "inst_mag", 3, 0),
            ("Instrumental Sig:", "inst_sig", 3, 2),
            ("Filename / Point ID:", "filename", 4, 0),
        ]

        self.col_combos: Dict[str, ttk.Combobox] = {}
        for label_text, key, r, c in fields:
            ttk.Label(col_frame, text=label_text).grid(row=r, column=c, sticky="w", padx=4, pady=2)
            cb = ttk.Combobox(col_frame, textvariable=self.col_vars[key], state="readonly", width=16)
            cb.grid(row=r, column=c + 1, sticky="w", padx=4, pady=2)
            self.col_combos[key] = cb

        # Time Conversion Frame
        time_frame = ttk.LabelFrame(bottom_box, text="Time System & Defaults", padding=6)
        time_frame.pack(side=RIGHT, fill=BOTH, padx=(4, 0))

        ttk.Label(time_frame, text="Time Format:").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        time_format_combo = ttk.Combobox(
            time_frame,
            textvariable=self.time_format_var,
            values=[
                "julian_date",
                "mjd",
                "iso_utc",
                "relative_days",
                "relative_hours",
                "relative_minutes",
                "relative_seconds",
            ],
            state="readonly",
            width=16,
        )
        time_format_combo.grid(row=0, column=1, sticky="w", padx=4, pady=2)

        ttk.Label(time_frame, text="Zero Epoch (JD):").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(time_frame, textvariable=self.zero_time_var, width=16).grid(row=1, column=1, sticky="w", padx=4, pady=2)

        ttk.Label(time_frame, text="Default Band:").grid(row=2, column=0, sticky="w", padx=4, pady=2)
        ttk.Combobox(
            time_frame,
            textvariable=self.default_band_var,
            values=["U", "B", "V", "R", "I", "g", "r", "i", "z", "Rc", "Ic"],
            state="normal",
            width=14,
        ).grid(row=2, column=1, sticky="w", padx=4, pady=2)

        ttk.Label(time_frame, text="Rejected Marker:").grid(row=3, column=0, sticky="w", padx=4, pady=2)
        ttk.Entry(time_frame, textvariable=self.rejected_val_var, width=16).grid(row=3, column=1, sticky="w", padx=4, pady=2)

        # Action Buttons (Cancel / OK)
        btn_frame = ttk.Frame(main_frame)
        btn_frame.pack(fill=X, pady=(4, 0))

        ttk.Button(btn_frame, text="Cancel", command=self.top.destroy, bootstyle="secondary").pack(side=RIGHT, padx=4)
        ttk.Button(btn_frame, text="Load Data", command=self._on_ok, bootstyle="success").pack(side=RIGHT, padx=4)

    def _on_delim_changed(self):
        val = self.delimiter_var.get()
        if val == "custom":
            self.custom_delim_entry.configure(state="normal")
        else:
            self.custom_delim_entry.configure(state="disabled")
        self._update_preview()

    def _get_active_delimiter(self) -> str:
        d = self.delimiter_var.get()
        if d == "\\t (tab)":
            return "\t"
        if d == "custom":
            return self.custom_delim_var.get() or "whitespace"
        return d

    def _update_preview(self):
        """Re-read file header and populate Treeview preview and mapping dropdowns."""
        delim = self._get_active_delimiter()
        sep = r"\s+" if delim == "whitespace" else delim
        skip = self.skip_rows_var.get()
        comment = self.comment_char_var.get() or None
        has_hdr = self.has_header_var.get()

        try:
            df_preview = pd.read_csv(
                self.filepath,
                sep=sep,
                skiprows=skip,
                comment=comment,
                header=0 if has_hdr else None,
                nrows=20,
                engine="python",
            )
        except Exception as e:
            # Clear treeview and show error in first item
            self.tree.delete(*self.tree.get_children())
            self.tree["columns"] = ["Error"]
            self.tree.heading("Error", text="Read Error")
            self.tree.insert("", "end", values=[f"Failed to parse preview: {e}"])
            return

        # Configure columns
        raw_cols = [str(c) for c in df_preview.columns]
        self.preview_cols = raw_cols

        # Format column identifiers for dropdowns: "[0] col_name" or "[0] Col 0"
        choices = ["[None]"]
        for idx, col in enumerate(raw_cols):
            label = f"[{idx}] {col}" if has_hdr else f"[{idx}] Column {idx}"
            choices.append(label)

        # Update Treeview columns
        self.tree.delete(*self.tree.get_children())
        self.tree["columns"] = raw_cols
        self.tree["show"] = "headings"
        for idx, col in enumerate(raw_cols):
            heading_text = f"[{idx}] {col}" if has_hdr else f"[{idx}]"
            self.tree.heading(col, text=heading_text)
            self.tree.column(col, width=100, anchor="center")

        # Insert rows
        for _, row in df_preview.iterrows():
            vals = [str(v) for v in row.values]
            self.tree.insert("", "end", values=vals)

        # Update mapping comboboxes
        for key, cb in self.col_combos.items():
            cb["values"] = choices
            # Attempt auto-selection if empty or unselected
            cur = self.col_vars[key].get()
            if not cur or cur not in choices:
                matched = self._find_best_col_match(key, raw_cols, has_hdr)
                if matched is not None:
                    self.col_vars[key].set(choices[matched + 1])
                else:
                    self.col_vars[key].set("[None]")

    def _find_best_col_match(self, key: str, cols: List[str], has_hdr: bool) -> Optional[int]:
        """Attempt to guess the best column match by index or header name."""
        # By index default if no headers
        default_indices = {
            "julian_date": 0,
            "mag": 1,
            "sig": 2,
        }
        if not has_hdr:
            idx = default_indices.get(key)
            if idx is not None and idx < len(cols):
                return idx
            return None

        # Fuzzy match by column name
        import re
        patterns = {
            "julian_date": [r"^jd$", r"julian", r"^time$", r"^date$", r"^obstime$", r"^epoch$", r"^mjd$"],
            "mag": [r"^mag$", r"magnitude", r"^mag_cal$", r"^v$", r"^r$", r"^target_mag$"],
            "sig": [r"^sig$", r"^err$", r"^error$", r"mag_err", r"magerr", r"sigma", r"magunc", r"^unc$"],
            "inst_mag": [r"inst_mag", r"mag_inst", r"instrumental"],
            "inst_sig": [r"inst_sig", r"sig_inst", r"inst_err"],
            "band": [r"^band$", r"^filter$", r"^filt$", r"passband"],
            "rejected": [r"^rejected$", r"^rej$", r"^flag$"],
            "batch_number": [r"^batch$", r"batch_number", r"batch_id"],
            "filename": [r"^filename$", r"^file$", r"^fname$", r"^image$", r"^frame$", r"point_?id"],
        }
        pats = patterns.get(key, [])
        for idx, col in enumerate(cols):
            c_low = col.strip().lower()
            if any(re.search(pat, c_low) for pat in pats):
                return idx

        # Fallback to default indices
        idx = default_indices.get(key)
        if idx is not None and idx < len(cols):
            return idx
        return None

    def _refresh_presets_combo(self):
        presets = list(self.cfg.text_import.presets.keys())
        self.preset_combo["values"] = ["[None / Last Used]"] + presets

    def _on_preset_selected(self, event=None):
        name = self.preset_var.get()
        if name in ("[None / Last Used]", ""):
            return
        preset = self.cfg.text_import.presets.get(name)
        if not preset:
            return

        self.delimiter_var.set(preset.get("delimiter", "whitespace"))
        self.skip_rows_var.set(preset.get("skip_rows", 0))
        self.comment_char_var.set(preset.get("comment_char", "#"))
        self.has_header_var.set(preset.get("has_header", True))
        self.time_format_var.set(preset.get("time_format", "julian_date"))
        self.default_band_var.set(preset.get("default_band", "V"))

        self._update_preview()

        # Apply column indices or names from preset
        choices = self.col_combos["julian_date"]["values"]
        for key in self.col_vars.keys():
            saved_val = preset.get(f"col_{key}")
            if saved_val:
                # Find matching choice
                for choice in choices:
                    if choice.startswith(f"[{saved_val}]") or choice.endswith(f" {saved_val}"):
                        self.col_vars[key].set(choice)
                        break

    def _save_preset(self):
        # Dialog to prompt for preset name
        preset_name = ttk.dialogs.dialogs.Querybox.get_string(
            parent=self.top,
            prompt="Enter a name for this import preset:",
            title="Save Import Preset",
        )
        if not preset_name or not preset_name.strip():
            return
        preset_name = preset_name.strip()

        preset_data = {
            "delimiter": self.delimiter_var.get(),
            "skip_rows": self.skip_rows_var.get(),
            "comment_char": self.comment_char_var.get(),
            "has_header": self.has_header_var.get(),
            "time_format": self.time_format_var.get(),
            "default_band": self.default_band_var.get(),
        }
        for key, var in self.col_vars.items():
            val = var.get()
            # Extract index or name
            if val and val != "[None]":
                idx_str = val.split("]")[0].replace("[", "").strip()
                preset_data[f"col_{key}"] = idx_str

        self.cfg.text_import.presets[preset_name] = preset_data
        save_user_config(self.cfg)
        self._refresh_presets_combo()
        self.preset_var.set(preset_name)
        messagebox.showinfo("Preset Saved", f"Preset '{preset_name}' saved successfully.", parent=self.top)

    def _delete_preset(self):
        name = self.preset_var.get()
        if name in ("[None / Last Used]", ""):
            return
        if messagebox.askyesno("Delete Preset", f"Delete preset '{name}'?", parent=self.top):
            if name in self.cfg.text_import.presets:
                del self.cfg.text_import.presets[name]
                save_user_config(self.cfg)
                self._refresh_presets_combo()
                self.preset_var.set("[None / Last Used]")

    def _on_ok(self):
        """Validate selections and generate the final TextImportConfig."""
        # Time and Mag are required
        time_choice = self.col_vars["julian_date"].get()
        mag_choice = self.col_vars["mag"].get()

        if not time_choice or time_choice == "[None]":
            messagebox.showerror("Missing Column", "Please map a column for Time / Julian Date.", parent=self.top)
            return
        if not mag_choice or mag_choice == "[None]":
            messagebox.showerror("Missing Column", "Please map a column for Magnitude.", parent=self.top)
            return

        col_mapping: Dict[str, Union[int, str]] = {}
        for key, var in self.col_vars.items():
            choice = var.get()
            if choice and choice != "[None]":
                # Parse 0-based index from "[N] Name"
                idx = int(choice.split("]")[0].replace("[", "").strip())
                if self.has_header_var.get() and idx < len(self.preview_cols):
                    # Use column name if header exists, otherwise int index
                    col_mapping[key] = self.preview_cols[idx]
                else:
                    col_mapping[key] = idx

        # Parse zero time if specified
        zero_val: Optional[float] = None
        zt_str = self.zero_time_var.get().strip()
        if zt_str:
            try:
                zero_val = float(zt_str)
            except ValueError:
                messagebox.showerror("Invalid Epoch", "Zero Epoch must be a valid float Julian Date.", parent=self.top)
                return

        cfg = TextImportConfig(
            delimiter=self._get_active_delimiter(),
            skip_rows=self.skip_rows_var.get(),
            comment_char=self.comment_char_var.get(),
            has_header=self.has_header_var.get(),
            col_mapping=col_mapping,
            time_format=self.time_format_var.get(),
            zero_time=zero_val,
            default_band=self.default_band_var.get() or "V",
            rejected_value=self.rejected_val_var.get().strip() or None,
        )

        # Persist as last_used in user config
        last_used = {
            "delimiter": self.delimiter_var.get(),
            "skip_rows": self.skip_rows_var.get(),
            "comment_char": self.comment_char_var.get(),
            "has_header": self.has_header_var.get(),
            "time_format": self.time_format_var.get(),
            "default_band": self.default_band_var.get(),
        }
        for key, var in self.col_vars.items():
            val = var.get()
            if val and val != "[None]":
                last_used[f"col_{key}"] = val.split("]")[0].replace("[", "").strip()

        self.cfg.text_import.last_used = last_used
        try:
            save_user_config(self.cfg)
        except Exception:
            pass

        self.result_config = cfg
        self.top.destroy()
