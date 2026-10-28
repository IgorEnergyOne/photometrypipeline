# -*- coding: utf-8 -*-
"""
Modular light-curve data loaders for pp_interactive_src.

Provides an extensible strategy pattern (BaseLoader) and registry (LoaderRegistry)
for importing lightcurve observations from diverse file formats:
- CSVLoader: standard pipeline CSV files.
- AtlasLoader: standard ATLAS (.ATL) submission files and raw text ATLAS files.
- CustomTextLoader: arbitrary delimited ASCII tables with custom column mappings.
"""
from __future__ import annotations

import os
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from astropy.time import Time

if TYPE_CHECKING:
    from .models import FitsContext


@dataclass
class LightCurveLoadResult:
    """Encapsulates the normalized lightcurve DataFrame and associated metadata."""
    df: pd.DataFrame
    metadata: Dict[str, Any] = field(default_factory=dict)
    fits_context: Optional[FitsContext] = None
    loader_name: str = ""


@dataclass
class TextImportConfig:
    """Configuration options for parsing custom delimited text files."""
    delimiter: str = "whitespace"  # "whitespace", ",", "\t", ";", "|", or custom string
    skip_rows: int = 0
    comment_char: str = "#"
    has_header: bool = True
    # Column mapping: canonical field name -> column name (str) or 0-based column index (int)
    # Required keys: 'julian_date' (or 'time'), 'mag'
    # Optional keys: 'sig', 'inst_mag', 'inst_sig', 'band', 'rejected', 'sextractor_flags', 'batch_number', 'filename'
    col_mapping: Dict[str, Union[int, str]] = field(default_factory=dict)
    time_format: str = "julian_date"  # "julian_date", "mjd", "iso_utc", "relative_days", "relative_hours", "relative_minutes", "relative_seconds"
    zero_time: Optional[float] = None  # Base Julian Date for relative time modes
    default_band: str = "V"
    rejected_value: Optional[str] = "!"  # String indicating rejected point, or None for truthy check


class BaseLoader(ABC):
    """Abstract base class for all lightcurve file loaders."""

    name: str = "base"
    display_name: str = "Base Loader"
    filetypes: List[Tuple[str, str]] = [("All Files", "*.*")]

    @abstractmethod
    def can_load(self, filepath: str) -> bool:
        """Return True if this loader can handle the specified file."""
        pass

    @abstractmethod
    def load(self, filepath: str, **kwargs) -> LightCurveLoadResult:
        """Load the file and return a normalized LightCurveLoadResult."""
        pass


class CSVLoader(BaseLoader):
    """Loader for standard photometrypipeline and generic CSV lightcurve files."""

    name: str = "csv"
    display_name: str = "Pipeline CSV (*.csv)"
    filetypes: List[Tuple[str, str]] = [("CSV Files", "*.csv *.CSV"), ("All Files", "*.*")]

    def can_load(self, filepath: str) -> bool:
        """Check if file has .csv extension and valid CSV structure."""
        if not os.path.isfile(filepath):
            return False
        if filepath.lower().endswith(".csv"):
            return True
        # Try inspecting first line
        try:
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                first_line = f.readline().strip()
                return "," in first_line
        except Exception:
            return False

    def load(self, filepath: str, **kwargs) -> LightCurveLoadResult:
        """Read CSV and normalize common column aliases to canonical names."""
        df = pd.read_csv(filepath)

        # Standard alias normalization
        alias_map = {
            'jd': 'julian_date',
            'jdate': 'julian_date',
            'time': 'julian_date',
            'filter': 'band',
            'filt': 'band',
            'mag_err': 'sig',
            'sigma': 'sig',
            'error': 'sig',
            'batch': 'batch_number',
            'flag': 'sextractor_flags',
            'flags': 'sextractor_flags',
            'fname': 'filename',
            'file': 'filename',
            'image': 'filename',
        }
        for alias, canonical in alias_map.items():
            if alias in df.columns and canonical not in df.columns:
                df.rename(columns={alias: canonical}, inplace=True)

        has_point_filenames = False
        if 'filename' in df.columns:
            valid_names = df['filename'].dropna().astype(str).str.strip()
            if len(valid_names) > 0:
                base = os.path.basename(filepath)
                stem = Path(filepath).stem
                if valid_names.nunique() > 1:
                    has_point_filenames = True
                elif valid_names.iloc[0] not in (base, stem, "lightcurve.csv"):
                    has_point_filenames = True

        metadata = {
            "source_file": filepath,
            "format": "csv",
            "has_point_filenames": has_point_filenames,
        }
        return LightCurveLoadResult(df=df, metadata=metadata, loader_name=self.name)


class AtlasLoader(BaseLoader):
    """Loader for ATLAS photometric lightcurves.

    Handles:
    1. Official pipeline submission format (.ATL): 18-line header, DATA: block,
       rejection flags ('!'), reduc_time relative to ZERO TIME.
    2. Raw converted/text ATLAS tables (.txt): ZERO TIME header and obsTime, mag, magUnc columns.
    """

    name: str = "atlas"
    display_name: str = "ATLAS Lightcurve (*.ATL, *.txt)"
    filetypes: List[Tuple[str, str]] = [("ATLAS Files", "*.ATL *.atl *.txt *.TXT"), ("All Files", "*.*")]

    def can_load(self, filepath: str) -> bool:
        """Detect ATLAS file by extension or signature headers."""
        if not os.path.isfile(filepath):
            return False
        ext = os.path.splitext(filepath)[1].lower()
        if ext in (".atl",):
            return True

        # Sniff first 30 lines for ATLAS headers
        try:
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                head = [f.readline() for _ in range(35)]
            text = "".join(head)
            if "ZERO TIME" in text:
                return True
            if "OBJECT." in text and ("DATA:" in text or "END OF OBJECT" in text):
                return True
        except Exception:
            return False
        return False

    def load(self, filepath: str, **kwargs) -> LightCurveLoadResult:
        """Parse an ATLAS lightcurve file."""
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()

        # Check whether this is a full .ATL format (has DATA: section)
        # or a raw text ATLAS format with ZERO TIME and whitespace table
        has_data_section = any(line.strip().startswith("DATA:") for line in lines)

        if has_data_section:
            return self._load_full_atl(filepath, lines)
        else:
            return self._load_raw_atlas_txt(filepath, lines)

    def _load_full_atl(self, filepath: str, lines: List[str]) -> LightCurveLoadResult:
        """Parse a full .ATL file with 18 header rows and DATA: section."""
        headers: Dict[str, str] = {}
        data_start_idx = -1

        for idx, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith("DATA:"):
                data_start_idx = idx + 1
                break
            if ":" in line:
                parts = line.split(":", 1)
                key = parts[0].replace(".", "").strip().upper()
                val = parts[1].strip()
                headers[key] = val

        if data_start_idx == -1:
            raise ValueError(f"No DATA: section found in ATLAS file {filepath}")

        # Extract zero_time
        zero_time = 0.0
        zero_time_str = headers.get("ZERO TIME", "")
        if zero_time_str:
            match = re.search(r"(\d+\.\d+)", zero_time_str)
            if match:
                zero_time = float(match.group(1))

        # Extract band/filter
        band = "V"
        columns_header = headers.get("COLUMNS", "")
        if columns_header:
            # Columns header format is often "#{band}.-f" (e.g. "#Rc.-f" or "#R-.")
            col_match = re.search(r"#([A-Za-z0-9_]+)", columns_header)
            if col_match:
                band_candidate = col_match.group(1)
                # Strip trailing dots or dashes
                band = band_candidate.split(".")[0].split("-")[0]
        if not band or band == "V":
            if "FILTER" in headers:
                band = headers["FILTER"].strip()
            elif "PHOT. SYSTEM" in headers:
                # E.g. "Johnson-Cousins"
                pass

        # Parse data rows
        rows = []
        for line in lines[data_start_idx:]:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("===") or stripped.startswith("END OF OBJECT"):
                break

            # Check rejection flag
            is_rejected = False
            content = stripped
            if content.startswith("!"):
                is_rejected = True
                content = content[1:].strip()

            tokens = content.split()
            if not tokens:
                continue

            try:
                # Expected columns: reduc_time, mag, inst_sig, sig, sextractor_flags
                # Or variable: reduc_time, mag, sig
                reduc_time = float(tokens[0])
                mag = float(tokens[1])
                inst_sig = float(tokens[2]) if len(tokens) > 2 else 0.0
                sig = float(tokens[3]) if len(tokens) > 3 else (inst_sig if len(tokens) == 3 else 0.0)
                flags = int(tokens[4]) if len(tokens) > 4 else 0

                jd = zero_time + reduc_time
                rows.append({
                    "julian_date": jd,
                    "reduc_time": reduc_time,
                    "mag": mag,
                    "inst_sig": inst_sig,
                    "sig": sig,
                    "sextractor_flags": flags,
                    "rejected": is_rejected,
                    "band": band,
                    "filename": os.path.basename(filepath),
                })
            except (ValueError, IndexError):
                continue

        df = pd.DataFrame(rows)
        fits_context = self._build_fits_context(filepath, headers)

        metadata = {
            "source_file": filepath,
            "format": "atlas_atl",
            "headers": headers,
            "zero_time": zero_time,
            "band": band,
            "has_point_filenames": False,
        }
        return LightCurveLoadResult(
            df=df,
            metadata=metadata,
            fits_context=fits_context,
            loader_name=self.name,
        )

    def _load_raw_atlas_txt(self, filepath: str, lines: List[str]) -> LightCurveLoadResult:
        """Parse raw ATLAS text files with ZERO TIME header and whitespace-delimited table."""
        zero_time = 0.0
        data_start_idx = 0
        headers: Dict[str, str] = {}

        for idx, line in enumerate(lines):
            stripped = line.strip()
            if "ZERO TIME" in stripped:
                match = re.search(r"(\d+\.\d+)", stripped)
                if match:
                    zero_time = float(match.group(1))
                headers["ZERO TIME"] = stripped
            elif stripped.startswith("OBJECT"):
                headers["OBJECT"] = stripped.split(":", 1)[-1].strip()
            elif stripped.startswith("#"):
                data_start_idx = idx + 1
            elif stripped and not stripped.startswith("#") and data_start_idx == 0:
                # First non-header, non-comment line
                data_start_idx = idx

        # Try to infer filter band from filename (e.g. 20250922UTC_119_2025FA22_0092_Rc.txt -> 'Rc')
        band = "V"
        filename_stem = Path(filepath).stem
        band_match = re.search(r"_([A-Za-z]{1,2}(?:c)?)$", filename_stem)
        if band_match:
            band = band_match.group(1)

        rows = []
        for line in lines[data_start_idx:]:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            is_rejected = False
            content = stripped
            if content.startswith("!"):
                is_rejected = True
                content = content[1:].strip()

            tokens = content.split()
            if len(tokens) < 2:
                continue

            try:
                time_val = float(tokens[0])
                mag_val = float(tokens[1])
                sig_val = float(tokens[2]) if len(tokens) > 2 else 0.0

                # If time_val is a day fraction (< 100), add zero_time
                if zero_time > 0.0 and time_val < 2000000.0:
                    jd = zero_time + time_val
                else:
                    jd = time_val

                rows.append({
                    "julian_date": jd,
                    "mag": mag_val,
                    "sig": sig_val,
                    "inst_sig": sig_val,
                    "band": band,
                    "rejected": is_rejected,
                    "sextractor_flags": 0,
                    "filename": os.path.basename(filepath),
                })
            except (ValueError, IndexError):
                continue

        df = pd.DataFrame(rows)
        fits_context = self._build_fits_context(filepath, headers)

        metadata = {
            "source_file": filepath,
            "format": "atlas_txt",
            "headers": headers,
            "zero_time": zero_time,
            "band": band,
            "has_point_filenames": False,
        }
        return LightCurveLoadResult(
            df=df,
            metadata=metadata,
            fits_context=fits_context,
            loader_name=self.name,
        )

    def _build_fits_context(self, filepath: str, headers: Dict[str, str]) -> Optional[FitsContext]:
        """Create a FitsContext populated with metadata parsed from ATLAS headers."""
        target = headers.get("OBJECT")
        obs_site = headers.get("OBSERVING SITE", "")
        obs_code = None
        if obs_site:
            code_match = re.search(r"code\s+([A-Za-z0-9]+)", obs_site, re.IGNORECASE)
            if code_match:
                obs_code = code_match.group(1)

        if not target and not obs_code:
            # Fallback target from filename
            name_stem = Path(filepath).stem
            # Typical format: 20250922UTC_119_2025FA22_0092_Rc
            stem_parts = name_stem.split("_")
            if len(stem_parts) >= 3:
                target = stem_parts[2]
            if len(stem_parts) >= 2 and stem_parts[1].isdigit():
                obs_code = stem_parts[1]

        try:
            from .models import FitsContext
            ctx = FitsContext.__new__(FitsContext)
            ctx.filepath = ""
            ctx.header = None
            ctx.obsparam = None
            ctx.target_object = target or "Unknown"
            ctx.observatory_code = obs_code or "500"
            return ctx
        except Exception:
            return None


class CustomTextLoader(BaseLoader):
    """Loader for arbitrary delimited text files with customizable column mapping and time formats."""

    name: str = "custom_text"
    display_name: str = "Custom Delimited Text (*.txt, *.dat, *.*)"
    filetypes: List[Tuple[str, str]] = [
        ("Text / Data Files", "*.txt *.TXT *.dat *.DAT *.csv *.CSV"),
        ("All Files", "*.*"),
    ]

    def can_load(self, filepath: str) -> bool:
        """Custom text loader can attempt to parse any readable text file."""
        return os.path.isfile(filepath)

    def load(self, filepath: str, config: Optional[TextImportConfig] = None, **kwargs) -> LightCurveLoadResult:
        """Parse text file using the provided or auto-detected TextImportConfig."""
        if config is None:
            config = self.auto_detect_config(filepath)

        df_raw = self._read_raw_dataframe(filepath, config)
        df_norm = self._normalize_columns(df_raw, config, filepath)

        has_point_filenames = (
            config.col_mapping.get("filename") is not None
            and "filename" in df_norm.columns
            and df_norm["filename"].dropna().astype(str).str.strip().ne("").any()
        )

        metadata = {
            "source_file": filepath,
            "format": "custom_text",
            "import_config": config,
            "has_point_filenames": has_point_filenames,
        }
        return LightCurveLoadResult(df=df_norm, metadata=metadata, loader_name=self.name)

    def _read_raw_dataframe(self, filepath: str, config: TextImportConfig) -> pd.DataFrame:
        """Read file into a raw pandas DataFrame according to delimiter and header settings."""
        sep = r"\s+" if config.delimiter == "whitespace" else config.delimiter
        header = 0 if config.has_header else None

        return pd.read_csv(
            filepath,
            sep=sep,
            skiprows=config.skip_rows,
            comment=config.comment_char if config.comment_char else None,
            header=header,
            engine="python",
        )

    def _normalize_columns(
        self,
        df_raw: pd.DataFrame,
        config: TextImportConfig,
        filepath: str,
    ) -> pd.DataFrame:
        """Extract and convert mapped columns into canonical LightCurveData columns."""
        df_out = pd.DataFrame()

        def _get_series(col_ref: Union[int, str, None]) -> Optional[pd.Series]:
            if col_ref is None:
                return None
            if isinstance(col_ref, int):
                if 0 <= col_ref < df_raw.shape[1]:
                    return df_raw.iloc[:, col_ref]
            elif isinstance(col_ref, str):
                if col_ref in df_raw.columns:
                    return df_raw[col_ref]
                # Try finding case-insensitive match
                for col in df_raw.columns:
                    if str(col).strip().lower() == col_ref.strip().lower():
                        return df_raw[col]
            return None

        # 1. Time / Julian Date
        time_ref = config.col_mapping.get("julian_date", config.col_mapping.get("time"))
        time_s = _get_series(time_ref)
        if time_s is None:
            # Fallback to column 0 if unmapped
            time_s = df_raw.iloc[:, 0]

        df_out["julian_date"] = self._convert_time_series(time_s, config.time_format, config.zero_time)

        # 2. Magnitude (required)
        mag_ref = config.col_mapping.get("mag")
        mag_s = _get_series(mag_ref)
        if mag_s is None:
            # Fallback to column 1 if unmapped
            mag_s = df_raw.iloc[:, 1] if df_raw.shape[1] > 1 else pd.Series(0.0, index=df_raw.index)
        df_out["mag"] = pd.to_numeric(mag_s, errors="coerce")

        # 3. Uncertainty / Sigma
        sig_ref = config.col_mapping.get("sig")
        sig_s = _get_series(sig_ref)
        if sig_s is not None:
            df_out["sig"] = pd.to_numeric(sig_s, errors="coerce").fillna(0.0)
        else:
            df_out["sig"] = 0.0

        # 4. Instrumental mag / sig
        inst_mag_s = _get_series(config.col_mapping.get("inst_mag"))
        if inst_mag_s is not None:
            df_out["inst_mag"] = pd.to_numeric(inst_mag_s, errors="coerce")
        else:
            df_out["inst_mag"] = df_out["mag"]

        inst_sig_s = _get_series(config.col_mapping.get("inst_sig"))
        if inst_sig_s is not None:
            df_out["inst_sig"] = pd.to_numeric(inst_sig_s, errors="coerce").fillna(0.0)
        else:
            df_out["inst_sig"] = df_out["sig"]

        # 5. Band
        band_s = _get_series(config.col_mapping.get("band"))
        if band_s is not None:
            df_out["band"] = band_s.astype(str).str.strip()
        else:
            df_out["band"] = config.default_band

        # 6. Rejection flag
        rej_s = _get_series(config.col_mapping.get("rejected"))
        if rej_s is not None:
            if config.rejected_value:
                df_out["rejected"] = rej_s.astype(str).str.strip().str.contains(
                    re.escape(config.rejected_value), regex=True, na=False
                )
            else:
                # Truthy conversion
                df_out["rejected"] = rej_s.astype(str).str.strip().str.lower().isin(
                    ["true", "1", "yes", "t", "y", "!"]
                )
        else:
            df_out["rejected"] = False

        # 7. Sextractor flags
        flags_s = _get_series(config.col_mapping.get("sextractor_flags"))
        if flags_s is not None:
            df_out["sextractor_flags"] = pd.to_numeric(flags_s, errors="coerce").fillna(0).astype(int)
        else:
            df_out["sextractor_flags"] = 0

        # 8. Batch number
        batch_s = _get_series(config.col_mapping.get("batch_number"))
        if batch_s is not None:
            df_out["batch_number"] = batch_s
        else:
            df_out["batch_number"] = None

        # Filename / Point ID
        fname_s = _get_series(config.col_mapping.get("filename"))
        if fname_s is not None:
            df_out["filename"] = fname_s.astype(str).str.strip()
        else:
            df_out["filename"] = os.path.basename(filepath)

        return df_out

    def _convert_time_series(
        self,
        s: pd.Series,
        time_format: str,
        zero_time: Optional[float] = None,
    ) -> pd.Series:
        """Convert a pandas Series into Julian Date float values based on time_format."""
        zero = zero_time if zero_time is not None else 0.0

        if time_format in ("julian_date", "jd"):
            return pd.to_numeric(s, errors="coerce")

        if time_format in ("mjd",):
            numeric = pd.to_numeric(s, errors="coerce")
            return numeric + 2400000.5

        if time_format in ("iso_utc", "isot", "datetime"):
            # Use astropy Time or pandas to_datetime
            try:
                times = Time(s.astype(str).tolist(), scale="utc")
                return pd.Series(times.jd, index=s.index)
            except Exception:
                dt = pd.to_datetime(s, utc=True)
                # Convert datetime to Julian Date: 2440587.5 + unix_epoch_seconds / 86400
                jd = 2440587.5 + dt.astype("int64") / 1e9 / 86400.0
                return pd.Series(jd, index=s.index)

        if time_format in ("relative_days", "days"):
            numeric = pd.to_numeric(s, errors="coerce")
            return zero + numeric

        if time_format in ("relative_hours", "hours"):
            numeric = pd.to_numeric(s, errors="coerce")
            return zero + (numeric / 24.0)

        if time_format in ("relative_minutes", "minutes"):
            numeric = pd.to_numeric(s, errors="coerce")
            return zero + (numeric / 1440.0)

        if time_format in ("relative_seconds", "seconds"):
            numeric = pd.to_numeric(s, errors="coerce")
            return zero + (numeric / 86400.0)

        # Fallback
        return pd.to_numeric(s, errors="coerce")

    @classmethod
    def auto_detect_config(cls, filepath: str) -> TextImportConfig:
        """Inspect file header and columns to propose a TextImportConfig with fuzzy matching."""
        cfg = TextImportConfig()

        try:
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                sample_lines = [f.readline() for _ in range(25)]
        except Exception:
            return cfg

        # Find skip rows and delimiter
        header_idx = -1
        delimiter = "whitespace"

        for idx, line in enumerate(sample_lines):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                # Potential commented header (e.g. #obsTime mag magUnc)
                tokens = stripped[1:].strip().split()
                if len(tokens) >= 2 and any(t.lower() in ("obstime", "time", "jd", "mag") for t in tokens):
                    cfg.skip_rows = idx + 1
                    cfg.has_header = False  # we will use column indices
                    header_tokens = tokens
                    header_idx = idx
                    break
                continue

            # First non-comment row: test delimiters
            comma_count = stripped.count(",")
            tab_count = stripped.count("\t")
            semi_count = stripped.count(";")

            if comma_count >= 1:
                delimiter = ","
            elif tab_count >= 1:
                delimiter = "\t"
            elif semi_count >= 1:
                delimiter = ";"
            else:
                delimiter = "whitespace"

            header_idx = idx
            break

        cfg.delimiter = delimiter

        # Try reading preview to identify column names
        try:
            preview_df = pd.read_csv(
                filepath,
                sep=r"\s+" if delimiter == "whitespace" else delimiter,
                skiprows=cfg.skip_rows,
                comment="#",
                nrows=5,
                engine="python",
            )
            col_names = [str(c).strip().lower() for c in preview_df.columns]
        except Exception:
            col_names = []

        mapping: Dict[str, Union[int, str]] = {}
        fuzzy_patterns = {
            "julian_date": [r"^jd$", r"julian", r"^time$", r"^date$", r"^obstime$", r"^epoch$", r"^mjd$"],
            "mag": [r"^mag$", r"magnitude", r"^mag_cal$", r"^v$", r"^r$", r"^target_mag$"],
            "sig": [r"^sig$", r"^err$", r"^error$", r"mag_err", r"magerr", r"sigma", r"magunc", r"^unc$"],
            "inst_mag": [r"inst_mag", r"mag_inst", r"instrumental"],
            "inst_sig": [r"inst_sig", r"sig_inst", r"inst_err"],
            "band": [r"^band$", r"^filter$", r"^filt$", r"passband"],
            "rejected": [r"^rejected$", r"^rej$", r"^flag$"],
            "batch_number": [r"^batch$", r"batch_number", r"batch_id"],
        }

        for canonical, patterns in fuzzy_patterns.items():
            for idx, col in enumerate(col_names):
                if any(re.search(pat, col) for pat in patterns):
                    mapping[canonical] = preview_df.columns[idx]
                    break

        cfg.col_mapping = mapping
        return cfg


class LoaderRegistry:
    """Registry managing available lightcurve file loaders and format dispatch."""

    def __init__(self):
        self._loaders: Dict[str, BaseLoader] = {}
        self.register(AtlasLoader())
        self.register(CSVLoader())
        self.register(CustomTextLoader())

    def register(self, loader: BaseLoader) -> None:
        """Register a new loader instance."""
        self._loaders[loader.name] = loader

    def get(self, name: str) -> Optional[BaseLoader]:
        """Get a registered loader by unique name."""
        return self._loaders.get(name)

    def all_loaders(self) -> List[BaseLoader]:
        """Return list of all registered loaders."""
        return list(self._loaders.values())

    def detect_loader(self, filepath: str) -> BaseLoader:
        """Automatically choose the most appropriate loader for a file.

        Priority order:
        1. AtlasLoader (handles .ATL and ZERO TIME headers)
        2. CSVLoader (handles standard CSV)
        3. CustomTextLoader (fallback for general text files)
        """
        atlas = self.get("atlas")
        if atlas and atlas.can_load(filepath):
            return atlas

        csv = self.get("csv")
        if csv and csv.can_load(filepath):
            return csv

        custom = self.get("custom_text")
        if custom:
            return custom

        raise ValueError(f"No suitable loader registered for file: {filepath}")

    def load(
        self,
        filepath: str,
        loader_name: Optional[str] = None,
        **kwargs,
    ) -> LightCurveLoadResult:
        """Load a file using either an explicit loader or auto-detection."""
        if loader_name:
            loader = self.get(loader_name)
            if not loader:
                raise ValueError(f"Unknown loader name '{loader_name}'. Registered: {list(self._loaders.keys())}")
        else:
            loader = self.detect_loader(filepath)

        return loader.load(filepath, **kwargs)


# Global default registry instance
loaders_registry = LoaderRegistry()
