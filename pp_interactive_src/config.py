# -*- coding: utf-8 -*-
"""
Configuration system for pp_interactive_src.

Two configuration layers are merged at startup:

1. ``default_config.toml``  - bundled with the package, **never** modified.
2. ``~/.config/pp_interactive/config.toml``  - user overrides (created on
   first save from the App Settings dialog).

Only keys present in the user file are overridden; everything else falls
back to the defaults.

Usage::

    from .config import get_config, load_config, save_user_config, AppConfig

    cfg = get_config()          # returns cached singleton
    save_user_config(cfg)       # persist changes to the user file
"""
from __future__ import annotations

import copy
import os
import sys
import warnings
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # type: ignore[no-redef]

try:
    import tomli_w  # type: ignore
    _HAS_TOMLIW = True
except ImportError:
    _HAS_TOMLIW = False

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_DEFAULT_TOML: Path = Path(__file__).parent / "default_config.toml"


def get_user_config_path() -> Path:
    """Return the platform-appropriate user config file path."""
    if sys.platform == "win32":
        base = Path(os.environ.get("APPDATA", Path.home()))
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return base / "pp_interactive" / "config.toml"


# ---------------------------------------------------------------------------
# Dataclasses - one per TOML section
# ---------------------------------------------------------------------------

@dataclass
class AppSection:
    debug: bool = False
    theme: str = "flatly"


@dataclass
class WindowSection:
    width: int = 1400
    height: int = 720


@dataclass
class DefaultsSection:
    mode: str = "target"
    time_mode: str = "minutes"
    show_rejected: bool = True
    errorbar_type: str = "calibrated"
    show_legend: bool = True
    show_grid: bool = True
    period: float = 4.0
    period_step: float = 0.02
    phase_max: float = 1.2
    offset_step: float = 0.05
    color_legend_loc: str = "upper left"


@dataclass
class PlotSection:
    figsize: List[float] = field(default_factory=lambda: [8.5, 5.2])
    marker_style: str = "o"
    marker_size: float = 4.0
    errorbar_capsize: float = 2.0
    errorbar_capthick: float = 1.0
    errorbar_linewidth: float = 1.0
    ghost_alpha: float = 0.45
    ghost_edge_width: float = 1.2
    grid_alpha: float = 0.15
    legend_fontsize: str = "small"
    label_fontsize: int = 9
    phase_wrap_color: str = "gray"
    phase_wrap_linewidth: float = 0.8
    phase_wrap_alpha: float = 0.6


@dataclass
class ColorsSection:
    flagged: str = "orange"
    rejected: str = "red"
    flagged_use_filter: bool = False
    default_palette: List[str] = field(default_factory=lambda: [
        "red", "orange", "olive", "green", "blue",
        "purple", "brown", "pink", "gray", "cyan",
    ])
    band: Dict[str, str] = field(default_factory=lambda: {
        "U": "indigo", "B": "royalblue", "V": "limegreen",
        "R": "#7c3150", "I": "dimgray",
        "g": "#348034", "r": "#944d4d", "i": "#59327d", "z": "#753427",
        "G": "green", "BP": "blue", "RP": "red",
    })


@dataclass
class AppConfig:
    app: AppSection = field(default_factory=AppSection)
    window: WindowSection = field(default_factory=WindowSection)
    defaults: DefaultsSection = field(default_factory=DefaultsSection)
    plot: PlotSection = field(default_factory=PlotSection)
    colors: ColorsSection = field(default_factory=ColorsSection)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* into *base* (non-destructive)."""
    result = copy.deepcopy(base)
    for key, val in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = copy.deepcopy(val)
    return result


def _dict_to_config(d: dict) -> AppConfig:
    """Populate an :class:`AppConfig` from a parsed TOML dict."""
    cfg = AppConfig()
    if "app" in d:
        for k, v in d["app"].items():
            if hasattr(cfg.app, k):
                setattr(cfg.app, k, v)
    if "window" in d:
        for k, v in d["window"].items():
            if hasattr(cfg.window, k):
                setattr(cfg.window, k, v)
    if "defaults" in d:
        for k, v in d["defaults"].items():
            if hasattr(cfg.defaults, k):
                setattr(cfg.defaults, k, v)
    if "plot" in d:
        for k, v in d["plot"].items():
            if hasattr(cfg.plot, k):
                setattr(cfg.plot, k, v)
    if "colors" in d:
        colors_d = d["colors"]
        for k, v in colors_d.items():
            if k == "band":
                if isinstance(v, dict):
                    cfg.colors.band.update(v)
            elif hasattr(cfg.colors, k):
                setattr(cfg.colors, k, v)
    return cfg


def _config_to_toml_str(cfg: AppConfig) -> str:
    """Serialize *cfg* to a TOML-formatted string (fallback when tomli_w absent)."""
    d = asdict(cfg)

    def _val(v: Any) -> str:
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, str):
            return f'"{v}"'
        if isinstance(v, list):
            items = ", ".join(_val(x) for x in v)
            return f"[{items}]"
        return str(v)

    lines: List[str] = []

    def _write_section(section_name: str, section_dict: dict) -> None:
        lines.append(f"\n[{section_name}]")
        sub: Dict[str, dict] = {}
        for k, v in section_dict.items():
            if isinstance(v, dict):
                sub[k] = v
            else:
                lines.append(f"{k} = {_val(v)}")
        for sub_key, sub_val in sub.items():
            _write_section(f"{section_name}.{sub_key}", sub_val)

    for sec_key, sec_val in d.items():
        if isinstance(sec_val, dict):
            _write_section(sec_key, sec_val)

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Singleton
# ---------------------------------------------------------------------------

_config: Optional[AppConfig] = None


def load_config() -> AppConfig:
    """Load & cache the merged configuration (default + user overrides)."""
    global _config

    # 1. Read bundled defaults
    with open(_DEFAULT_TOML, "rb") as fh:
        base_dict = tomllib.load(fh)

    # 2. Deep-merge user overrides (if present)
    user_path = get_user_config_path()
    if user_path.exists():
        try:
            with open(user_path, "rb") as fh:
                user_dict = tomllib.load(fh)
            base_dict = _deep_merge(base_dict, user_dict)
        except Exception as exc:
            warnings.warn(f"pp_interactive: could not read user config {user_path}: {exc}")

    _config = _dict_to_config(base_dict)
    return _config


def get_config() -> AppConfig:
    """Return the cached :class:`AppConfig` (calls :func:`load_config` on first access)."""
    global _config
    if _config is None:
        _config = load_config()
    return _config


def save_user_config(cfg: AppConfig) -> None:
    """Persist *cfg* to the user config file and update the in-memory singleton."""
    global _config
    user_path = get_user_config_path()
    user_path.parent.mkdir(parents=True, exist_ok=True)

    d = asdict(cfg)
    if _HAS_TOMLIW:
        with open(user_path, "wb") as fh:
            tomli_w.dump(d, fh)
    else:
        with open(user_path, "w", encoding="utf-8") as fh:
            fh.write(_config_to_toml_str(cfg))

    _config = cfg

