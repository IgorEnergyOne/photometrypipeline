from __future__ import annotations

import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Union


@dataclass
class AtlasRecord:
    """Represents a single observation data line."""
    raw_time: float
    time: float
    raw_mag: float
    mag: float
    mag_err: float
    flag: str = ""
    extra_cols: List[Union[float, str]] = field(default_factory=list)


@dataclass
class AtlasBlock:
    """Represents an observation block containing metadata and lightcurve points."""
    headers: Dict[str, Any]
    records: List[AtlasRecord] = field(default_factory=list)

    @property
    def object_name(self) -> Optional[str]:
        return self.headers.get("OBJECT")

    @property
    def zero_time(self) -> float:
        return self.headers.get("PARSED_ZERO_TIME", 0.0)

    @property
    def zero_mag(self) -> float:
        return self.headers.get("PARSED_ZERO_MAG", 0.0)

    @property
    def unit_of_time(self) -> str:
        return self.headers.get("UNIT OF TIME", "1 day")

    def to_dict(self) -> Dict[str, List[Any]]:
        """Exports data records to a columnar dictionary."""
        return {
            "time": [r.time for r in self.records],
            "raw_time": [r.raw_time for r in self.records],
            "mag": [r.mag for r in self.records],
            "raw_mag": [r.raw_mag for r in self.records],
            "mag_err": [r.mag_err for r in self.records],
            "flag": [r.flag for r in self.records],
            "extra_cols": [r.extra_cols for r in self.records],
        }

    def to_dataframe(self):
        """Converts block records to a pandas DataFrame (requires pandas)."""
        try:
            import pandas as pd
            df = pd.DataFrame(self.to_dict())
            return df
        except ImportError:
            raise ImportError("Install pandas (`pip install pandas`) to use to_dataframe().")


class BaseAtlasParser(ABC):
    """
    Abstract base class providing structural parsing pipeline.
    Extend this class to support different Atlas file dialects.
    """

    def parse_file(self, filepath: Union[str, Path]) -> List[AtlasBlock]:
        """Parse an Atlas file from a filesystem path."""
        with open(filepath, "r", encoding="utf-8") as f:
            return self.parse_lines(f)

    def parse_string(self, text: str) -> List[AtlasBlock]:
        """Parse Atlas contents from a string."""
        return self.parse_lines(text.splitlines())

    def parse_lines(self, lines: Union[Iterator[str], List[str]]) -> List[AtlasBlock]:
        """Core template method iterating over lines and assembling blocks."""
        blocks: List[AtlasBlock] = []
        current_header: Dict[str, Any] = {}
        current_records: List[AtlasRecord] = []
        in_data_section = False

        for raw_line in lines:
            line = raw_line.strip()

            if not line:
                continue

            # Check for block separators or boundaries
            if self.is_block_delimiter(line):
                if current_header or current_records:
                    block = self.create_block(current_header, current_records)
                    blocks.append(self.post_process_block(block))
                    current_header = {}
                    current_records = []
                    in_data_section = False
                continue

            # Check if entering data section
            if self.is_data_start(line):
                in_data_section = True
                continue

            if in_data_section:
                record = self.parse_data_row(line, current_header)
                if record is not None:
                    current_records.append(record)
            else:
                self.parse_header_line(line, current_header)

        # Flush final block if file does not end with a delimiter
        if current_header or current_records:
            block = self.create_block(current_header, current_records)
            blocks.append(self.post_process_block(block))

        return blocks

    @abstractmethod
    def is_block_delimiter(self, line: str) -> bool:
        """Return True if line separates distinct blocks."""
        pass

    @abstractmethod
    def is_data_start(self, line: str) -> bool:
        """Return True if line indicates the start of the DATA table."""
        pass

    @abstractmethod
    def parse_header_line(self, line: str, header_dict: Dict[str, Any]) -> None:
        """Parse a single header line and populate header_dict."""
        pass

    @abstractmethod
    def parse_data_row(self, line: str, header_dict: Dict[str, Any]) -> Optional[AtlasRecord]:
        """Parse a single measurement row into an AtlasRecord."""
        pass

    def create_block(self, header: Dict[str, Any], records: List[AtlasRecord]) -> AtlasBlock:
        """Factory hook for creating an AtlasBlock instance."""
        return AtlasBlock(headers=header, records=records)

    def post_process_block(self, block: AtlasBlock) -> AtlasBlock:
        """Optional hook to run validation or transforms after block is assembled."""
        return block


class DefaultAtlasParser(BaseAtlasParser):
    """Default parser implementation for standard Asteroid Photometric Atlas files."""

    FLOAT_PATTERN = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
    DELIMITER_PATTERN = re.compile(r"^(?:={3,}-+={3,}|END OF OBJECT)\s*$")

    def is_block_delimiter(self, line: str) -> bool:
        return bool(self.DELIMITER_PATTERN.match(line))

    def is_data_start(self, line: str) -> bool:
        return line.upper().startswith("DATA:")

    def parse_header_line(self, line: str, header_dict: Dict[str, Any]) -> None:
        if ":" not in line:
            return

        key_part, val_part = line.split(":", 1)
        key = key_part.rstrip(". \t")
        val = val_part.strip()

        # Extract numeric zero points for quick calculations
        if key == "ZERO TIME":
            match = self.FLOAT_PATTERN.search(val)
            header_dict["PARSED_ZERO_TIME"] = float(match.group(0)) if match else 0.0
        elif key == "ZERO MAG":
            match = self.FLOAT_PATTERN.search(val)
            header_dict["PARSED_ZERO_MAG"] = float(match.group(0)) if match else 0.0

        # Handle repeated keys (e.g. INFORMATION) as lists
        if key in header_dict:
            if isinstance(header_dict[key], list):
                header_dict[key].append(val)
            else:
                header_dict[key] = [header_dict[key], val]
        else:
            header_dict[key] = val

    def parse_data_row(self, line: str, header_dict: Dict[str, Any]) -> Optional[AtlasRecord]:
        if line.startswith("#"):
            return None

        flag = ""
        # Handle observation quality / exclusion flags at row start (e.g., '!')
        if line[0] in "!*?":
            flag = line[0]
            line = line[1:].strip()

        tokens = line.split()
        if len(tokens) < 3:
            return None

        try:
            raw_time = float(tokens[0])
            raw_mag = float(tokens[1])
            mag_err = float(tokens[2])
        except ValueError:
            return None

        zero_time = header_dict.get("PARSED_ZERO_TIME", 0.0)
        zero_mag = header_dict.get("PARSED_ZERO_MAG", 0.0)

        # Baseline calibrations
        abs_time = zero_time + raw_time
        abs_mag = zero_mag + raw_mag

        extra_cols: List[Union[float, str]] = []
        for tok in tokens[3:]:
            try:
                extra_cols.append(float(tok))
            except ValueError:
                extra_cols.append(tok)

        return AtlasRecord(
            raw_time=raw_time,
            time=abs_time,
            raw_mag=raw_mag,
            mag=abs_mag,
            mag_err=mag_err,
            flag=flag,
            extra_cols=extra_cols,
        )