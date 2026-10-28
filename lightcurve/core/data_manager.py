import pandas as pd
from typing import Optional, Dict, Any

class DataManager:
    """Manages data loading and manipulation for lightcurve visualization."""
    
    def __init__(self):
        self.df: Optional[pd.DataFrame] = None
        self.filename: Optional[str] = None
    
    def load(self, filepath: str) -> None:
        """Load data from a CSV file."""
        self.filename = filepath
        self.df = pd.read_csv(filepath)
        
        # Ensure required columns exist
        if 'rejected' not in self.df.columns:
            self.df['rejected'] = False
        if 'sextractor_flags' not in self.df.columns:
            self.df['sextractor_flags'] = 0
    
    def save(self, filepath: Optional[str] = None) -> None:
        """Save data to a CSV file."""
        if self.df is not None:
            self.df.to_csv(filepath or self.filename, index=False)
    
    def toggle_rejection(self, index: int) -> None:
        """Toggle the rejection status of a data point."""
        if self.df is not None and 0 <= index < len(self.df):
            self.df.at[index, 'rejected'] = not self.df.at[index, 'rejected']
    
    def get_target_name(self) -> str:
        """Get the target name from the data."""
        if self.df is not None and not self.df.empty:
            try:
                return str(self.df['target'].iloc[0])
            except KeyError:
                return 'Lightcurve'
        return 'Lightcurve'
