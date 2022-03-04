from pathlib import Path

def list_dirfiles(dirpath: Path):
    files = list(Path('.').glob('*.'))
    return files