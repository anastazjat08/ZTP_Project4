import sys
from pathlib import Path

# Dodanie folder scripts do PYTHONPATH
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))