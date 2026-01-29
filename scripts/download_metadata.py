from pathlib import Path
from load_data import *
import pandas as pd

# Parametry z pliku Snakefile
output_file = Path(snakemake.output[0])
# Tworzenie katalogu, jeśli nie istnieje
output_file.parent.mkdir(parents=True, exist_ok=True)

# Pobieranie i zapis metadanych
metadata_df = load_metadata()
if metadata_df is None:
    raise ValueError("Nie udało się pobrać metadanych.")
else:
    metadata_df.to_pickle(output_file)