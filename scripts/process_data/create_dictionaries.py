import pandas as pd
import pickle
from load_data import *
from pathlib import Path

# Parametry z pliku Snakefile
metadata = snakemake.input[0]
output_path = Path(snakemake.output[0])
output_path.parent.mkdir(parents=True, exist_ok=True)

metadata_df = pd.read_pickle(metadata)
# Wyciąganie starych kodów i miejscowości (tworzenie słowników)
old_codes, cities, provinces = get_old_station_codes(metadata_df)

# Zapisanie trzech słowników do pliku
with open(output_path, "wb") as f:
    pickle.dump({"old_codes": old_codes, "cities": cities, "provinces": provinces}, f)