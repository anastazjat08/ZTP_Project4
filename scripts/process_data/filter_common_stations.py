from load_data import *
from pathlib import Path
import pickle
import pandas as pd

"""
Filtrowanie danych PM2.5 do wspólnych stacji pomiarowych i dodanie multiindeksu
"""

# Parametry z pliku Snakefile
dicts_path = Path(snakemake.input[0])
with open(dicts_path, "rb") as f:
    dicts = pickle.load(f)

old_codes = dicts["old_codes"]
cities = dicts["cities"]
provinces = dicts["provinces"]

common_stations_path = Path(snakemake.input[1])
clean_data_path = Path(snakemake.input[2])
output_filepath = Path(snakemake.output[0])
output_filepath.parent.mkdir(parents=True, exist_ok=True)

# Wczytanie listy wspólnych stacji
with open(common_stations_path, "rb") as f:
    common_stations = pickle.load(f)

# Wczytanie danych czystych dla danego roku
df_clean = pd.read_pickle(clean_data_path)

# Pierwsza kolumna to daty
date_col = df_clean.columns[0]

# Zostawiamy daty + tylko kolumny, które są w common_stations
cols_to_keep = [date_col] + [c for c in df_clean.columns[1:] if c in common_stations]

df_common = df_clean[cols_to_keep]

# Dodanie multiindeksu do DataFrame
df_multiindexed = add_multiindex_to_df(df_common, cities, provinces)
# Zapisanie przefiltrowanych danych
df_multiindexed.to_pickle(output_filepath)