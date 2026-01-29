from load_data import *
from pathlib import Path
import pickle
import pandas as pd

# Parametry z pliku Snakefile
dicts_path = Path(snakemake.input[0])
with open(dicts_path, "rb") as f:
    dicts = pickle.load(f)

old_codes = dicts["old_codes"]
cities = dicts["cities"]
provinces = dicts["provinces"]

dfs = [pd.read_pickle(Path(f)) for f in snakemake.input[1:]]

output_filepath = Path(snakemake.output[0])
output_filepath.parent.mkdir(parents=True, exist_ok=True)


# Łączenie danych PM2.5 z wielu lat w jeden DataFrame
df = merge_dataframes(dfs, cities, provinces)
df.to_pickle(output_filepath)