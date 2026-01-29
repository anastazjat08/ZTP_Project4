from load_data import *
import pandas as pd
import pickle
from pathlib import Path


# Parametry z pliku Snakefile
dfs = dfs = [pd.read_pickle(Path(f)) for f in snakemake.input]
output_filepath = Path(snakemake.output[0])
output_filepath.parent.mkdir(parents=True, exist_ok=True)

# Znalezienie wspólnych stacji we wszystkich DataFrame
common_stations_list = common_stations(dfs)
# Zapisanie listy wspólnych stacji do pliku
with open(output_filepath, "wb") as f:
    pickle.dump(common_stations_list, f)

print(len(common_stations_list), "common stations found.")