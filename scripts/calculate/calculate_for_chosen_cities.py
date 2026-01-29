from scripts.process_data.load_data import *
import pandas as pd
from pathlib import Path

# Parametry
cities = snakemake.params['cities']
# Wczytanie listy DataFrame'ów dla wybranych lat
dfs = [pd.read_csv(Path(f)) for f in snakemake.input]
output_filepath = Path(snakemake.output[0])
output_filepath.parent.mkdir(parents=True, exist_ok=True)

if len(dfs) > 1:
    concat_df = pd.concat(dfs, axis=0, ignore_index=True)
else:
    concat_df = dfs[0]
# Filtrowanie danych dla wybranych miast
filtered_df = get_cities_years(concat_df, cities)

# Zapisywanie wyników do pliku CSV
filtered_df.to_csv(output_filepath, index=False)