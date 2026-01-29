from load_data import *
from pathlib import Path
import pickle

# Parametry z pliku Snakefile
dicts_path = Path(snakemake.input[0])
with open(dicts_path, "rb") as f:
    dicts = pickle.load(f)

old_codes = dicts["old_codes"]
cities = dicts["cities"]
provinces = dicts["provinces"]

df = pd.read_pickle(snakemake.input[1])

output_filepath = Path(snakemake.output[0])



# Czyszczenie danych PM2.5 dla danego roku
df_cleaned = clean_pm25_data(df)
df_new_codes = replace_old_codes(df_cleaned, old_codes)
df_corrected_dates = correct_dates(df_new_codes)

# Tworzenie katalogu, jeśli nie istnieje
output_filepath.parent.mkdir(parents=True, exist_ok=True)
df_corrected_dates.to_pickle(output_filepath)