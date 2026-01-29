from calculations import *
import pandas as pd
from pathlib import Path

# Parametry z pliku Snakefile
df = pd.read_pickle(Path(snakemake.input[0]))
output_filepath = Path(snakemake.output[0])
# Tworzenie katalogu, jeśli nie istnieje
output_filepath.parent.mkdir(parents=True, exist_ok=True)

result = calculate_days_exceeding_limit(df, limit=15)
result.to_csv(output_filepath) # może lepiej na pickle, zachowuje strukture