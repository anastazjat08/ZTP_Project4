from calculations import *
import pandas as pd
from pathlib import Path

# Parametry z pliku Snakefile
df = pd.read_pickle(Path(snakemake.input[0]))
output_filepath = Path(snakemake.output[0])
# Tworzenie katalogu, jeśli nie istnieje
output_filepath.parent.mkdir(parents=True, exist_ok=True)

result = calculate_days_exceeding_limit_by_province(df, limit=15)
result = result.T.reset_index()
result.columns = ['Województwa', 'Dni z przekroczeniem średniej normy w roku']
result.to_csv(output_filepath, index = False) # może lepiej na pickle, zachowuje strukture