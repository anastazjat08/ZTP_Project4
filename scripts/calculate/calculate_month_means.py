from calculations import *
import pandas as pd
from pathlib import Path

# Calculating on merged data

# Parametry z pliku Snakefile
input_filepath = Path(snakemake.input[0])
output_filepath = Path(snakemake.output[0])
# Tworzenie katalogu, jeśli nie istnieje
output_filepath.parent.mkdir(parents=True, exist_ok=True)

# Obliczanie miesięcznych średnich wartości PM2.5 dla każdej stacji
df = pd.read_pickle(input_filepath)
# Liczenie średnich miesięcznych na stacje
month_means = calculate_station_monthly_averages(df)
# Liczenie średnich miesięcznych na miasta
result = calculate_city_monthly_averages(month_means)
result["Rok"] = int(snakemake.wildcards.year)
result.to_csv(output_filepath)