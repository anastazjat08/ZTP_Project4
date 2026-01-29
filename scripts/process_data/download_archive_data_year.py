from pathlib import Path
import pandas as pd
from load_data import *

# Parametry z pliku Snakefile
year = snakemake.params['year']
gios_archive_url = snakemake.params['gios_archive_url']
gios_id = snakemake.params['gios_id']
filename = snakemake.params['filename']
output_filepath = Path(snakemake.output[0])

# Pobranie danych archiwalnych GIOŚ dla danego roku i zapis do pliku
df = download_gios_archive(year, gios_archive_url, gios_id, filename)

# Tworzenie katalogu, jeśli nie istnieje
output_filepath.parent.mkdir(parents=True, exist_ok=True)
# Zapisywanie danych do pliku
df.to_pickle(output_filepath)
print(f"Zapisano dane PM2.5 dla roku {year} do {output_filepath}")