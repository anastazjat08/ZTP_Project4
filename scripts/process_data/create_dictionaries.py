import pandas as pd
import pickle
from pathlib import Path

'''
Tworzenie słowników z metadanych do ujednolicenia surowych danych PM2.5
'''
def get_old_station_codes(metadata_df):
    """ Wyciąga stare kody stacji z metadanych

    Args:
        metadata_df (pd.DataFrame): dane metadanych GIOS

    Returns:
        tuple: krotka składająca sie z 3 słowników:
            1) słownik mapujący stare kody stacji na nowe
            2) słownik mapujący kody stacji na nazwy miejscowości
            3) słownik mapujący kody stacji na nazwy wojewódstw

    """
    metadata_filtered = metadata_df[metadata_df["Stary Kod stacji"].notna()]
    old_codes = {}
    for k, row in metadata_filtered.iterrows():
        old = row['Stary Kod stacji']
        new = row['Kod stacji']
        if isinstance(old, str):
            for code in old.split(','): # w jednej komórce metadanych może być kilka starych kodów rozdzielonych przecinkiem
                old_codes[code.strip()] = new
    cities = dict(zip(metadata_df["Kod stacji"], metadata_df["Miejscowość"]))
    provinces = dict(zip(metadata_df["Kod stacji"], metadata_df["Województwo"]))
    return old_codes, cities, provinces


def main():
    # Parametry z pliku Snakefile
    metadata = snakemake.input[0]
    output_path = Path(snakemake.output[0])

    metadata_df = pd.read_pickle(metadata)
    # Wyciąganie starych kodów i miejscowości (tworzenie słowników)
    old_codes, cities, provinces = get_old_station_codes(metadata_df)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Zapisanie trzech słowników do pliku
    with open(output_path, "wb") as f:
        pickle.dump({"old_codes": old_codes, "cities": cities, "provinces": provinces}, f)

if __name__ == "__main__":
    main()