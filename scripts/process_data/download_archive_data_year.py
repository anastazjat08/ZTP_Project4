from pathlib import Path
import pandas as pd
import io
import requests
import zipfile

'''
Pobieranie surowych danych - pomiarów PM2.5 z określonych lat
'''
def download_gios_archive(year, gios_archive_url, gios_id, filename):
    """ Ściąganie podanego archiwum GIOS i wczytanie pliku z danymi PM2.5 do DataFrame
    Args:
        year (int): rok
        gios_archive_url (str): URL do archiwum GIOS
        gios_id (str): ID archiwum GIOS
        filename (str): nazwa pliku do pobrania

    Returns:
        pd.DataFrame: dane PM2.5 dla podanego roku
    """
    # Pobranie archiwum ZIP do pamięci
    url = f"{gios_archive_url}{gios_id}"
    response = requests.get(url)
    response.raise_for_status()  # jeśli błąd HTTP, zatrzymaj
    df = pd.DataFrame()
    
    # Otwórz zip w pamięci
    with zipfile.ZipFile(io.BytesIO(response.content)) as z:
        # znajdź właściwy plik z PM2.5
        if not filename:
            print(f"Błąd: nie znaleziono {filename}.")
        else:
            # wczytaj plik do pandas
            with z.open(filename) as f:
                try:
                    df = pd.read_excel(f, header=None)
                except Exception as e:
                    print(f"Błąd przy wczytywaniu {year}: {e}")
    return df

def main():
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
    # Zapisywanie danych do pliku (pickle zachowuje strukture)
    df.to_pickle(output_filepath)
    print(f"Zapisano dane PM2.5 dla roku {year} do {output_filepath}")

if __name__ == "__main__":
    main()