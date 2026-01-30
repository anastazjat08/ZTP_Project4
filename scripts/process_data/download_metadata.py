from pathlib import Path
import pandas as pd
import requests
from bs4 import BeautifulSoup
from io import BytesIO

'''
Pobieranie metadanych z GIOŚ
'''
def load_metadata():
    """ Wyszukuje najnowszy plik metadanych GIOS na stronie archiwum,
        pobiera go i zwraca jako DataFrame.

    Returns:
        pd.DataFrame: dane metadanych GIOS
    """
    
    archive_url = "https://powietrze.gios.gov.pl/pjp/archives"

    try:
        r = requests.get(archive_url)
        r.raise_for_status()
    except Exception as e:
        print(f"Błąd pobierania strony archiwum: {e}")
        return None

    soup = BeautifulSoup(r.text, "html.parser")

    # Linki z 'downloadFile/...'
    links = soup.find_all("a", href=True)
    candidates = []

    for a in links:
        href = a["href"]
        text = a.get_text(strip=True).lower()

        # warunek: tekst zawiera metadane itp.
        if "meta" in text and "downloadFile" in href:
            candidates.append((text, href))

    if not candidates:
        print("Nie znaleziono pliku metadanych!")
        return None

    text, href = candidates[0]

    file_url = "https://powietrze.gios.gov.pl" + href

    try:
        r = requests.get(file_url)
        r.raise_for_status()
    except Exception as e:
        print(f"Błąd pobierania pliku metadanych: {e}")
        return None

    try:
        df = pd.read_excel(BytesIO(r.content), header=0)
        df = df.rename(columns={'Stary Kod stacji \n(o ile inny od aktualnego)': 'Stary Kod stacji'})
    except Exception as e:
        print(f"Błąd odczytu pliku metadanych: {e}")
        return None
    return df

def main():
    # Parametry z pliku Snakefile
    output_file = Path(snakemake.output[0])
    # Tworzenie katalogu, jeśli nie istnieje
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Pobieranie i zapis metadanych
    metadata_df = load_metadata()
    if metadata_df is None:
        raise ValueError("Nie udało się pobrać metadanych.")
    else:
        metadata_df.to_pickle(output_file)

if __name__ == "__main__":
    main()