from pathlib import Path
import pickle
import pandas as pd
import re

'''
Czyszczenie surowych dany:
 - zostawienie tylko potrzebnych wierszy
 - stare kody stacji na nowy
 - poprawienie dat
'''
def clean_pm25_data(df):
    """Czyści Dataframe z danymi PM2.5

    Args:
        df (pd.DataFrame): DataFrame z danymi PM2.5
    Returns:
        dict: słownik z oczyszczonymi DataFrame dla każdego roku
    """

    cleaned_df = df.copy()
    date_format = re.compile(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}')

    # Zostawiamy tylko wiersze z potrzebnymi danymi
    mask = (cleaned_df.iloc[:, 0].astype(str).str.match(date_format) |
            (cleaned_df.iloc[:, 0] == 'Kod stacji'))
    
    cleaned_df = cleaned_df[mask].reset_index(drop=True)

    # Ustawienie wiersza gdzie jest 'Kod stacji' jako nagłówki kolumn
    id = cleaned_df[cleaned_df.iloc[:, 0] == 'Kod stacji'].index[0]
    cleaned_df.columns = cleaned_df.loc[id].tolist()
    cleaned_df = cleaned_df.drop(index=id).reset_index(drop=True)

    # Przemianowanie kolumny z datami i zmiana na format datetime
    cleaned_df = cleaned_df.rename(columns={'Kod stacji': 'Data'})
    cleaned_df['Data'] = pd.to_datetime(cleaned_df['Data'])

    # Zamiana przecinków na kropki (jeśli plik używa przecinków jako separatora dziesiętnego, np. 2018)
    cleaned_df = cleaned_df.replace(',', '.', regex=True).infer_objects(copy=False)

    # Konwersja kolumn stacji na typ float
    for col in cleaned_df.columns:
        if col != 'Data':
            cleaned_df[col] = pd.to_numeric(cleaned_df[col], errors='coerce')

    return cleaned_df

def replace_old_codes(df, old_codes):
    """Zamienia stare kody stacji na nowe w Dataframe

    Args:
        df (pd.DataFrame): DataFrame z danymi PM2.5
        old_codes (dict): słownik mapujący stare kody stacji na nowe

    Returns:
        pd.DataFrame: DataFrame z zamienionymi kodami stacji
    """

    changed_df = df.copy()
    stations = changed_df.columns.tolist()
    changes = 0
    sample_changes = []

    for station in stations[1:]:
        if station in old_codes:
            new_code = old_codes[station]
            if len(sample_changes) < 5:
                sample_changes.append((station, new_code))
            
            stations[stations.index(station)] = new_code
            changes += 1

    changed_df.columns = stations

    return changed_df

def correct_dates(df):
    """Poprawia daty

    Args:
        df (pd.DataFrame): DataFrame z danymi PM2.5
    Returns:
        pd.DataFrame: DataFrame z poprawionymi datami
    """
    changed_df = df.copy()

    cutoff = pd.Timedelta(seconds=59)
    mask_midnight = changed_df['Data'].dt.time <= (pd.Timestamp("00:00:00") + cutoff).time()

    if mask_midnight.any():
        # bierzemy pierwszy
        idx = mask_midnight.idxmax()#bierzemy true
        example_before = changed_df.loc[idx, 'Data']
    else:
        example_before = None

    # korekta
    changed_df.loc[mask_midnight, 'Data'] = (changed_df.loc[mask_midnight, 'Data'].dt.normalize() - pd.Timedelta(seconds=1))


    return changed_df

def add_multiindex(df, cities, provinces):
    """Dodaje multiindeks

    Args:
        dfs (dict): słownik z DataFrame dla każdego roku
        cities (dict): słownik mapujący kody stacji na nazwy miejscowości
        provinces (dict): słownik mapujący kody stacji na nazwy wojewódstw

    Returns:
        pd.DataFrame: DataFrame z multindeksem
    """
    
    # Zamiana na MultiIndex
    new_columns = []
    for col in df.columns:
        if col == "Data":
            new_columns.append(("Data", "", ""))  # np. zostaw "Data" jako kolumnę dat
        else:
            miejscowosc = cities.get(col, "Nieznana")  # default jeśli brak w metadanych
            wojewodstwo = provinces.get(col, "Nieznane") # default jeśli brak w metadanych
            new_columns.append((wojewodstwo,miejscowosc,col))

    df.columns = pd.MultiIndex.from_tuples(new_columns,names=["Wojewodztwo", "Miejscowosc", "Stacja"])

    # Konwersja kolumn do odpowiednich typów
    cols_to_convert = df.columns[1:]
    df[cols_to_convert] = df[cols_to_convert].apply(pd.to_numeric, errors="coerce")

    return df


def main():
    # Parametry z pliku Snakefile
    dicts_path = Path(snakemake.input.dicts)
    with open(dicts_path, "rb") as f:
        dicts = pickle.load(f)

    old_codes = dicts["old_codes"]
    cities = dicts["cities"]
    provinces = dicts["provinces"]

    df = pd.read_pickle(snakemake.input.raw_data)
    output_filepath = Path(snakemake.output[0])

    # Czyszczenie danych PM2.5 dla danego roku
    df_cleaned = clean_pm25_data(df)
    df_new_codes = replace_old_codes(df_cleaned, old_codes)
    df_corrected_dates = correct_dates(df_new_codes)
    df_multiindex = add_multiindex(df_corrected_dates, cities, provinces)

    # Tworzenie katalogu, jeśli nie istnieje
    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    df_multiindex.to_pickle(output_filepath)

if __name__ == '__main__':
    main()  