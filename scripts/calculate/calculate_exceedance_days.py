import pandas as pd
from pathlib import Path
'''
Liczenie dni z przekroczeniem PM2.5 dla każdego województwa, dla każdego roku
'''

def calculate_daily_station_averages(df):
    """
    Oblicza dzienne średnie wartości PM2.5 dla każdej stacji w każdym roku

    Args:
        df (pd.DataFrame): DataFrame z danymi PM2.5, gdzie kolumny to kody stacji, a indeks to daty.

    Returns:
        pd.DataFrame: DataFrame z dziennymi średnimi wartościami PM2.5.
    """
    df_copy = df.copy()
    # Obliczanie średnich dziennych stężeń na stacje
    daily_means = (
        df_copy.groupby(df_copy["Data"].dt.floor("D")).mean(numeric_only=True)
    )
    return daily_means

def calculate_days_exceeding_limit_by_province(df, limit=15):
    """
    Oblicza liczbę dni w roku, kiedy średnia dzienna wartość PM2.5
    przekracza określony limit w danym województwie
    (jeśli przynajmniej jedna stacja w województwie przekroczyła limit).

    Args:
        df (pd.DataFrame): DataFrame z danymi PM2.5.
        limit (float): Limit przekroczenia PM2.5 w µg/m^3. Domyślnie 15.

    Returns:
        pd.DataFrame: DataFrame z liczbą dni przekroczeń dla każdego województwa i roku.
    """

    df_copy = df.copy()
    # Obliczanie średnich dziennych
    daily_means = calculate_daily_station_averages(df_copy)
    # Sprawdzenie przekroczeń dla każdej stacji
    exceeded = daily_means > limit
    # Sprawdzam czy w danym dniu było przekroczenie w województwie
    exceeded_by_province = (exceeded.groupby(axis=1, level="Wojewodztwo").any())
    # Zliczanie dni w poszczególnych latach
    result = (exceeded_by_province.groupby(exceeded_by_province.index.year).sum())
    return result

def main():
    # Parametry z pliku Snakefile
    df = pd.read_pickle(Path(snakemake.input[0]))
    output_filepath = Path(snakemake.output[0])

    result = calculate_days_exceeding_limit_by_province(df, limit=15)
    result = result.T.reset_index() # Zamiana na długi format - do raportu
    result.columns = ['Województwa', 'Dni z przekroczeniem średniej normy w roku']

    # Tworzenie katalogu, jeśli nie istnieje
    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_filepath, index = False) # może lepiej na pickle, zachowuje strukture

if __name__ == "__main__":
    main()