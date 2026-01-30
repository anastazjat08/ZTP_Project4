import pandas as pd
from pathlib import Path
'''
Obliczanie miesięcznych średnich wartości PM2.5 dla miasta
'''

def calculate_station_monthly_averages(df):
    """
    Oblicza miesięczne średnie wartości PM2.5 dla każdej stacji w każdym roku
    
    Args:
        df (pd.DataFrame): DataFrame z danymi PM2.5, gdzie kolumny to kody stacji, a indeks to daty.
        
    Returns:
        pd.DataFrame: DataFrame z miesięcznymi średnimi wartościami PM2.5.
    """
    df_copy = df.copy()
    months_means = (
        df_copy.groupby([df_copy["Data"].dt.year, df_copy["Data"].dt.month]).mean(numeric_only=True)
    )

    # Ustawienie nazw indeksów
    months_means.index.names = ["Rok", "Miesiąc"]

    return months_means

def calculate_city_monthly_averages(df):
    """
    Oblicza miesięczne średnie wartości PM2.5 dla każdego miasta w każdym roku
    
    Args:
        df (pd.DataFrame): DataFrame ze średnimi dla każdej stacji.
        
    Returns:
        pd.DataFrame: DataFrame z miesięcznymi średnimi wartościami PM2.5 dla miejscowości.
    """
    df_copy = df.copy()
    city_month_means = df_copy.T.groupby(level="Miejscowosc").mean().T

    return city_month_means


def main():
    # Parametry z pliku Snakefile
    input_filepath = Path(snakemake.input[0])
    output_filepath = Path(snakemake.output[0])

    df = pd.read_pickle(input_filepath)
    # Liczenie średnich miesięcznych na stacje
    month_means = calculate_station_monthly_averages(df)
    # Liczenie średnich miesięcznych na miasta
    result = calculate_city_monthly_averages(month_means)

    # Tworzenie katalogu, jeśli nie istnieje
    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_filepath, index=False)

if __name__ == "__main__":
    main()