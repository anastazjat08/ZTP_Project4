import pandas as pd
from pathlib import Path

'''
Wyciąganie średnich miesięcznych dla wybranych miast
'''

def get_cities_years(df, cities):
    """Zwraca Dataframe z danymi dla podanych miast i lat
    Args:
        df (pd.DataFrame): DataFrame z danymi PM2.5 dla wybranych lat.
        cities (list): lista nazw miast do wybrania.
    Returns:
        pd.DataFrame: DataFrame z danymi dla podanych miast i lat.
    """
    # Kolumny do zachowania: Rok, Miesiąc i wybrane miasta
    columns_to_keep = ["Rok", "Miesiąc"] + cities
    
    # Czy wszystkie istnieją w DF
    columns_to_keep = [col for col in columns_to_keep if col in df.columns]
    
    # Zostawiamy tylko wybrane kolumny
    result_df = df[columns_to_keep].copy()
    
    return result_df


def main():
    # Parametry
    cities = snakemake.params['cities']
    df = pd.read_csv(Path(snakemake.input[0]))
    output_filepath = Path(snakemake.output[0])

    # Filtrowanie danych dla wybranych miast
    filtered_df = get_cities_years(df, cities)

    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    # Zapisywanie wyników do pliku CSV
    filtered_df.to_csv(output_filepath, index=False)

if __name__ == "__main__":
    main()