import argparse
import yaml
import pandas as pd
from Bio import Entrez
from pathlib import Path
import ssl

'''
Skrypt do pobierania artykułów z PubMed na podstawie roku publikacji i zapytania wyszukiwania.
Metadane artykułów są zapisywane do plików CSV.
'''

ssl._create_default_https_context = ssl._create_unverified_context


def parse_args():
    '''
    Parsuje argumenty wiersza poleceń
    '''
    parser = argparse.ArgumentParser(description="Przeszukaj PubMed i pobierz artykuły z danego roku")
    parser.add_argument("--year", type=int, required=True)
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--pubmed_data", type=str, required=True)
    parser.add_argument("--summary_by_year", type=str, required=True)
    parser.add_argument("--top_journals", type=str, required=True)
    parser.add_argument("--summary_by_month", type=str, required=True)
    return parser.parse_args()

def load_config(config_path):
    '''
    Wczytuje plik konfiguracyjny YAML
    
    Args:
        config_path (str): Ścieżka do pliku konfiguracyjnego YAML
    Returns:
        dict: Słownik z danymi konfiguracyjnymi
    '''
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
    
def fetch_pubmed_articles(year, query, email, max_results):
    '''
    Pobiera identyfikatory artykułów z PubMed dla danego roku i zapytania
    
    Args:
        year (int): Rok publikacji artykułów
        query (str): Zapytanie wyszukiwania
        email (str): Adres email użytkownika
        max_results (int): Maksymalna liczba wyników do pobrania
    Returns:
        list: Lista identyfikatorów artykułów PubMed
    '''
    Entrez.email = email
    handle = Entrez.esearch(
        db="pubmed",
        term=f'{query} AND ("{year}/01/01"[Date - Publication] : "{year}/12/31"[Date - Publication])',
        retmax=max_results)
    
    results = Entrez.read(handle)
    handle.close()

    return results["IdList"]

def parse_date(pub_date):
    '''
    Parsuje datę publikacji z formatu PubMed na słownik z kluczami Year, Month, Day
    
    Args:
        pub_date (dict): Słownik z datą publikacji w formacie PubMed
    Returns:
        dict: Słownik z kluczami Year, Month, Day
    '''
    if not pub_date:
        return {"Year": None, "Month": None, "Day": None}

    year = pub_date.get("Year")
    month = pub_date.get("Month")
    day = pub_date.get("Day")

    # Miesiąc jako liczba
    # month_map = {
    #     "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4,
    #     "May": 5, "Jun": 6, "Jul": 7, "Aug": 8,
    #     "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12
    # }
    # month_map = {
    #     "Jan": 'Sty', "Feb": 'Lut', "Mar": 'Mar', "Apr": 'Kwi',
    #     "May": 'Maj', "Jun": 'Czer', "Jul": 'Lip', "Aug": 'Sier',
    #     "Sep": , "Oct": 10, "Nov": 11, "Dec": 12
    # }
    # month_num = month_map.get(month) if month else None

    return {
        "Year": int(year) if year else None,
        "Month": month if month else None,
        "Day": int(day) if day else None
    }

def fetch_pubmed_metadata(id_list, email):
    '''
    Pobiera metadane artykułów z PubMed na podstawie listy identyfikatorów
    
    Args:
        id_list (list): Lista identyfikatorów artykułów PubMed
        email (str): Adres email użytkownika
    Returns:
        pd.DataFrame: DataFrame z metadanymi artykułów
    '''
    if not id_list:
        return pd.DataFrame()  # Zwraca pustą listę, jeśli nie ma identyfikatorów
    else:
        # Pobieranie metadanych
        Entrez.email = email
        ids = ",".join(id_list)
        handle = Entrez.efetch(
            db="pubmed",
            id=ids,
            rettype="full", retmode="xml")
        records = Entrez.read(handle)

        metadata = []
        # Przetwarzanie każdego artykułu, wyciąganie metadanych
        for article in records.get("PubmedArticle", []):
            medline = article.get("MedlineCitation", {})
            article_data = medline.get("Article", {})

            pmid = medline.get("PMID", "")

            authors_list = article_data.get("AuthorList", [])
            article_title = article_data.get("ArticleTitle", "")

            journal_title = article_data.get("Journal", {}).get("Title", "")
            pub_date = article_data.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})

            authors_str = ", ".join(
                [
                    f"{author.get('LastName', '')} {author.get('Initials', '')}"
                    for author in authors_list
                    if 'LastName' in author and 'Initials' in author
                ]
            )
            pub_date = parse_date(pub_date)

            metadata.append({
                "PMID": pmid,
                "Tytuł": article_title,
                "Autorzy": authors_str,
                "Czasopismo": journal_title,
                "Rok": pub_date["Year"],
                "Miesiąc": pub_date["Month"],
                "Dzień": pub_date["Day"]
            })
        handle.close()
        df = pd.DataFrame(metadata)
        for col in ["Rok","Dzień"]:
            df[col] = df[col].astype("Int64")

        return df
    
def save_metadata_to_csv(df, output_path):
    '''
    Zapisuje metadane artykułów do pliku CSV
    
    Args:
        df (pd.DataFrame): DataFrame z metadanymi artykułów
        output_path (str): Ścieżka do pliku wyjściowego CSV
    '''
    df.to_csv(output_path, index=False)

def aggregation(df, year):
    
    # Podsumowanie lat publikacji
    summary_by_year = (
        df.groupby('Rok')
        .size()
        .reset_index(name='Liczba artykułów')
        .sort_values(by='Rok')
    )

    # Podsumowanie czasopism
    top_journals = (
        df.groupby('Czasopismo')
        .size()
        .reset_index(name='Liczba artykułów')
        .sort_values(by='Liczba artykułów', ascending=False)
        .head(10)
    )

    summary_by_month = df[(df['Rok'] == year) & (df['Miesiąc'].notna())]
    # Grupowanie po miesiącach tylko w danym roku
    summary_by_month = (
    summary_by_month.groupby('Miesiąc')
               .size()
               .reset_index(name='Liczba artykułów')
               .sort_values('Miesiąc')
    )


    return summary_by_year, top_journals, summary_by_month


def main():
    # Parsowanie argumentów i wczytywanie konfiguracji
    args = parse_args()
    config = load_config(args.config)
    query = config['pubmed']['query']
    email = config['pubmed']['email']
    max_results = config['pubmed'].get('max_results', 500)
    pubmed_data_path = Path(args.pubmed_data)
    summary_by_year_path = Path(args.summary_by_year)
    top_journals_path = Path(args.top_journals)
    summary_by_month_path = Path(args.summary_by_month)

    # Pobieranie identyfikatorów artykułów z PubMed
    print(f"Pobieranie artykułów z PubMed dla roku {args.year} i zapytania '{query}'...")
    id_list = fetch_pubmed_articles(args.year, query, email, max_results)
    print(f"Znaleziono {len(id_list)} artykułów.")

    # Pobieranie metadanych artykułów
    print("Pobieranie metadanych artykułów...")
    metadata_df = fetch_pubmed_metadata(id_list, email)


    # Zapisywanie metadanych do pliku CSV
    pubmed_data_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Zapisywanie metadanych do pliku {pubmed_data_path}...")
    save_metadata_to_csv(metadata_df, pubmed_data_path)
    print("Gotowe!")

    # Agregacja i podsumowanie danych
    print("Tworzenie podsumowań...")
    summary_by_year, top_journals, summary_by_month = aggregation(metadata_df, args.year)

    summary_by_year.to_csv(summary_by_year_path, index=False)
    print("Zapisano podsumowanie lat publikacji.")
    
    top_journals.to_csv(top_journals_path, index=False)
    print("Zapisano podsumowanie najpopularniejszych czasopism.")
    
    summary_by_month.to_csv(summary_by_month_path, index=False)
    print("Zapisano podsumowanie miesięczne publikacji.")



if __name__ == "__main__":
    main()
    