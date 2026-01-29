import argparse
import yaml
import os
import pandas as pd
from Bio import Entrez
from pathlib import Path
import ssl

ssl._create_default_https_context = ssl._create_unverified_context


def parse_args():
    '''
    Parsuje argumenty wiersza poleceń
    '''
    parser = argparse.ArgumentParser(description="Przeszukaj PubMed i pobierz artykuły z danego roku")
    parser.add_argument("--year", type=int, required=True, help="Rok opublikowania artykułów do pobrania")
    parser.add_argument("--config", type=str, required=True, help="Ścieżka do pliku konfiguracyjnego YAML z danymi użytkownika")
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
        term=f"{query} AND {year}[PDAT]",
        retmax=max_results)
    
    results = Entrez.read(handle)
    handle.close()

    return results["IdList"]

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

            authors_str = ",".join(
                [
                    f"{author.get('LastName', '')} {author.get('Initials', '')}"
                    for author in authors_list
                    if 'LastName' in author and 'Initials' in author
                ]
            )

            metadata.append({
                "PMID": pmid,
                "Title": article_title,
                "Authors": authors_str,
                "Journal": journal_title,
                "PublicationDate": pub_date
            })
        handle.close()
        df = pd.DataFrame(metadata)

        return df
    
def save_metadata_to_csv(df, output_path):
    '''
    Zapisuje metadane artykułów do pliku CSV
    
    Args:
        df (pd.DataFrame): DataFrame z metadanymi artykułów
        output_path (str): Ścieżka do pliku wyjściowego CSV
    '''
    df.to_csv(output_path, index=False)

def summary(df, year):
    pass



def main():
    # Parsowanie argumentów i wczytywanie konfiguracji
    args = parse_args()
    config = load_config(args.config)
    query = config['pubmed']['query']
    email = config['pubmed']['email']
    max_results = config['pubmed'].get('max_results', 500)

    # Pobieranie identyfikatorów artykułów z PubMed
    print(f"Pobieranie artykułów z PubMed dla roku {args.year} i zapytania '{query}'...")
    id_list = fetch_pubmed_articles(args.year, query, email, max_results)
    print(f"Znaleziono {len(id_list)} artykułów.")

    # Pobieranie metadanych artykułów
    print("Pobieranie metadanych artykułów...")
    metadata_df = fetch_pubmed_metadata(id_list, email)


    # Zapisywanie metadanych do pliku CSV
    output_path = Path(f"results/literature/{args.year}/pubmed_data.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Zapisywanie metadanych do pliku {output_path}...")
    save_metadata_to_csv(metadata_df, output_path)
    print("Gotowe!")


if __name__ == "__main__":
    main()
    