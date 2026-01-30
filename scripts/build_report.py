from pathlib import Path
import pandas as pd

'''
Budowanie raportu dla zebranych danych PM2.5 i PubMed
'''


def build_report(years, pm25_files, pubmed_files, output_file: Path):

    output_file.parent.mkdir(parents=True, exist_ok=True)

    lines = ["## Raport - Analiza danych PM2.5 i PubMed\n"]

    for year, pm25_file, pubmed_list in zip(years, pm25_files, pubmed_files):
        lines.append(f"### Rok {year}\n")

        # Wczytanie danych PM2.5
        exceedance_df = pd.read_csv(pm25_file)
        lines.append("#### Liczba dni przekraczających normę PM2.5 w danych województwach\n")
        lines.append(exceedance_df.to_markdown(index=False))

        # Wczytanie danych PubMed
        pubmed_data = pd.read_csv(Path(pubmed_list[0])) # Zachowana kolejność plików
        summary_by_year = pd.read_csv(Path(pubmed_list[1]))
        top_journals = pd.read_csv(Path(pubmed_list[2]))
        summary_by_month = pd.read_csv(Path(pubmed_list[3]))
        
        lines.append("\n#### Dane PubMed\n")
        # 1. Liczba artykułów
        lines.append(f"Liczba znalezionych artykułów w PubMed: {len(pubmed_data)}\n") 
        # 2. Rozkład publikacji w latach
        lines.append("Podsumowanie artykułów według lat:\n")
        lines.append(summary_by_year.to_markdown(index=False))
        # 3. Rozkład miesięczny publikacji
        lines.append('\nRozkład miesięczny publikacji w danym roku:\n')
        lines.append(summary_by_month.to_markdown(index=False))
        # 4. Czasopisma z największą liczbą publikacji w interesującym roku
        lines.append("\nCzasopisma, które opublikowały najwięcej artykułów:\n")
        lines.append(top_journals.to_markdown(index=False))
        # 5. Przykłądowe tytuły
        lines.append("\nPrzykładowe tytuły artykułów:\n")
        for title in pubmed_data['Tytuł'].head(5):
            lines.append(f"- {title}")
        lines.append("\n")

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))

def main():
    # Parsowanie argumentów i wczytywanie konfiguracji
    years = snakemake.params.years
    pm25_files = [Path(file) for file in snakemake.input.pm25_files]
    pubmed_files = snakemake.input.pubmed_files
    pubmed_files_grouped = [pubmed_files[i:i+4] for i in range(0, len(pubmed_files), 4)]
    report = Path(snakemake.output.report)

    build_report(years, pm25_files, pubmed_files_grouped, report)
    

if __name__ == "__main__":
    main()
