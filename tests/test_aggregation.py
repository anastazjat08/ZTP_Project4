import pandas as pd
from scripts.literature.pubmed_fetch import aggregation


def test_aggregation():

    df = pd.DataFrame({
        "PMID": ["1", "2", "3", "4"],
        "Tytuł": ["Art1", "Art2", "Art3", "Art4"],
        "Autorzy": ["Nowak A", "Nowak B", "Nowak C", "Nowak D"],
        "Rok": [2023, 2023, 2022, 2023],
        "Miesiąc": ["Jan", "Feb", "Jan", "Jan"],
        "Dzień": [1, 2, 3, 4],
        "Czasopismo": ["Nature", "Nature", "Science", "Cell"]
    })

    summary_by_year, top_journals, summary_by_month = aggregation(df, 2023)
    
    # Sprawdzanie liczby publikacji w 2023
    count_2023 = summary_by_year.loc[
        summary_by_year['Rok'] == 2023, 'Liczba artykułów'
    ].iloc[0]

    assert count_2023 == 3

    # Sprawdzanie top czasopism
    assert top_journals.iloc[0]["Czasopismo"] == "Nature"
    assert top_journals.iloc[0]["Liczba artykułów"] == 2

    # Sprawdzenie agregacji miesięcznych dla 2023
    jan_count = summary_by_month.loc[
        summary_by_month["Miesiąc"] == "Jan", "Liczba artykułów"
    ].iloc[0]

    assert jan_count == 2

def test_aggregation_none():
    df = pd.DataFrame(columns = ["PMID", "Tytuł", "Autorzy", "Rok", "Miesiąc", "Dzień", "Czasopismo"])

    summary_by_year, top_journals, summary_by_month = aggregation(df, 2023)

    assert summary_by_year.empty
    assert top_journals.empty
    assert summary_by_month.empty
