configfile: "config/task4.yaml"

YEARS = config["years"]
INTERESTED_YEARS = config["interested_years"]

# -----------------------------
# Wszystkie wymagane pliki
# -----------------------------

rule all:
    input:
        r"data\metadata.pkl",
        expand(r"data\raw\pm25_data_{year}.pkl", year=YEARS),
        r"data\dictionaries.pkl",
        expand(r"data\clean\pm25_data_clean_{year}.pkl", year=YEARS),
        expand(r"results\pm25\{year}\month_means.csv", year=YEARS),
        expand(r"results\pm25\{year}\chosen_cities_month_means_{year}.csv", year=INTERESTED_YEARS),
        expand(r"results\pm25\{year}\exceedance_days.csv", year=YEARS), 
        expand(r"results\literature\{year}\pubmed_data.csv", year=YEARS),
        expand(r"results\literature\{year}\summary_by_year.csv", year=YEARS),
        expand(r"results\literature\{year}\top_journals.csv", year=YEARS),
        expand(r"results\literature\{year}\summary_by_month.csv", year=YEARS),
        r"results\Raport.md"
        

# -----------------------------
# Reguły przetwarzania danych PM2.5
# -----------------------------

# 1. Pobranie metadanych
rule pm25_metadata:
    output:
        r"data\metadata.pkl"
    script:
        r"scripts\process_data\download_metadata.py"

# 2. Pobranie surowych danych PM2.5 dla każdego roku
rule pm25_archive_download:
    output:
        r"data\raw\pm25_data_{year}.pkl"
    params:
        year = lambda wildcards: int(wildcards.year),
        gios_archive_url = config["gios"]["archive_url"],
        gios_id = lambda wildcards: config["gios"]["ids"][int(wildcards.year)],
        filename = lambda wildcards: config["gios"]["pm25_filenames"][int(wildcards.year)]
    script:
        r"scripts\process_data\download_archive_data_year.py"

# 3. Utworzenie słowników kodów
rule pm25_dictionaries:
    input:
        r"data\metadata.pkl"
    output:
        r"data\dictionaries.pkl"
    script:
        r"scripts\process_data\create_dictionaries.py"

# 4. Czyszczenie danych PM2.5 dla każdego roku
rule pm25_clean:
    input:
        dicts = r"data\dictionaries.pkl",
        raw_data = r"data\raw\pm25_data_{year}.pkl"
    output:
        r"data\clean\pm25_data_clean_{year}.pkl"
    script:
        r"scripts\process_data\clean_data.py"


# -----------------------------
# Wyliczenia
# -----------------------------

# Liczenie średnich miesięcznych PM2.5 dla każdego miasta
rule pm25_month_means:
    input:
        r"data\clean\pm25_data_clean_{year}.pkl"
    output:
        r"results\pm25\{year}\month_means.csv"
    script:
        r"scripts\calculate\calculate_month_means.py"

# Wyciąganie średnich miesięcznych PM2.5 dla wybranych miast i lat
rule pm25_means_for_chosen_cities:
    input:
        r"results\pm25\{year}\month_means.csv"
    output:
        r"results\pm25\{year}\chosen_cities_month_means_{year}.csv"
    params:
        cities = config["cities"]
    script:
        r"scripts\calculate\calculate_for_chosen_cities.py"

# Liczenie dni z przekroczeniami PM2.5 dla każdego województwa, dla każdego roku
rule pm25_exceedance_days:
    input:
        r"data\clean\pm25_data_clean_{year}.pkl"
    output:
        r"results\pm25\{year}\exceedance_days.csv"
    script:
        r"scripts\calculate\calculate_exceedance_days.py"

# ---------------------------
# Pubmed
# ---------------------------

# Pobieranie danych z PubMed i agregacja wyników
rule pubmed_download:
    input:
        config = r"config\task4.yaml"
    output:
        pubmed_data = r"results\literature\{year}\pubmed_data.csv",
        summary_by_year = r"results\literature\{year}\summary_by_year.csv",
        top_journals = r"results\literature\{year}\top_journals.csv",
        summary_by_month = r"results\literature\{year}\summary_by_month.csv"
    params:
        year = lambda wildcards: int(wildcards.year)
    shell:
        """
        python scripts\literature\pubmed_fetch.py \
            --year {params.year} \
            --config {input.config} \
            --pubmed_data {output.pubmed_data} \
            --summary_by_year {output.summary_by_year} \
            --top_journals {output.top_journals} \
            --summary_by_month {output.summary_by_month}
        """
  
#-----------------------------
# Budowanie raportu
# ----------------------------

rule report:
    input:
        pm25_files = expand(
            r"results\pm25\{year}\exceedance_days.csv", year=YEARS),
        pubmed_files = expand(
            r"results\literature\{year}\{file}", year=YEARS,
            file=["pubmed_data.csv", "summary_by_year.csv", "top_journals.csv", "summary_by_month.csv"])
    params:
        years = YEARS
    output:
        report = r"results\Raport.md"
    script:
        r"scripts\build_report.py"