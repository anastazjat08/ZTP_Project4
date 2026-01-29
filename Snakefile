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
        r"data\pm25_data_merged.pkl",
        r"data\common_stations.pkl",
        expand(r"data\clean\pm25_data_common_{year}.pkl", year=YEARS),
        expand(r"results\pm25\{year}\month_means.csv", year=YEARS),
        r"results\pm25\chosen_cities_month_means.csv",
        expand(r"results\pm25\{year}\exceedance_days.csv", year=YEARS), 
        expand(r"results\literature\{year}\pubmed_data.csv", year=YEARS)
        

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
        r"data\dictionaries.pkl",
        r"data\raw\pm25_data_{year}.pkl"
    output:
        r"data\clean\pm25_data_clean_{year}.pkl"
    script:
        r"scripts\process_data\clean_data.py"

# 5. Znalezienie wspólnych stacji pomiarowych we wszystkich latach
rule pm25_find_common:
    input:
        expand(r"data\clean\pm25_data_clean_{year}.pkl", year=YEARS)
    output:
        r"data\common_stations.pkl"
    script:
        r"scripts\process_data\find_common_stations.py"

# 6. Filtrowanie danych PM2.5 do wspólnych stacji pomiarowych
rule pm25_filter_common:
    input:
        r"data\dictionaries.pkl",
        r"data\common_stations.pkl",
        r"data\clean\pm25_data_clean_{year}.pkl"
    output:
        r"data\clean\pm25_data_common_{year}.pkl"
    script:
        r"scripts\process_data\filter_common_stations.py"

#   5. Scalanie danych PM2.5 ze wszystkich lat
# rule pm25_merge:
#     input:
#         r"data\dictionaries.pkl",
#         expand(r"data\clean\pm25_data_clean_{year}.pkl", year=YEARS)
#     output:
#         r"data\pm25_data_merged.pkl"
#     script:
#         r"scripts\merge_data.py"

# -----------------------------
# Wyliczenia
# -----------------------------

# Liczenie średnich miesięcznych PM2.5 dla każdego miasta
rule pm25_month_means:
    input:
        r"data\clean\pm25_data_common_{year}.pkl"
    output:
        r"results\pm25\{year}\month_means.csv"
    script:
        r"scripts\calculate\calculate_month_means.py"

# Wyciąganie średnich miesięcznych PM2.5 dla wybranych miast i lat
rule pm25_means_for_chosen_cities:
    input:
        expand(r"results\pm25\{year}\month_means.csv", year=INTERESTED_YEARS)
    output:
        r"results\pm25\chosen_cities_month_means.csv"
    params:
        cities = config["cities"]
    script:
        r"scripts\calculate\calculate_for_chosen_cities.py"

# Liczenie dni z przekroczeniami PM2.5 dla każdego miasta, dla każdego roku
rule pm25_exceedance_days:
    input:
        r"data\clean\pm25_data_common_{year}.pkl"
    output:
        r"results\pm25\{year}\exceedance_days.csv"
    script:
        r"scripts\calculate\calculate_exceedance_days.py"

# ---------------------------
# Pubmed
# ---------------------------

rule pubmed_download:
    input:
        config = r"config\task4.yaml"
    output:
        r"results\literature\{year}\pubmed_data.csv"
    params:
        year = lambda wildcards: int(wildcards.year)
    shell:
        """
        python scripts\literature\pubmed_fetch.py \
            --year {params.year} \
            --config {input.config}
        """
  