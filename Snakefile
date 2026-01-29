configfile: "config/task4.yaml"

YEARS = config["years"]

# -----------------------------
# Wszystkie wymagane pliki
# -----------------------------

rule all:
    input:
        r"data\metadata.pkl",
        expand(r"data\raw\pm25_data_{year}.pkl", year=YEARS),
        r"data\dictionaries.pkl",
        expand(r"data\clean\pm25_data_clean_{year}.pkl", year=YEARS),
        r"data\pm25_data_merged.pkl"

# -----------------------------
# Reguły przetwarzania danych
#   1. Pobranie metadanych
#   2. Pobranie surowych danych PM2.5 dla każdego roku
#   3. Utworzenie słowników kodów
#   4. Czyszczenie danych PM2.5 dla każdego roku
#   5. Scalanie danych PM2.5 ze wszystkich lat
# -----------------------------

rule pm25_metadata:
    output:
        r"data\metadata.pkl"
    script:
        r"scripts\download_metadata.py"


rule pm25_archive_download:
    output:
        r"data\raw\pm25_data_{year}.pkl"
    params:
        year = lambda wildcards: int(wildcards.year),
        gios_archive_url = config["gios"]["archive_url"],
        gios_id = lambda wildcards: config["gios"]["ids"][int(wildcards.year)],
        filename = lambda wildcards: config["gios"]["pm25_filenames"][int(wildcards.year)]
    script:
        r"scripts\download_archive_data_year.py"
        
rule pm25_dictionaries:
    input:
        r"data\metadata.pkl"
    output:
        r"data\dictionaries.pkl"
    script:
        r"scripts\create_dictionaries.py"

rule pm25_clean:
    input:
        r"data\dictionaries.pkl",
        r"data\raw\pm25_data_{year}.pkl"
    output:
        r"data\clean\pm25_data_clean_{year}.pkl"
    script:
        r"scripts\clean_data.py"

rule pm25_merge:
    input:
        r"data\dictionaries.pkl",
        expand(r"data\clean\pm25_data_clean_{year}.pkl", year=YEARS)
    output:
        r"data\pm25_data_merged.pkl"
    script:
        r"scripts\merge_data.py"