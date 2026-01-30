## Projekt 4 - Analiza danych PM2.5 i PubMed
Ten pipeline łączy analizę stężeń PM2.5 w wybranych latach z automatycznym przeglądem literatury naukowej (PubMed) dla wskazanych lat. Wyniki z pomiarów powietrza i metadane artykułów są agregowane w raport końcowy. Pipeline jest **incremental** - Snakemake dba o to, aby nie przeliczać ponownie danych dla lat, które już zostały przetworzone, co pozwala na szybkie uruchamianie i łatwe uzupełnianie wyników o nowe lata.

### Struktura katalogów
```
config/     - pliki konfiguracyjne
data/       - surowe pliki pobrane z GIOŚ oraz przetworzone i ujednolicone
scripts/    - skrypty do PM2.5 i PubMed
results/    - wyniki dla kalkulacji PM2.5 oraz przeszukiwania literatury PubMed
```


### Uruchomienie
W terminalu:
```
snakemake --j <liczba_wątków>
```

### Scenariusze uruchomienia
#### Pierwsze wywołanie
1. Ustawienie w pliku config `years: [2021, 2024]`
2. Uruchomienie snakemake
3. Pipeline liczy:
    - PM2.5 dla 2021 i 2024
    - PubMed dla 2021 i 2024
    - raport dla {2021, 2024}
#### Drugie wywołanie
1. Zmiana w pliku config na `years: [2019, 2024]`
2. Ponowne uruchomienie snakemake
3. Pipeline liczy tylko brakujące kroki:
    - policzyć PM2.5 dla 2019 (2024 ma zostać pominięty),
    - pobrać/analitykę PubMed dla 2019 (2024 ma zostać pominięty),
    - wygenerować nowy raport dla {2019, 2024}.

### Weryfikacja poprawności uruchomienia
Jeżeli wszystkie pliki już istnieją i nie zmieniły się inputy (config lub skrypt), Snakemake wyświetli w logu komunikat `Nothing to be done`.  
Jeśli brakuje któregoś pliku, w logu zostaną wylistowane odpowiednie `jobs` wraz z podaną przyczyną ich wykonania (`reason`) - czyli bralujących plików, które będą liczone.
