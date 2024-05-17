# PubMed-FAST-Scrape

PubMed-FAST-Scrape is a Python package intended for Biomedical NLP scientists. It is designed for fast and efficient scraping of PubMed article metadata by hybridizing packet-based lookups with Bio.Entrez and PubMed-direct parsing for bulk scraping of articles x200 faster than otherwise. It enables researchers and data scientists to easily gather articles based on a field of interest, year range, and minimum number of citations.

## Features

- Fast scraping of PubMed abstracts & article metadata.
- Filter articles by field of interest, year range, and minimum citations.
- Easy integration into data analysis pipelines.

## Installation

Install PubMed-FAST-Scrape using pip:

```bash
pip install git+https://github.com/jimnoneill/pubmed-fast-scrape.git
```

## Usage
PubMed-FAST-Scrape can be used as a command-line tool or imported into your Python scripts.

# Command Line Interface
To use PubMed-FAST-Scrape from the command line:

```bash
pubmed-fast-scrape --field "Cancer Research" --start_year 2010 --end_year 2020
```
# In Python Scripts
```python
from pubmed_fast_scrape.scraper import PubMedScraper

Note: currently API keys are not supported
scraper = PubMedScraper(email='email_not_required@makesitfaster.com')
results = scraper.scrape('Cancer Research', (2023, 2024), 1) # field, year range, n-min citations

results.head()

````

## Development
To contribute to PubMed-FAST-Scrape, clone the repository and create a new branch for your feature or bug fix.

```bash
git clone https://github.com/yourusername/pubmed-fast-scrape.git
cd pubmed-fast-scrape
git checkout -b your-feature-branch

```

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
- The PubMed API for providing access to their invaluable database of articles.
- BioPython and BeautifulSoup for making data extraction easier.


## Disclaimer
This tool is intended for academic and research purposes. Please ensure you adhere to PubMed's terms of use when using this scraper.

For more information, issues, or questions about the PubMed-FAST-Scrape, please visit the GitHub repository.
