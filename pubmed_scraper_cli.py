import argparse
from datetime import datetime
from pubmed_fast_scrape.scraper import PubMedScraper
import pandas
import os
import requests
import time

def main():
    current_year = datetime.now().year

    parser = argparse.ArgumentParser(description="PubMed Article Scraper based on field of interest, year range, and minimum citations.")
    parser.add_argument("--email", help="User's email for Entrez queries", required=False)
    parser.add_argument("--api_key", help="PubMed API key for enhanced rate limit", default=False, required=False)
    parser.add_argument("--field", type=str, default="cancer", help="Field of interest for the PubMed search. Default is 'cancer'.")
    parser.add_argument("--start_year", type=int, default=current_year, help=f"Start year for the search range. Default is the current year ({current_year}).")
    parser.add_argument("--end_year", type=int, default=current_year, help=f"End year for the search range. Default is the current year ({current_year}).")
    parser.add_argument("--min_citations", type=int, default=1, help="Minimum number of citations for articles to include. Default is 1.")

    args = parser.parse_args()
    start = time.time()
    # Initialize the scraper with email and api_key if provided
    scraper = PubMedScraper(email=args.email, api_key=args.api_key)

    # Perform the scrape based on the provided arguments
    results = scraper.scrape_articles(args.field, (args.start_year, args.end_year), args.min_citations)

    # Create the data/ directory if it doesn't exist
    os.makedirs("data", exist_ok=True)

    # Format the filename with the search field and current datetime for uniqueness
    filename = f"data/{args.field}_articles_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"

    # Save the DataFrame to a TSV file
    results.to_csv(filename, sep="\t", index=False)
    stop = time.time()

    print(f"Results saved to {filename}. Took {start-stop} seconds")

if __name__ == "__main__":
    main()
