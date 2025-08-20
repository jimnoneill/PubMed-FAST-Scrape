import argparse
from datetime import datetime
from pubmed_fast_scrape.scraper import PubMedScraper
import os
import time

def main():
    current_year = datetime.now().year

    parser = argparse.ArgumentParser(description="PubMed Article Scraper with improved pagination and PDF download capabilities.")
    parser.add_argument("--email", help="User's email for Entrez queries", required=False, default="your.email@example.com")
    parser.add_argument("--api_key", help="PubMed API key for enhanced rate limit", default="", required=False)
    parser.add_argument("--topic", type=str, required=True, help="Topic of interest for the PubMed search (required).")
    parser.add_argument("--keywords", type=str, default=None, help="Optional additional keywords separated by commas in quotes (e.g., 'treatment,therapy,drug')")
    parser.add_argument("--start_year", type=int, default=current_year, help=f"Start year for the search range. Default is the current year ({current_year}).")
    parser.add_argument("--end_year", type=int, default=current_year, help=f"End year for the search range. Default is the current year ({current_year}).")
    parser.add_argument("--min_citations", type=int, default=0, help="Minimum number of citations for articles to include. Default is 0.")
    parser.add_argument("--download_pdfs", action="store_true", default=False, help="Enable PDF download for articles without abstracts. Default is False.")

    args = parser.parse_args()
    start = time.time()
    
    # Initialize the scraper with email and api_key if provided
    scraper = PubMedScraper(email=args.email, api_key=args.api_key)

    print(f"Starting PubMed scrape for topic: '{args.topic}'")
    if args.keywords:
        print(f"Additional keywords: {args.keywords}")
    print(f"Year range: {args.start_year}-{args.end_year}")
    print(f"PDF download: {'Enabled' if args.download_pdfs else 'Disabled'}")
    print("-" * 50)

    # Perform the scrape based on the provided arguments
    results = scraper.scrape_articles(
        field=args.topic,
        year_range=(args.start_year, args.end_year),
        keywords=args.keywords,
        min_citations=args.min_citations,
        download_pdfs=args.download_pdfs
    )

    # Create the data/ directory if it doesn't exist
    os.makedirs("data", exist_ok=True)

    # Format the filename with the search field and current datetime for uniqueness
    topic_filename = args.topic.replace(" ", "_")
    filename = f"data/{topic_filename}_articles_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"

    # Save the DataFrame to a TSV file
    results.to_csv(filename, sep="\t", index=False)
    stop = time.time()

    print("-" * 50)
    print(f"Scraping completed!")
    print(f"Found {len(results)} articles")
    print(f"Results saved to {filename}")
    print(f"Time taken: {stop-start:.2f} seconds")

if __name__ == "__main__":
    main()
