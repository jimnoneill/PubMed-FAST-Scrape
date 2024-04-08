import unittest
from unittest.mock import patch
from pubmed_fast_scrape.scraper import PubMedScraper
import pandas as pd

class TestPubMedScraper(unittest.TestCase):

    def setUp(self):
        self.scraper = PubMedScraper(api_key="your_api_key", email="your_email@example.com")

    def test_convert_search_terms_single(self):
        search_terms = ['antineoplastic']
        expected = "%28antineoplastic%29"
        result = self.scraper.convert_search_terms(search_terms)
        self.assertEqual(result, expected)

    def test_convert_search_terms_combined(self):
        search_terms = ['antineoplastic', ('prophylactic', 'cancer')]
        expected = "%28antineoplastic%29+OR+%28prophylactic+AND+cancer%29"
        result = self.scraper.convert_search_terms(search_terms)
        self.assertEqual(result, expected)

    @patch('pubmed_fast_scrape.scraper.PubMedScraper.search_mesh_terms')
    def test_search_mesh_terms(self, mock_search_mesh_terms):
        mock_search_mesh_terms.return_value = ['Cancer', 'Neoplasms']
        topic = "cancer"
        result = self.scraper.search_mesh_terms(topic)
        self.assertEqual(result, ['Cancer', 'Neoplasms'])

    @patch('pubmed_fast_scrape.scraper.PubMedScraper.scrape_articles')
    def test_scrape_articles(self, mock_scrape_articles):
        mock_scrape_articles.return_value = pd.DataFrame({
            'PMID': ['123456', '789012'],
            'Abstract': ['Abstract 1', 'Abstract 2']
        })
        search_terms = ['antineoplastic']
        year_range = (2020, 2021)
        min_citations = 50
        results = self.scraper.scrape_articles(search_terms, year_range, min_citations)
        self.assertEqual(len(results), 2)
        self.assertEqual(list(results['PMID']), ['123456', '789012'])

if __name__ == '__main__':
    unittest.main()
