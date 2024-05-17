import requests
from bs4 import BeautifulSoup
import sys
import time
import pandas as pd
from Bio import Entrez
from collections import Counter
import urllib.parse
from datetime import datetime
from geonamescache import GeonamesCache
import pgeocode
import re

class PubMedScraper:
    def __init__(self, api_key='', email="your@email.com"):
        self.api_key = api_key
        self.email = email
        if self.email:
            from Bio import Entrez
            Entrez.email = self.email

    def convert_search_terms(self, search_terms):
        """Converts a list of search terms into search-term URL format."""
        converted_terms = []

        for term in search_terms:
            if isinstance(term, tuple):
                converted_term = "+AND+".join(term)
                converted_terms.append("(%s)" % converted_term)
            else:
                converted_terms.append("(%s)" % term)

        encoded_terms = [urllib.parse.quote(term) for term in converted_terms]
        search_term_string = "+OR+".join(encoded_terms)
        return search_term_string

    def search_mesh_terms(self, topic, top_n=10):
        """Grab top 10 MeSH terms for a topic from top 100 articles and filter them by relevance."""
        from Bio import Entrez
        from datetime import datetime
        from collections import Counter

        # Get the current date
        current_date = datetime.now()

        # Format the date in the specified format
        formatted_date = current_date.strftime("%Y/%m/%d")

        handle = Entrez.esearch(db="pubmed", term=topic, retmax=100, mindate="2022/01/01", maxdate=formatted_date)
        record = Entrez.read(handle)
        handle.close()
        pubmed_ids = record["IdList"]

        mesh_terms = []
        for pubmed_id in pubmed_ids:
            handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            if "PubmedArticle" in records:
                for article in records["PubmedArticle"]:
                    mesh_headings = article.get("MedlineCitation", {}).get("MeshHeadingList", [])
                    for mesh_heading in mesh_headings:
                        descriptor_name = mesh_heading["DescriptorName"]
                        mesh_terms.append(str(descriptor_name))

        # Count the frequency of each MeSH term
        term_counts = Counter(mesh_terms)

        # Get the top N most common MeSH terms
        most_common_terms = term_counts.most_common(top_n)

        # Extract just the terms from the most common terms
        top_mesh_terms = [term for term, count in most_common_terms]

        # Add the topic to the list of MeSH terms
        top_mesh_terms.append(topic)

        return top_mesh_terms

    def scrape_articles(self, field, year_range, min_citations, api_key='', email='your@email.com'):
        """Scrapes articles based on search terms and year range."""
        import requests
        from bs4 import BeautifulSoup
        import sys
        import time
        import pandas as pd
        from collections import Counter
        import urllib.parse
        from datetime import datetime
        from geonamescache import GeonamesCache
        import pgeocode
        import re

        search_terms = self.search_mesh_terms(field)
        field = field.replace(" ", "_")
        search_term_string = self.convert_search_terms(search_terms)
        articles = []

        for year in range(int(year_range[0]), int(year_range[1]) + 1):
            if self.api_key:
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&rettype=count&term={search_term_string}+AND+{year}[pdat]&api_key={self.api_key}"
            else:
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&rettype=count&term={search_term_string}+AND+{year}[pdat]&usehistory=y"

            for i in range(1, 10):
                try:
                    page = requests.get(url, timeout=3)
                    if page.status_code == 200:
                        break
                except Exception as e:
                    sys.stderr.write('Got error when requesting URL "' + url + '": ' + str(e) + '\n')
                    if i == 9:
                        raise e
                    time.sleep(3 * (i - 1))

            time.sleep(5)
            soup = BeautifulSoup(page.content, 'lxml')
            pmidCount = int(soup.find('count').get_text())
            numIter = pmidCount // 100000

            for idx_pmid in range(numIter + 1):
                url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=100000&retstart={idx_pmid * 100000}&term={search_term_string}+AND+{year}[pdat]&usehistory=y&api_key={api_key}'

                for i in range(1, 16):
                    try:
                        page = requests.get(url, timeout=3)
                        if page.status_code == 200:
                            break
                        elif page.status_code == 429:
                            time.sleep(3 * i)
                    except Exception as e:
                        sys.stderr.write('Got error when requesting URL "' + url + '": ' + str(e) + '\n')
                        if i == 15:
                            continue
                        time.sleep(3 * (i - 1))

                soup = BeautifulSoup(page.content, 'lxml')
                pmidVector = [s.get_text() for s in soup.find_all('id')]

                numIter_2 = len(pmidVector) // 250
                for idx_abs in range(numIter_2 + 1):
                    pmidSubvector = pmidVector[(idx_abs * 250): (idx_abs + 1) * 250]
                    pmid_packet = ','.join(pmidSubvector)
                    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid_packet}&retmode=XML&rettype=abstract&api_key={api_key}'

                    for i in range(1, 16):
                        try:
                            page = requests.get(url, timeout=3)
                            if page.status_code == 200:
                                break
                            elif page.status_code == 429:
                                time.sleep(3 * i)
                        except Exception as e:
                            sys.stderr.write('Got error when requesting URL "' + url + '": ' + str(e) + '\n')
                            if i == 15:
                                continue
                            time.sleep(3 * (i - 1))

                    soup = BeautifulSoup(page.content, 'lxml')
                    pubmedArticles = soup.find_all('pubmedarticle')

                    class PubMedArticle:
                        gc = GeonamesCache()
                        nomi = pgeocode.Nominatim('us')
                        nomi_ca = pgeocode.Nominatim('ca')

                        def __init__(self, article):
                            self.pmid = self.get_text_or_empty(article, 'pmid')
                            self.article_title = self.get_text_or_empty(article, 'articletitle')
                            self.abstract = self.get_text_or_empty(article, 'abstracttext')
                            self.pub_date_year, self.pub_date_month = self.get_pub_date(article)
                            self.journal = self.get_text_or_empty(article, 'title')
                            self.keywords = self.get_joined_text(article, 'keyword')
                            self.mesh_terms = self.get_joined_text(article, 'descriptorname')
                            self.affiliations = self.get_joined_text(article, 'affiliation')
                            self.authors = self.get_authors(article)
                            self.funding_sources = self.get_funding_sources(article)
                            self.location = self.get_location(article)
                            self.email = self.get_email(article)

                        @staticmethod
                        def get_text_or_empty(article, tag):
                            element = article.find(tag)
                            return element.get_text() if element else ''

                        @staticmethod
                        def get_pub_date(article):
                            pub_date = article.find('pubdate')
                            if pub_date:
                                year = int(pub_date.year.text) if pub_date.year else ''
                                month = pub_date.month.text if pub_date.month else ''
                                return year, month
                            return '', ''

                        @staticmethod
                        def get_joined_text(article, tag):
                            elements = article.find_all(tag)
                            return ";".join([el.get_text() for el in elements]) if elements else ''

                        @staticmethod
                        def get_authors(article):
                            aulist = article.find_all('authorlist')
                            if aulist:
                                authors = ';'.join([author.split('</forename>')[0].title() for author in str(aulist).replace('</lastname><forename>', ' ').split('<lastname>')][1:])
                                return authors
                            return ''

                        @staticmethod
                        def get_funding_sources(article):
                            grant_list = article.find_all('grant')
                            if grant_list:
                                funding_sources = ";".join(list(set([grant.agency.text for grant in grant_list if grant.agency])))
                                return funding_sources
                            return ''

                        @staticmethod
                        def get_location(article):
                            affiliation = PubMedArticle.get_joined_text(article, 'affiliation')
                            if not affiliation:
                                return ''

                            countries = {country['name']: country for country in PubMedArticle.gc.get_countries().values()}
                            locations = []

                            found_country = None
                            for country_name in countries.keys():
                                if country_name.lower() in affiliation.lower():
                                    found_country = country_name
                                    break

                            city_state_regex = re.compile(r'([A-Za-z\s]+),\s*([A-Za-z\s]+)')
                            matches = city_state_regex.findall(affiliation)

                            if matches:
                                for match in matches:
                                    city, state_province = match
                                    location_parts = [city.strip(), state_province.strip()]

                                    if found_country:
                                        location_parts.append(found_country)

                                    location_string = ', '.join(location_parts)
                                    locations.append(location_string)

                            return ";".join(list(set(locations)))

                        @staticmethod
                        def get_email(article):
                            affiliation = PubMedArticle.get_joined_text(article, 'affiliation')
                            if not affiliation:
                                return ''

                            email_regex = r'[\w\.-]+@[\w\.-]+\.\w+'
                            email_matches = re.findall(email_regex, affiliation)

                            return ";".join(list(set(email_matches)))

                    for idxa, article in enumerate(pubmedArticles):
                        pubmed_article = PubMedArticle(article)
                        if pubmed_article.pmid:
                            articles.append(pubmed_article)
                        else:
                            sys.stderr.write(f"Article missing PMID: {article}\n")

        articles_dict_list = [{
            'PMID': article.pmid,
            'Abstract': article.abstract,
            'ArticleTitle': article.article_title,
            'PubDate': article.pub_date_year,
            'PubMonth': article.pub_date_month,
            'Author': article.authors,
            'KeyWords': article.keywords,
            'MeSH': article.mesh_terms,
            'Journals': article.journal,
            'Affiliations': article.affiliations,
            'FundingSources': article.funding_sources,
            'Location': article.location,
            'Email': article.email
        } for article in articles if article.pmid]  # Ensure only articles with PMID are included

        all_articles_df = pd.DataFrame(articles_dict_list)
        all_articles_df.to_csv(f'{field}_articles.csv', index=False, sep="\t")

        return all_articles_df


# Example usage
# api_key = 'your_api_key_here'
# scraper = PubMedScraper(api_key)
# search_terms = ['antineoplastic', 'anticarcinogenic', ('prophylactic','cancer')]
#year_range = (2020, 2023)
# articles = scraper.scrape_articles(search_terms, year_range)
''' citation stuff
        #import requests
        #from requests.adapters import HTTPAdapter
        #from requests.packages.urllib3.util.retry import Retry
        #from bs4 import BeautifulSoup
        #import numpy as np

        #def requests_retry_session(retries=3, backoff_factor=0.3, status_forcelist=(500, 502, 504)):
            #session = requests.Session()
            #retry = Retry(
                #total=retries,
                #read=retries,
                #connect=retries,
                #backoff_factor=backoff_factor,
                #status_forcelist=status_forcelist,
            #)
            #adapter = HTTPAdapter(max_retries=retry)
            #session.mount('http://', adapter)
            #session.mount('https://', adapter)
            #return session
        #def get_citation_counts(pmids_batch, session):
            #pmid_str = '&id='.join([str(i) + "+" for i in pmids_batch])

            #url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&linkname=pubmed_pubmed_citedin&id={pmid_str}&tool=citation_tool&email=jnoneill@ucdavis.edu"
            #try:
                #response = session.get(url, timeout=30)
                #if response.status_code == 200:
                    #soup = BeautifulSoup(response.content, 'lxml')
                    #citation_counts = {}

                    #citations = []
                    #for id in soup.findAll('linkset'):
                        #cites = len(id.findAll('id')) - 1 if id.findAll('id') else 0
                        #citations.append(cites)

                    #citation_counts = dict(zip(pmids_batch,citations))#id_Citations.tolist()))
                    #return citation_counts
                #else:
                    #print(f"Received status code {response.status_code} for PMIDs: {pmids_batch}")
                    #return {pmid: 0 for pmid in pmids_batch}
            #except Exception as e:
                #print(f"Error requesting URL '{url}': {e}")
                #return {pmid: 0 for pmid in pmids_batch}

        ## Creating batches of 200 PMIDs each for the citation counts
        #batches = [PMID_List[i:i + 200] for i in range(0, len(PMID_List), 200)]

        ## Using ThreadPoolExecutor to parallelize the requests
        #session = requests_retry_session()
        #citation_data = {}
        #with ThreadPoolExecutor(max_workers=5) as executor:
            #future_to_pmids = {executor.submit(get_citation_counts, batch, session): batch for batch in batches}
            #for future in as_completed(future_to_pmids):
                #citation_data = {**citation_data, **future.result()}


        ## Merge the citations into the original articles DataFrame
        #for pmid, citation_count in citation_data.items():
            #all_articles_df.loc[all_articles_df['PMID'] == pmid, 'Citations'] = citation_count

        # Ensure there are no missing values in 'Citations'
'''
