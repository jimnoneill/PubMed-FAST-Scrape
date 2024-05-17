import requests
from bs4 import BeautifulSoup
import sys
import time
import pandas as pd
from Bio import Entrez
from collections import Counter
from datetime import datetime
from geonamescache import GeonamesCache
import pgeocode
import re
import urllib

class PubMedScraper:

    def __init__(self, api_key='', email="your@email.com"):
        self.api_key = api_key
        self.email = email
        if self.email:
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
        current_date = datetime.now()
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

        term_counts = Counter(mesh_terms)
        most_common_terms = term_counts.most_common(top_n)
        top_mesh_terms = [term for term, count in most_common_terms]
        top_mesh_terms.append(topic)

        return top_mesh_terms

    def scrape_articles(self, field, year_range, min_citations, api_key='', email='your@email.com'):
        import requests
        from bs4 import BeautifulSoup
        import sys
        import time
        import pandas as pd
        from Bio import Entrez
        from collections import Counter
        from datetime import datetime
        from geonamescache import GeonamesCache
        import pgeocode
        import re
        import urllib

        search_terms = self.search_mesh_terms(field)
        search_term_string = self.convert_search_terms(search_terms)
        field = field.replace(" ", "_")

        columns = ['PMID', 'Abstract', 'ArticleTitle', 'PubDate', 'PubMonth', 'Author', 'KeyWords', 'MeSH', 'Journals', 'Affiliations', 'FundingSources', 'Location', 'Email']
        all_articles_df = pd.DataFrame(columns=columns)

        for year in range(int(year_range[0]), int(year_range[1]) + 1):
            articles = []

            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&rettype=count&term={search_term_string}+AND+{year}[pdat]&usehistory=y'#&api_key={self.api_key}'
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
                url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=100000&retstart={idx_pmid * 100000}&term={search_term_string}+AND+{year}[pdat]&usehistory=y'#&api_key={self.api_key}'
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
                    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid_packet}&retmode=XML&rettype=abstract'#&api_key={self.api_key}'

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
                } for article in articles]

                articles_df = pd.DataFrame(articles_dict_list)
                all_articles_df = pd.concat([all_articles_df, articles_df], ignore_index=True)

        #def get_citation_count(pmid):
            #url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&linkname=pubmed_pubmed_citedin&id={pmid}"
            #for i in range(1, 16):
                #try:
                    #response = requests.get(url, timeout=3)
                    #if response.status_code == 200:
                        #break
                    #elif response.status_code == 429:
                        #time.sleep(3 * i)
                #except Exception as e:
                    #sys.stderr.write('Got error when requesting URL "' + url + '": ' + str(e) + '\n')
                    #if i == 15:
                        #continue
                    #time.sleep(3 * (i - 1))

            #soup = BeautifulSoup(response.content, 'xml')
            #link_set_dbs = soup.find_all('LinkSetDb')
            #if link_set_dbs:
                #citation_ids = link_set_dbs[0].find_all('Id')
                #return len(citation_ids)
            #else:
                #return 0

        #all_articles_df['Citations'] = all_articles_df['PMID'].apply(get_citation_count)
        all_articles_df.to_csv(f'{field}_articles_{year_range[0]}-{year_range[1]}.csv', index=False, sep="\t")
        return all_articles_df
