import requests
from bs4 import BeautifulSoup
import sys
import time
import pandas as pd
from Bio import Entrez
from collections import Counter
from datetime import datetime
import pgeocode
import re
import urllib
import os
import asyncio
import aiohttp
from aiohttp import ClientSession
from asyncio_throttle import Throttler
import math

class PubMedArticle:
    nomi = pgeocode.Nominatim('us')
    nomi_ca = pgeocode.Nominatim('ca')

    def __init__(self, article, api_key='', email='your.email@example.com'):
        self.api_key = api_key
        self.email = email
        self.pmid = self.get_text_or_empty(article, 'PMID')
        self.article_title = self.get_text_or_empty(article, 'ArticleTitle')
        self.abstract = self.get_text_or_empty(article, 'AbstractText')
        self.pub_date_year, self.pub_date_month = self.get_pub_date(article)
        self.journal = self.get_text_or_empty(article, 'Title')
        self.keywords = self.get_joined_text(article, 'Keyword')
        self.mesh_terms = self.get_joined_text(article, 'DescriptorName')
        self.affiliations = self.get_joined_text(article, 'Affiliation')
        self.funding_sources = self.get_funding_sources(article)
        self.location = self.get_location(article)
        self.email_addresses = self.get_email(article)
        self.authors_info = self.get_authors(article)

    @staticmethod
    def get_text_or_empty(article, tag):
        element = article.find(tag)
        return element.get_text() if element else ''

    @staticmethod
    def get_pub_date(article):
        pub_date = article.find('PubDate')
        if pub_date:
            year = pub_date.find('Year')
            month = pub_date.find('Month')
            year_text = year.text if year else ''
            month_text = month.text if month else ''
            return year_text, month_text
        return '', ''

    @staticmethod
    def get_joined_text(article, tag):
        elements = article.find_all(tag)
        return ";".join([el.get_text() for el in elements]) if elements else ''

    @staticmethod
    def get_authors(article):
        authors_info = []
        for author in article.find_all('Author'):
            last_name = author.find('LastName')
            fore_name = author.find('ForeName')
            initials = author.find('Initials')
            orcid_id = ''
            identifier = author.find('Identifier', {'Source': 'ORCID'})
            if identifier:
                orcid_id = identifier.text
            name_parts = []
            if fore_name:
                name_parts.append(fore_name.text)
            elif initials:
                name_parts.append(initials.text)
            if last_name:
                name_parts.append(last_name.text)
            name = ' '.join(name_parts)
            if name:
                authors_info.append({'Name': name.strip(), 'ORCID': orcid_id.strip()})
        return authors_info

    @staticmethod
    def get_funding_sources(article):
        grant_list = article.find_all('Grant')
        if grant_list:
            funding_sources = ";".join(list(set([grant.find('Agency').text for grant in grant_list if grant.find('Agency')])))
            return funding_sources
        return ''

    @staticmethod
    def get_location(article):
        affiliation = PubMedArticle.get_joined_text(article, 'Affiliation')
        if not affiliation:
            return ''

        countries = ['USA', 'United States', 'Canada']  # Add more countries as needed
        locations = []

        for country in countries:
            if country.lower() in affiliation.lower():
                locations.append(country)

        return ";".join(list(set(locations)))

    @staticmethod
    def get_email(article):
        affiliation = PubMedArticle.get_joined_text(article, 'Affiliation')
        if not affiliation:
            return ''

        email_regex = r'[\w\.-]+@[\w\.-]+\.\w+'
        email_matches = re.findall(email_regex, affiliation)

        return ";".join(list(set(email_matches)))

    async def download_pdf(self, save_dir, session, throttler):
        """Asynchronously download the PDF of the PMC article using the PMC ID."""
        pmc_id = await self.get_pmc_id(session, throttler)
        if not pmc_id:
            # PDF not available
            return

        pmc_id_str = f'PMC{pmc_id}'
        url = f'https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id_str}/pdf/'
        headers = {
            'User-Agent': 'Mozilla/5.0',
            'Accept': 'application/pdf',
            'Accept-Language': 'en-US,en;q=0.9',
            'Connection': 'keep-alive',
            'email': self.email
        }
        pdf_path = os.path.join(save_dir, f'{self.pmid}.pdf')
        # Check if PDF already exists
        if os.path.exists(pdf_path):
            return

        async with throttler:
            for i in range(1, 6):
                try:
                    async with session.get(url, headers=headers, allow_redirects=True, timeout=10) as response:
                        if response.status == 200:
                            content = await response.read()
                            with open(pdf_path, 'wb') as f:
                                f.write(content)
                            return
                        elif response.status == 429:
                            await asyncio.sleep(3 * i)
                        else:
                            return
                except Exception as e:
                    if i == 5:
                        return
                    await asyncio.sleep(3 * (i - 1))
            else:
                pass

    async def get_pmc_id(self, session, throttler):
        """Asynchronously get the PMC ID associated with the PMID, if available."""
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?' \
              f'dbfrom=pubmed&db=pmc&id={self.pmid}'
        if self.api_key:
            url += f'&api_key={self.api_key}'
        headers = {
            'email': self.email
        }
        async with throttler:
            for i in range(1, 6):
                try:
                    async with session.get(url, headers=headers, allow_redirects=True, timeout=10) as response:
                        if response.status == 200:
                            text = await response.text()
                            soup = BeautifulSoup(text, 'lxml-xml')
                            pmc_ids = []
                            for linksetdb in soup.find_all('LinkSetDb'):
                                linkname = linksetdb.LinkName
                                if linkname and linkname.text == 'pubmed_pmc':
                                    for link in linksetdb.find_all('Link'):
                                        pmc_id = link.Id.text
                                        pmc_ids.append(pmc_id)
                            if pmc_ids:
                                return pmc_ids[0]  # Return the first PMC ID
                            else:
                                return None
                        elif response.status == 429:
                            await asyncio.sleep(3 * i)
                        else:
                            return None
                except Exception as e:
                    if i == 5:
                        return None
                    await asyncio.sleep(3 * (i - 1))
            else:
                return None

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

    def search_mesh_terms(self, topic, top_n=5):
        """Grab top MeSH terms for a topic from top 100 articles and filter them by relevance."""
        current_date = datetime.now()
        formatted_date = current_date.strftime("%Y/%m/%d")

        handle = Entrez.esearch(db="pubmed", term=topic, retmax=100, mindate="2022/01/01", maxdate=formatted_date)
        record = Entrez.read(handle, validate=False)  # Add validate=False
        handle.close()
        pubmed_ids = record["IdList"]

        mesh_terms = []
        for pubmed_id in pubmed_ids:
            handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
            records = Entrez.read(handle, validate=False)  # Add validate=False
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

        return top_mesh_terms

    def create_refined_search_terms(self, main_topic, mesh_terms):
        """Creates refined search terms by combining the main topic with MeSH terms using AND."""
        refined_terms = []
        for term in mesh_terms:
            refined_terms.append((main_topic, term))
        return refined_terms

    def gather_pmids_time_slice(self, term, start_date, end_date):
        """
        Recursively (or iteratively) gather all PMIDs from 'start_date' to 'end_date'
        even if count > 9999, by subdividing the date range into multiple slices
        proportional to the result count.

        start_date, end_date: "YYYY/MM/DD" strings
        returns: list of PMID strings
        """

        # 1) minimal eSearch to get count
        count_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
            f"db=pubmed&term={term}&rettype=count"
            f"&mindate={start_date}&maxdate={end_date}"
            "&usehistory=n"
        )
        if self.api_key:
            count_url += f"&api_key={self.api_key}"

        resp = None
        for attempt in range(1, 6):
            try:
                resp = requests.get(count_url, timeout=30)
                if resp.status_code == 200:
                    break
                elif resp.status_code == 429:
                    wait_s = 3 * attempt
                    print(f"[time-slice count] Rate-limit, waiting {wait_s}s ...")
                    time.sleep(wait_s)
            except Exception as e:
                print(f"Error counting range {start_date}-{end_date}: {str(e)}")
                if attempt == 5:
                    break
                time.sleep(3 * attempt)

        if not resp or resp.status_code != 200:
            print(f"[time-slice] Could not get count for {start_date}-{end_date}")
            return []

        soup = BeautifulSoup(resp.content, 'lxml-xml')
        c_tag = soup.find('Count')
        if not c_tag:
            print(f"[time-slice] No count found for {start_date}-{end_date}")
            return []

        c = int(c_tag.text)
        print(f"Range {start_date}-{end_date} => count = {c}")

        if c == 0:
            return []

        # 2) if c <= 9999 => normal chunk retrieval
        if c <= 9999:
            pmid_list = []
            retmax = 9999  # Set to 9999 to match the limit
            start_offset = 0
            # retrieve in standard lumps
            while start_offset < c:
                chunk_url = (
                    f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                    f"db=pubmed&retmode=xml&term={term}"
                    f"&mindate={start_date}&maxdate={end_date}"
                    f"&retstart={start_offset}&retmax={retmax}"
                    "&usehistory=n"
                )
                if self.api_key:
                    chunk_url += f"&api_key={self.api_key}"

                pmids_this_chunk = []
                for attempt in range(1, 6):
                    try:
                        ch_resp = requests.get(chunk_url, timeout=30)
                        if ch_resp.status_code == 200:
                            s2 = BeautifulSoup(ch_resp.content, 'lxml-xml')
                            pmids_this_chunk = [x.text for x in s2.find_all('Id')]
                            break
                        elif ch_resp.status_code == 429:
                            w_s = 3 * attempt
                            print(f"[time-slice normal chunk] Rate-limit offset {start_offset}, wait {w_s}s.")
                            time.sleep(w_s)
                    except Exception as e:
                        print(f"[time-slice chunk] Error offset={start_offset}: {e}")
                        if attempt == 5:
                            break
                        time.sleep(3*attempt)

                if not pmids_this_chunk:
                    print(f"[time-slice chunk] No PMIDs found offset={start_offset}, breaking.")
                    break

                pmid_list.extend(pmids_this_chunk)
                got_count = len(pmids_this_chunk)
                print(f"[time-slice normal chunk] {start_date}-{end_date}: got {got_count} pmids at offset={start_offset}, total={len(pmid_list)}")
                if got_count < retmax:
                    break
                start_offset += got_count

            return pmid_list

        # 3) if c > 9999 => we must subdivide into intervals_needed
        # Let's do date slicing in equal intervals by day
        # so intervals_needed = ceil(c / 9999)
        intervals_needed = math.ceil(c / 9999)
        print(f"[time-slice] c={c} >9999 => subdividing into {intervals_needed} slices")

        sd = datetime.strptime(start_date, "%Y/%m/%d").date()
        ed = datetime.strptime(end_date, "%Y/%m/%d").date()
        total_days = (ed - sd).days + 1
        # each slice covers ~ total_days / intervals_needed days
        # watch rounding
        slice_days = math.ceil(total_days / intervals_needed)

        all_pmids = []
        slice_start = sd
        while slice_start <= ed:
            slice_end = slice_start + pd.Timedelta(days=slice_days - 1)
            if slice_end > ed:
                slice_end = ed
            slice_start_str = slice_start.strftime("%Y/%m/%d")
            slice_end_str = slice_end.strftime("%Y/%m/%d")

            print(f"[time-slice] Sub-slice => {slice_start_str} to {slice_end_str}")
            sub_pmids = self.gather_pmids_time_slice(term, slice_start_str, slice_end_str)
            all_pmids.extend(sub_pmids)

            # move to next slice
            slice_start = slice_end + pd.Timedelta(days=1)

        return all_pmids

    async def download_pdfs_async(self, articles, save_dir):
        # Limit to N requests per second (e.g., 10 requests per second with API key)
        max_requests_per_second = 9  # Set to 9 to stay within the limit
        throttler = Throttler(rate_limit=max_requests_per_second)
        connector = aiohttp.TCPConnector(limit=10)  # Limit the number of connections

        async with ClientSession(connector=connector) as session:
            tasks = []
            for article in articles:
                if not article.abstract.strip():
                    task = asyncio.ensure_future(article.download_pdf(save_dir, session, throttler))
                    tasks.append(task)
            await asyncio.gather(*tasks)

    def scrape_articles(self, field, year_range, search_term_string=None, keywords=None, min_citations=0, download_pdfs=False):
        if min_citations > 0:
            print(f"Warning: min_citations filtering is not yet implemented. All articles will be returned.")
        
        # Build search term
        if not search_term_string:
            # If keywords provided, use them instead of MeSH terms
            if keywords:
                # Parse keywords and create search string
                keyword_list = [k.strip() for k in keywords.split(',')]
                search_terms = [field] + keyword_list
                search_term_string = self.convert_search_terms(search_terms)
            else:
                # Use simple field search without MeSH refinement
                search_term_string = urllib.parse.quote(field)
        
        field_filename = field.replace(" ", "_")

        columns = ['PMID', 'Abstract', 'ArticleTitle', 'PubDate', 'PubMonth', 'Author', 'AuthorORCID', 'KeyWords',
                   'MeSH', 'Journals', 'Affiliations', 'FundingSources', 'Location', 'Email']
        all_articles_df = pd.DataFrame(columns=columns)

        for year in range(int(year_range[0]), int(year_range[1]) + 1):
            # Use time slice method for each year
            start_date = f"{year}/01/01"
            end_date = f"{year}/12/31"

            print(f"Getting PMIDs for year {year} using time-slice method...")
            pmids = self.gather_pmids_time_slice(search_term_string, start_date, end_date)

            if not pmids:
                print(f"No PMIDs found for year {year}")
                continue

            print(f"Retrieved {len(pmids)} PMIDs for year {year}")

            # Process PMIDs in batches
            articles = []

            # Process PMIDs in chunks of 250 (NCBI's recommended batch size)
            chunks = [pmids[i:i + 250] for i in range(0, len(pmids), 250)]
            for i, chunk in enumerate(chunks):
                if not chunk:
                    continue

                pmid_string = ','.join(chunk)
                url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' \
                      f'db=pubmed&id={pmid_string}&retmode=xml&rettype=abstract'
                if self.api_key:
                    url += f'&api_key={self.api_key}'

                for attempt in range(1, 6):
                    try:
                        print(f"Fetching details for chunk {i+1}/{len(chunks)} ({len(chunk)} PMIDs)")
                        response = requests.get(url, timeout=30)
                        if response.status_code == 200:
                            break
                        elif response.status_code == 429:
                            wait_time = 3 * attempt
                            print(f"Rate limit exceeded, waiting {wait_time}s...")
                            time.sleep(wait_time)
                        else:
                            print(f"Error status code: {response.status_code}")
                            if attempt == 5:
                                break
                            time.sleep(3 * attempt)
                    except Exception as e:
                        print(f"Error fetching chunk {i+1}: {e}")
                        if attempt == 5:
                                break
                        time.sleep(3 * attempt)

                if response.status_code != 200:
                    print(f"Failed to fetch details for chunk {i+1}, skipping...")
                    continue

                soup = BeautifulSoup(response.content, 'lxml-xml')
                pubmed_articles = soup.find_all('PubmedArticle')

                for article in pubmed_articles:
                    pubmed_article = PubMedArticle(article, api_key=self.api_key, email=self.email)
                    if pubmed_article.pmid:
                        articles.append(pubmed_article)

                # Avoid overwhelming the server
                time.sleep(1)

            print(f"Year {year}: Processing {len(articles)} articles.")

            # Download PDFs if requested
            if download_pdfs:
                # Prepare directory for PDFs
                save_dir = f'{field_filename}_pdfs/{year}/'
                os.makedirs(save_dir, exist_ok=True)

                # Download PDFs asynchronously for articles without abstracts
                asyncio.run(self.download_pdfs_async(articles, save_dir))

            # Convert articles to dictionary list for DataFrame
                articles_dict_list = [{
                    'PMID': article.pmid,
                    'Abstract': article.abstract,
                    'ArticleTitle': article.article_title,
                    'PubDate': article.pub_date_year,
                    'PubMonth': article.pub_date_month,
                    'Author': ';'.join([author['Name'] for author in article.authors_info]),
                    'AuthorORCID': ';'.join([author['ORCID'] for author in article.authors_info]),
                    'KeyWords': article.keywords,
                    'MeSH': article.mesh_terms,
                    'Journals': article.journal,
                    'Affiliations': article.affiliations,
                    'FundingSources': article.funding_sources,
                    'Location': article.location,
                'Email': article.email_addresses
                } for article in articles]

            # Add this year's articles to the master DataFrame
                articles_df = pd.DataFrame(articles_dict_list)
                all_articles_df = pd.concat([all_articles_df, articles_df], ignore_index=True)

        # Save the complete dataset
        output_file = f'{field_filename}_articles_{year_range[0]}-{year_range[1]}.tsv'
        all_articles_df.to_csv(output_file, index=False, sep="\t")
        print(f"Saved all articles to {output_file}")

        return all_articles_df
