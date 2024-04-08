from setuptools import setup, find_packages

setup(
    name='pubmed-fast-scrape',
    version='0.1.0',
    author='Jamey ONeill',
    author_email='joneilliii@sdsu.edu',
    description='A fast scraper for PubMed articles based on field of interest, year range, and minimum citations.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/pubmed-fast-scrape',  # Change to your GitHub repo URL
    packages=find_packages(),
    install_requires=[
        'requests',
        'biopython',
        'beautifulsoup4',  # The package name for bs4 is beautifulsoup4
        'pandas',
        'lxml',  # If you're parsing XML in BeautifulSoup, lxml is highly recommended
        # You might also need 'pgeocode' based on your code snippets
        # Include any other dependency not part of the standard library
    ],
    entry_points={
        'console_scripts': [
            'pubmed-fast-scrape=pubmed_fast_scrape.pubmed_scraper_cli:main',
        ]
    },
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.org/classifiers/
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',  # Choose the appropriate license
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.7',
    # Add any additional URLs that are relevant to your project
    project_urls={
        'Source': 'https://github.com/yourusername/pubmed-fast-scrape',
        # 'Documentation': 'URL for documentation (if available)'
    },
)
