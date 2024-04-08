import os
from datetime import datetime
import pandas as pd

def ensure_data_dir_exists():
    """
    Checks if the data directory exists in the project structure and creates it if not.
    """
    data_dir = os.path.join(os.getcwd(), 'data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    return data_dir

def save_results_to_tsv(dataframe, filename):
    """
    Saves a given pandas DataFrame to a TSV file within the data directory.

    Parameters:
    - dataframe: Pandas DataFrame containing the scraped data.
    - filename: The name of the file to save the data to.
    """
    data_dir = ensure_data_dir_exists()
    file_path = os.path.join(data_dir, filename)
    dataframe.to_csv(file_path, sep='\t', index=False)
    print(f"Data saved to {file_path}")

def current_year():
    """
    Returns the current year.
    """
    return datetime.now().year

def load_dataset(file_path):
    """
    Loads a dataset from a given TSV file path.

    Parameters:
    - file_path: Path to the TSV file.

    Returns:
    - Pandas DataFrame containing the loaded dataset.
    """
    return pd.read_csv(file_path, sep='\t')

