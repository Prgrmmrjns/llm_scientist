from Bio import Entrez
import time
from typing import List, Optional
from models import StudyMetadata
from study_storage import StudyStorage
import logging

def build_pubmed_query(topic, journal=None, date_range=None):
    """
    Build a PubMed query string from parameters.
    """
    # Base query with topic
    query_parts = [topic]
    
    # Optional journal-specific search
    if journal:
        if isinstance(journal, list):
            journal_query = ' OR '.join(f'"{j}"[Journal]' for j in journal)
            query_parts.append(f"({journal_query})")
        else:
            query_parts.append(f'"{journal}"[Journal]')
    
    # Optional date range
    if date_range:
        start_date, end_date = date_range
        query_parts.append(f"({start_date}:{end_date}[Date - Publication])")

    return " AND ".join(query_parts)

def fetch_pubmed_details(pmid: str, email: str) -> Optional[StudyMetadata]:
    """
    Fetch detailed information for a PubMed ID and return as StudyMetadata.
    """
    Entrez.email = email
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
        record = handle.read()
        handle.close()
        
        # Initialize metadata fields
        title = None
        abstract = None
        journal = None
        year = None
        authors = []
        doi = None
        
        # Parse the MEDLINE format
        current_field = None
        for line in record.split('\n'):
            if line.startswith('TI  - '):
                title = line[6:].strip()
            elif line.startswith('AB  - '):
                abstract = line[6:].strip()
            elif line.startswith('JT  - '):
                journal = line[6:].strip()
            elif line.startswith('DP  - '):
                # Extract year from date
                year_str = line[6:].strip()[:4]
                try:
                    year = int(year_str)
                except ValueError:
                    year = None
            elif line.startswith('AU  - '):
                authors.append(line[6:].strip())
            elif line.startswith('LID - '):
                # Look for DOI in article ID
                if '[doi]' in line:
                    doi = line[6:].strip().replace(' [doi]', '')
        
        # Create URL (this is a placeholder - actual URL would depend on the journal)
        url = f"https://doi.org/{doi}" if doi else None
        
        # Create and return StudyMetadata object
        return StudyMetadata(
            pmid=pmid,
            title=title or "No title available",
            authors=authors,
            journal=journal,
            year=year,
            doi=doi,
            abstract=abstract,
            url=url,
            download_status="not_attempted"
        )
        
    except Exception as e:
        logging.error(f"Error fetching details for PMID {pmid}: {str(e)}")
        return None

def search_literature(system_params) -> List[StudyMetadata]:
    """
    Search PubMed for relevant literature and store results.
    """
    print(f"Searching PubMed for literature on: {system_params.topic}")
    
    # Initialize storage
    storage = StudyStorage()
    
    # Build the query
    query = build_pubmed_query(
        system_params.topic,
        system_params.journals,
        system_params.date_range
    )
    
    # Set up Entrez
    Entrez.email = system_params.email
    
    try:
        # Search PubMed
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=system_params.max_results
        )
        record = Entrez.read(handle)
        handle.close()
        
        pmids = record["IdList"]
        print(f"Found {len(pmids)} papers matching the search criteria")
        
        # Fetch details for each paper
        studies = []
        for pmid in pmids:
            # Add delay to respect NCBI's rate limits
            time.sleep(0.34)  # ~3 requests per second
            
            # Check if we already have this study
            existing_study = storage.get_study(pmid)
            if existing_study:
                studies.append(existing_study)
                continue
            
            # Fetch new study details
            study_metadata = fetch_pubmed_details(pmid, system_params.email)
            if study_metadata:
                # Store the study
                storage.add_study(study_metadata)
                studies.append(study_metadata)
        
        print(f"Successfully retrieved details for {len(studies)} papers")
        
        # Try to download PDFs
        print("Attempting to download available PDFs...")
        storage.bulk_download_pdfs()
        
        return studies
        
    except Exception as e:
        logging.error(f"Error during PubMed search: {str(e)}")
        return []
