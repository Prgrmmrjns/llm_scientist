from Bio import Entrez
import time

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

def fetch_pubmed_details(pmid, email):
    """
    Fetch detailed information for a PubMed ID.
    """
    Entrez.email = email
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
        record = handle.read()
        handle.close()
        
        # Parse the MEDLINE format
        title = None
        abstract = None
        journal = None
        year = None
        
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
        
        return {
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "pmid": pmid
        }
    except Exception as e:
        print(f"Error fetching details for PMID {pmid}: {str(e)}")
        return None

def search_literature(system_params):
    """
    Search PubMed for relevant literature based on system parameters.
    """
    print(f"Searching PubMed for literature on: {system_params.topic}")
    
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
            paper_details = fetch_pubmed_details(pmid, system_params.email)
            if paper_details:
                studies.append(paper_details)
        
        print(f"Successfully retrieved details for {len(studies)} papers")
        return studies
        
    except Exception as e:
        print(f"Error during PubMed search: {str(e)}")
