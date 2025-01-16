import os
import json
import requests
from datetime import datetime
from pathlib import Path
from typing import List, Dict
from models import StudyMetadata
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('study_downloads.log'),
        logging.StreamHandler()
    ]
)

class StudyStorage:
    def __init__(self, base_dir: str = "systematic_review_data"):
        """Initialize storage directories."""
        self.base_dir = Path(base_dir)
        self.articles_dir = self.base_dir / "articles"
        self.metadata_file = self.base_dir / "study_metadata.json"
        self.text_dir = self.base_dir / "text_content"  # New directory for text content
        
        # Create directories if they don't exist
        self.articles_dir.mkdir(parents=True, exist_ok=True)
        self.text_dir.mkdir(parents=True, exist_ok=True)  # Create text directory
        
        # Initialize or load metadata
        self.studies: Dict[str, StudyMetadata] = {}
        if self.metadata_file.exists():
            self._load_metadata()
    
    def _load_metadata(self):
        """Load study metadata from JSON file."""
        try:
            with open(self.metadata_file, 'r') as f:
                data = json.load(f)
                self.studies = {
                    pmid: StudyMetadata(**study_data)
                    for pmid, study_data in data.items()
                }
        except Exception as e:
            logging.error(f"Error loading metadata: {str(e)}")
            self.studies = {}
    
    def _save_metadata(self):
        """Save study metadata to JSON file."""
        try:
            with open(self.metadata_file, 'w') as f:
                json.dump(
                    {pmid: study.dict() for pmid, study in self.studies.items()},
                    f,
                    indent=2,
                    default=str  # Handle datetime serialization
                )
        except Exception as e:
            logging.error(f"Error saving metadata: {str(e)}")
    
    def store_study_text(self, study: StudyMetadata):
        """Store study content as text file."""
        try:
            text_path = self.text_dir / f"{study.pmid}.txt"
            
            with open(text_path, 'w', encoding='utf-8') as f:
                # Write basic metadata
                f.write(f"Title: {study.title}\n")
                f.write(f"Authors: {', '.join(study.authors)}\n")
                f.write(f"Journal: {study.journal or 'N/A'}\n")
                f.write(f"Year: {study.year or 'N/A'}\n")
                f.write(f"DOI: {study.doi or 'N/A'}\n")
                f.write(f"PMID: {study.pmid}\n")
                f.write("\n")
                
                # Write abstract
                f.write("Abstract:\n")
                f.write(study.abstract or "No abstract available")
                f.write("\n\n")
                
                # Write relevance assessment
                f.write("Relevance Assessment:\n")
                f.write(f"Score: {study.relevance_score or 'N/A'}\n")
                f.write(f"Explanation: {study.inclusion_explanation or 'N/A'}\n")
                f.write("\n")
                
                # Write PDF status
                f.write("PDF Status:\n")
                f.write(f"Download Status: {study.download_status}\n")
                if study.download_date:
                    f.write(f"Download Date: {study.download_date}\n")
                if study.pdf_path:
                    f.write(f"PDF Path: {study.pdf_path}\n")
            
            logging.info(f"Stored text content for study {study.pmid}")
            return True
            
        except Exception as e:
            logging.error(f"Error storing text content for study {study.pmid}: {str(e)}")
            return False
    
    def add_study(self, study: StudyMetadata):
        """Add or update a study in the storage."""
        self.studies[study.pmid] = study
        self._save_metadata()
        # Also store text content
        self.store_study_text(study)
    
    def get_study(self, pmid: str) -> StudyMetadata:
        """Get a study by its PubMed ID."""
        return self.studies.get(pmid)
    
    def get_all_studies(self) -> List[StudyMetadata]:
        """Get all stored studies."""
        return list(self.studies.values())
    
    def download_pdf(self, pmid: str, url: str) -> bool:
        """
        Download PDF for a study.
        
        Args:
            pmid: PubMed ID of the study
            url: URL to download the PDF from
            
        Returns:
            bool: True if download was successful, False otherwise
        """
        if pmid not in self.studies:
            logging.error(f"Study {pmid} not found in metadata")
            return False
        
        study = self.studies[pmid]
        pdf_path = self.articles_dir / f"{pmid}.pdf"
        
        try:
            # Attempt to download the PDF
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            # Check if it's actually a PDF
            content_type = response.headers.get('content-type', '').lower()
            if 'pdf' not in content_type:
                raise ValueError(f"URL does not point to a PDF (content-type: {content_type})")
            
            # Save the PDF
            with open(pdf_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            # Update metadata
            study.pdf_path = str(pdf_path)
            study.download_status = "success"
            study.download_date = datetime.now()
            self._save_metadata()
            
            logging.info(f"Successfully downloaded PDF for study {pmid}")
            return True
            
        except Exception as e:
            # Update metadata with failure
            study.download_status = f"failed: {str(e)}"
            study.download_date = datetime.now()
            self._save_metadata()
            
            logging.error(f"Error downloading PDF for study {pmid}: {str(e)}")
            return False
    
    def bulk_download_pdfs(self):
        """Attempt to download PDFs for all studies with URLs."""
        for study in self.studies.values():
            if study.url and study.download_status == "not_attempted":
                logging.info(f"Attempting to download PDF for study {study.pmid}")
                self.download_pdf(study.pmid, study.url) 