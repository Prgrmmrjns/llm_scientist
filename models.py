from pydantic import BaseModel, Field
from typing import List, Optional
from datetime import datetime

class IntroductionContent(BaseModel):
    """Pydantic model for the introduction content returned by the LLM."""
    title: str = Field(..., description="The main title for the systematic review")
    subtitle: str = Field(None, description="Optional subtitle for the systematic review")
    key_points: List[str] = Field(..., description="List of key introduction points")
    research_questions: List[str] = Field(..., description="List of main research questions to address")
    
    class Config:
        json_schema_extra = {
            "example": {
                "title": "The Impact of Machine Learning in Healthcare: A Systematic Review",
                "subtitle": "Trends, Applications, and Future Directions",
                "key_points": [
                    "Machine learning is transforming healthcare delivery and outcomes",
                    "Recent advances show promising results in diagnosis and treatment",
                    "Integration challenges remain a significant concern"
                ],
                "research_questions": [
                    "How effective are ML models in clinical decision support?",
                    "What are the main implementation challenges?",
                    "What are the ethical considerations?"
                ]
            }
        } 

class StudyMetadata(BaseModel):
    """Pydantic model for storing study metadata."""
    pmid: str = Field(..., description="PubMed ID of the study")
    title: str = Field(..., description="Title of the study")
    authors: List[str] = Field(default_factory=list, description="List of authors")
    journal: Optional[str] = Field(None, description="Journal name")
    year: Optional[int] = Field(None, description="Publication year")
    doi: Optional[str] = Field(None, description="Digital Object Identifier")
    abstract: Optional[str] = Field(None, description="Study abstract")
    url: Optional[str] = Field(None, description="URL to the full text")
    pdf_path: Optional[str] = Field(None, description="Local path to stored PDF")
    download_status: str = Field(default="not_attempted", description="Status of PDF download attempt")
    download_date: Optional[datetime] = Field(None, description="Date when PDF was downloaded")
    relevance_score: Optional[int] = Field(None, description="Relevance score from LLM assessment")
    inclusion_explanation: Optional[str] = Field(None, description="Explanation for the relevance score and inclusion decision")
    
    class Config:
        json_schema_extra = {
            "example": {
                "pmid": "12345678",
                "title": "Example Study Title",
                "authors": ["Author One", "Author Two"],
                "journal": "Example Journal",
                "year": 2023,
                "doi": "10.1234/example.doi",
                "abstract": "Example abstract text",
                "url": "https://example.com/paper",
                "pdf_path": "articles/12345678.pdf",
                "download_status": "success",
                "download_date": "2024-02-21T10:30:00",
                "relevance_score": 85,
                "inclusion_explanation": "Highly relevant due to direct application of methods"
            }
        } 