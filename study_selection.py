import os
from dotenv import load_dotenv
from openai import OpenAI
from models import StudyMetadata
from typing import Tuple, List
import logging

load_dotenv()

client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

def assess_relevance(study: StudyMetadata, system_params) -> Tuple[int, str]:
    """
    Use LLM to assess the relevance of a study based on its abstract.
    
    Args:
        study: StudyMetadata object containing study information
        system_params: System parameters including the topic
        
    Returns:
        Tuple of (score, explanation)
    """
    prompt = f"""
    Please assess the relevance of this study for a systematic review on {system_params.topic}.
    
    Title: {study.title}
    Abstract: {study.abstract if study.abstract else 'No abstract available'}
    
    Rate the relevance on a scale of 0-100 and provide a brief explanation.
    Return your response in the format: "Score: [number]|[explanation]"
    Example: "Score: 85|Highly relevant due to direct application of LLMs"
    """
    
    response = client.chat.completions.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "system", "content": "You are a research assistant helping with systematic review study selection. Always respond with the exact format: 'Score: [0-100]|[explanation]'"},
            {"role": "user", "content": prompt}
        ],
        temperature=0.3
    )
    
    result = response.choices[0].message.content.strip()
    try:
        # Extract score and explanation
        score_part, explanation = result.split('|', 1)
        score = int(score_part.replace('Score:', '').strip())
        return score, explanation.strip()
    except (ValueError, AttributeError) as e:
        logging.warning(f"Warning: Could not parse LLM response properly: {result}")
        # Return a default score if parsing fails
        return 70, "Score parsing failed, using default score"

def filter_studies(search_results: List[StudyMetadata], system_params) -> List[StudyMetadata]:
    """
    Filter studies based on LLM relevance assessment.
    
    Args:
        search_results: List of StudyMetadata objects
        system_params: System parameters including minimum score threshold
        
    Returns:
        List of filtered and ranked StudyMetadata objects
    """
    filtered_studies = []
    
    for study in search_results:
        score, explanation = assess_relevance(study, system_params)
        # Update the study object with the relevance information
        study.relevance_score = score
        study.inclusion_explanation = explanation
        
        if score >= system_params.min_relevance_score:
            filtered_studies.append(study)
    
    # Sort by relevance score
    filtered_studies.sort(key=lambda x: x.relevance_score or 0, reverse=True)
    return filtered_studies
