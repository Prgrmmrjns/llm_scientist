import os
from dotenv import load_dotenv
from openai import OpenAI
import json
from models import IntroductionContent

load_dotenv()

client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

def generate_introduction_content(topic: str) -> IntroductionContent:
    """
    Generate introduction content for the systematic review using OpenAI API.
    
    Args:
        topic (str): The main topic for the systematic review
        
    Returns:
        IntroductionContent: Structured content for the introduction
    """
    prompt = f"""Generate a structured introduction for a systematic review on the topic: {topic}
    
    Please provide:
    1. A compelling main title (should be concise and academic, ending with ": A Systematic Review")
    2. An optional subtitle (should complement the title and provide additional context)
    3. 4-5 key introduction points that:
       - Highlight the importance and context of this review
       - Explain current challenges or gaps
       - Discuss potential impacts
       - Each point should be a complete sentence
       - Avoid using special characters or mathematical symbols
    4. 2-3 specific research questions that:
       - Are clear and focused
       - Can be answered through literature review
       - Start with words like "What", "How", or "To what extent"
    
    Format your response as a JSON object with the following structure:
    {{
        "title": "string",
        "subtitle": "string",
        "key_points": ["string"],
        "research_questions": ["string"]
    }}
    
    Make the content academic and professional in tone, avoiding any special characters or formatting that might cause issues in LaTeX."""

    try:
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {
                    "role": "system", 
                    "content": "You are a research assistant helping to structure systematic reviews. Provide clear, academic content that will work well in LaTeX documents. Avoid special characters and ensure proper formatting."
                },
                {"role": "user", "content": prompt}
            ],
            response_format={"type": "json_object"},
            temperature=0.7
        )
        
        # Parse the response and validate with Pydantic
        content_dict = json.loads(response.choices[0].message.content)
        introduction_content = IntroductionContent(**content_dict)
        
        return introduction_content
    
    except Exception as e:
        print(f"Error generating introduction content: {str(e)}")
        # Return a default structure if there's an error
        return IntroductionContent(
            title=f"A Systematic Review of {topic}",
            subtitle="Current State, Challenges, and Future Directions",
            key_points=[
                "This field has shown significant growth and development in recent years.",
                "Current research indicates promising applications and potential benefits.",
                "Several challenges and limitations need to be addressed.",
                "A comprehensive review is needed to guide future research efforts."
            ],
            research_questions=[
                "What are the current trends and developments in this field?",
                "What are the main challenges and potential solutions?",
                "How can future research address existing limitations?"
            ]
        ) 