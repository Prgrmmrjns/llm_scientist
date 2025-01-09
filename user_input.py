def collect_user_input():
    """
    Collect user preferences and constraints from CLI or GUI.
    Return a dictionary with the topic, journals, date range, etc.
    """
    user_params = {
        "topic": "Machine Learning in Healthcare",
        "journals": ["Nature", "Science", "Lancet", "BMJ"],
        "date_range": ("2020-01-01", "2023-12-31"),
        "min_relevance_score": 70,
        "email": "jonascw@web.de",  # Replace with your email
        "max_results": 10  # Number of papers to fetch from PubMed
    }
    return user_params