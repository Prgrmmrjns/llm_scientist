class SystemParams:
    """
    System parameters for the systematic review process.
    """
    def __init__(self, user_params):
        self.topic = user_params.get('topic', '')
        self.journals = user_params.get('journals', [])
        self.date_range = user_params.get('date_range', None)
        self.min_relevance_score = user_params.get('min_relevance_score', 70)
        self.email = user_params.get('email', 'your_email@example.com')  # Required for PubMed API
        self.max_results = user_params.get('max_results', 10)  # Number of papers to fetch
        
        # LaTeX template settings
        self.latex_settings = {
            'document_class': 'article',
            'font_size': '11pt',
            'margin': '2.5cm',
            'line_spacing': '1.5'
        }
