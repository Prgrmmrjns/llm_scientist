from user_input import collect_user_input
from params import SystemParams
from literature_search import search_literature
from study_selection import filter_studies
from reporting import generate_latex_report, generate_plots
from content_generation import generate_introduction_content

def run_systematic_review():
    """
    Main function to run the systematic review process.
    """
    # 1. Gather user preferences
    user_params = collect_user_input()
    system_params = SystemParams(user_params)

    # 2. Generate introduction content
    print("\nGenerating introduction content...")
    intro_content = generate_introduction_content(system_params.topic)
    print(f"\nTitle: {intro_content.title}")
    print(f"Subtitle: {intro_content.subtitle}")
    print("\nKey Points:")
    for point in intro_content.key_points:
        print(f"- {point}")
    print("\nResearch Questions:")
    for question in intro_content.research_questions:
        print(f"- {question}")

    # 3. Fetch relevant literature
    print("\nSearching for relevant literature...")
    search_results = search_literature(system_params)

    # 4. Filter and rank studies
    print("\nAnalyzing and filtering studies...")
    final_studies = filter_studies(search_results, system_params)

    # 5. Generate outputs
    print(f"\nFinal set of {len(final_studies)} studies included in the systematic review:")
    for study in final_studies:
        print(f"- {study.title} (Score: {study.relevance_score or 'N/A'})")

    # Generate a LaTeX report
    print("\nGenerating LaTeX report...")
    generate_latex_report(final_studies, intro_content=intro_content)
    
    # Generate plots
    print("\nGenerating visualization plots...")
    generate_plots(final_studies)

if __name__ == "__main__":
    run_systematic_review()