from user_input import collect_user_input
from params import SystemParams
from literature_search import search_literature
from study_selection import filter_studies
from reporting import generate_latex_report, generate_plots

def run_systematic_review():
    """
    Main function to run the systematic review process.
    """
    # 1. Gather user preferences
    user_params = collect_user_input()
    system_params = SystemParams(user_params)

    # 2. Fetch relevant literature
    search_results = search_literature(system_params)

    # 3. Filter and rank studies
    final_studies = filter_studies(search_results, system_params)

    # 4. Generate outputs
    print(f"\nFinal set of {len(final_studies)} studies included in the systematic review:")
    for study in final_studies:
        print(f"- {study['title']} (Score: {study.get('relevance_score', 'N/A')})")

    # Generate a LaTeX report
    generate_latex_report(final_studies)
    
    # Generate plots
    generate_plots(final_studies)

if __name__ == "__main__":
    run_systematic_review()