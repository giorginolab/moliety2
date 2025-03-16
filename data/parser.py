# Your task is to parse the attached pseudo-html file into a data frame.
# The file contains a list of SMARTS rules and descriptions, organized
# as follows. The HTML tagging is not well-formed, so run tests until
# the pandas dataframe looks right.

# Rules are organised hierarchically:

# H2 tags starting with a number are main topics
# H2 without a number introduce subtopics (ignore badly formatted html)
# H3 are sub-sub-topics
# DT has rule name
# DD rule contents (smarts)
# if further DDs are encountered, they are comments
# if several patterns appear in the DD separated by BR, split them



from bs4 import BeautifulSoup
import pandas as pd

# read the pseudo-html file
with open("smarts_examples.ff.txt", "r", encoding="utf-8") as f:
    html_content = f.read()

# use BeautifulSoup with the built-in parser (which can handle messy HTML)
soup = BeautifulSoup(html_content, "html.parser")

# variables to keep track of the current hierarchy
current_main_topic = None
current_subtopic = None
current_sub_sub_topic = None

rows = []

# iterate over all tags in document order
for element in soup.recursiveChildGenerator():
    if hasattr(element, "name"):
        # Process H2 tags
        if element.name == "h2":
            text = element.get_text(strip=True, separator=' ')
            # if H2 text starts with a digit, it is a main topic
            if text and text[0].isdigit():
                current_main_topic = text
                current_subtopic = None
                current_sub_sub_topic = None
            else:
                # H2 that does not start with a number: treat as subtopic
                current_subtopic = text
                current_sub_sub_topic = None

        # Process H3 tags as sub-sub-topics
        elif element.name == "h3":
            current_sub_sub_topic = element.get_text(strip=True, separator=' ')

        # Process DT tags which contain the rule name
        elif element.name == "dt":
            rule_name = element.get_text(strip=True, separator=' ')
            # Find the next DD sibling
            dd = element.find_next_sibling("dd")
            if dd:
                # If the dd has <br> tags, split on them to get separate patterns.
                if dd.find("br"):
                    # decode the inner HTML and split on <br>, ignoring newlines
                    raw_html = dd.decode_contents().replace('\n', '')
                    pattern_list = [part.strip() for part in raw_html.split("<br/>") if part.strip()]
                    # For multiple patterns, append index to rule name
                    if len(pattern_list) > 1:
                        rule_names = [f"{rule_name} ({i+1})" for i in range(len(pattern_list))]
                    else:
                        rule_names = [rule_name]
                else:
                    pattern_list = [dd.get_text(strip=True, separator=' ')]
                    rule_names = [rule_name]
                
                # Look for an additional DD that might contain comments and normalize whitespace
                comment_dd = dd.find_next_sibling("dd")
                if comment_dd:
                    # Replace newlines with spaces and normalize whitespace
                    comments = ' '.join(comment_dd.get_text(separator=' ').split())
                else:
                    comments = ""

                # Create a row for each pattern with its corresponding rule name
                for pattern, name in zip(pattern_list, rule_names):
                    rows.append({
                        "Main Topic": current_main_topic,
                        "Subtopic": current_subtopic,
                        "Sub-sub-topic": current_sub_sub_topic,
                        "Rule Name": name,
                        "Pattern": pattern,
                        "Comments": comments
                    })

# build the DataFrame
df = pd.DataFrame(rows)
# display a sample of the dataframe
df.to_csv("smarts_examples_parsed.csv", index=False)
df.to_clipboard()

