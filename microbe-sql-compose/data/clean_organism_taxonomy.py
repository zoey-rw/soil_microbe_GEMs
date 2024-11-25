import pandas as pd

# Load the CSV file
input_file = "organism_taxonomy.csv"
output_file = "cleaned_organism_taxonomy.csv"

# Read the CSV into a DataFrame
df = pd.read_csv(input_file)

# Replace "NA" with None (interpreted as NULL in PostgreSQL)
df.replace("NA", None, inplace=True)

# Save the cleaned DataFrame to a new CSV file
df.to_csv(output_file, index=False)

print(f"Cleaned CSV saved to {output_file}")
