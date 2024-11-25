import pandas as pd

# Input and output file paths
input_file = "species_abundance_filt.csv"
output_file = "cleaned_species_abundance_filt.csv"

# Load the CSV into a pandas DataFrame
df = pd.read_csv(input_file)

# Replace "NA" or any invalid placeholder with None (interpreted as NULL in PostgreSQL)
df.replace("NA", None, inplace=True)

# Save the cleaned DataFrame to a new CSV file
df.to_csv(output_file, index=False)

print(f"Cleaned CSV saved to {output_file}")
