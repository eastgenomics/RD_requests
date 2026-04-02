import pandas as pd

# Read the CSV files
export_df = pd.read_csv(
  '/home/katherine/Downloads/250603_01_wgs_rd_nuh.csv',
  usecols=['Test Order Priority','Order Request Date Time', 'WGS Referral ID']
)
current_df = pd.read_excel(
  '/home/katherine/Downloads/NUH-Ready for analysis via Epic.xlsx',
  usecols=['NGIS ID', 'Date Case reported at NUH', 'Download link']
)

# Clean the ID fields (remove any whitespace, or spaces in R numbers)
export_df['WGS Referral ID'] = export_df['WGS Referral ID'].str.strip().replace(" ", "")
current_df['NGIS ID'] = current_df['NGIS ID'].str.strip().replace(" ", "")

# Create a set of IDs from current.csv that have either a date interpretion started or a download link
already_done = current_df.loc[
    (current_df['Date Case reported at NUH'].notna()) |  # Date is not blank
    (current_df['Download link'].notna())    # Link is not blank
, 'NGIS ID'].to_list()

# Filter Clarity export to remove these
filtered_records = export_df[~export_df['WGS Referral ID'].isin(already_done)]

# Save the result to a new CSV
filtered_records.to_csv('filtered_results.csv', index=False, header=['Priority', 'Date', 'R Number'])
