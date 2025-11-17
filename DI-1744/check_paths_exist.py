import os
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Check if files and their directories exist from a CSV column.")
    parser.add_argument('--csv', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('--rows', type=int, default=None, help='Number of rows to check (default: all).')
    parser.add_argument('--output', type=str, default=None, help='Optional output CSV path.')

    args = parser.parse_args()

    # Load CSV
    df = pd.read_csv(args.csv)

    # Ensure 'path' column exists
    if 'path' not in df.columns:
        raise ValueError("CSV does not contain a 'path' column.")

    # Subset rows if requested
    if isinstance(args.rows, int):
        df = df.head(args.rows)

    # File existence check
    df['file_exists'] = df['path'].apply(lambda p: os.path.isfile(p))

    # Directory existence check
    df['directory_exists'] = df['path'].apply(lambda p: os.path.isdir(os.path.dirname(p)))

    # Print results
    print(df[['path', 'file_exists', 'directory_exists']])

    # Save output if specified
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\nResults saved to {args.output}")

if __name__ == '__main__':
    main()