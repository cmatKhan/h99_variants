import sys
import pandas as pd

valid_chr_names = [
    'CP003820.1','CP003821.1','CP003822.1','CP003823.1','CP003824.1',
    'CP003825.1','CP003826.1','CP003827.1','CP003828.1','CP003829.1',
    'CP003830.1','CP003831.1','CP003832.1','CP003833.2','CP003834.1'
]

seqlengths = [
    2291499, 1621675, 1575141, 1084805, 1814975, 1422463, 1399503,
    1398693, 1186808, 1059964, 1561994, 774062, 756744, 942867, 24919
]

seqname_to_seqlength = dict(zip(valid_chr_names, seqlengths))

def validate_bed(file_path):
    try:
        # Load the BED file
        bed = pd.read_csv(file_path, sep='\t', header=None, comment='#')

        # Check for minimum required columns (chrom, start, end, ., 0, .)
        if bed.shape[1] != 6:
            raise ValueError("BED6 file must have exactly 6 columns")

        # Check if start and end columns are integers
        if not pd.api.types.is_integer_dtype(bed[1]) or not pd.api.types.is_integer_dtype(bed[2]):
            raise ValueError("Start and end columns must be integers")

        # Check if start < end
        if not (bed[1] < bed[2]).all():
            raise ValueError("Start positions must be less than the end positions")

        # Check if chromosome names are valid
        if not bed[0].isin(valid_chr_names).all():
            raise ValueError("Invalid chromosome names found in the BED file")

        # Check if start values are >= 0 and <= corresponding seqlengths
        for idx, row in bed.iterrows():
            chrom = row[0]
            start = row[1]
            end = row[2]
            if start < 0 or start > seqname_to_seqlength[chrom]:
                raise ValueError(f"Start position out of bounds for chromosome {chrom} at row {idx + 1}")
            if end < 0 or end > seqname_to_seqlength[chrom]:
                raise ValueError(f"End position out of bounds for chromosome {chrom} at row {idx + 1}")

        # Check that columns 3, 4, and 5 have the values ".", "0", and "." respectively
        if not (bed[3] == '.').all():
            raise ValueError("Column 3 must have the value '.'")
        if not (bed[4] == 0).all():
            raise ValueError("Column 4 must have the value '0'")
        if not (bed[5] == '.').all():
            raise ValueError("Column 5 must have the value '.'")

        print("BED file is valid.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python validate_bed.py <path_to_bed_file>")
        sys.exit(1)

    bed_file_path = sys.argv[1]
    validate_bed(bed_file_path)

