import string
from Bio import SeqIO
from Bio.Seq import Seq
from argparse import RawTextHelpFormatter
import argparse
import os

def remove_gaps_from_a3m(input_file, output_file):
    # Read the A3M file
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # Parse the A3M file using SeqIO
        for record in SeqIO.parse(infile, "fasta"):
            if "consensus" not in record.description:
            # Remove gaps (dashes) from the sequence
                ungapped_sequence = str(record.seq).replace("-", "")
            # Update the record's sequence with the ungapped sequence
                record.seq = Seq(ungapped_sequence)
            # Write the ungapped sequence to the output file in FASTA format
                SeqIO.write(record, outfile, "fasta")


def get_input():
    parser = argparse.ArgumentParser(
        description="filter_fastas.py",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--indir", action="store", help="Input dir with subdirs aa..zz containing a3m format MSAs from enVhogs."
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        help="Directory to write the output to.",
        default=os.path.join(os.getcwd(), "output/"),
    )

    args = parser.parse_args()

    return args

def main():
    args = get_input()
    
    os.makedirs(args.outdir, exist_ok=True) 

    subdirs = [f"{first}{second}" for first in string.ascii_lowercase for second in string.ascii_lowercase]

# Create the subdirectories if they don't already exist
    for subdir in subdirs:
        os.makedirs(os.path.join(args.outdir, subdir), exist_ok=True)

    for subdir in subdirs:
        # test done file
        print(f"parsing all a3m alignments in {subdir}")
        touch_file=os.path.join(args.outdir,f"{subdir}.done")
        if os.path.isfile(touch_file) is False:
            for f in os.listdir(os.path.join(args.indir,subdir)):
                if os.path.isfile(os.path.join(args.indir, subdir, f)) :
                    inpath=os.path.join(args.indir, subdir, f)
                    outpath=os.path.join(args.outdir, subdir, f"{f}.fasta")
                    remove_gaps_from_a3m(inpath, outpath)
            # touch the file once done
            with open(touch_file, 'w') as file:
                file.write("done")
        else:
            print(f"{subdir}.done exists. Moving on")

    
        







if __name__ == "__main__":
    main()



# Replace 'input.a3m' with the path to your A3M file and 'output.fasta' with the desired output file path
#remove_gaps_from_a3m('a3ms/1000', 'output.fasta')

