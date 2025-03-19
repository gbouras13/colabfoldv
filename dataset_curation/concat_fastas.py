import string
from Bio import SeqIO
from Bio.Seq import Seq
from argparse import RawTextHelpFormatter
import argparse
import os



def get_input():
    parser = argparse.ArgumentParser(
        description="concat_fastas.py",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--indir", action="store", help="Input dir with subdirs aa..zz containing fasta format proteins from enVhogs."
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
    output_file_path = os.path.join(args.outdir, "envhogs_proteins_combined.fasta")
    subdirs = [f"{first}{second}" for first in string.ascii_lowercase for second in string.ascii_lowercase]

    with open(output_file_path, 'w') as output_file:
	
        for subdir in subdirs:
            # test done file
            print(f"concatenating all fastas in {subdir}")
            for f in os.listdir(os.path.join(args.indir,subdir)):
                if f.endswith('.fasta') :
                    inpath=os.path.join(args.indir, subdir, f)
                    with open(inpath, 'r') as fasta_file:
                        # Read the content of the FASTA file
                        content = fasta_file.read()

                        # Write the content to the output file
                        output_file.write(content)
            

    
        







if __name__ == "__main__":
    main()



# Replace 'input.a3m' with the path to your A3M file and 'output.fasta' with the desired output file path
#remove_gaps_from_a3m('a3ms/1000', 'output.fasta')

