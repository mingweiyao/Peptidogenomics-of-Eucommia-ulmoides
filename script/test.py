# import os
# from Bio import SeqIO

# def test_transcript(directory):
#     for f in os.listdir(directory):
#         if f.endswith('.gff3'):
#             with open(os.path.join(directory, f), 'r') as file:
#                 for line in file:
#                     if 'transcript' in line:
#                         print(f"{f}")
#                         break
#         elif f.endswith('.fasta'):
#             for record in SeqIO.parse(os.path.join(directory, f), "fasta"):
#                 if 'transcript' in record.description:
#                     print(f"{f}")
#                     break
# def main():
#     gff3_dictorty = r"G:\Eu_peptido\20251018 imeta\file\00raw\GFF3"
#     fasta_directory = r"G:\Eu_peptido\20251018 imeta\file\00raw\fasta"
#     test_transcript(gff3_dictorty)
#     test_transcript(fasta_directory)

# if __name__ == "__main__":
#     main()

# from Bio import SeqIO
# count = 0
# for rec in SeqIO.parse(r"G:\Eu_peptido\20251018 imeta\file\00raw\output\file2_nonredundant_coding_transcript_pep.fasta", "fasta"):
#     count += 1
# print(count)


count = 0
with open(r"G:\Eu_peptido\20251018 imeta\file\00raw\output\file1_nonredundant_coding_transcripts.gff3") as f:
    for line in f:
        if "gene" in line:
            count+=1
print(count)