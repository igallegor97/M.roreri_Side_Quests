import pandas as pd
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import os

# Config
GFF_FILE = "data/MrorC26.groups.gff3"
FASTA_FILE = "data/Mror_C26.groups.fasta"
DB_FILE = "data/genome.db"
DEG_FILE = "features/deg_dataset_v2.csv"
FASTA_OUT = "promoters/promoters.fasta"
CSV_OUT = "promoters/promoters.csv"
PROMOTER_UP = 1000
PROMOTER_DOWN = 100

# Build or load GFF database
if not os.path.exists(DB_FILE):
    print("Creating GFF database...")
    gffutils.create_db(GFF_FILE, DB_FILE, force=True, keep_order=True,
                       merge_strategy='merge', sort_attribute_values=True)
db = gffutils.FeatureDB(DB_FILE, keep_order=True)

# Load genome FASTA 
print("Loading genome FASTA...")
genome = SeqIO.to_dict(SeqIO.parse(FASTA_FILE, "fasta"))

# Load DEG annotations
print("Loading DEG labels...")
deg_df = pd.read_csv(DEG_FILE)
deg_labels = dict(zip(deg_df["gene_id"], deg_df["DEG"]))

# Extract promoter sequences
print("Extracting promoter regions...")
promoter_records = []
metadata = []

for gene in tqdm(db.features_of_type('gene')):
    gene_id = gene.id
    contig = gene.seqid
    strand = gene.strand

    if contig not in genome:
        continue

    seq_len = len(genome[contig])
    chrom_seq = genome[contig].seq

    # Define promoter region depending on strand
    if strand == "+":
        start = max(gene.start - PROMOTER_UP, 1)
        end = min(gene.start + PROMOTER_DOWN, seq_len)
        promoter_seq = chrom_seq[start-1:end]  # python is 0-based
    elif strand == "-":
        start = max(gene.end - PROMOTER_DOWN, 1)
        end = min(gene.end + PROMOTER_UP, seq_len)
        promoter_seq = chrom_seq[start-1:end].reverse_complement()
    else:
        continue  # skip if strand info is missing

    if len(promoter_seq) != (PROMOTER_UP + PROMOTER_DOWN):
        continue  # discard incomplete promoters

    # Get DEG label (default = 0)
    deg_label = deg_labels.get(gene_id, 0)

    # Save as FASTA record
    record = SeqRecord(promoter_seq, id=gene_id, description="")
    promoter_records.append(record)

    # Save metadata
    metadata.append({
        "gene_id": gene_id,
        "contig": contig,
        "strand": strand,
        "promoter_seq": str(promoter_seq),
        "DEG": deg_label
    })

# Write outputs
print(f"Writing {len(promoter_records)} promoters to FASTA and CSV...")
os.makedirs("promoters", exist_ok=True)
SeqIO.write(promoter_records, FASTA_OUT, "fasta")
pd.DataFrame(metadata).to_csv(CSV_OUT, index=False)
print("Done.")
