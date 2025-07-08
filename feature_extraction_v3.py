import gffutils
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd
from tqdm import tqdm
import os

# CONFIG
GFF_FILE = "/home/isagallegor/M.roreri_Side_Quests/data/MrorC26.groups.gff3"
FASTA_FILE = "/home/isagallegor/M.roreri_Side_Quests/data/MrorC26.groups.fasta"
OUTPUT = "features/gene_features_v3.csv"
DB_FILE = "/home/isagallegor/M.roreri_Side_Quests/data/genome.db"

# Create database from GFF if not exists
if not os.path.exists(DB_FILE):
    print("Creating database from GFF...")
    gffutils.create_db(GFF_FILE, DB_FILE, force=True, keep_order=True,
                       merge_strategy='merge', sort_attribute_values=True)
db = gffutils.FeatureDB(DB_FILE, keep_order=True)

# Load genome sequence
print("Loading genome FASTA...")
seq_dict = SeqIO.to_dict(SeqIO.parse(FASTA_FILE, "fasta"))

# Extract features from genes
print("Extracting gene features...")
genes = list(db.features_of_type('gene'))
features = []

for gene in tqdm(genes):
    gene_id = gene.id
    contig = gene.seqid
    start = gene.start
    end = gene.end
    strand = gene.strand
    gene_length = end - start + 1
    start_pos_norm = start / len(seq_dict[contig]) if contig in seq_dict else None

    # Look for mRNA children if exists
    mrnas = list(db.children(gene, featuretype='mRNA', level=1))
    if mrnas:
        target = mrnas[0]
    else:
        target = gene

    # CDS 
    cds_list = list(db.children(target, featuretype='CDS', level=1))
    is_coding = int(len(cds_list) > 0)
    cds_length = sum(cds.end - cds.start + 1 for cds in cds_list)

    # Exons
    exons = list(db.children(target, featuretype='exon', level=1))
    exon_count = len(exons)

    # Sequence and GC content
    try:
        seq = seq_dict[contig].seq[start-1:end]
        gc_content = gc_fraction(seq) * 100  # percent
    except:
        gc_content = None

    features.append({
        "gene_id": gene_id,
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "gene_length": gene_length,
        "start_position_norm": start_pos_norm,
        "exon_count": exon_count,
        "cds_length": cds_length,
        "gc_content": gc_content,
        "is_coding": is_coding
    })

# Save
df = pd.DataFrame(features)
df.to_csv(OUTPUT, index=False)
print(f"Features saved to {OUTPUT}")
