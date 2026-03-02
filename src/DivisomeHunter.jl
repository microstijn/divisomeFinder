module DivisomeHunter

export hunt_divisome, DivisomeCriteria, ProteinCandidate

using DataFrames

include("Types.jl")
using .Types

include("AccessionMapper.jl")
include("DataParser.jl")
include("UniProtAPI.jl")
include("Classifier.jl")

using .DataParser
using .UniProtAPI
using .Classifier

"""
    hunt_divisome(tradis_csv_path::String, fasta_path::String; criteria::Any=DivisomeCriteria(), from_db::String="UniProtKB_AC-ID") -> Tuple{Vector{ProteinCandidate}, DataFrame}

Main entry point for DivisomeHunter.jl. Automates the discovery of novel cell division machinery
candidates in non-canonical bacteria by applying rigid biophysical and structural filters
to essential genes defined by TraDIS.

# Arguments
- `tradis_csv_path::String`: Path to the CSV file containing essentiality data (e.g., TraDIS output).
- `fasta_path::String`: Path to the Proteome FASTA file mapping locus tags to UniProt Accessions.
- `criteria::Any`: Optional configuration struct defining the bounds of the search
  (mass thresholds, keywords). Defaults to `DivisomeCriteria()`.
- `from_db::String`: Optional UniProt mapping database source. Defaults to `"UniProtKB_AC-ID"`.
  If parsing GenBank/RefSeq FASTA files with `[protein_id=...]`, use `"EMBL-GenBank-DDBJ_CDS"`.

# Returns
- A tuple containing:
  1. A `Vector{ProteinCandidate}` array of the highly-typed structs representing each hit.
  2. A `DataFrame` projection of the candidates for easy CSV export and analysis.

# Biological Overview
1.  Identifies genes explicitly marked as "essential" or "essential#" in the provided TraDIS dataset.
2.  Parses the FASTA file to establish a link between the gene `Locus_tag` and its `UniProt Accession`.
3.  Utilizes the UniProt REST API (asynchronously) to download massive GZIP-compressed metadata entries for the candidates.
4.  Applies the `DivisomeCriteria` filters to evaluate each candidate based on mass, sequence length, topology, and InterPro domains.
"""
function hunt_divisome(tradis_csv_path::String, fasta_path::String; criteria::Any=DivisomeCriteria(), from_db::String="UniProtKB_AC-ID")
    # 1. Parse essential genes
    println("🧬 Step 1: Parsing Essential Locus Tags...")
    essential_loci = DataParser.parse_essential_genes(tradis_csv_path)
    println("✅ Found $(length(essential_loci)) essential locus tags.")

    # 2. Map Locus to Accession via FASTA
    println("\n🧬 Step 2: Mapping Locus Tags to UniProt Accessions...")
    accessions, locus_to_acc = DataParser.map_locus_to_accession(fasta_path, essential_loci)
    println("✅ Mapped $(length(accessions)) valid Accessions.")

    # 3. Fetch UniProt JSON Data
    println("\n🧬 Step 3: Fetching Data from UniProt...")
    json_data = UniProtAPI.fetch_uniprot_data(accessions; from_db=from_db)

    # 4. Parse JSON and Classify
    println("\n🧬 Step 4: Classifying Candidates based on Criteria...")
    candidates = Classifier.classify_candidates(json_data, locus_to_acc, criteria)

    # Convert to DataFrame for easier inspection/export
    df_out = DataFrame([
        (
            Locus_tag = c.locus_tag,
            Accession = c.accession,
            Candidate_Type = c.candidate_type,
            Protein_Name = c.protein_name,
            Mass_kDa = c.mass_kda,
            Length_AA = c.length_aa,
            Transmembrane = c.transmembrane,
            Signal_Peptide = c.signal_peptide,
            Domains = c.domains
        ) for c in candidates
    ])

    println("🎉 Divisome Hunt Complete! Found $(length(candidates)) candidates.")
    return candidates, df_out
end

end # module
