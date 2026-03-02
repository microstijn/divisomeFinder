using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using divisomeFinder


fa = raw"D:\divisome\data\fasta\CP001744.1.fa.txt"
moesm = raw"D:\divisome\data\essentiality\moesm.csv"

    println("🧬 Step 1: Parsing Essential Locus Tags...")
    using divisomeFinder.DataParser
    essential_loci = DataParser.parse_essential_genes(moesm)
    println("✅ Found $(length(essential_loci)) essential locus tags.")

    # 2. Map Locus to Accession via FASTA
    println("\n🧬 Step 2: Mapping Locus Tags to UniProt Accessions...")
    accessions, locus_to_acc = DataParser.map_locus_to_accession(fa, essential_loci)
    println("✅ Mapped $(length(accessions)) valid Accessions.")

    # 3. Fetch UniProt JSON Data
    println("\n🧬 Step 3: Fetching Data from UniProt...")
     using divisomeFinder.UniProtAPI
    json_data = UniProtAPI.fetch_uniprot_data(accessions)

    # 4. Parse JSON and Classify
    println("\n🧬 Step 4: Classifying Candidates based on Criteria...")
    using divisomeFinder.Classifier
    candidates = Classifier.classify_candidates(json_data, locus_to_acc, DivisomeCriteria())
using DataFrames
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

df_out[:, [:Candidate_Type, :Signal_Peptide, :Transmembrane, :Domains]]
using CSV
CSV.write(raw"D:\divisome\data\out\out.csv", df_out)