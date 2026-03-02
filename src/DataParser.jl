module DataParser

export parse_essential_genes, map_locus_to_accession

using CSV
using DataFrames

"""
    parse_essential_genes(tradis_csv_path::String) -> Set{String}

Reads a TraDIS CSV file and extracts the locus tags of essential genes.
Rows are considered essential if `Essentiality` is exactly "essential" or "essential#".
"""
function parse_essential_genes(tradis_csv_path::String)
    df = CSV.read(tradis_csv_path, DataFrame, missingstring=["", "NA"])
    essential_df = filter(row -> !ismissing(row.Essentiality) &&
                                (row.Essentiality == "essential" || row.Essentiality == "essential#"), df)

    # Ensure locus tags are converted to a strict Set{String}
    essential_loci = Set{String}()
    for locus in essential_df.Locus_tag
        if !ismissing(locus)
            push!(essential_loci, String(locus))
        end
    end

    return essential_loci
end

"""
    map_locus_to_accession(fasta_path::String, essential_loci::Set{String}) -> Tuple{Vector{String}, Dict{String, String}}

Parses a FASTA file containing headers mapped to locus tags and UniProt accessions.
Extracts the accessions for genes found in the `essential_loci` set.

# Returns
- A tuple containing:
  1. A `Vector{String}` of UniProt accessions.
  2. A `Dict{String, String}` mapping accessions back to their original locus tags.
"""
function map_locus_to_accession(fasta_path::String, essential_loci::Set{String})
    accessions = String[]
    locus_to_acc = Dict{String, String}()

    open(fasta_path, "r") do io
        for line in eachline(io)
            if startswith(line, ">")
                parts = split(strip(line[2:end]), "|")
                if length(parts) >= 2
                    acc = parts[2]
                    for locus in essential_loci
                        if occursin(locus, line)
                            push!(accessions, String(acc))
                            locus_to_acc[String(acc)] = locus
                            break
                        end
                    end
                end
            end
        end
    end

    return accessions, locus_to_acc
end

end
