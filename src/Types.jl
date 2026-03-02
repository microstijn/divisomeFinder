module Types

export DivisomeCriteria, ProteinCandidate

"""
    DivisomeCriteria

A configuration struct containing the biological constraints required to identify
divisome candidates.

# Fields
- `scaffold_min_mass_kda::Float64`: The minimum mass for a structural scaffold protein (default: 70.0).
  Scaffolds need to be large to form division rings or cross-link peptidoglycan.
- `motor_min_mass_kda::Float64`: The minimum mass for a motor protein (default: 35.0).
- `motor_max_mass_kda::Float64`: The maximum mass for a motor protein (default: 65.0).
  Motor proteins typically fall into this specific size bracket.
- `scaffold_keywords::Vector{String}`: Keywords for structural domains (e.g., Ig-like, invasin, DUFs).
- `motor_keywords::Vector{String}`: Keywords for nucleotide-binding and force-generating domains (e.g., ATPase, GTPase, AAA+, P-loop).
"""
Base.@kwdef struct DivisomeCriteria
    scaffold_min_mass_kda::Float64 = 70.0
    motor_min_mass_kda::Float64 = 35.0
    motor_max_mass_kda::Float64 = 65.0
    scaffold_keywords::Vector{String} = ["ig-like", "big_", "invasin", "intimin", "fibronectin", "propeller", "duf", "hypothetical"]
    motor_keywords::Vector{String} = ["gtpase", "atpase", "p-loop", "aaa+"]
end

"""
    ProteinCandidate

Represents a biological protein candidate evaluated for its potential role in the divisome.

# Fields
- `locus_tag::String`: The unique gene identifier (e.g., "Plim_0001").
- `accession::String`: The UniProt Accession ID.
- `candidate_type::String`: The classification result ("Scaffold", "Motor", or "None").
- `protein_name::String`: The descriptive name of the protein.
- `mass_kda::Float64`: Molecular weight in kilodaltons.
- `length_aa::Int`: Sequence length in amino acids.
- `transmembrane::Bool`: True if it possesses transmembrane helices.
- `signal_peptide::Bool`: True if it possesses a signal peptide for secretion/periplasmic localization.
- `domains::String`: A comma-separated string of identified structural or functional domains.
"""
Base.@kwdef struct ProteinCandidate
    locus_tag::String
    accession::String
    candidate_type::String
    protein_name::String
    mass_kda::Float64
    length_aa::Int
    transmembrane::Bool
    signal_peptide::Bool
    domains::String
end

end
