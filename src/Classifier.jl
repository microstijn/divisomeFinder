module Classifier

export classify_candidates, evaluate_candidate

using ..Types

"""
    extract_protein_name(entry::Dict) -> String

Extracts the recommended or submitted fullName from the UniProt entry description.
"""
function extract_protein_name(entry)
    p_desc = get(entry, :proteinDescription, nothing)
    name = "Unknown Protein"
    if !isnothing(p_desc)
        if haskey(p_desc, :recommendedName) && haskey(p_desc.recommendedName, :fullName)
            name = p_desc.recommendedName.fullName.value
        elseif haskey(p_desc, :submissionNames) && !isempty(p_desc.submissionNames)
            first_name = p_desc.submissionNames[1]
            if haskey(first_name, :fullName)
                name = first_name.fullName.value
            end
        end
    end
    return name
end

"""
    extract_topology(entry::Dict) -> Tuple{Bool, Bool}

Determines if the protein has Transmembrane and/or Signal features.
Returns `(is_membrane, has_signal)`.
"""
function extract_topology(entry)
    is_membrane = false
    has_signal = false
    if haskey(entry, :features)
        for feat in entry.features
            feat_type = get(feat, :type, "")
            if feat_type == "Transmembrane"
                is_membrane = true
            elseif feat_type == "Signal"
                has_signal = true
            end
        end
    end
    return is_membrane, has_signal
end

"""
    extract_interpro_domains(entry::Dict) -> String

Parses the cross-references for InterPro, Pfam, or SUPFAM properties
and extracts a unique, comma-separated string of domain names.
"""
function extract_interpro_domains(entry)
    domains = String[]
    if haskey(entry, :uniProtKBCrossReferences)
        for xref in entry.uniProtKBCrossReferences
            if get(xref, :database, "") in ["InterPro", "Pfam", "SUPFAM"]
                if haskey(xref, :properties)
                    prop = findfirst(p -> p.key == "EntryName", xref.properties)
                    !isnothing(prop) && push!(domains, xref.properties[prop].value)
                end
            end
        end
    end
    return join(unique(domains), ", ")
end

"""
    evaluate_candidate(
        mass_kda::Float64,
        is_membrane::Bool,
        has_signal::Bool,
        domain_str::String,
        name::String,
        criteria::Any
    ) -> String

Evaluates biological constraints (based on size, structural topology, and active domains)
against a defined `DivisomeCriteria`.
Returns "Scaffold", "Motor", or "None".

# Biological Rationale
- **Scaffold**: The division ring requires massive structural adhesins (≥ `scaffold_min_mass_kda`) packed with specific folds (Ig-like, invasin, etc.). They must localize to the membrane/periplasm (Transmembrane/Signal).
- **Motor**: The force-generating hub. Lacks structural domains but is defined by nucleotide-binding capabilities (ATPase, GTPase) and a tightly conserved mass window (`motor_min_mass_kda` to `motor_max_mass_kda`).
"""
function evaluate_candidate(mass_kda::Float64, is_membrane::Bool, has_signal::Bool, domain_str::String, name::String, criteria::Any)
    domain_str_lower = lowercase(domain_str)
    name_lower = lowercase(name)
    
    is_scaffold = any(kw -> occursin(kw, domain_str_lower) || occursin(kw, name_lower), criteria.scaffold_keywords) && 
                  (mass_kda >= criteria.scaffold_min_mass_kda) && (is_membrane || has_signal)
                  
    is_motor = any(kw -> occursin(kw, domain_str_lower) || occursin(kw, name_lower), criteria.motor_keywords) && 
               (criteria.motor_min_mass_kda <= mass_kda <= criteria.motor_max_mass_kda)

    return is_scaffold ? "Scaffold" : (is_motor ? "Motor" : "None")
end

"""
    classify_candidates(json_data::Any, locus_to_acc::Dict{String, String}, criteria::Any) -> Vector{Types.ProteinCandidate}

Iterates through the raw JSON results from the UniProt API, extracts protein properties,
and evaluates each entry against the defined criteria. Returns a vector of valid `ProteinCandidate`s.
"""
function classify_candidates(json_data, locus_to_acc::Dict{String, String}, criteria::Any)
    records = Types.ProteinCandidate[]
    
    isnothing(json_data) && return records
    !haskey(json_data, :results) && return records

    for item in json_data.results
        entry = get(item, :to, nothing)
        isnothing(entry) && continue
        
        acc = String(item.from)
        locus_tag = get(locus_to_acc, acc, "Unknown")

        # Name Extraction
        name = extract_protein_name(entry)

        # Physics & Topology
        seq_info = get(entry, :sequence, Dict())
        len_aa = get(seq_info, :length, 0)
        mass_kda = get(seq_info, :molWeight, 0) / 1000.0

        is_membrane, has_signal = extract_topology(entry)

        # InterPro Domains
        domain_str = extract_interpro_domains(entry)
        
        # Classification Logic
        candidate_type = evaluate_candidate(mass_kda, is_membrane, has_signal, domain_str, name, criteria)

        if candidate_type != "None"
            push!(records, Types.ProteinCandidate(
                locus_tag = locus_tag,
                accession = acc,
                candidate_type = candidate_type,
                protein_name = name,
                mass_kda = round(mass_kda, digits=1),
                length_aa = len_aa,
                transmembrane = is_membrane,
                signal_peptide = has_signal,
                domains = domain_str
            ))
        end
    end
    
    # Sort for best presentation: Scaffolds first, then Motors, ordered by Mass
    sort!(records, by = x -> (x.candidate_type, x.mass_kda), rev=true)
    
    return records
end

end
