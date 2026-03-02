module UniProtAPI

export fetch_uniprot_data

using ..AccessionMapper

"""
    fetch_uniprot_data(accessions::Vector{String}; from_db::String="UniProtKB_AC-ID")

Fetches a stream of UniProt JSON data for a given list of accessions.
This uses the custom `AccessionMapper` functions under the hood.

# Process:
1. Submits an asynchronous mapping job to the UniProt REST API.
2. Waits until the job is "FINISHED".
3. Downloads and parses the resulting GZIP-compressed JSON stream.

# Returns
- A `JSON3.Object` containing the raw results.
"""
function fetch_uniprot_data(accessions::Vector{String}; from_db::String="EMBL-GenBank-DDBJ_CDS")
    # Make sure we remove duplicates
    unique_accs = unique(accessions)

    if isempty(unique_accs)
        println("⚠️  Warning: No accessions provided to UniProt API.")
        return nothing
    end

    println("🚀 Submitting job to UniProt for $(length(unique_accs)) accessions (from_db=$from_db)...")
    job_id = AccessionMapper.submit_job(unique_accs, from_db, "UniProtKB")
    
    results_url = ""
    try
        results_url = AccessionMapper.wait_for_job(job_id)
    catch e
        error("❌ Failed while waiting for UniProt job $job_id: $e")
    end

    json_data = nothing
    try
        json_data = AccessionMapper.download_json_stream(results_url)
    catch e
        error("❌ Failed to download/parse JSON from $results_url: $e")
    end

    println("✅ Data successfully retrieved.")
    return json_data
end

end
