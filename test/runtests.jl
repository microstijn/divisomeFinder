using Test

include("../src/divisomeFinder.jl")

using .divisomeFinder
using .divisomeFinder.Types
using .divisomeFinder.Classifier

@testset "DivisomeHunter Tests" begin
    @testset "Types" begin
        # Test defaults
        criteria = Types.DivisomeCriteria()
        @test criteria.scaffold_min_mass_kda == 70.0
        @test criteria.motor_min_mass_kda == 35.0
        @test criteria.motor_max_mass_kda == 65.0
        
        # Test custom initialization
        custom_criteria = Types.DivisomeCriteria(scaffold_min_mass_kda=100.0)
        @test custom_criteria.scaffold_min_mass_kda == 100.0
        
        # Test ProteinCandidate struct
        candidate = Types.ProteinCandidate(
            locus_tag="Gene_01",
            accession="A0A123",
            candidate_type="Scaffold",
            protein_name="Big Adhesin",
            mass_kda=150.5,
            length_aa=1200,
            transmembrane=true,
            signal_peptide=false,
            domains="Ig-like, Invasin"
        )
        @test candidate.locus_tag == "Gene_01"
        @test candidate.mass_kda == 150.5
    end

    @testset "Classifier Logic" begin
        criteria = Types.DivisomeCriteria()
        
        # 1. Valid Scaffold
        @test Classifier.evaluate_candidate(
            120.0, true, false, "ig-like, invasin", "Some Protein", criteria
        ) == "Scaffold"

        # 2. Invalid Scaffold (too small)
        @test Classifier.evaluate_candidate(
            60.0, true, false, "ig-like", "Some Protein", criteria
        ) == "None"

        # 3. Invalid Scaffold (no membrane localization)
        @test Classifier.evaluate_candidate(
            120.0, false, false, "ig-like", "Some Protein", criteria
        ) == "None"

        # 4. Valid Motor
        @test Classifier.evaluate_candidate(
            50.0, false, false, "AAA+ ATPase", "FtsZ-like", criteria
        ) == "Motor"

        # 5. Invalid Motor (too small)
        @test Classifier.evaluate_candidate(
            20.0, false, false, "AAA+ ATPase", "Small Protein", criteria
        ) == "None"

        # 6. Invalid Motor (too large)
        @test Classifier.evaluate_candidate(
            100.0, false, false, "AAA+ ATPase", "Big Protein", criteria
        ) == "None"
    end
end
