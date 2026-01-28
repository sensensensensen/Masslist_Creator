using Test
using Masslist_Creator

# ==============================================================================
# Helper: Access Internal Modules
# ==============================================================================
# The functions we want to test (parse_formula, calculate_neutral_mass) are 
# located inside: Masslist_Creator -> Masslist_final -> ChemUtils.
# Since ChemUtils is not exported, we access it via the fully qualified path.
const CU = Masslist_Creator.Masslist_final.ChemUtils

@testset "Masslist_Creator Unit Tests" begin

    # ==========================================================================
    # 1. Testing Formula Parsing
    # ==========================================================================
    @testset "ChemUtils: parse_formula" begin
        
        # Test Case 1: Standard simple molecule (Water)
        # Expected: Dictionary with H=2, O=1
        f1 = CU.parse_formula("H2O")
        @test f1["H"] == 2
        @test f1["O"] == 1

        # Test Case 2: Molecule with dual-letter elements (Salt)
        # [cite_start]Verifies logic: isuppercase(f[i]) && islowercase(f[i+1]) [cite: 20]
        f2 = CU.parse_formula("NaCl")
        @test f2["Na"] == 1
        @test f2["Cl"] == 1
        @test !haskey(f2, "N") # Ensure 'Na' wasn't split into N and a

        # Test Case 3: Isotope handling with parentheses
        # [cite_start]Verifies logic: C(13) is converted to internal key "C13" [cite: 18, 26]
        f3 = CU.parse_formula("C(13)H4")
        @test haskey(f3, "C13")
        @test f3["C13"] == 1
        @test f3["H"] == 4
        
        # Test Case 4: Complex heavy isotope
        # [cite_start]Verifies logic: S(34) -> S34 conversion [cite: 29]
        f4 = CU.parse_formula("H2S(34)O4")
        @test f4["S34"] == 1
        @test f4["O"] == 4

        # Test Case 5: Single atom without number (Methane)
        # [cite_start]Verifies logic: Default count is 1 if no digit follows [cite: 25]
        f5 = CU.parse_formula("CH4")
        @test f5["C"] == 1
        @test f5["H"] == 4
    end

    # ==========================================================================
    # 2. Testing Mass Calculation
    # ==========================================================================
    @testset "ChemUtils: calculate_neutral_mass" begin
        
        # [cite_start]Reference atomic masses from source code [cite: 2]
        # H = 1.0078250, O = 15.9949146, C13 = 13.0033548
        H_mass = 1.0078250
        O_mass = 15.9949146
        C13_mass = 13.0033548

        # Test Case 1: Basic calculation (H2O)
        # Using isapprox (≈) to handle floating point errors
        expected_water = 2 * H_mass + O_mass
        @test CU.calculate_neutral_mass("H2O") ≈ expected_water atol=1e-5

        # Test Case 2: Handling Tuple input
        # [cite_start]Verifies logic: val = f_input isa Tuple ? f_input[1] : f_input [cite: 30]
        # This protects against inputs coming from DataFrame rows or iteration
        @test CU.calculate_neutral_mass(("H2O",)) ≈ expected_water atol=1e-5

        # Test Case 3: Empty string handling
        # [cite_start]Verifies logic: if isempty(f_str); return 0.0; end [cite: 31]
        @test CU.calculate_neutral_mass("") == 0.0

        # Test Case 4: Isotope Mass Calculation
        # Verifies that "C(13)" correctly fetches mass for "C13"
        @test CU.calculate_neutral_mass("C(13)") ≈ C13_mass atol=1e-5
    end

end