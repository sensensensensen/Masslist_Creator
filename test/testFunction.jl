cd(raw"D:\BaiduSyncdisk\zhensen_working\Work in Europe\Tools\Masslist_Creator")

#
using DataFrames, CSV
using Printf


#####################################################################################################
"""
    extract_mcm_column(filepath::String; col_name::Symbol=:Formula)

Reads an MCM mechanism export file (TSV format), automatically handles citation headers,
and extracts a specific column as a vector of strings.

# Arguments
- `filepath::String`: The path to the MCM .tsv file.
- `col_name::Symbol`: (Optional) The column to extract. Defaults to `:Formula`.
                      Other common options from MCM files include `:Name`, `:Mass`, `:Smiles`, `:Inchi`.

# Returns
- `Vector{String}`: A list of data from the requested column.
"""
function extract_mcm_column(filepath::String; col_name::Symbol=:Formula)
    # 1. Validate file existence
    if !isfile(filepath)
        error("Error: File not found at path \"$filepath\"")
    end

    try
        # 2. Read the data using CSV.jl
        # [cite_start]- delim='\t': MCM files are tab-separated [cite: 1]
        # [cite_start]- comment="*": MCM files utilize '*' for citation headers and comments [cite: 1]
        # - header=1: The first non-comment line contains the column names (Name, Smiles, Formula, etc.)
        df = CSV.read(filepath, DataFrame; delim='\t', comment="*", header=1)

        # 3. Verify the requested column exists in the DataFrame
        if !hasproperty(df, col_name)
            available_cols = names(df)
            error("Error: Column :$col_name not found in file.\nAvailable columns are: $available_cols")
        end

        # 4. Extract the column and convert to String
        # Using string.() ensures type stability even if extracting numeric columns like :Mass
        return string.(df[!, col_name])

    catch e
        # Catch and display file reading errors
        error("An error occurred while reading the MCM file: $e")
    end
end



# ================= Usage Examples =================
filename = raw"D:\BaiduSyncdisk\zhensen_working\Work in Europe\Tools\Masslist_Creator\mcm_export_species.tsv"

# Example 1: Default usage (extracts Formulas)
formulas = extract_mcm_column(filename)
# CSV.write(raw"D:\BaiduSyncdisk\zhensen_working\Work in Europe\Tools\Masslist_Creator\MCM_AP_compounds.csv", DataFrame(Formula = formulas))
# println("--- First 5 Formulas ---")
# println(first(formulas, 5))

# # Example 2: Extracting other columns (e.g., Mass or SMILES)
# masses = extract_mcm_column(filename, col_name=:Mass)
# println("\n--- First 5 Masses ---")
# println(first(masses, 5))




#####################################################################################################
using CSV
using DataFrames
using Printf

# ================= 1. Configuration =================

# Remove 'const' to avoid redefinition warnings
ISOTOPE_MASSES = Dict(
    "C"     => 12.0,
    "C(13)" => 13.00335,
    "H"     => 1.00783,
    "H+"    => 1.007276,
    "N"     => 14.00307,
    "O"     => 15.99492,
    "O(18)" => 17.99916,
    "S"     => 31.97207,
    "S(34)" => 33.967867,
    "I"     => 126.90448,
    "Si"    => 27.9763
)

OUTPUT_COLUMNS = ["C", "C(13)", "H", "H+", "N", "O", "O(18)", "S", "S(34)", "I", "Si"]

ION_COMPOSITION = Dict(
    "H+" => Dict("H+" => 1),
    "NH4+" => Dict("N" => 1, "H" => 3, "H+" => 1)
)

# ================= 2. Helper Functions =================

"""
    parse_formula_to_atoms(formula::AbstractString)

Parses formula string into atom dictionary.
Accepts AbstractString to handle SubString/String15 types.
"""
function parse_formula_to_atoms(formula::AbstractString)
    atoms = Dict{String, Int}()
    # Convert to String for regex safety
    regex = r"([A-Z][a-z]?)(\d*)"
    for m in eachmatch(regex, String(formula)) 
        elem = m.captures[1]
        count_str = m.captures[2]
        count = isempty(count_str) ? 1 : parse(Int, count_str)
        atoms[elem] = get(atoms, elem, 0) + count
    end
    return atoms
end

"""
    standardize_input(input_data)

Converts input (File or Vector) to a standard DataFrame.
"""
function standardize_input(input_data::AbstractString)
    if !isfile(input_data)
        error("Error: Input file '$input_data' not found.")
    end
    df = CSV.read(input_data, DataFrame)
    cols = names(df)
    f_idx = findfirst(c -> occursin("formula", lowercase(string(c))), cols)
    n_idx = findfirst(c -> occursin("name", lowercase(string(c))), cols)
    
    if isnothing(f_idx)
        error("Error: Input file must contain a 'Formula' column.")
    end
    
    std_df = DataFrame()
    # Force conversion to standard String
    std_df.Formula = string.(df[!, f_idx])
    std_df.Name = isnothing(n_idx) ? std_df.Formula : string.(df[!, n_idx])
    return std_df
end

function standardize_input(input_data::Vector{<:AbstractString})
    return DataFrame(Formula = input_data, Name = input_data)
end

# ================= 3. Main Generator Function =================

"""
    generate_masslist(input_data; ion_mode::AbstractString="NH4+")

Generates a mass list DataFrame AND a formatted file string.

# Arguments
- `input_data`: File path (String) OR Vector of formulas.
- `ion_mode`: "NH4+" (Default) or "H+".

# Returns
- `(DataFrame, String)`: 
    1. The DataFrame with calculated data.
    2. The full string content formatted strictly for Tofware (with Headers).
"""
function generate_masslist(input_data; ion_mode::AbstractString="NH4+")
    
    # 1. Validation
    if !haskey(ION_COMPOSITION, ion_mode)
        error("Error: Invalid ion_mode '$ion_mode'. Supported: $(keys(ION_COMPOSITION))")
    end
    
    # 2. Input Processing
    df_input = standardize_input(input_data)
    ion_atoms_to_add = ION_COMPOSITION[ion_mode]
    
    # 3. Calculation Lists
    data_cols = Dict{String, Vector{Float64}}()
    for col in OUTPUT_COLUMNS
        data_cols[col] = Float64[]
    end
    result_masses = Float64[]
    result_names = String[]

    # 4. Process Data
    for row in eachrow(df_input)
        # A. Parse
        current_atoms = parse_formula_to_atoms(row.Formula)
        
        # B. Ionize
        for (elem, count) in ion_atoms_to_add
            current_atoms[elem] = get(current_atoms, elem, 0) + count
        end
        
        # C. Calculate
        total_mass = 0.0
        for col in OUTPUT_COLUMNS
            count = float(get(current_atoms, col, 0))
            push!(data_cols[col], count)
            if haskey(ISOTOPE_MASSES, col)
                total_mass += count * ISOTOPE_MASSES[col]
            end
        end
        
        push!(result_masses, total_mass)
        push!(result_names, "$(row.Name) $ion_mode")
    end
    
    # 5. Build DataFrame Output (Variable 1)
    df_out = DataFrame()
    for col in OUTPUT_COLUMNS
        df_out[!, col] = data_cols[col]
    end
    df_out.Mass = result_masses
    df_out.Name = result_names
    
    # 6. Build Formatted String Output (Variable 2)
    # Using IOBuffer to efficiently build the large string
    buf = IOBuffer()
    
    # -- Header Section --
    println(buf, "# Elements: ")
    println(buf, join(OUTPUT_COLUMNS, "\t"))
    println(buf, join([string(ISOTOPE_MASSES[c]) for c in OUTPUT_COLUMNS], "\t"))
    println(buf, "") # Empty line
    println(buf, "#Masses: ")
    println(buf, join(OUTPUT_COLUMNS, "\t") * "\tMass\tName")
    
    # -- Data Section --
    # Iterate through the generated DataFrame to ensure consistency
    for row in eachrow(df_out)
        # Convert row values to tab-separated string
        # Tip: Use explicit string formatting if specific precision is needed, e.g. @sprintf
        vals = [string(row[col]) for col in OUTPUT_COLUMNS]
        line = join(vals, "\t") * "\t$(row.Mass)\t$(row.Name)"
        println(buf, line)
    end
    
    file_content_str = String(take!(buf))
    
    return df_out, file_content_str
end

# ================= 4. Usage Examples =================

# file_path = "MCM_AP_compounds.csv"
# # Run function
# df_result, file_string = generate_masslist(file_path, ion_mode="NH4+")

# # 1. Use the DataFrame
# println("First 5 rows of DataFrame:")
# println(first(df_result, 5))

# # 2. Use the File String (e.g., write to disk)
# write("output_masslist.csv", file_string)
# println("File string length: ", length(file_string))



###########################################################################################################




















