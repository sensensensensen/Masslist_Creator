module MassSpecApp

using GLMakie
using NativeFileDialog
using DataFrames
using CSV
using XLSX
using HDF5
using Statistics
using Dates

# Export main function so it can be called at the end
export main

# ==============================================================================
# 1. Chemistry Calculation Module
# ==============================================================================
module ChemUtils
    # Comprehensive Isotope Dictionary
    const ATOMIC_MASS = Dict(
        "C"     => 12.00000,
        "C(13)" => 13.00335,
        "H"     => 1.00783,
        "H+"    => 1.007276,
        "N"     => 14.00307,
        "O"     => 15.99491,
        "O(18)" => 17.99916,
        "S"     => 31.97207,
        "S(34)" => 33.967867,
        "Cl"    => 34.96885,
        "Na"    => 22.98977, 
        "K"     => 39.0983,
        "I"     => 126.90448,
        "Si"    => 27.97693,
        "F"     => 18.99840
    )

    # Common Adduct Definitions
    const ADDUCTS = Dict(
        "[M+H]+"   => (mass=1.007276,  label="+ H⁺"),
        "[M+NH4]+" => (mass=18.034374, label="+ NH₄⁺"),
        "[M+Na]+"  => (mass=22.989218, label="+ Na⁺"),
        "[M+K]+"   => (mass=38.963158, label="+ K⁺"),
        "[M-H]-"   => (mass=-1.007276, label="- H⁺")
    )

    """
    Parse chemical formula and calculate Neutral Mass.
    Supports AbstractString to handle CSV.jl's optimized string types.
    """
    function calculate_mw(formula::AbstractString)
        mass = 0.0
        formula_str = String(formula) # Convert to standard String
        
        # Regex match: (Element Symbol)(Optional Count)
        # E.g.: "C" "10", "H" "", "O" "2"
        for m in eachmatch(r"([A-Z][a-z]?(\(\d+\))?)(\d*)", formula_str)
            element = m.captures[1]
            count_str = m.captures[3]
            count = isempty(count_str) ? 1 : parse(Int, count_str)
            
            if haskey(ATOMIC_MASS, element)
                mass += ATOMIC_MASS[element] * count
            end
        end
        return mass
    end
end

using .ChemUtils

# ==============================================================================
# 2. Global State Management (Observables)
# ==============================================================================

mutable struct AppState
    # Data Storage
    mass_list_df::DataFrame         # Stores Formula, NeutralMass
    
    # Plotting Observables (Signals)
    obs_spec_mz::Observable{Vector{Float64}}
    obs_spec_int::Observable{Vector{Float64}}
    obs_targets_mz::Observable{Vector{Float64}}    # Currently visible target lines (m/z)
    
    # Control Parameters
    obs_ion_mode::Observable{String}
    obs_threshold::Observable{Float64}
    obs_yscale_str::Observable{String}             # Tracks "Linear" or "Log10"
    
    # Status Info
    obs_status::Observable{String}
    function AppState()
        new(
            DataFrame(),
            Observable([NaN]), # Init with NaN to prevent plot errors
            Observable([NaN]),
            Observable([NaN]),
            Observable("[M+H]+"),
            Observable(0.0),
            Observable("Linear"), # Default Y-Scale
            Observable("Ready. Please load data.")
        )
    end
end

# Use a const REF or recreated instance inside module
const STATE = AppState()

# ==============================================================================
# 3. Core Algorithms (Cleaning, Loading, Filtering)
# ==============================================================================

"""
Clean negative values.
Replaces all intensity values < 0 with NaN.
"""
function clean_negative_values(data::AbstractVector)
    return [(!ismissing(x) && x >= 0) ? Float64(x) : NaN for x in data]
end

"""
Load HDF5 Spectrum Logic
"""
function load_hdf5_spectrum(path::String)
    try
        h5open(path, "r") do file
            if !haskey(file, "MassAxis")
                STATE.obs_status[] = "Error: MassAxis not found in HDF5."
                return
            end
            
            mz = read(file["MassAxis"])
            
            local intens
            if haskey(file, "AvgSpectrum")
                intens = read(file["AvgSpectrum"])
            elseif haskey(file, "SumSpecs")
                temp = read(file["SumSpecs"])
                intens = vec(mean(temp, dims=2))
            else
                STATE.obs_status[] = "Error: No Spectrum data found."
                return
            end
            
            STATE.obs_spec_mz[] = Float64.(mz)
            STATE.obs_spec_int[] = clean_negative_values(intens)
            STATE.obs_status[] = "HDF5 Loaded: $(basename(path))"
            
            update_visible_targets()
        end
    catch e
        STATE.obs_status[] = "Error loading HDF5: $e"
        println(e)
    end
end

"""
Load Mass List Logic (MCM, CSV, XLSX)
Includes Type Handling and smart column detection.
"""
function load_mass_list(path::String)
    try
        df = DataFrame()
        
        if endswith(lowercase(path), ".xlsx")
            xf = XLSX.readxlsx(path)
            sname = XLSX.sheetnames(xf)[1]
            df = DataFrame(XLSX.gettable(xf[sname]))
        else
            # CSV/MCM: Use stringtype=String to force standard strings
            df = CSV.read(path, DataFrame; comment="#", silencewarnings=true, normalizenames=true, stringtype=String)
        end
        
        col_names = names(df)
        target_col = nothing

        # 1. Strategy A: Explicit header "Formula"
        f_idx = findfirst(c -> occursin("formula", lowercase(string(c))), col_names)
        
        if f_idx !== nothing
            target_col = col_names[f_idx]
        else
            # 2. Strategy B: Smart Content Scan
            best_score = 0
            best_col = nothing
            
            for col in col_names
                sample_rows = collect(skipmissing(df[!, col]))
                if isempty(sample_rows); continue; end
                
                check_count = min(5, length(sample_rows))
                score = 0
                
                for i in 1:check_count
                    val = string(sample_rows[i])
                    # Rule: Contains 'C' or 'H', has Digit, NO '=', '/', '('
                    if !all(isdigit, val) && 
                       (occursin(r"[C|H]", val) && occursin(r"\d", val)) &&
                       !occursin(r"[=/]", val) 
                        score += 1
                    end
                end
                
                if score > best_score
                    best_score = score
                    best_col = col
                end
            end
            
            if best_col !== nothing
                target_col = best_col
                STATE.obs_status[] = "Auto-detected formula column: $target_col"
            else
                target_col = col_names[1]
                STATE.obs_status[] = "Warning: Formula column not found. Using $(target_col)."
            end
        end
        
        # Process Data
        raw_values = df[!, target_col]
        formulas = [strip(string(v)) for v in raw_values if !ismissing(v)]
        
        try
            neutral_masses = ChemUtils.calculate_mw.(formulas)
            
            STATE.mass_list_df = DataFrame(
                Formula = formulas,
                NeutralMass = neutral_masses,
                Label = ["" for _ in formulas],
                MZ = zeros(length(formulas))
            )
            
            recalculate_mz()
            STATE.obs_status[] = "Mass List Loaded: $(length(formulas)) entries."
        catch calc_err
             STATE.obs_status[] = "Error calculating MW. Check format."
             println(calc_err)
        end
        
    catch e
        STATE.obs_status[] = "Error loading Mass List: $e"
        println(e)
    end
end

"""
Recalculate theoretical m/z based on selected Ion Mode
"""
function recalculate_mz()
    if isempty(STATE.mass_list_df)
        return
    end
    
    mode = STATE.obs_ion_mode[]
    adduct = ChemUtils.ADDUCTS[mode]
    
    STATE.mass_list_df.MZ = STATE.mass_list_df.NeutralMass .+ adduct.mass
    STATE.mass_list_df.Label = [
        "$(r.Formula) $(adduct.label)" for r in eachrow(STATE.mass_list_df)
    ]
    
    update_visible_targets()
end

"""
Smart Filter: Update visible red lines based on Noise Threshold.
"""
function update_visible_targets()
    # 1. If list is empty, show NaN
    if isempty(STATE.mass_list_df)
        STATE.obs_targets_mz[] = [NaN]
        return
    end
    
    raw_targets = STATE.mass_list_df.MZ
    spec_mz = STATE.obs_spec_mz[]
    spec_int = STATE.obs_spec_int[]
    thresh = STATE.obs_threshold[]
    
    # 2. If no spectrum loaded or threshold is 0, show all
    if length(spec_mz) <= 1 || isnan(spec_mz[1]) || thresh <= 0
        STATE.obs_targets_mz[] = raw_targets
        return
    end
    
    # 3. Apply Threshold Filter
    visible = Float64[]
    tolerance = 0.05 # Da window
    
    for t in raw_targets
        # Find index range in spectrum
        idx_start = searchsortedfirst(spec_mz, t - tolerance)
        idx_end   = searchsortedfirst(spec_mz, t + tolerance)
        
        if idx_start < length(spec_mz) && idx_end > 1 && idx_start < idx_end
            # Extract window data (ignoring NaNs)
            window = filter(!isnan, spec_int[idx_start:idx_end])
            if !isempty(window)
                # Check if max peak in window exceeds threshold
                if maximum(window) >= thresh
                    push!(visible, t)
                end
            end
        end
    end
    
    STATE.obs_targets_mz[] = isempty(visible) ? [NaN] : visible
end

# ==============================================================================
# 4. Main GUI (GLMakie)
# ==============================================================================

function main()
    # Create Main Window
    fig = Figure(size = (1400, 800), fontsize=18)
    
    # --- Layout Definitions ---
    sidebar = fig[1, 1] = GridLayout()
    plot_area = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 1, Fixed(300))
    
    # --- Plot Area Logic ---
    
    title_str = lift((m, s) -> "Spectrum ($m) - $s", STATE.obs_ion_mode, STATE.obs_status)
    
    # Y-Axis Scale Logic (Typed Observable)
    yscale_func = Observable{Function}(identity)
    on(STATE.obs_yscale_str) do s
        yscale_func[] = (s == "Log10") ? log10 : identity
    end

    ax = Axis(plot_area[1, 1],
        title = title_str,
        xlabel = "m/z",
        ylabel = "Intensity",
        yscale = yscale_func,
        titlesize = 16
    )
    
    lines!(ax, STATE.obs_spec_mz, STATE.obs_spec_int, color=:black, linewidth=1.0)
    vlines!(ax, STATE.obs_targets_mz, color=(:red, 0.6), linestyle=:dash, linewidth=1.5)
    
    # Tooltip
    tooltip_label = Label(plot_area[1, 1], "", 
        tellwidth=false, tellheight=false, 
        halign=:right, valign=:top, padding=(10,10,10,10), 
        color=:blue, font=:bold, fontsize=20
    )
    
    on(events(ax).mouseposition) do _
        targets = STATE.obs_targets_mz[]
        if isempty(targets) || isnan(targets[1]); return; end
        
        mouse_p = mouseposition(ax)
        x_m = mouse_p[1]
        
        min_dist = 0.2
        nearest = nothing
        
        for t in targets
            d = abs(t - x_m)
            if d < min_dist
                min_dist = d
                nearest = t
            end
        end
        
        if nearest !== nothing
            row_idx = findfirst(x -> abs(x - nearest) < 0.001, STATE.mass_list_df.MZ)
            if row_idx !== nothing
                lbl = STATE.mass_list_df[row_idx, :Label]
                tooltip_label.text = "$lbl\nm/z: $(round(nearest, digits=4))"
                tooltip_label.visible = true
            end
        else
            tooltip_label.visible = false
        end
    end

    # --- Sidebar Controls ---
    Label(sidebar[1, 1], "Controls", fontsize=24, font=:bold)
    
    btn_h5 = Button(sidebar[2, 1], label="1. Load HDF5 Spectrum", buttoncolor=:orange)
    btn_list = Button(sidebar[3, 1], label="2. Load Mass List\n(MCM / CSV / Excel)", buttoncolor=:lightblue)
    
    Label(sidebar[4, 1], "3. Ion Mode:", halign=:left)
    menu_ion = Menu(sidebar[5, 1], options=collect(keys(ChemUtils.ADDUCTS)), default="[M+H]+")
    
    Label(sidebar[6, 1], "4. Y-Axis Scale:", halign=:left)
    menu_scale = Menu(sidebar[7, 1], options=["Linear", "Log10"], default="Linear")
    
    # --- UPDATED: Textbox + Button for Threshold ---
    Label(sidebar[8, 1], "5. Noise Threshold:", halign=:left)
    
    # Textbox for input (Stored as string)
    tb_thresh = Textbox(sidebar[9, 1], placeholder="0.0", width=nothing)
    
    # Button to apply
    btn_apply_thresh = Button(sidebar[10, 1], label="Apply Filter", buttoncolor=:lightgreen)
    
    btn_reset = Button(sidebar[11, 1], label="Reset View / Auto Limit", buttoncolor=:lightgray)
    
    # --- Interactions ---
    on(btn_h5.clicks) do _
        path = pick_file()
        if !isempty(path)
            load_hdf5_spectrum(path)
            autolimits!(ax)
        end
    end
    
    on(btn_list.clicks) do _
        path = pick_file()
        if !isempty(path)
            load_mass_list(path)
        end
    end
    
    on(menu_ion.selection) do val
        STATE.obs_ion_mode[] = val
        recalculate_mz()
    end
    
    on(menu_scale.selection) do val
        STATE.obs_yscale_str[] = val
        autolimits!(ax)
    end
    
    # --- UPDATED Interaction: Apply Threshold from Textbox ---
    on(btn_apply_thresh.clicks) do _
        input_str = tb_thresh.stored_string[]
        # Try to parse the input string to Float64
        val = tryparse(Float64, input_str)
        
        if val !== nothing
            STATE.obs_threshold[] = val
            update_visible_targets()
            STATE.obs_status[] = "Threshold set to: $val"
        else
            STATE.obs_status[] = "Invalid threshold input! Please enter a number."
        end
    end
    
    on(btn_reset.clicks) do _
        autolimits!(ax)
    end
    
    rowgap!(sidebar, 10)
    Label(sidebar[12, 1], "") # Filler
    
    display(fig)
end

end # End of Module

# ==============================================================================
# Execute App
# ==============================================================================
using .MassSpecApp
MassSpecApp.main()

