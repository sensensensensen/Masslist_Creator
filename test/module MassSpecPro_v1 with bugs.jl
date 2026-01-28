module MassSpecPro

using GLMakie
using NativeFileDialog
using DataFrames
using CSV
using HDF5
using Statistics
using LsqFit
using Combinatorics
using Interpolations
using Printf
using LinearAlgebra

export main

# ==============================================================================
# 1. Chemistry & Physics Module
# ==============================================================================
module ChemUtils
    using Combinatorics

    # -- Exact Masses (SIS Based) --
    const ATOMIC_MASS = Dict(
        "C"     => 12.0000000,
        "H"     => 1.0078250,
        "N"     => 14.0030740,
        "O"     => 15.9949146,
        "S"     => 31.9720710,
        "C13"   => 13.0033548,
        "N15"   => 15.0001089,
        "O18"   => 17.9991604,
        "S33"   => 32.9714585,
        "S34"   => 33.9678668
    )

    const ISOTOPE_ABUNDANCE = Dict(
        "C" => [("C", 0.9893), ("C13", 0.0107)],
        "N" => [("N", 0.99632), ("N15", 0.00368)],
        "O" => [("O", 0.99757), ("O18", 0.00205)],
        "S" => [("S", 0.9493), ("S33", 0.0076), ("S34", 0.0429)]
    )

    # -- Ion Modes (Fixed NamedTuple Syntax) --
    # [FIX] Added 'mass=' and 'priority=' to ALL entries
    const ION_MODES = Dict(
        "H+" => [
            (name="[M+H]+",       mass=1.007276, priority=1),
            (name="[M+H3O]+",     mass=19.01784, priority=2),
            (name="[M+(H2O)2H]+", mass=37.02840, priority=3),
            (name="[M+NH4]+",     mass=18.03383, priority=4)
        ],
        "NH4+" => [
            (name="[M+NH4]+",      mass=18.03383, priority=1),
            (name="[M+NH4(NH3)]+", mass=35.06038, priority=2),
            (name="[M+NH4(H2O)]+", mass=36.04439, priority=3),
            (name="[M+H]+",        mass=1.007276, priority=4)
        ]
    )

    struct FormulaCandidate
        comp::Dict{String, Int}
        mass::Float64
        formula_str::String
    end

    """
    [ENHANCED PARSER]
    Handles standard elements and isotopes (e.g., C13) safely.
    """
    function parse_formula(f::AbstractString)
        d = Dict{String, Int}()
        i = 1
        len = length(f)
        
        while i <= len
            # 1. Check for Two-Letter Elements or Isotopes
            if i + 1 <= len && isuppercase(f[i]) && (islowercase(f[i+1]) || isdigit(f[i+1]))
                # Check for Isotope format "C13" (Upper + Digit + Digit)
                if i + 2 <= len && isuppercase(f[i]) && all(isdigit, f[i+1:i+2])
                    elem = f[i:i+2]
                    i += 3
                elseif islowercase(f[i+1])
                    elem = f[i:i+1]
                    i += 2
                else
                    elem = string(f[i])
                    i += 1
                end
            # 2. Single Letter Elements
            elseif isuppercase(f[i])
                elem = string(f[i])
                i += 1
            else
                i += 1
                continue
            end
            
            # 3. Parse Count
            start_num = i
            while i <= len && isdigit(f[i])
                i += 1
            end
            
            count = 1
            if i > start_num
                count = parse(Int, f[start_num:i-1])
            end
            
            d[elem] = get(d, elem, 0) + count
        end
        return d
    end

    """
    [CRITICAL FIX] Calculate Neutral Mass with Defensive Unpacking.
    Accepts Any to handle Tuples/SubStrings.
    """
    function calculate_neutral_mass(f_input)
        # Defense: Unpack tuple if wrapped
        val = f_input isa Tuple ? f_input[1] : f_input
        f_str = strip(string(val))
        
        atoms = parse_formula(f_str)
        mass = 0.0
        for (el, count) in atoms
            if haskey(ATOMIC_MASS, el)
                mass += ATOMIC_MASS[el] * count
            end
        end
        return mass
    end

    function generate_formula_db(min_m, max_m)
        db = FormulaCandidate[]
        mC, mH, mN, mO, mS = ATOMIC_MASS["C"], ATOMIC_MASS["H"], ATOMIC_MASS["N"], ATOMIC_MASS["O"], ATOMIC_MASS["S"]

        for nC in 1:40
            mass_C = nC * mC
            if mass_C > max_m; break; end
            for nH in 1:80
                mass_CH = mass_C + nH * mH
                if mass_CH > max_m; break; end
                if !(0.3 <= nH / nC <= 2.5); continue; end

                for nN in 0:5
                    mass_CHN = mass_CH + nN * mN
                    if mass_CHN > max_m; break; end
                    if (nC - 0.5*nH + 0.5*nN + 1) < 0; continue; end

                    for nO in 0:30
                        mass_CHNO = mass_CHN + nO * mO
                        if mass_CHNO > max_m; break; end
                        if !(0.0 <= nO / nC <= 3.0); continue; end

                        for nS in 0:2
                            mass_final = mass_CHNO + nS * mS
                            if mass_final > max_m; break; end
                            if mass_final < min_m; continue; end
                            
                            f_str = "C$nC" * "H$nH"
                            if nO > 0; f_str *= "O$nO"; end
                            if nN > 0; f_str *= "N$nN"; end
                            if nS > 0; f_str *= "S$nS"; end
                            
                            push!(db, FormulaCandidate(Dict("C"=>nC, "H"=>nH, "O"=>nO, "N"=>nN, "S"=>nS), mass_final, f_str))
                        end
                    end
                end
            end
        end
        sort!(db, by = x -> x.mass)
        return db
    end

    function calculate_isotopes(formula_comp::Dict{String, Int})
        peaks = Dict{Float64, Float64}()
        peaks[0.0] = 1.0
        threshold = 0.001

        nC = get(formula_comp, "C", 0)
        if nC > 0
            p13 = ISOTOPE_ABUNDANCE["C"][2][2]
            for i in 1:min(nC, 5) 
                prob = binomial(nC, i) * (p13^i) * ((1-p13)^(nC-i))
                if prob > threshold
                    mass_shift = i * (ATOMIC_MASS["C13"] - ATOMIC_MASS["C"])
                    peaks[mass_shift] = get(peaks, mass_shift, 0.0) + prob
                end
            end
        end
        
        nS = get(formula_comp, "S", 0)
        if nS > 0
             p34 = ISOTOPE_ABUNDANCE["S"][3][2]
             if p34 * nS > threshold
                 mass_shift = ATOMIC_MASS["S34"] - ATOMIC_MASS["S"]
                 peaks[mass_shift] = get(peaks, mass_shift, 0.0) + p34 * nS
             end
        end
        
        sorted_peaks = sort(collect(peaks), by=x->x[1])
        max_int = maximum(x->x[2], sorted_peaks)
        return [(p[1], p[2]/max_int) for p in sorted_peaks]
    end
end

# ==============================================================================
# 2. State & Data Structures
# ==============================================================================
mutable struct AppState
    hdf5_path::String
    mass_axis::Any          
    spectrum_int::Any       
    
    peak_shape_centers::Vector{Float64}
    peak_shape_data::Matrix{Float64}
    peak_shape_axis::Vector{Float64}
    
    mass_list_df::DataFrame             
    formula_db::Vector{ChemUtils.FormulaCandidate} 
    fitted_peaks::DataFrame             
    
    obs_plot_mz::Observable{Vector{Float64}}
    obs_plot_int::Observable{Vector{Float64}}
    obs_status::Observable{String}
    obs_progress::Observable{Float64}
    is_data_loaded::Observable{Bool}
    obs_log_history::Observable{String}
    
    param_min_mass::Observable{Float64}
    param_max_mass::Observable{Float64}
    param_threshold::Observable{Float64}
    param_ion_mode::Observable{String}

    function AppState()
        new("", [], [], Float64[], Matrix{Float64}(undef, 0,0), Float64[],
            DataFrame(), [], DataFrame(),
            Observable(Float64[]), Observable(Float64[]), 
            Observable("Ready"), Observable(0.0), Observable(false),
            Observable("[Log Started]\nWaiting for data..."), 
            Observable(100.0), Observable(500.0), Observable(10.0), Observable("NH4+"))
    end
end

const STATE = AppState()

# ==============================================================================
# 3. Core Algorithms (Fit & Search)
# ==============================================================================

function clean_negative_values(data::AbstractVector)
    return [(!ismissing(x) && x >= 0) ? Float64(x) : NaN for x in data]
end

function get_peak_model(mz_center)
    if !isempty(STATE.peak_shape_centers)
        idx = searchsortedfirst(STATE.peak_shape_centers, mz_center)
        if idx > length(STATE.peak_shape_centers); idx = length(STATE.peak_shape_centers); end
        if idx < 1; idx = 1; end
        shape_profile = STATE.peak_shape_data[:, idx]
        shape_x = STATE.peak_shape_axis 
        itp = LinearInterpolation(shape_x, shape_profile, extrapolation_bc=0.0)
        return (x, p) -> p[1] .* itp.(x .- p[2])
    end
    R = 10000.0
    sigma = mz_center / (R * 2.355)
    return (x, p) -> p[1] .* exp.(-0.5 .* ((x .- p[2]) ./ sigma).^2)
end

function fit_known_species(mz_chunk, int_chunk, known_formulas, ion_mode)
    if isempty(known_formulas); return [], zeros(length(mz_chunk)); end
    relevant_species = []
    chunk_min, chunk_max = mz_chunk[1], mz_chunk[end]
    
    if !haskey(ChemUtils.ION_MODES, ion_mode)
        error("Unknown ion mode: $ion_mode")
    end
    adducts = ChemUtils.ION_MODES[ion_mode]
    base_adduct = adducts[1] 
    
    # [FIX] Safe access to NamedTuple
    adduct_mass = base_adduct.mass

    for row in eachrow(known_formulas)
        ion_mass = row.NeutralMass + adduct_mass
        if chunk_min - 2.0 < ion_mass < chunk_max + 1.0
            push!(relevant_species, (row, ion_mass))
        end
    end
    
    if isempty(relevant_species); return [], zeros(length(mz_chunk)); end

    n_points = length(mz_chunk)
    n_species = length(relevant_species)
    A = zeros(n_points, n_species)
    
    for (j, (spec, center_mass)) in enumerate(relevant_species)
        atoms = ChemUtils.parse_formula(spec.Formula) 
        isotopes = ChemUtils.calculate_isotopes(atoms)
        profile = zeros(n_points)
        model_func = get_peak_model(center_mass)
        for (mass_shift, rel_abund) in isotopes
            iso_center = center_mass + mass_shift
            if chunk_min < iso_center < chunk_max
                mask = abs.(mz_chunk .- iso_center) .< 0.5 
                if any(mask)
                    profile[mask] .+= model_func(mz_chunk[mask], [rel_abund, iso_center])
                end
            end
        end
        A[:, j] = profile
    end
    
    try
        fit = curve_fit((x, p) -> A * p, mz_chunk, int_chunk, zeros(n_species), lower=zeros(n_species))
        coeffs = fit.param
        fitted_signal = A * coeffs
        
        results = DataFrame(Mass=Float64[], Formula=String[], Adduct=String[], Intensity=Float64[], Type=String[])
        threshold = STATE.param_threshold[]
        
        for (j, val) in enumerate(coeffs)
            if val > threshold 
                spec, m = relevant_species[j]
                push!(results, (m, spec.Formula, base_adduct.name, val, "Known"))
            end
        end
        return results, fitted_signal
    catch
        return [], zeros(length(mz_chunk))
    end
end

function fit_unknowns(mz_chunk, residual_chunk, ion_mode)
    found_peaks = DataFrame(Mass=Float64[], Formula=String[], Adduct=String[], Intensity=Float64[], Type=String[])
    threshold = STATE.param_threshold[]
    if maximum(residual_chunk) < threshold; return found_peaks; end
    
    current_resid = copy(residual_chunk)
    max_peaks = 5 
    for i in 1:max_peaks
        val, idx = findmax(current_resid)
        if val < threshold; break; end
        center_mz = mz_chunk[idx]
        model = get_peak_model(center_mz)
        best_match = nothing
        adduct_used = ""
        
        for adduct in ChemUtils.ION_MODES[ion_mode]
            neutral_mass = center_mz - adduct.mass
            tol = neutral_mass * 50e-6
            r_start = searchsortedfirst(STATE.formula_db, neutral_mass - tol, by=x->x.mass)
            r_end = searchsortedlast(STATE.formula_db, neutral_mass + tol, by=x->x.mass)
            if r_start <= r_end
                candidates = STATE.formula_db[r_start:r_end]
                sort!(candidates, by=x->abs(x.mass - neutral_mass))
                best_match = candidates[1]
                adduct_used = adduct.name
                break 
            end
        end
        if best_match !== nothing
            push!(found_peaks, (center_mz, best_match.formula_str, adduct_used, val, "Unknown"))
            simulated_peak = model(mz_chunk, [val, center_mz])
            current_resid .-= simulated_peak
        else
            mask = abs.(mz_chunk .- center_mz) .< 0.1
            current_resid[mask] .= 0
        end
    end
    return found_peaks
end

function run_analysis_task()
    try
        STATE.obs_status[] = "Initializing..."
        min_mz = STATE.param_min_mass[]
        max_mz = STATE.param_max_mass[]
        mode = STATE.param_ion_mode[]
        has_masslist = !isempty(STATE.mass_list_df)
        
        if isempty(STATE.mass_axis)
            STATE.obs_status[] = "Error: No HDF5 loaded."
            return
        end
        
        if isempty(STATE.formula_db)
            STATE.obs_status[] = "Generating Formula DB..."
            STATE.formula_db = ChemUtils.generate_formula_db(min_mz-50, max_mz) 
        end
        
        STATE.obs_status[] = "Preparing Data..."
        idx_start = searchsortedfirst(STATE.mass_axis, min_mz)
        idx_end = searchsortedlast(STATE.mass_axis, max_mz)
        
        if idx_start > idx_end
            STATE.obs_status[] = "Error: Invalid Mass Range."
            return
        end
        
        sub_mz = STATE.mass_axis[idx_start:idx_end]
        sub_int = STATE.spectrum_int[idx_start:idx_end]
        final_results = DataFrame(Mass=Float64[], Formula=String[], Adduct=String[], Intensity=Float64[], Type=String[])
        
        chunk_size = 1.0
        current_mz = min_mz
        total_steps = (max_mz - min_mz)
        
        known_list = has_masslist ? STATE.mass_list_df : DataFrame()
        STATE.obs_status[] = "Analyzing..."
        
        while current_mz < max_mz
            c_start = searchsortedfirst(sub_mz, current_mz)
            c_end = searchsortedlast(sub_mz, current_mz + chunk_size)
            if c_start < c_end
                m_chunk = sub_mz[c_start:c_end]
                i_chunk = sub_int[c_start:c_end]
                
                if length(m_chunk) > 5
                    if has_masslist
                        res_known, sig_known = fit_known_species(m_chunk, i_chunk, known_list, mode)
                        append!(final_results, res_known)
                        resid_chunk = i_chunk .- sig_known
                    else
                        resid_chunk = i_chunk
                    end
                    res_unknown = fit_unknowns(m_chunk, resid_chunk, mode)
                    append!(final_results, res_unknown)
                end
            end
            current_mz += chunk_size
            progress_val = (current_mz - min_mz) / total_steps
            STATE.obs_progress[] = progress_val
            yield()
        end
        
        STATE.fitted_peaks = final_results
        STATE.obs_status[] = "Done! Found $(nrow(final_results)) peaks."
        STATE.obs_progress[] = 1.0
        STATE.is_data_loaded[] = true
        
        show_results_window(sub_mz, sub_int, final_results)
    catch e
        err_msg = sprint(showerror, e)
        println("Analysis Crash: $err_msg")
        STATE.obs_status[] = "Error: See Console"
    end
end

function show_results_window(mz, int, results)
    f = Figure(size=(1000, 600))
    ax = Axis(f[1,1], title="Fit Results Overview", xlabel="m/z", ylabel="Intensity", titlesize=20)
    step = max(1, length(mz) ÷ 10000)
    lines!(ax, mz[1:step:end], clean_negative_values(int[1:step:end]), color=:black, label="Raw Data")
    if !isempty(results)
        mask = (results.Mass .>= mz[1]) .& (results.Mass .<= mz[end])
        if any(mask)
            sub_res = results[mask, :]
            stem!(ax, sub_res.Mass, sub_res.Intensity, color=:red, label="Fitted")
        end
    end
    axislegend(ax)
    display(GLMakie.Screen(), f)
end

# ==============================================================================
# 5. GUI Construction
# ==============================================================================

function main()
    fig = Figure(size = (1400, 900), fontsize=16)

    sidebar = GridLayout(fig[1, 1], width=300)
    ax = Axis(fig[1, 2], 
        title="Spectrum View", 
        xlabel="m/z", 
        ylabel="Counts", 
        titlesize=22,
        width=nothing, 
        height=nothing
    )
    lines!(ax, STATE.obs_plot_mz, STATE.obs_plot_int, color=:black)

    colsize!(fig.layout, 1, Fixed(300))
    colgap!(fig.layout, 20)

    # --- Sidebar ---
    Label(sidebar[1, 1], "MassSpec Pro", font=:bold, fontsize=26, color=:blue, halign=:left)
    Label(sidebar[2, 1], "Controls", font=:bold, fontsize=20, halign=:left)
    
    btn_load_h5 = Button(sidebar[3, 1], label="Load HDF5", buttoncolor=:orange, width=nothing)
    btn_load_list = Button(sidebar[4, 1], label="Load List", buttoncolor=:lightblue, width=nothing)
    
    Label(sidebar[5, 1], "Plot Controls:", halign=:left, font=:bold)
    
    ctrl_grid = GridLayout(sidebar[6, 1])
    btn_reset = Button(ctrl_grid[1, 1], label="Reset View", buttoncolor=:lightgray, width=nothing)
    
    scale_grid = GridLayout(sidebar[7, 1])
    Label(scale_grid[1, 1], "Y-Scale:", halign=:left)
    menu_scale = Menu(scale_grid[1, 2], options=["Linear", "Log10"], default="Linear", width=nothing)
    
    Label(sidebar[8, 1], "Mass Range:", halign=:left, font=:bold)
    range_grid = GridLayout(sidebar[9, 1])
    tb_min = Textbox(range_grid[1, 1], placeholder="100", stored_string="100", width=nothing)
    tb_max = Textbox(range_grid[1, 2], placeholder="500", stored_string="500", width=nothing)
    colsize!(range_grid, 1, Relative(0.5))
    
    Label(sidebar[10, 1], "Ion Mode:", halign=:left, font=:bold)
    menu_ion = Menu(sidebar[11, 1], options=["NH4+", "H+"], default="NH4+", width=nothing)
    
    Label(sidebar[12, 1], "Threshold:", halign=:left, font=:bold)
    tb_thresh = Textbox(sidebar[13, 1], placeholder="10", stored_string="10", width=nothing)
    
    btn_run = Button(sidebar[14, 1], label="Run Analysis", buttoncolor=:green, height=50, font=:bold, fontsize=18, width=nothing)
    slider_prog = Slider(sidebar[15, 1], range=0:0.01:1, startvalue=0, width=nothing)
    btn_export = Button(sidebar[16, 1], label="Export Results", buttoncolor=:lightgray, width=nothing)
    
    Label(sidebar[17, 1], lift(s->"Status:\n$s", STATE.obs_status), fontsize=14, color=:gray, width=nothing, word_wrap=true)
    
    # Boxed Log
    Label(sidebar[18, 1], "Log History:", halign=:left, font=:bold, fontsize=12)
    log_box_layout = GridLayout(sidebar[19, 1], height=150)
    Box(log_box_layout[1, 1], color=(:whitesmoke, 1.0), strokewidth=1, strokecolor=:gray)
    Label(log_box_layout[1, 1], STATE.obs_log_history, fontsize=11, color=:black, halign=:left, valign=:top, padding=(8, 8, 8, 8), word_wrap=true)
    
    Box(sidebar[20, 1], color=:transparent)
    rowsize!(sidebar, 20, Auto())
    rowgap!(sidebar, 10)

    # --- Interactions ---
    
    on(btn_load_h5.clicks) do _
        path = pick_file()
        if !isempty(path)
            try
                fid = h5open(path, "r")
                STATE.mass_axis = HDF5.readmmap(fid["MassAxis"])
                if haskey(fid, "AvgSpectrum")
                    STATE.spectrum_int = HDF5.readmmap(fid["AvgSpectrum"])
                elseif haskey(fid, "SumSpecs")
                    raw = read(fid["SumSpecs"]); STATE.spectrum_int = vec(mean(raw, dims=2))
                else
                    STATE.obs_status[] = "No Spectrum"; return
                end
                if haskey(fid, "MassDepPeakshape")
                    try
                        STATE.peak_shape_centers = read(fid["MassDepPeakshapeCenterMasses"])
                        STATE.peak_shape_data = read(fid["MassDepPeakshape"])
                        n_bins = size(STATE.peak_shape_data, 1)
                        STATE.peak_shape_axis = range(-0.5, 0.5, length=n_bins)
                    catch; end
                end
                STATE.hdf5_path = path
                step = max(1, length(STATE.mass_axis) ÷ 20000)
                STATE.obs_plot_mz[] = STATE.mass_axis[1:step:end]
                STATE.obs_plot_int[] = clean_negative_values(STATE.spectrum_int[1:step:end])
                STATE.is_data_loaded[] = true
                fname = basename(path)
                STATE.obs_status[] = "Loaded HDF5: $fname"
                STATE.obs_log_history[] *= "\n[HDF5] $fname"
                autolimits!(ax)
            catch e; STATE.obs_status[] = "Load Error"; println(e); end
        end
    end
    
    on(btn_load_list.clicks) do _
        path = pick_file()
        if !isempty(path)
            try
                # 1. Skip Metadata
                header_row = 1
                lines = readlines(path)
                for (i, line) in enumerate(lines)
                    if occursin("formula", lowercase(line))
                        header_row = i
                        break
                    end
                end
                
                # 2. Force Types = String
                is_tsv = endswith(lowercase(path), ".tsv")
                df = is_tsv ? CSV.read(path, DataFrame; delim='\t', header=header_row, types=String, silencewarnings=true) :
                              CSV.read(path, DataFrame; header=header_row, types=String, silencewarnings=true)
                
                target_col = nothing
                col_names = names(df)
                
                # 3. Detect Column
                f_idx = findfirst(c -> occursin("formula", lowercase(string(c))), col_names)
                if f_idx !== nothing; target_col = col_names[f_idx]
                else
                    best_score = 0
                    for col in col_names
                        sample = collect(skipmissing(df[!, col]))
                        if isempty(sample); continue; end
                        first_val = string(sample[1])
                        if occursin("InChI", first_val) || occursin("=", first_val) || occursin("/", first_val); continue; end
                        score = 0
                        for i in 1:min(5, length(sample))
                            val = string(sample[i])
                            if occursin(r"[C|H]", val) && occursin(r"\d", val) && !occursin(r"[=/]", val); score += 1; end
                        end
                        if score > best_score; best_score = score; target_col = col; end
                    end
                    if target_col === nothing; target_col = col_names[1]; end
                end
                
                nm_list = Float64[]
                raw_formulas = String[]
                col_data = df[!, target_col]
                
                # 4. Strict Tuple Unpacking
                for i in 1:length(col_data)
                    raw_val = col_data[i]
                    if ismissing(raw_val); continue; end
                    
                    # [SAFE UNPACK]
                    val_str = raw_val isa Tuple ? string(raw_val[1]) : string(raw_val)
                    f_str = strip(val_str)
                    
                    if isempty(f_str); continue; end
                    if lowercase(f_str) == "formula"; continue; end 
                    # Blacklist check
                    if occursin("=", f_str) || occursin("/", f_str); continue; end
                    
                    mw = ChemUtils.calculate_neutral_mass(f_str)
                    if mw > 0
                        push!(raw_formulas, f_str)
                        push!(nm_list, mw)
                    end
                end
                STATE.mass_list_df = DataFrame(Formula=raw_formulas, NeutralMass=nm_list)
                fname = basename(path)
                count = length(nm_list)
                STATE.obs_status[] = "List Loaded: $fname ($count items)"
                STATE.obs_log_history[] *= "\n[List] $fname ($count items)"
            catch e; STATE.obs_status[] = "List Error"; println("List Error: $e"); end
        end
    end
    
    on(menu_scale.selection) do val
        if STATE.is_data_loaded[]
            ax.yscale = (val == "Log10") ? log10 : identity
            autolimits!(ax)
        end
    end
    
    on(btn_reset.clicks) do _
        autolimits!(ax)
    end
    
    on(btn_run.clicks) do _
        STATE.param_min_mass[] = parse(Float64, tb_min.stored_string[])
        STATE.param_max_mass[] = parse(Float64, tb_max.stored_string[])
        STATE.param_threshold[] = parse(Float64, tb_thresh.stored_string[])
        STATE.param_ion_mode[] = menu_ion.selection[]
        STATE.obs_progress[] = 0.0
        @async begin; run_analysis_task(); end
    end
    
    on(btn_export.clicks) do _
        if isempty(STATE.fitted_peaks); return; end
        path = save_file()
        if !isempty(path)
            if !endswith(path, ".csv"); path *= ".csv"; end
            CSV.write(path, STATE.fitted_peaks, delim='\t')
            STATE.obs_status[] = "Saved: $(basename(path))"
            STATE.obs_log_history[] *= "\n[Export] Saved"
        end
    end
    
    on(STATE.obs_progress) do p
        set_close_to!(slider_prog, p)
    end

    display(fig)
end

end # module

using .MassSpecPro
MassSpecPro.main()
