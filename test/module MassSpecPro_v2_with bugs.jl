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
# 1. Chemistry & Physics Module (High Robustness)
# ==============================================================================
module ChemUtils
    using Combinatorics

    # -- Exact Masses (SIS Based) --
    # O17 and S36 are EXCLUDED as requested
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

    # -- Adduct Definitions (Strictly Typed) --
    # NO SODIUM (Na+) allowed
    # Using Vector of NamedTuples to prevent type ambiguity
    const ION_MODES = Dict{String, Vector{NamedTuple{(:name, :mass, :priority), Tuple{String, Float64, Int64}}}}(
        "H+" => [
            (name="[M+H]+",        mass=1.007276,  priority=1),
            (name="[M+H3O]+",      mass=19.01784,  priority=2),
            (name="[M+(H2O)2H]+",  mass=37.02840,  priority=3),
            (name="[M+NH4]+",      mass=18.03383,  priority=4)  # Background
        ],
        "NH4+" => [
            (name="[M+NH4]+",         mass=18.03383, priority=1),
            (name="[M+NH4(NH3)]+",    mass=35.06038, priority=2),
            (name="[M+NH4(H2O)]+",    mass=36.04439, priority=3),
            (name="[M+H]+",           mass=1.007276, priority=4)
        ]
    )

    struct FormulaCandidate
        comp::Dict{String, Int}
        mass::Float64
        formula_str::String
    end

    """
    [ROBUST PARSER] Handles complex formulas including isotopes e.g. C10H16N(15)O
    """
    function parse_formula(f::AbstractString)
        d = Dict{String, Int}()
        i = 1
        len = length(f)
        
        while i <= len
            # 1. Detect Element
            elem = ""
            if i + 4 <= len && f[i+1] == '(' # e.g. C(13)
                end_paren = findnext(')', f, i)
                if end_paren !== nothing
                    elem = f[i:end_paren] 
                    i = end_paren + 1
                else
                    i += 1; continue 
                end
            elseif i + 1 <= len && isuppercase(f[i]) && islowercase(f[i+1])
                elem = f[i:i+1]
                i += 2
            elseif isuppercase(f[i])
                elem = string(f[i])
                i += 1
            else
                i += 1
                continue 
            end
            
            # 2. Detect Count
            start_num = i
            while i <= len && isdigit(f[i])
                i += 1
            end
            
            count = 1
            if i > start_num
                count = parse(Int, f[start_num:i-1])
            end
            
            # Map isotopes to internal keys
            if elem == "C(13)"; elem = "C13"; end
            if elem == "N(15)"; elem = "N15"; end
            if elem == "O(18)"; elem = "O18"; end
            if elem == "S(34)"; elem = "S34"; end
            
            d[elem] = get(d, elem, 0) + count
        end
        return d
    end

    """
    Calculate Neutral Mass safely.
    """
    function calculate_neutral_mass(f_input)
        val = f_input isa Tuple ? f_input[1] : f_input
        f_str = strip(string(val))
        if isempty(f_str); return 0.0; end
        
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

        # Optimized loop C:1-40, H:1-80 etc.
        for nC in 1:40
            mass_C = nC * mC
            if mass_C > max_m; break; end
            for nH in 1:80
                mass_CH = mass_C + nH * mH
                if mass_CH > max_m; break; end
                
                # H/C Check: 0.3 - 2.5 (Aromatics to Alkanes)
                hc = nH / nC
                if !(0.3 <= hc <= 2.5); continue; end

                for nN in 0:5
                    mass_CHN = mass_CH + nN * mN
                    if mass_CHN > max_m; break; end
                    
                    # DBE Check
                    if (nC - 0.5*nH + 0.5*nN + 1) < 0; continue; end

                    for nO in 0:30
                        mass_CHNO = mass_CHN + nO * mO
                        if mass_CHNO > max_m; break; end
                        
                        # O/C Check: 0.0 - 3.0 (Allows Isoprene TMB etc)
                        if !(0.0 <= nO / nC <= 3.0); continue; end

                        for nS in 0:2
                            mass_final = mass_CHNO + nS * mS
                            if mass_final > max_m; break; end
                            if mass_final < min_m; continue; end
                            
                            f_str = "C$nC" * "H$nH"
                            if nO > 0; f_str *= "O$nO"; end
                            if nN > 0; f_str *= "N$nN"; end
                            if nS > 0; f_str *= "S$nS"; end
                            
                            push!(db, FormulaCandidate(Dict("C"=>nC,"H"=>nH,"O"=>nO,"N"=>nN,"S"=>nS), mass_final, f_str))
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
    # Data Sources
    hdf5_path::String
    # Use Any for mass_axis/spectrum to allow memory mapped arrays
    mass_axis::Any         
    spectrum_int::Any      
    
    # HDF5 Peak Shape Data
    peak_shape_centers::Vector{Float64}
    peak_shape_data::Matrix{Float64}
    peak_shape_axis::Vector{Float64}
    
    # User Input Data
    mass_list_df::DataFrame            
    formula_db::Vector{ChemUtils.FormulaCandidate} 
    
    # Process Results
    peak_table::DataFrame       # Step 1: Pure Peaks
    assignment_table::DataFrame # Step 2: Formulas
    
    # Plot Observables
    obs_plot_mz::Observable{Vector{Float64}}
    obs_plot_int::Observable{Vector{Float64}}
    obs_fit_mz::Observable{Vector{Float64}}
    obs_fit_int::Observable{Vector{Float64}}
    obs_resid_int::Observable{Vector{Float64}}
    
    # Labels (x, y, text)
    obs_labels::Observable{Vector{Tuple{Float64, Float64, String}}} 
    
    # Status & Params
    obs_status::Observable{String}
    obs_log_history::Observable{String}
    obs_progress::Observable{Float64}
    is_data_loaded::Observable{Bool}
    is_peaks_fitted::Observable{Bool}
    
    param_min_mass::Observable{Float64}
    param_max_mass::Observable{Float64}
    param_threshold::Observable{Float64}
    param_ion_mode::Observable{String}
    param_yscale_str::Observable{String}

    function AppState()
        new("", [], [], Float64[], Matrix{Float64}(undef, 0,0), Float64[],
            DataFrame(), [], DataFrame(), DataFrame(),
            Observable(Float64[]), Observable(Float64[]), # Raw
            Observable(Float64[]), Observable(Float64[]), # Fit
            Observable(Float64[]), # Resid
            Observable(Tuple{Float64, Float64, String}[]),
            Observable("Ready"), Observable("[Log Started]"), Observable(0.0), Observable(false), Observable(false),
            Observable(100.0), Observable(500.0), Observable(10.0), Observable("NH4+"), Observable("Linear"))
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
    # 1. Try HDF5 Shape Interpolation
    if !isempty(STATE.peak_shape_centers)
        idx = searchsortedfirst(STATE.peak_shape_centers, mz_center)
        if idx > length(STATE.peak_shape_centers); idx = length(STATE.peak_shape_centers); end
        if idx < 1; idx = 1; end
        shape_profile = STATE.peak_shape_data[:, idx]
        shape_x = STATE.peak_shape_axis 
        itp = LinearInterpolation(shape_x, shape_profile, extrapolation_bc=0.0)
        return (x, p) -> p[1] .* itp.(x .- p[2])
    end
    # 2. Fallback Gaussian
    R = 10000.0
    sigma = mz_center / (R * 2.355)
    return (x, p) -> p[1] .* exp.(-0.5 .* ((x .- p[2]) ./ sigma).^2)
end

# --- Step 1: Peak Fitting (Chunked + Joint Opt) ---

function fit_peaks_in_chunk(mz_chunk, i_chunk, known_formulas, ion_mode, threshold)
    chunk_min, chunk_max = mz_chunk[1], mz_chunk[end]
    peaks_found = DataFrame(Mass=Float64[], Intensity=Float64[], Type=String[], RefFormula=String[])
    
    # A. Targeted Fit (Knowns)
    relevant_species = []
    if !isempty(known_formulas) && haskey(ChemUtils.ION_MODES, ion_mode)
        # [ROBUSTNESS FIX] Ensure we get the NamedTuple correctly
        adduct_list = ChemUtils.ION_MODES[ion_mode]
        base_adduct = adduct_list[1] # First priority adduct
        
        for row in eachrow(known_formulas)
            ion_mass = row.NeutralMass + base_adduct.mass
            if chunk_min - 2.0 < ion_mass < chunk_max + 1.0
                push!(relevant_species, (row, ion_mass))
            end
        end
    end

    n_points = length(mz_chunk)
    n_known = length(relevant_species)
    A_known = zeros(n_points, n_known)
    
    for (j, (spec, center_mass)) in enumerate(relevant_species)
        atoms = ChemUtils.parse_formula(spec.Formula) 
        isotopes = ChemUtils.calculate_isotopes(atoms)
        model_func = get_peak_model(center_mass)
        for (mass_shift, rel_abund) in isotopes
            iso_center = center_mass + mass_shift
            if chunk_min < iso_center < chunk_max
                mask = abs.(mz_chunk .- iso_center) .< 0.5 
                if any(mask); A_known[mask, j] .+= model_func(mz_chunk[mask], [rel_abund, iso_center]); end
            end
        end
    end

    coeffs_known = zeros(n_known)
    if n_known > 0
        try
            # NNLS via LsqFit with bounds
            fit = curve_fit((x, p) -> A_known * p, mz_chunk, i_chunk, zeros(n_known), lower=zeros(n_known))
            coeffs_known = fit.param
        catch; end
    end
    
    fitted_signal = A_known * coeffs_known
    
    # Filter knowns by threshold
    for (j, val) in enumerate(coeffs_known)
        if val > threshold
            spec, m = relevant_species[j]
            push!(peaks_found, (m, val, "Known", spec.Formula))
        end
    end

    # B. Unknown Search (Residual)
    resid = i_chunk .- fitted_signal
    max_unknowns = 5
    
    for _ in 1:max_unknowns
        val, idx = findmax(resid)
        if val < threshold; break; end
        center_mz = mz_chunk[idx]
        
        model = get_peak_model(center_mz)
        # Greedy subtract
        resid .-= model(mz_chunk, [val, center_mz])
        push!(peaks_found, (center_mz, val, "Unknown", ""))
    end

    # C. Joint Optimization (Refine all heights)
    n_total = nrow(peaks_found)
    if n_total > 0
        A_joint = zeros(n_points, n_total)
        for (j, row) in enumerate(eachrow(peaks_found))
            if row.Type == "Known"
                atoms = ChemUtils.parse_formula(row.RefFormula)
                isotopes = ChemUtils.calculate_isotopes(atoms)
                model_func = get_peak_model(row.Mass)
                for (mass_shift, rel_abund) in isotopes
                    iso_center = row.Mass + mass_shift
                    if chunk_min < iso_center < chunk_max
                        mask = abs.(mz_chunk .- iso_center) .< 0.5 
                        if any(mask); A_joint[mask, j] .+= model_func(mz_chunk[mask], [rel_abund, iso_center]); end
                    end
                end
            else
                model_func = get_peak_model(row.Mass)
                mask = abs.(mz_chunk .- row.Mass) .< 0.5
                if any(mask); A_joint[mask, j] .+= model_func(mz_chunk[mask], [1.0, row.Mass]); end
            end
        end
        
        try
            fit_j = curve_fit((x, p) -> A_joint * p, mz_chunk, i_chunk, peaks_found.Intensity, lower=zeros(n_total))
            final_coeffs = fit_j.param
            final_signal = A_joint * final_coeffs
            peaks_found.Intensity = final_coeffs
            return peaks_found, final_signal
        catch
            return peaks_found, fitted_signal
        end
    end

    return peaks_found, fitted_signal
end

function task_fit_peaks()
    try
        STATE.obs_status[] = "Fitting Peaks..."
        STATE.is_peaks_fitted[] = false
        
        min_mz, max_mz = STATE.param_min_mass[], STATE.param_max_mass[]
        threshold = STATE.param_threshold[]
        ion_mode = STATE.param_ion_mode[]
        
        # Prepare Data View
        idx_start = searchsortedfirst(STATE.mass_axis, min_mz)
        idx_end = searchsortedlast(STATE.mass_axis, max_mz)
        if idx_start > idx_end; STATE.obs_status[]="Range Error"; return; end
        
        sub_mz = STATE.mass_axis[idx_start:idx_end]
        sub_int = STATE.spectrum_int[idx_start:idx_end]
        
        all_peaks = DataFrame(Mass=Float64[], Intensity=Float64[], Type=String[], RefFormula=String[])
        total_fit_curve = zeros(length(sub_mz))
        
        # Chunking (1.0 Th)
        chunk_size = 1.0
        current_mz = min_mz
        total_steps = max_mz - min_mz
        
        known_list = !isempty(STATE.mass_list_df) ? STATE.mass_list_df : DataFrame()
        
        while current_mz < max_mz
            c_start = searchsortedfirst(sub_mz, current_mz)
            c_end = searchsortedlast(sub_mz, current_mz + chunk_size)
            if c_start < c_end
                m_c = sub_mz[c_start:c_end]
                i_c = sub_int[c_start:c_end]
                
                # Only fit if data exists above noise
                if length(m_c) > 5 && maximum(i_c) > threshold * 0.5
                    p_df, sig = fit_peaks_in_chunk(m_c, i_c, known_list, ion_mode, threshold)
                    append!(all_peaks, p_df)
                    total_fit_curve[c_start:c_end] .= sig
                else
                    total_fit_curve[c_start:c_end] .= 0.0
                end
            end
            current_mz += chunk_size
            STATE.obs_progress[] = (current_mz - min_mz) / total_steps
            yield()
        end
        
        # Calculate Residual
        residual = sub_int .- total_fit_curve
        
        # Update Observables
        STATE.peak_table = all_peaks
        STATE.obs_fit_mz[] = sub_mz
        STATE.obs_fit_int[] = total_fit_curve
        STATE.obs_resid_int[] = residual
        
        STATE.obs_status[] = "Fit Done. $(nrow(all_peaks)) peaks."
        STATE.obs_log_history[] *= "\n[Fit] Found $(nrow(all_peaks)) peaks."
        STATE.is_peaks_fitted[] = true
        STATE.obs_progress[] = 1.0
        
    catch e
        STATE.obs_status[] = "Fit Error: $e"
        println(e)
    end
end

# --- Step 2: Formula Assignment ---

function task_assign_formulas()
    if !STATE.is_peaks_fitted[]
        STATE.obs_status[] = "Fit Peaks First!"
        return
    end
    try
        STATE.obs_status[] = "Assigning Formulas..."
        peaks = STATE.peak_table
        ion_mode = STATE.param_ion_mode[]
        
        # DB Gen
        if isempty(STATE.formula_db)
            min_m, max_m = STATE.param_min_mass[], STATE.param_max_mass[]
            STATE.formula_db = ChemUtils.generate_formula_db(min_m-50, max_m)
        end
        
        results = copy(peaks)
        results.AssignedFormula = fill("", nrow(results))
        results.ErrorPPM = fill(NaN, nrow(results))
        labels_to_plot = Tuple{Float64, Float64, String}[]
        
        for i in 1:nrow(results)
            row = results[i, :]
            
            # Knowns
            if row.Type == "Known"
                results[i, :AssignedFormula] = row.RefFormula
                push!(labels_to_plot, (row.Mass, row.Intensity, row.RefFormula))
                continue
            end
            
            # Unknowns: Back-trace
            best_match = nothing
            # Robust Loop over Adducts
            adduct_list = ChemUtils.ION_MODES[ion_mode]
            
            for adduct in adduct_list
                neutral_mass = row.Mass - adduct.mass
                tol = neutral_mass * 50e-6 # 50 ppm
                
                # Binary Search
                r_start = searchsortedfirst(STATE.formula_db, neutral_mass - tol, by=x->x.mass)
                r_end = searchsortedlast(STATE.formula_db, neutral_mass + tol, by=x->x.mass)
                
                if r_start <= r_end
                    candidates = STATE.formula_db[r_start:r_end]
                    sort!(candidates, by=x->abs(x.mass - neutral_mass))
                    best_match = candidates[1]
                    break
                end
            end
            
            if best_match !== nothing
                results[i, :AssignedFormula] = best_match.formula_str
                # Calc PPM error
                err = abs(row.Mass - (best_match.mass + (row.Mass-best_match.mass))) / row.Mass * 1e6
                results[i, :ErrorPPM] = err
                push!(labels_to_plot, (row.Mass, row.Intensity, best_match.formula_str))
            else
                results[i, :AssignedFormula] = "?"
            end
        end
        
        STATE.assignment_table = results
        STATE.obs_labels[] = labels_to_plot
        STATE.obs_status[] = "Assignments Complete."
        STATE.obs_log_history[] *= "\n[Assign] Matches found."
        
    catch e
        STATE.obs_status[] = "Assign Error: $e"
        println(e)
    end
end

# ==============================================================================
# 5. GUI Construction (Strict Layout + Robust)
# ==============================================================================

function main()
    # 1. Define Figure & Main Grid
    fig = Figure(size = (1400, 850), fontsize=18)
    
    # Layout: Fixed 300px sidebar Left, Plot Right
    sidebar = fig[1, 1] = GridLayout()
    plot_area = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 1, Fixed(300))
    
    # ------------------ WIDGETS (Strict Styling) ------------------
    # Using width=260 to ensure padding fits within 300px col
    
    # 1. Header
    Label(sidebar[1, 1], "MassSpec Pro", font=:bold, fontsize=24, color=:blue, halign=:left, tellwidth=false)
    Label(sidebar[2, 1], "Data Loading", font=:bold, fontsize=18, halign=:left, tellwidth=false)
    
    btn_load_h5 = Button(sidebar[3, 1], label="1. Load HDF5", buttoncolor=:orange, halign=:left, width=260, tellwidth=false)
    btn_load_list = Button(sidebar[4, 1], label="2. Load MassList", buttoncolor=:lightblue, halign=:left, width=260, tellwidth=false)
    
    # 2. Parameters
    Label(sidebar[5, 1], "Parameters", font=:bold, fontsize=18, halign=:left, tellwidth=false)
    
    # View Controls: Reset + Log Scale
    Label(sidebar[6, 1], "View Control:", halign=:left, fontsize=14, tellwidth=false)
    btn_reset = Button(sidebar[7, 1], label="Reset View", buttoncolor=:lightgray, halign=:left, width=260, tellwidth=false)
    
    # Nested Grid for Scale to keep aligned
    Label(sidebar[8, 1], "Y-Axis Scale:", halign=:left, fontsize=14, tellwidth=false)
    menu_scale = Menu(sidebar[9, 1], options=["Linear", "Log10"], default="Linear", width=260, halign=:left, tellwidth=false)
    
    # Mass Range (Stacked)
    Label(sidebar[10, 1], "Mass Range:", halign=:left, fontsize=14, tellwidth=false)
    # Stacked vertically for robustness
    tb_min = Textbox(sidebar[11, 1], placeholder="Min (100)", stored_string="100", width=260, halign=:left, tellwidth=false)
    tb_max = Textbox(sidebar[12, 1], placeholder="Max (500)", stored_string="500", width=260, halign=:left, tellwidth=false)
    
    Label(sidebar[13, 1], "Ion Mode:", halign=:left, fontsize=14, tellwidth=false)
    menu_ion = Menu(sidebar[14, 1], options=["NH4+", "H+"], default="NH4+", halign=:left, width=260, tellwidth=false)
    
    Label(sidebar[15, 1], "Noise Threshold:", halign=:left, fontsize=14, tellwidth=false)
    tb_thresh = Textbox(sidebar[16, 1], placeholder="10", stored_string="10", halign=:left, width=260, tellwidth=false)
    
    # Analysis Buttons
    Label(sidebar[17, 1], "Analysis", font=:bold, fontsize=18, halign=:left, tellwidth=false)
    btn_fit = Button(sidebar[18, 1], label="A. Fit Peaks", buttoncolor=:green, height=40, halign=:left, width=260, tellwidth=false)
    btn_assign = Button(sidebar[19, 1], label="B. Assign Formulas", buttoncolor=:teal, height=40, halign=:left, width=260, tellwidth=false)
    
    slider_prog = Slider(sidebar[20, 1], range=0:0.01:1, startvalue=0, halign=:left, width=260, tellwidth=false)
    btn_export = Button(sidebar[21, 1], label="Export Results", buttoncolor=:lightgray, halign=:left, width=260, tellwidth=false)
    
    # Log Box (Robust)
    Label(sidebar[22, 1], "Log:", font=:bold, fontsize=14, halign=:left, tellwidth=false)
    log_area = GridLayout(sidebar[23, 1], height=120, tellwidth=false, halign=:left)
    Box(log_area[1, 1], color=(:whitesmoke, 1.0), strokewidth=1, strokecolor=:gray)
    lbl_log = Label(log_area[1, 1], STATE.obs_log_history, fontsize=10, halign=:left, valign=:top, padding=(5,5,5,5), word_wrap=true)
    colsize!(log_area, 1, Fixed(260)) # Force box width
    
    # Status Line
    Label(sidebar[24, 1], lift(s->"Status: $s", STATE.obs_status), fontsize=12, color=:gray, halign=:left, tellwidth=false)
    
    rowgap!(sidebar, 8)
    Label(sidebar[25, 1], "", tellheight=true) # Spacer
    
    # ------------------ PLOT AREA ------------------
    # Define AXIS before binding events to it
    yscale_func = Observable{Function}(identity)
    ax = Axis(plot_area[1, 1], title="Spectrum Analysis", xlabel="m/z", ylabel="Counts", yscale=yscale_func, titlesize=20)
    
    # Plot Objects
    lines!(ax, STATE.obs_plot_mz, STATE.obs_plot_int, color=:black, linewidth=1.0, label="Raw")
    lines!(ax, STATE.obs_fit_mz, STATE.obs_fit_int, color=(:red, 0.8), linewidth=1.5, label="Fit")
    band!(ax, STATE.obs_fit_mz, 0.0, STATE.obs_resid_int, color=(:gray, 0.3), label="Resid")
    text!(ax, 
        lift(d->[x[1] for x in d], STATE.obs_labels),
        lift(d->[x[2] for x in d], STATE.obs_labels),
        text = lift(d->[x[3] for x in d], STATE.obs_labels),
        align = (:center, :bottom), fontsize=12, color=:blue
    )
    axislegend(ax)

    # ------------------ LOGIC BINDING ------------------
    
    on(btn_reset.clicks) do _
        autolimits!(ax)
    end
    
    on(menu_scale.selection) do val
        yscale_func[] = (val == "Log10") ? log10 : identity
        autolimits!(ax)
    end
    
    on(btn_load_h5.clicks) do _
        path = pick_file()
        if !isempty(path)
            try
                fid = h5open(path, "r")
                STATE.mass_axis = HDF5.readmmap(fid["MassAxis"])
                if haskey(fid, "AvgSpectrum"); STATE.spectrum_int = HDF5.readmmap(fid["AvgSpectrum"])
                elseif haskey(fid, "SumSpecs"); raw = read(fid["SumSpecs"]); STATE.spectrum_int = vec(mean(raw, dims=2))
                else; STATE.obs_status[] = "No Spectrum"; return; end
                
                # Robust HDF5 Peakshape Check
                if haskey(fid, "MassDepPeakshape")
                    try
                        STATE.peak_shape_centers = read(fid["MassDepPeakshapeCenterMasses"])
                        STATE.peak_shape_data = read(fid["MassDepPeakshape"])
                        n = size(STATE.peak_shape_data, 1); STATE.peak_shape_axis = range(-0.5, 0.5, length=n)
                        STATE.obs_log_history[] *= "\n[PeakShape] Loaded."
                    catch e; STATE.obs_log_history[] *= "\n[PeakShape] Error: $e"; end
                else
                    STATE.obs_log_history[] *= "\n[PeakShape] Not Found (Using Gaussian)."
                end
                
                STATE.hdf5_path = path
                step = max(1, length(STATE.mass_axis) ÷ 20000)
                STATE.obs_plot_mz[] = STATE.mass_axis[1:step:end]
                STATE.obs_plot_int[] = clean_negative_values(STATE.spectrum_int[1:step:end])
                STATE.is_data_loaded[] = true
                STATE.obs_status[] = "HDF5 Loaded"
                autolimits!(ax)
            catch e; STATE.obs_status[] = "HDF5 Error"; println(e); end
        end
    end
    
    on(btn_load_list.clicks) do _
        path = pick_file()
        if !isempty(path)
            try
                # Intelligent CSV Reader (Restored)
                lines = readlines(path); header_row=1
                for (i,l) in enumerate(lines); if occursin("formula", lowercase(l)); header_row=i; break; end; end
                
                df = CSV.read(path, DataFrame; header=header_row, types=String, silencewarnings=true)
                col = "Formula" in names(df) ? "Formula" : names(df)[1] 
                
                fs=String[]; nm=Float64[]
                for f in df[!, col]; 
                    s=strip(string(f)); 
                    if !isempty(s) && !occursin("=",s); push!(fs,s); push!(nm, ChemUtils.calculate_neutral_mass(s)); end
                end
                STATE.mass_list_df = DataFrame(Formula=fs, NeutralMass=nm)
                count = length(nm)
                STATE.obs_status[] = "List Loaded ($count)"
                STATE.obs_log_history[] *= "\n[List] $count formulas."
            catch e; STATE.obs_status[] = "List Error"; println(e); end
        end
    end
    
    on(btn_fit.clicks) do _
        STATE.param_min_mass[] = parse(Float64, tb_min.stored_string[])
        STATE.param_max_mass[] = parse(Float64, tb_max.stored_string[])
        STATE.param_threshold[] = parse(Float64, tb_thresh.stored_string[])
        STATE.param_ion_mode[] = menu_ion.selection[]
        empty!(STATE.obs_labels[])
        STATE.obs_progress[] = 0.0
        @async task_fit_peaks()
    end
    
    on(btn_assign.clicks) do _
        @async task_assign_formulas()
    end
    
    on(btn_export.clicks) do _
        if isempty(STATE.peak_table); return; end
        path = save_file()
        if !isempty(path)
            if !endswith(path, ".csv"); path *= ".csv"; end
            CSV.write(path, isempty(STATE.assignment_table) ? STATE.peak_table : STATE.assignment_table, delim='\t')
            STATE.obs_status[] = "Saved"
        end
    end
    
    on(STATE.obs_progress) do p; set_close_to!(slider_prog, p); end

    display(fig)
end

end # module

using .MassSpecPro
MassSpecPro.main()


