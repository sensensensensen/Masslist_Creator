module Masslist_Pro_inorganic

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

    const ATOMIC_MASS = Dict(
        "C"     => 12.0000000,
        "H"     => 1.0078250,
        "N"     => 14.0030740,
        "O"     => 15.9949146,
        "S"     => 31.9720710,
        "Si"    => 27.9769265, 
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
        "S" => [("S", 0.9493), ("S33", 0.0076), ("S34", 0.0429)],
        "Si"=> [("Si", 0.922), ("Si29", 0.047), ("Si30", 0.031)]
    )

    # [TRACK 1] Standard Adducts
    const ION_MODES = Dict{String, Vector{NamedTuple{(:name, :mass, :priority), Tuple{String, Float64, Int64}}}}(
        "H+" => [
            (name="[M+H]+",        mass=1.007276,  priority=1),
            (name="[M+H3O]+",      mass=19.01784,  priority=2),
            (name="[M+(H2O)2H]+",  mass=37.02840,  priority=3),
            (name="[M+NH4]+",      mass=18.03383,  priority=4) 
        ],
        "NH4+" => [
            (name="[M+NH4]+",         mass=18.03383, priority=1),
            (name="[M+NH4(NH3)]+",    mass=35.06038, priority=2),
            (name="[M+NH4(H2O)]+",    mass=36.04439, priority=3),
            (name="[M+H]+",           mass=1.007276, priority=4)
        ]
    )

    # [TRACK 2] Preset Background Ions
    const BACKGROUND_TARGETS = Dict{String, Vector{NamedTuple{(:name, :mass), Tuple{String, Float64}}}}(
        "H+" => [
            (name="H3O+",          mass=19.017841),
            (name="N2H+",          mass=29.013424),
            (name="H3O(H2O)+",     mass=37.028406),
            (name="H3O(H2O)2+",    mass=55.038971),
            (name="H3O(H2O)3+",    mass=73.049536)
        ],
        "NH4+" => [
            (name="NH4+",              mass=18.033826),
            (name="H3O+",              mass=19.017841),
            (name="N2H+",              mass=29.013424),
            (name="NH4(NH3)+",         mass=35.060375),
            (name="NH4(H2O)+",         mass=36.044391),
            (name="H3O(H2O)+",         mass=37.028406),
            (name="NH4(NH3)2+",        mass=52.086924),
            (name="NH4(NH3)(H2O)+",    mass=53.070940),
            (name="NH4(H2O)2+",        mass=54.054956),
            (name="H3O(H2O)2+",        mass=55.038971),
            (name="NH4(NH3)3+",        mass=69.113473),
            (name="NH4(NH3)2(H2O)+",   mass=70.097489),
            (name="H3O(H2O)3+",        mass=73.049536)
        ]
    )

    # [EXPORT HELPERS] Offsets for Decomposition
    const ADDUCT_OFFSETS = Dict(
        "[M+H]+"         => Dict("H"=>0),
        "[M+H3O]+"       => Dict("H"=>2, "O"=>1),
        "[M+(H2O)2H]+"   => Dict("H"=>4, "O"=>2),
        "[M+NH4]+"       => Dict("N"=>1, "H"=>3),
        "[M+NH4(NH3)]+"  => Dict("N"=>2, "H"=>6),
        "[M+NH4(H2O)]+"  => Dict("N"=>1, "H"=>5, "O"=>1)
    )
    
    const BG_DECOMP = Dict(
        "NH4+"           => Dict("N"=>1, "H"=>3),
        "H3O+"           => Dict("H"=>2, "O"=>1),
        "N2H+"           => Dict("N"=>2),
        "H3O(H2O)+"      => Dict("H"=>4, "O"=>2),
        "H3O(H2O)2+"     => Dict("H"=>6, "O"=>3),
        "H3O(H2O)3+"     => Dict("H"=>8, "O"=>4),
        "NH4(NH3)+"      => Dict("N"=>2, "H"=>6),
        "NH4(H2O)+"      => Dict("N"=>1, "H"=>5, "O"=>1),
        "NH4(NH3)2+"     => Dict("N"=>3, "H"=>9),
        "NH4(NH3)(H2O)+" => Dict("N"=>2, "H"=>8, "O"=>1),
        "NH4(H2O)2+"     => Dict("N"=>1, "H"=>7, "O"=>2),
        "NH4(NH3)3+"     => Dict("N"=>4, "H"=>12),
        "NH4(NH3)2(H2O)+"=> Dict("N"=>3, "H"=>11, "O"=>1)
    )

    struct FormulaCandidate
        comp::Dict{String, Int}
        mass::Float64
        formula_str::String
    end

    function parse_formula(f::AbstractString)
        d = Dict{String, Int}()
        i = 1
        len = length(f)
        while i <= len
            elem = ""
            if i + 4 <= len && f[i+1] == '('
                end_paren = findnext(')', f, i)
                if end_paren !== nothing; elem = f[i:end_paren]; i = end_paren + 1
                else; i += 1; continue; end
            elseif i + 1 <= len && isuppercase(f[i]) && islowercase(f[i+1])
                elem = f[i:i+1]; i += 2
            elseif isuppercase(f[i])
                elem = string(f[i]); i += 1
            else; i += 1; continue; end
            
            start_num = i
            while i <= len && isdigit(f[i]); i += 1; end
            count = (i > start_num) ? parse(Int, f[start_num:i-1]) : 1
            
            if elem == "C(13)"; elem = "C13"; end
            if elem == "N(15)"; elem = "N15"; end
            if elem == "O(18)"; elem = "O18"; end
            if elem == "S(34)"; elem = "S34"; end
            d[elem] = get(d, elem, 0) + count
        end
        return d
    end

    function calculate_neutral_mass(f_input)
        val = f_input isa Tuple ? f_input[1] : f_input
        f_str = strip(string(val))
        if isempty(f_str); return 0.0; end
        atoms = parse_formula(f_str)
        mass = 0.0
        for (el, count) in atoms
            if haskey(ATOMIC_MASS, el); mass += ATOMIC_MASS[el] * count; end
        end
        return mass
    end

    function generate_formula_db(min_m, max_m)
        db = FormulaCandidate[]
        mC, mH, mN, mO, mS, mSi = ATOMIC_MASS["C"], ATOMIC_MASS["H"], ATOMIC_MASS["N"], ATOMIC_MASS["O"], ATOMIC_MASS["S"], ATOMIC_MASS["Si"]
        
        # [MODIFIED] Allow C=0 for inorganic sulfur species (e.g., H2SO4, SO2)
        for nC in 0:40 
            mass_C = nC * mC
            if mass_C > max_m; break; end
            
            # [MODIFIED] Allow H=0 for oxides like SO2
            for nH in 0:80
                mass_CH = mass_C + nH * mH
                if mass_CH > max_m; break; end
                
                # [CONSTRAINT 1] Empty molecule check
                if nC == 0 && nH == 0 && 0 == 0; # Will check other atoms inside
                    # Placeholder, will rely on inner loops
                end 

                # [CONSTRAINT 2] Organic Rules (Strict)
                if nC > 0
                    hc = nH / nC
                    if !(0.3 <= hc <= 2.5); continue; end
                end
                # [CONSTRAINT 3] Inorganic Rules (Safety)
                if nC == 0
                    if nH > 12; continue; end # Avoid huge H clusters
                end

                for nN in 0:5
                    mass_CHN = mass_CH + nN * mN
                    if mass_CHN > max_m; break; end
                    
                    if nC > 0
                        if (nC - 0.5*nH + 0.5*nN + 1) < 0; continue; end
                        if !(0.0 <= nN/nC <= 1.0); continue; end
                    end

                    for nO in 0:30
                        mass_CHNO = mass_CHN + nO * mO
                        if mass_CHNO > max_m; break; end
                        
                        # [CONSTRAINT 4] Organic O/C Ratio
                        if nC > 0
                            if !(0.0 <= nO / nC <= 3.0); continue; end
                        end

                        for nS in 0:2
                            mass_CHNOS = mass_CHNO + nS * mS
                            if mass_CHNOS > max_m; break; end
                            
                            # [CONSTRAINT 5] Inorganic Safety: Must contain S or N
                            # If C=0, we strictly require S or N to avoid pure H-O clusters (background noise)
                            if nC == 0
                                if nS == 0 && nN == 0
                                    continue 
                                end
                            end
                            
                            for nSi in 0:2
                                mass_final = mass_CHNOS + nSi * mSi
                                if mass_final > max_m; break; end
                                if mass_final < min_m; continue; end

                                # Construct Formula String
                                f_str = ""
                                if nC > 0; f_str *= "C$nC"; end
                                if nH > 0; f_str *= "H$nH"; end
                                if nO > 0; f_str *= "O$nO"; end
                                if nN > 0; f_str *= "N$nN"; end
                                if nS > 0; f_str *= "S$nS"; end
                                if nSi > 0; f_str *= "Si$nSi"; end
                                
                                if !isempty(f_str)
                                    push!(db, FormulaCandidate(Dict("C"=>nC,"H"=>nH,"O"=>nO,"N"=>nN,"S"=>nS,"Si"=>nSi), mass_final, f_str))
                                end
                            end
                        end
                    end
                end
            end
        end
        sort!(db, by = x -> x.mass)
        return db
    end

    function calculate_isotopes(formula_comp::Dict{String, Int})
        peaks = Dict{Float64, Float64}(0.0 => 1.0)
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
    
    peak_table::DataFrame       
    assignment_table::DataFrame 
    
    obs_plot_mz::Observable{Vector{Float64}}
    obs_plot_int::Observable{Vector{Float64}}
    obs_fit_mz::Observable{Vector{Float64}}
    obs_fit_int::Observable{Vector{Float64}}
    obs_resid_int::Observable{Vector{Float64}}
    
    obs_peak_locs::Observable{Vector{Float64}}
    obs_labels::Observable{Vector{Tuple{Float64, Float64, String}}} 
    
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
            Observable(Float64[]), Observable(Float64[]),
            Observable(Float64[]), Observable(Float64[]),
            Observable(Float64[]),
            Observable([NaN]), 
            Observable(Tuple{Float64, Float64, String}[]),
            Observable("Ready"), Observable("[Log Started]"), Observable(0.0), Observable(false), Observable(false),
            Observable(17.0), Observable(400.0), Observable(10.0), Observable("NH4+"), Observable("Linear"))
    end
end

const STATE = AppState()

# ==============================================================================
# 3. Core Algorithms
# ==============================================================================

function clean_plot_data(data::AbstractVector)
    return [(!ismissing(x) && x > 0) ? Float64(x) : NaN for x in data]
end

function clean_analysis_data(data::AbstractVector)
    return [(!ismissing(x) && x > 0) ? Float64(x) : 0.0 for x in data]
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

function fit_peaks_in_chunk(mz_chunk, i_chunk, known_formulas, ion_mode, threshold)
    chunk_min, chunk_max = mz_chunk[1], mz_chunk[end]
    peaks_found = DataFrame(Mass=Float64[], Intensity=Float64[], Type=String[], RefFormula=String[])
    
    relevant_species = []
    if !isempty(known_formulas) && haskey(ChemUtils.ION_MODES, ion_mode)
        base_adduct = ChemUtils.ION_MODES[ion_mode][1]
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
            fit = curve_fit((x, p) -> A_known * p, mz_chunk, i_chunk, zeros(n_known), lower=zeros(n_known))
            coeffs_known = fit.param
        catch; end
    end
    
    fitted_signal = A_known * coeffs_known
    for (j, val) in enumerate(coeffs_known)
        if val > threshold
            spec, m = relevant_species[j]
            push!(peaks_found, (m, val, "Known", spec.Formula))
        end
    end

    resid = i_chunk .- fitted_signal
    for _ in 1:5
        val, idx = findmax(resid)
        if val < threshold; break; end
        center_mz = mz_chunk[idx]
        model = get_peak_model(center_mz)
        resid .-= model(mz_chunk, [val, center_mz])
        push!(peaks_found, (center_mz, val, "Unknown", ""))
    end

    n_total = nrow(peaks_found)
    if n_total > 0
        A_joint = zeros(n_points, n_total)
        for (j, row) in enumerate(eachrow(peaks_found))
            model_func = get_peak_model(row.Mass)
            mask = abs.(mz_chunk .- row.Mass) .< 0.5
            if any(mask); A_joint[mask, j] .+= model_func(mz_chunk[mask], [1.0, row.Mass]); end
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
        STATE.obs_peak_locs[] = [NaN]
        
        min_mz, max_mz = STATE.param_min_mass[], STATE.param_max_mass[]
        threshold = STATE.param_threshold[]
        ion_mode = STATE.param_ion_mode[]
        
        idx_start = searchsortedfirst(STATE.mass_axis, min_mz)
        idx_end = searchsortedlast(STATE.mass_axis, max_mz)
        if idx_start > idx_end; STATE.obs_status[]="Range Error"; return; end
        
        sub_mz = STATE.mass_axis[idx_start:idx_end]
        sub_int = clean_analysis_data(STATE.spectrum_int[idx_start:idx_end])
        
        all_peaks = DataFrame(Mass=Float64[], Intensity=Float64[], Type=String[], RefFormula=String[])
        total_fit_curve = zeros(length(sub_mz))
        
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
        
        residual = sub_int .- total_fit_curve
        STATE.peak_table = all_peaks
        STATE.obs_fit_mz[] = sub_mz
        STATE.obs_fit_int[] = clean_plot_data(total_fit_curve)
        STATE.obs_resid_int[] = clean_plot_data(residual)
        
        if !isempty(all_peaks)
            STATE.obs_peak_locs[] = all_peaks.Mass
        else
            STATE.obs_peak_locs[] = [NaN]
        end
        
        STATE.obs_status[] = "Fit Done. $(nrow(all_peaks)) peaks."
        STATE.obs_log_history[] *= "\n[Fit] Found $(nrow(all_peaks)) peaks."
        STATE.is_peaks_fitted[] = true
        STATE.obs_progress[] = 1.0
        
    catch e
        STATE.obs_status[] = "Fit Error: $e"
        println(e)
        showerror(stdout, e, catch_backtrace())
    end
end

function task_assign_formulas()
    if !STATE.is_peaks_fitted[]; STATE.obs_status[] = "Fit Peaks First!"; return; end
    try
        STATE.obs_status[] = "Assigning Formulas..."
        peaks = STATE.peak_table
        ion_mode = STATE.param_ion_mode[]
        
        if isempty(STATE.formula_db)
            min_m, max_m = STATE.param_min_mass[], STATE.param_max_mass[]
            STATE.formula_db = ChemUtils.generate_formula_db(min_m-50, max_m)
        end
        
        results = copy(peaks)
        results.AssignedFormula = fill("", nrow(results))
        results.ErrorPPM = fill(NaN, nrow(results))
        results.AdductType = fill("", nrow(results))
        labels_to_plot = Tuple{Float64, Float64, String}[]
        
        bg_list = get(ChemUtils.BACKGROUND_TARGETS, ion_mode, [])
        adduct_list = ChemUtils.ION_MODES[ion_mode]
        db_len = length(STATE.formula_db)
        
        function compare_mixed(a, b)
            val_a = a isa ChemUtils.FormulaCandidate ? a.mass : a
            val_b = b isa ChemUtils.FormulaCandidate ? b.mass : b
            return val_a < val_b
        end

        for i in 1:nrow(results)
            row = results[i, :]
            
            if row.Type == "Known"
                base_adduct = adduct_list[1]
                label_str = "$(row.RefFormula) $(base_adduct.name)"
                results[i, :AssignedFormula] = row.RefFormula
                results[i, :AdductType] = base_adduct.name
                push!(labels_to_plot, (row.Mass, row.Intensity, label_str))
                continue
            end
            
            match_found = false
            
            # --- TRACK 2: Background Match First ---
            best_bg_match = nothing
            best_bg_ppm = 50.0 
            for bg in bg_list
                theoretical = bg.mass
                err = abs(row.Mass - theoretical) / theoretical * 1e6
                if err < best_bg_ppm
                    best_bg_ppm = err
                    best_bg_match = bg
                end
            end
            
            if best_bg_match !== nothing
                results[i, :AssignedFormula] = best_bg_match.name
                results[i, :ErrorPPM] = best_bg_ppm
                results[i, :AdductType] = "Background"
                push!(labels_to_plot, (row.Mass, row.Intensity, best_bg_match.name))
                match_found = true
            end
            
            # --- TRACK 1: DB Search ---
            if !match_found
                for adduct in adduct_list
                    adduct_mass = Float64(adduct.mass)
                    target_neutral_mass = row.Mass - adduct_mass
                    tol = target_neutral_mass * 50e-6 
                    
                    r_start = searchsortedfirst(STATE.formula_db, target_neutral_mass - tol, lt=compare_mixed)
                    r_end = searchsortedlast(STATE.formula_db, target_neutral_mass + tol, lt=compare_mixed)
                    
                    if r_start <= r_end && r_start <= db_len
                        best_local_match = nothing
                        best_local_ppm = 50.0
                        for k in r_start:min(r_end, db_len)
                            cand = STATE.formula_db[k]
                            theoretical_ion = cand.mass + adduct_mass
                            err = abs(row.Mass - theoretical_ion) / theoretical_ion * 1e6
                            if err < best_local_ppm
                                best_local_ppm = err
                                best_local_match = cand
                            end
                        end
                        if best_local_match !== nothing
                            label_str = "$(best_local_match.formula_str) $(adduct.name)"
                            results[i, :AssignedFormula] = best_local_match.formula_str
                            results[i, :ErrorPPM] = best_local_ppm
                            results[i, :AdductType] = adduct.name
                            push!(labels_to_plot, (row.Mass, row.Intensity, label_str))
                            match_found = true
                            break 
                        end
                    end
                end
            end
            
            if !match_found
                results[i, :AssignedFormula] = "?"
            end
        end
        
        # [NEW] Peak Cleaning (De-splitting)
        sort!(results, :Mass)
        rows_to_keep = trues(nrow(results))
        merge_window = 0.02 # Da
        
        i = 1
        while i < nrow(results)
            cluster_idxs = [i]
            j = i + 1
            while j <= nrow(results) && (results.Mass[j] - results.Mass[j-1] < merge_window)
                push!(cluster_idxs, j)
                j += 1
            end
            
            if length(cluster_idxs) > 1
                has_formula = [results.AssignedFormula[idx] != "?" for idx in cluster_idxs]
                if any(has_formula)
                    valid_idxs = cluster_idxs[has_formula]
                    best_idx = valid_idxs[argmax(results.Intensity[valid_idxs])]
                else
                    best_idx = cluster_idxs[argmax(results.Intensity[cluster_idxs])]
                end
                
                for idx in cluster_idxs
                    if idx != best_idx
                        rows_to_keep[idx] = false
                    end
                end
            end
            i = j 
        end
        
        results = results[rows_to_keep, :]
        
        empty!(labels_to_plot)
        final_locs = Float64[]
        for row in eachrow(results)
            push!(final_locs, row.Mass)
            if row.AssignedFormula != "?"
                lbl = row.AdductType == "Background" ? row.AssignedFormula : "$(row.AssignedFormula) $(row.AdductType)"
                push!(labels_to_plot, (row.Mass, row.Intensity, lbl))
            end
        end
        
        STATE.assignment_table = results
        STATE.obs_labels[] = labels_to_plot
        STATE.obs_peak_locs[] = isempty(final_locs) ? [NaN] : final_locs
        
        STATE.obs_status[] = "Assignments Complete (Cleaned)."
        STATE.obs_log_history[] *= "\n[Assign] Matches found & cleaned."
        
    catch e
        STATE.obs_status[] = "Assign Error: $e"
        println(e)
        showerror(stdout, e, catch_backtrace())
    end
end

# ==============================================================================
# 5. GUI Construction & Export
# ==============================================================================

function main()
    fig = Figure(size = (1400, 850), fontsize=16)
    
    sidebar = fig[1, 1] = GridLayout(alignmode=Outside(5))
    plot_area = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 1, Fixed(260))
    WIDGET_WIDTH = 240
    
    # ------------------ WIDGETS ------------------
    Label(sidebar[1, 1], "Masslist Pro", font=:bold, fontsize=20, color=:blue, halign=:left, width=WIDGET_WIDTH)
    Label(sidebar[2, 1], "Data Loading", font=:bold, fontsize=16, halign=:left, width=WIDGET_WIDTH)
    btn_load_h5 = Button(sidebar[3, 1], label="1. Load HDF5", buttoncolor=:orange, halign=:left, width=WIDGET_WIDTH)
    btn_load_list = Button(sidebar[4, 1], label="2. Load MassList", buttoncolor=:lightblue, halign=:left, width=WIDGET_WIDTH)
    Label(sidebar[5, 1], "Parameters", font=:bold, fontsize=16, halign=:left, width=WIDGET_WIDTH)
    Label(sidebar[6, 1], "View Control:", halign=:left, fontsize=14, width=WIDGET_WIDTH)
    btn_reset = Button(sidebar[7, 1], label="Reset View", buttoncolor=:lightgray, halign=:left, width=WIDGET_WIDTH)
    Label(sidebar[8, 1], "Y-Axis Scale:", halign=:left, fontsize=14, width=WIDGET_WIDTH)
    menu_scale = Menu(sidebar[9, 1], options=["Linear", "Log10"], default="Linear", width=WIDGET_WIDTH, halign=:left)
    Label(sidebar[10, 1], "Mass Range:", halign=:left, fontsize=14, width=WIDGET_WIDTH)
    tb_min = Textbox(sidebar[11, 1], placeholder="Min (17)", stored_string="17", width=WIDGET_WIDTH, halign=:left)
    tb_max = Textbox(sidebar[12, 1], placeholder="Max (400)", stored_string="400", width=WIDGET_WIDTH, halign=:left)
    Label(sidebar[13, 1], "Ion Mode:", halign=:left, fontsize=14, width=WIDGET_WIDTH)
    menu_ion = Menu(sidebar[14, 1], options=["NH4+", "H+"], default="NH4+", halign=:left, width=WIDGET_WIDTH)
    Label(sidebar[15, 1], "Noise Threshold:", halign=:left, fontsize=14, width=WIDGET_WIDTH)
    tb_thresh = Textbox(sidebar[16, 1], placeholder="10", stored_string="10", halign=:left, width=WIDGET_WIDTH)
    Label(sidebar[17, 1], "Analysis", font=:bold, fontsize=16, halign=:left, width=WIDGET_WIDTH)
    btn_fit = Button(sidebar[18, 1], label="A. Fit Peaks", buttoncolor=:green, height=40, halign=:left, width=WIDGET_WIDTH)
    btn_assign = Button(sidebar[19, 1], label="B. Assign Formulas", buttoncolor=:teal, height=40, halign=:left, width=WIDGET_WIDTH)
    slider_prog = Slider(sidebar[20, 1], range=0:0.01:1, startvalue=0, halign=:left, width=WIDGET_WIDTH)
    btn_export = Button(sidebar[21, 1], label="Export Masslist", buttoncolor=:lightgray, halign=:left, width=WIDGET_WIDTH)
    Label(sidebar[22, 1], "Log:", font=:bold, fontsize=14, halign=:left, width=WIDGET_WIDTH)
    log_area = GridLayout(sidebar[23, 1], height=120, halign=:left, width=WIDGET_WIDTH)
    Box(log_area[1, 1], color=(:whitesmoke, 1.0), strokewidth=1, strokecolor=:gray)
    lbl_log = Label(log_area[1, 1], STATE.obs_log_history, fontsize=10, halign=:left, valign=:top, padding=(5,5,5,5), word_wrap=true)
    colsize!(log_area, 1, Fixed(WIDGET_WIDTH))
    Label(sidebar[24, 1], lift(s->"Status: $s", STATE.obs_status), fontsize=12, color=:gray, halign=:left, width=WIDGET_WIDTH)
    rowgap!(sidebar, 8)
    Label(sidebar[25, 1], "", tellheight=true) 

    # ------------------ PLOT AREA ------------------
    yscale_func = Observable{Function}(identity)
    ax = Axis(plot_area[1, 1], title="Spectrum Analysis", xlabel="m/z", ylabel="Counts", yscale=yscale_func, titlesize=20)
    
    vlines!(ax, STATE.obs_peak_locs, color=(:green, 0.5), linestyle=:dash, linewidth=1.0, label="Peak Centers")
    lines!(ax, STATE.obs_plot_mz, STATE.obs_plot_int, color=:black, linewidth=1.0, label="Raw")
    lines!(ax, STATE.obs_fit_mz, STATE.obs_fit_int, color=(:red, 0.8), linewidth=1.5, label="Fit")
    band!(ax, STATE.obs_fit_mz, 0.0, STATE.obs_resid_int, color=(:gray, 0.3), label="Resid")
    
    text!(ax, 
        lift(d->[x[1] for x in d], STATE.obs_labels),
        lift(d->[x[2] for x in d], STATE.obs_labels),
        text = lift(d->[x[3] for x in d], STATE.obs_labels),
        rotation = pi/2,          
        align = (:left, :center), 
        offset = (5, 0),          
        fontsize=12, color=:blue
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
                STATE.obs_plot_int[] = clean_plot_data(STATE.spectrum_int[1:step:end])
                STATE.is_data_loaded[] = true
                STATE.obs_status[] = "HDF5 Loaded"
                autolimits!(ax)
            catch e; STATE.obs_status[] = "HDF5 Error"; println(e); end
        end
    end
    
    # [NEW] Enhanced MassList Loader
    on(btn_load_list.clicks) do _
        path = pick_file()
        if !isempty(path)
            try
                lines = readlines(path)
                
                # Check for Tofware Format
                is_tofware = false
                data_start_idx = 0
                name_col_idx = 0
                
                for (i, line) in enumerate(lines)
                    if startswith(strip(line), "#Masses:")
                        is_tofware = true
                        if i + 1 <= length(lines)
                            headers = split(lines[i+1], '\t')
                            name_col_idx = findfirst(x -> strip(x) == "Name", headers)
                            data_start_idx = i + 2 
                        end
                        break
                    end
                end
                
                fs = String[]
                nm = Float64[]
                
                if is_tofware && name_col_idx !== nothing
                    # Tofware Parsing
                    for i in data_start_idx:length(lines)
                        line = strip(lines[i])
                        if isempty(line); continue; end
                        cols = split(line, '\t')
                        
                        if length(cols) >= name_col_idx
                            name_val = strip(string(cols[name_col_idx]))
                            if isempty(name_val); continue; end
                            
                            # Parse Name: Split by space
                            parts = split(name_val, ' ')
                            if !isempty(parts)
                                candidate = strip(parts[1])
                                
                                # Filter: Skip pure ions (ending with + or -)
                                if endswith(candidate, "+") || endswith(candidate, "-")
                                    continue
                                end
                                
                                # Parse valid neutral formula
                                mass = ChemUtils.calculate_neutral_mass(candidate)
                                if mass > 0
                                    push!(fs, candidate)
                                    push!(nm, mass)
                                end
                            end
                        end
                    end
                    STATE.obs_status[] = "Tofware List Loaded ($(length(fs)))"
                    STATE.obs_log_history[] *= "\n[List] Tofware format: $(length(fs)) formulas."
                    
                else
                    # Legacy CSV Parsing
                    lines = readlines(path); header_row=1
                    for (i,l) in enumerate(lines); if occursin("formula", lowercase(l)); header_row=i; break; end; end
                    df = CSV.read(path, DataFrame; header=header_row, types=String, silencewarnings=true)
                    col = "Formula" in names(df) ? "Formula" : names(df)[1] 
                    for f in df[!, col]; 
                        s=strip(string(f)); 
                        if !isempty(s) && !occursin("=",s); push!(fs,s); push!(nm, ChemUtils.calculate_neutral_mass(s)); end
                    end
                    STATE.obs_status[] = "CSV List Loaded ($(length(fs)))"
                    STATE.obs_log_history[] *= "\n[List] CSV format: $(length(fs)) formulas."
                end
                
                STATE.mass_list_df = DataFrame(Formula=fs, NeutralMass=nm)
                
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
        if isempty(STATE.assignment_table); return; end
        path = save_file()
        if !isempty(path)
            if !endswith(path, ".csv"); path *= ".csv"; end
            
            written_rows = Set{String}()
            
            open(path, "w") do io
                println(io, "# Elements:")
                println(io, "C\tC(13)\tH\tH+\tN\tO\tO(18)\tS\tS(34)\tI\tSi")
                println(io, "")
                println(io, "12.0\t13.00335\t1.00783\t1.007276\t14.00307\t15.99492\t17.99916\t31.97207\t33.967867\t126.90448\t27.9769265")
                println(io, "")
                println(io, "#Masses:")
                println(io, "C\tC(13)\tH\tH+\tN\tO\tO(18)\tS\tS(34)\tI\tSi\tMass\tName")

                for row in eachrow(STATE.assignment_table)
                    c=0; c13=0; h=0; h_plus=0; n=0; o=0; o18=0; s=0; s34=0; i_el=0; si=0
                    name_str = ""
                    
                    if row.AssignedFormula == "?" || ismissing(row.AssignedFormula)
                        name_str = ""
                    else
                        offsets = Dict{String, Int}()
                        if row.AdductType == "Background"
                            offsets = get(ChemUtils.BG_DECOMP, row.AssignedFormula, Dict{String,Int}())
                            name_str = row.AssignedFormula
                            h_plus = 1
                            
                            c  = get(offsets, "C", 0)
                            h  = get(offsets, "H", 0)
                            n  = get(offsets, "N", 0)
                            o  = get(offsets, "O", 0)
                            s  = get(offsets, "S", 0)
                            si = get(offsets, "Si", 0)
                        else
                            base_atoms = ChemUtils.parse_formula(row.AssignedFormula)
                            add_atoms = get(ChemUtils.ADDUCT_OFFSETS, row.AdductType, Dict{String,Int}())
                            c = get(base_atoms, "C", 0) + get(add_atoms, "C", 0)
                            h = get(base_atoms, "H", 0) + get(add_atoms, "H", 0)
                            n = get(base_atoms, "N", 0) + get(add_atoms, "N", 0)
                            o = get(base_atoms, "O", 0) + get(add_atoms, "O", 0)
                            s = get(base_atoms, "S", 0) + get(add_atoms, "S", 0)
                            si = get(base_atoms, "Si", 0) + get(add_atoms, "Si", 0)
                            c13 = get(base_atoms, "C13", 0)
                            o18 = get(base_atoms, "O18", 0)
                            s34 = get(base_atoms, "S34", 0)
                            h_plus = 1 
                            name_str = "$(row.AssignedFormula) $(row.AdductType)"
                        end
                    end
                    
                    row_str = "$c\t$c13\t$h\t$h_plus\t$n\t$o\t$o18\t$s\t$s34\t$i_el\t$si\t$(row.Mass)\t$name_str"
                    
                    if !(row_str in written_rows)
                        println(io, row_str)
                        push!(written_rows, row_str)
                    end
                end
            end
            STATE.obs_status[] = "Exported Tofware Format"
        end
    end
    
    on(STATE.obs_progress) do p; set_close_to!(slider_prog, p); end

    display(fig)
end

end # module

using .Masslist_Pro_inorganic
Masslist_Pro_inorganic.main()
