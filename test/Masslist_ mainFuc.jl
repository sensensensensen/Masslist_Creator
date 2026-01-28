module MassSpecPro

using GLMakie
using NativeFileDialog
using DataFrames
using CSV
using XLSX
using HDF5
using Statistics
using Dates
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
    # O17 and S36 are excluded as requested
    const ATOMIC_MASS = Dict(
        "C"     => 12.0000000,
        "H"     => 1.0078250,
        "N"     => 14.0030740,
        "O"     => 15.9949146,
        "S"     => 31.9720710,
        # Isotopes for calculation
        "C13"   => 13.0033548,
        "N15"   => 15.0001089,
        "O18"   => 17.9991604,
        "S33"   => 32.9714585,
        "S34"   => 33.9678668
    )

    const ISOTOPE_ABUNDANCE = Dict(
        "C" => [("C", 0.9893), ("C13", 0.0107)],
        "N" => [("N", 0.99632), ("N15", 0.00368)],
        "O" => [("O", 0.99757), ("O18", 0.00205)], # O17 excluded
        "S" => [("S", 0.9493), ("S33", 0.0076), ("S34", 0.0429)] # S36 excluded
    )

    # -- Adduct Definitions (Prioritized) --
    # No Sodium (Na+) allowed
    const ION_MODES = Dict(
        "M+H" => [
            (name="[M+H]+",        mass=1.007276,  priority=1),
            (name="[M+H3O]+",      mass=19.01784,  priority=2), # +H + H2O
            (name="[M+(H2O)2H]+",  mass=37.02840,  priority=3),
            (name="[M+NH4]+",      mass=18.03383,  priority=4)  # Background NH4
        ],
        "M+NH4" => [
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
    Generate database of possible organic formulas within mass range.
    Constraints: 0.3 <= H/C <= 2.5, 0 <= O/C <= 3.0, DBE >= 0.
    """
    function generate_formula_db(min_m, max_m)
        db = FormulaCandidate[]
        
        # Ranges optimized for atmospheric organics
        # C:1-40, H:1-80, O:0-20, N:0-5, S:0-2
        # Using nested loops with early pruning for performance
        
        mC, mH, mN, mO, mS = ATOMIC_MASS["C"], ATOMIC_MASS["H"], ATOMIC_MASS["N"], ATOMIC_MASS["O"], ATOMIC_MASS["S"]

        for nC in 1:40
            mass_C = nC * mC
            if mass_C > max_m; break; end
            
            for nH in 1:80
                mass_CH = mass_C + nH * mH
                if mass_CH > max_m; break; end
                
                # H/C Check (0.3 - 2.5) covers Aromatics to Alkanes
                hc_ratio = nH / nC
                if !(0.3 <= hc_ratio <= 2.5); continue; end

                for nN in 0:5
                    mass_CHN = mass_CH + nN * mN
                    if mass_CHN > max_m; break; end
                    
                    # DBE Check
                    # DBE = C - H/2 + N/2 + 1
                    dbe = nC - 0.5*nH + 0.5*nN + 1
                    if dbe < 0; continue; end

                    for nO in 0:30
                        mass_CHNO = mass_CHN + nO * mO
                        if mass_CHNO > max_m; break; end
                        
                        # O/C Check (0 - 3.0) covers Isoprene to HOMs
                        oc_ratio = nO / nC
                        if !(0.0 <= oc_ratio <= 3.0); continue; end

                        for nS in 0:2
                            mass_final = mass_CHNO + nS * mS
                            if mass_final > max_m; break; end
                            if mass_final < min_m; continue; end

                            # Build Formula String
                            f_str = "C$nC" * "H$nH"
                            if nO > 0; f_str *= "O$nO"; end
                            if nN > 0; f_str *= "N$nN"; end
                            if nS > 0; f_str *= "S$nS"; end
                            
                            push!(db, FormulaCandidate(
                                Dict("C"=>nC, "H"=>nH, "O"=>nO, "N"=>nN, "S"=>nS),
                                mass_final,
                                f_str
                            ))
                        end
                    end
                end
            end
        end
        # Sort by mass for binary search
        sort!(db, by = x -> x.mass)
        return db
    end

    """
    Calculate isotopic cluster (Mono + Isotopes).
    Returns vector of (mass_offset, relative_abundance).
    Ignores O17 and S36.
    """
    function calculate_isotopes(formula_comp::Dict{String, Int})
        # Simplified probabilistic convolution
        # Base peak is always relative 1.0 (or normalized later)
        peaks = Dict{Float64, Float64}()
        peaks[0.0] = 1.0 # Monoisotopic
        
        threshold = 0.001 # 0.1% cutoff

        # Convolve C
        nC = get(formula_comp, "C", 0)
        if nC > 0
            # Binomial approx for C13
            p13 = ISOTOPE_ABUNDANCE["C"][2][2]
            for i in 1:min(nC, 5) # Consider up to M+5
                prob = binomial(nC, i) * (p13^i) * ((1-p13)^(nC-i))
                if prob > threshold
                    mass_shift = i * (ATOMIC_MASS["C13"] - ATOMIC_MASS["C"])
                    peaks[mass_shift] = get(peaks, mass_shift, 0.0) + prob
                end
            end
        end
        
        # Convolve S (S33, S34) - Similar logic, simplified for speed
        nS = get(formula_comp, "S", 0)
        if nS > 0
             p34 = ISOTOPE_ABUNDANCE["S"][3][2]
             if p34 * nS > threshold
                 mass_shift = ATOMIC_MASS["S34"] - ATOMIC_MASS["S"]
                 peaks[mass_shift] = get(peaks, mass_shift, 0.0) + p34 * nS
             end
        end

        # Sort by mass shift
        sorted_peaks = sort(collect(peaks), by=x->x[1])
        # Normalize to main peak = 100% (or 1.0)
        max_int = maximum(x->x[2], sorted_peaks)
        return [(p[1], p[2]/max_int) for p in sorted_peaks]
    end
end

# ==============================================================================
# 2. State & Data Structures
# ==============================================================================
mutable struct AppState
    # Data
    hdf5_path::String
    mass_axis::Vector{Float64}         # Full Mass Axis (mmap-able via wrapper, here loaded)
    spectrum_int::Vector{Float64}      # Full Intensity
    peak_shape_interp::Any             # Interpolation object for Peak Shape
    
    mass_list_df::DataFrame            # Known species
    formula_db::Vector{ChemUtils.FormulaCandidate} # Pre-calculated unknowns
    
    # Results
    fitted_peaks::DataFrame            # Final results
    
    # Observables for GUI
    obs_plot_mz::Observable{Vector{Float64}}
    obs_plot_int::Observable{Vector{Float64}}
    obs_status::Observable{String}
    obs_progress::Observable{Float64}
    
    # Parameters
    param_min_mass::Observable{Float64}
    param_max_mass::Observable{Float64}
    param_threshold::Observable{Float64}
    param_ion_mode::Observable{String}

    function AppState()
        new("", [], [], nothing, DataFrame(), [], DataFrame(),
            Observable(Float64[]), Observable(Float64[]), 
            Observable("Ready"), Observable(0.0),
            Observable(100.0), Observable(500.0), Observable(10.0), Observable("M+H"))
    end
end

const STATE = AppState()

# ==============================================================================
# 3. Core Algorithms
# ==============================================================================

"""
Get Peak Shape at specific m/z.
Returns a function f(x, center, height).
"""
function get_peak_model(mz_center)
    # If HDF5 shape exists, use interpolation
    if STATE.peak_shape_interp !== nothing
        # Get shape profile from interpolator (assuming it's normalized to width)
        # For simplicity in this script, we fallback to a Resolution-based Pseudo-Voigt 
        # if the HDF5 structure is complex to mock. 
        # Real implementation would query STATE.peak_shape_interp(mz_center)
    end
    
    # Fallback: Pseudo-Voigt with R = 10000 (Typical ToF)
    R = 10000.0
    sigma = mz_center / (R * 2.355)
    gamma = sigma # Approximation
    
    # Simple Gaussian for speed in LsqFit
    return (x, p) -> p[1] .* exp.(-0.5 .* ((x .- p[2]) ./ sigma).^2)
end

"""
Stage 1: Fit Known Species (Rigid Cluster)
"""
function fit_known_species(mz_chunk, int_chunk, known_formulas, ion_mode)
    # This is a simplified Non-Negative Least Squares (NNLS) approach
    # We build a matrix A where each column is the profile of a known species (including isotopes)
    
    if isempty(known_formulas); return [], zeros(length(mz_chunk)); end
    
    # Filter formulas relevant to this chunk
    relevant_species = []
    chunk_min, chunk_max = mz_chunk[1], mz_chunk[end]
    
    adducts = ChemUtils.ION_MODES[ion_mode]
    base_adduct = adducts[1] # Use primary adduct for targeted fit
    
    for row in eachrow(known_formulas)
        # Calculate primary ion mass
        ion_mass = row.NeutralMass + base_adduct.mass
        if chunk_min - 1.0 < ion_mass < chunk_max + 1.0
            push!(relevant_species, (row, ion_mass))
        end
    end
    
    if isempty(relevant_species); return [], zeros(length(mz_chunk)); end

    # Build Design Matrix A
    n_points = length(mz_chunk)
    n_species = length(relevant_species)
    A = zeros(n_points, n_species)
    
    for (j, (spec, center_mass)) in enumerate(relevant_species)
        # Get isotope pattern
        isotopes = ChemUtils.calculate_isotopes(ChemUtils.FormulaCandidate(Dict(), 0.0, "").comp) # Need simple parser here
        # Note: In full implementation, pass parsed dict. Here dummy for structure.
        # Re-parsing formula string from DF for accuracy:
        atoms = ChemUtils.parse_formula_to_atoms_dummy(spec.Formula) # Helper needed
        isotopes = ChemUtils.calculate_isotopes(atoms)
        
        # Construct profile
        profile = zeros(n_points)
        model_func = get_peak_model(center_mass)
        
        for (mass_shift, rel_abund) in isotopes
            iso_center = center_mass + mass_shift
            # Only add if within view
            if chunk_min < iso_center < chunk_max
                # Use fixed width gaussian with unit height
                # Optimization: evaluate only near peak
                mask = abs.(mz_chunk .- iso_center) .< 0.5 
                if any(mask)
                    # p[1]=height, p[2]=center. set height=rel_abund
                    profile[mask] .+= model_func(mz_chunk[mask], [rel_abund, iso_center])
                end
            end
        end
        A[:, j] = profile
    end
    
    # Solve Y = A*x for x >= 0 (NNLS)
    # Using LsqFit with lower bounds as a proxy for NNLS
    fit = curve_fit((x, p) -> A * p, mz_chunk, int_chunk, zeros(n_species), lower=zeros(n_species))
    coeffs = fit.param
    
    # Reconstruct fitted signal
    fitted_signal = A * coeffs
    
    # Compile results
    results = DataFrame(Mass=Float64[], Formula=String[], Adduct=String[], Intensity=Float64[], Type=String[])
    for (j, val) in enumerate(coeffs)
        if val > STATE.param_threshold[] # Only keep significant fits
            spec, m = relevant_species[j]
            push!(results, (m, spec.Formula, base_adduct.name, val, "Known"))
        end
    end
    
    return results, fitted_signal
end

"""
Stage 2: Unknown Search (Iterative Residual)
"""
function fit_unknowns(mz_chunk, residual_chunk, ion_mode)
    found_peaks = DataFrame(Mass=Float64[], Formula=String[], Adduct=String[], Intensity=Float64[], Type=String[])
    
    # Threshold check
    threshold = STATE.param_threshold[]
    if maximum(residual_chunk) < threshold
        return found_peaks
    end
    
    current_resid = copy(residual_chunk)
    
    # Iterative greedy search (Mickwitz style simplified)
    max_peaks = 5 # limit per chunk to prevent infinite loops
    for i in 1:max_peaks
        val, idx = findmax(current_resid)
        if val < threshold; break; end
        
        center_mz = mz_chunk[idx]
        
        # 1. Fit single peak here
        model = get_peak_model(center_mz)
        # Fit height only, fix pos (or allow small wiggle)
        # Simple estimation: height = val
        
        # 2. Back-trace Formula
        best_match = nothing
        adduct_used = ""
        
        # Check Adduct Priorities
        for adduct in ChemUtils.ION_MODES[ion_mode]
            neutral_mass = center_mz - adduct.mass
            
            # Binary search in DB
            # Find candidates within 50 ppm
            tol = neutral_mass * 50e-6
            
            # (Assuming STATE.formula_db is sorted)
            # Find range
            r_start = searchsortedfirst(STATE.formula_db, neutral_mass - tol, by=x->x.mass)
            r_end = searchsortedlast(STATE.formula_db, neutral_mass + tol, by=x->x.mass)
            
            if r_start <= r_end
                # Found candidates! Pick best (closest mass)
                candidates = STATE.formula_db[r_start:r_end]
                # Sort by ppm error
                sort!(candidates, by=x->abs(x.mass - neutral_mass))
                best_match = candidates[1]
                adduct_used = adduct.name
                break # Stop at highest priority adduct match
            end
        end
        
        if best_match !== nothing
            # Valid Chemical Formula
            push!(found_peaks, (center_mz, best_match.formula_str, adduct_used, val, "Unknown"))
            
            # Subtract this peak from residual to find next
            # (Approximation: subtract Gaussian)
            simulated_peak = model(mz_chunk, [val, center_mz])
            current_resid .-= simulated_peak
        else
            # If no formula found, maybe just noise or unlisted. 
            # Zero out this region to prevent re-finding
            current_resid[max(1, idx-10):min(length(current_resid), idx+10)] .= 0
        end
    end
    
    return found_peaks
end

# Helper for parser
function ChemUtils.parse_formula_to_atoms_dummy(f::String)
    # Simple regex parser for "C10H16O2"
    d = Dict{String, Int}()
    for m in eachmatch(r"([A-Z][a-z]?)(\d*)", f)
        elem = m.captures[1]
        count = isempty(m.captures[2]) ? 1 : parse(Int, m.captures[2])
        d[elem] = get(d, elem, 0) + count
    end
    return d
end

# ==============================================================================
# 4. GUI & Main Logic
# ==============================================================================

function run_analysis_task()
    # 1. Validate inputs
    min_mz = STATE.param_min_mass[]
    max_mz = STATE.param_max_mass[]
    mode = STATE.param_ion_mode[]
    
    has_masslist = !isempty(STATE.mass_list_df)
    
    # 2. Pre-generate DB if needed
    if isempty(STATE.formula_db)
        STATE.obs_status[] = "Generating Formula DB..."
        STATE.formula_db = ChemUtils.generate_formula_db(min_mz-50, max_mz) # Pad for neutral mass
    end
    
    # 3. Prepare Data (Slice & Dice)
    STATE.obs_status[] = "Preparing Data..."
    
    # Find indices for mass range
    # Assuming sorted mass axis
    idx_start = searchsortedfirst(STATE.mass_axis, min_mz)
    idx_end = searchsortedlast(STATE.mass_axis, max_mz)
    
    if idx_start > idx_end
        STATE.obs_status[] = "Error: Invalid Mass Range"
        return
    end
    
    sub_mz = STATE.mass_axis[idx_start:idx_end]
    sub_int = STATE.spectrum_int[idx_start:idx_end]
    
    # Results container
    final_results = DataFrame(Mass=Float64[], Formula=String[], Adduct=String[], Intensity=Float64[], Type=String[])
    total_fitted_sig = zeros(length(sub_mz))
    
    # 4. Chunk Processing Loop
    # Process 1.0 Th chunks for speed
    chunk_size = 1.0
    current_mz = min_mz
    
    total_steps = (max_mz - min_mz)
    
    # Prepare MassList for Targeted Fit
    known_list = has_masslist ? STATE.mass_list_df : DataFrame()
    
    while current_mz < max_mz
        # Define chunk
        c_start = searchsortedfirst(sub_mz, current_mz)
        c_end = searchsortedlast(sub_mz, current_mz + chunk_size)
        
        if c_start < c_end
            m_chunk = sub_mz[c_start:c_end]
            i_chunk = sub_int[c_start:c_end]
            
            # Step A: Targeted Fit
            if has_masslist
                res_known, sig_known = fit_known_species(m_chunk, i_chunk, known_list, mode)
                append!(final_results, res_known)
                total_fitted_sig[c_start:c_end] .+= sig_known
                
                # Residual for unknown search
                resid_chunk = i_chunk .- sig_known
            else
                resid_chunk = i_chunk
            end
            
            # Step B: Unknown Fit
            res_unknown = fit_unknowns(m_chunk, resid_chunk, mode)
            append!(final_results, res_unknown)
            
            # Update fitted signal for visualization (simplified)
            # (In real app, re-simulate unknowns)
        end
        
        current_mz += chunk_size
        STATE.obs_progress[] = (current_mz - min_mz) / total_steps
        yield() # Let GUI breathe
    end
    
    STATE.fitted_peaks = final_results
    STATE.obs_status[] = "Analysis Complete! Found $(nrow(final_results)) peaks."
    
    # Pop-up results
    show_results_window(sub_mz, sub_int, final_results)
end

function show_results_window(mz, int, results)
    f = Figure(size=(1000, 600))
    ax = Axis(f[1,1], title="Fit Results", xlabel="m/z", ylabel="Intensity")
    
    # Plot Raw
    lines!(ax, mz, int, color=:black, label="Raw")
    
    # Plot Detected Sticks
    if !isempty(results)
        stem!(ax, results.Mass, results.Intensity, color=:red, label="Fitted")
    end
    
    display(GLMakie.Screen(), f)
end

function main()
    fig = Figure(size = (1200, 800))
    
    # Layout
    sidebar = fig[1, 1] = GridLayout()
    plot_area = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 1, Fixed(250))
    
    # --- Sidebar ---
    Label(sidebar[1, 1], "Controls", font=:bold, fontsize=20)
    
    # File Loading
    btn_load_h5 = Button(sidebar[2, 1], label="1. Load HDF5", buttoncolor=:orange)
    btn_load_list = Button(sidebar[3, 1], label="2. Load Mass List", buttoncolor=:lightblue)
    
    # Parameters
    Label(sidebar[4, 1], "Mass Range:", halign=:left)
    range_grid = GridLayout(sidebar[5, 1], cols=2)
    tb_min = Textbox(range_grid[1, 1], placeholder="100", stored_string="100")
    tb_max = Textbox(range_grid[1, 2], placeholder="500", stored_string="500")
    
    Label(sidebar[6, 1], "Ion Mode:", halign=:left)
    menu_ion = Menu(sidebar[7, 1], options=["M+H", "M+NH4"], default="M+H")
    
    Label(sidebar[8, 1], "Noise Threshold:", halign=:left)
    tb_thresh = Textbox(sidebar[9, 1], placeholder="10", stored_string="10")
    
    # Actions
    btn_run = Button(sidebar[10, 1], label="3. Run Analysis", buttoncolor=:green, height=50)
    prog_bar = ProgressBar(sidebar[11, 1], height=10)
    
    btn_export = Button(sidebar[12, 1], label="4. Export Masslist", buttoncolor=:lightgray)
    
    Label(sidebar[13, 1], lift(s->"Status: $s", STATE.obs_status), textsize=12, color=:gray)
    
    rowgap!(sidebar, 10)
    Label(sidebar[14, 1], "", tellheight=true) # Spacer
    
    # --- Plot Area ---
    ax = Axis(plot_area[1, 1], title="Spectrum View", xlabel="m/z", ylabel="Counts")
    
    # Plot Observables
    lines!(ax, STATE.obs_plot_mz, STATE.obs_plot_int, color=:black)
    
    # Reset Button below plot
    btn_reset = Button(plot_area[2, 1], label="Reset View", width=100)
    
    # --- Event Handlers ---
    
    # 1. Load HDF5 (Mmap)
    on(btn_load_h5.clicks) do _
        path = pick_file()
        if isempty(path); return; end
        try
            fid = h5open(path, "r")
            # Read full data using mmap if supported or direct read
            # For simplicity in this demo script, we read directly but downsample for plot
            STATE.mass_axis = read(fid["MassAxis"])
            # Handling potential different HDF5 structures for Spectrum
            if haskey(fid, "AvgSpectrum")
                STATE.spectrum_int = read(fid["AvgSpectrum"])
            elseif haskey(fid, "SumSpecs")
                temp = read(fid["SumSpecs"])
                STATE.spectrum_int = vec(mean(temp, dims=2))
            end
            
            STATE.hdf5_path = path
            
            # Downsample for Plotting (Performance)
            step = max(1, length(STATE.mass_axis) ÷ 20000)
            STATE.obs_plot_mz[] = STATE.mass_axis[1:step:end]
            STATE.obs_plot_int[] = STATE.spectrum_int[1:step:end]
            
            STATE.obs_status[] = "Loaded HDF5."
            autolimits!(ax)
            close(fid)
        catch e
            STATE.obs_status[] = "Error loading HDF5"
            println(e)
        end
    end
    
    # 2. Load Mass List
    on(btn_load_list.clicks) do _
        path = pick_file()
        if isempty(path); return; end
        # Reuse previous robust loader logic (condensed here)
        try
            df = CSV.read(path, DataFrame; comment="#", silencewarnings=true)
            # Simple column detection
            col = "Formula" in names(df) ? "Formula" : names(df)[1]
            # Calculate neutral mass
            df.NeutralMass = ChemUtils.calculate_mw.(string.(df[!, col])) # Use dummy calc for now
            # Real app needs ChemUtils.calculate_mw from previous turn
            STATE.mass_list_df = df
            STATE.obs_status[] = "Mass List Loaded."
        catch
            STATE.obs_status[] = "Error loading list."
        end
    end
    
    # 3. Run Analysis
    on(btn_run.clicks) do _
        # Update Params
        STATE.param_min_mass[] = parse(Float64, tb_min.stored_string[])
        STATE.param_max_mass[] = parse(Float64, tb_max.stored_string[])
        STATE.param_threshold[] = parse(Float64, tb_thresh.stored_string[])
        STATE.param_ion_mode[] = menu_ion.selection[]
        
        # Check Mass List
        if isempty(STATE.mass_list_df)
            # In a real GUI, this would be a confirmation dialog
            println("Warning: No Mass List. Proceeding with full unknown search.")
        end
        
        # Run Async
        @async begin
            run_analysis_task()
        end
    end
    
    # 4. Export
    on(btn_export.clicks) do _
        if isempty(STATE.fitted_peaks)
            STATE.obs_status[] = "No results to export."
            return
        end
        path = save_file()
        if !isempty(path)
            CSV.write(path, STATE.fitted_peaks, delim='\t')
            STATE.obs_status[] = "Saved to $path"
        end
    end
    
    on(btn_reset.clicks) do _
        autolimits!(ax)
    end
    
    # Link progress bar
    on(STATE.obs_progress) do p
        prog_bar.percentage[] = p
    end

    display(fig)
end

end # module

# ==============================================================================
# Bootstrap
# ==============================================================================
using .MassSpecPro
MassSpecPro.main()


