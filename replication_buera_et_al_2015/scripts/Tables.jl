# Produce tables for the report 
# Content
# 1. Steady state calibration table 

# 1. Steady state calibration table 
function fnPrintCalibrationLaTeX(params, endo; variant::Symbol = :grid,
                                 path::String = joinpath("tables", "calibration.tex"))

    @assert variant in (:grid, :interp) "variant must be :grid or :interp"
    @unpack Nᶻ, Nᵃ, a⃗, ψ = params

    # A. Compute all moments (same logic as fnPrintCalibrationElements)
    ωᵃ = dropdims(sum(endo.g, dims=(3, 4)), dims=(3, 4))
    model_r         = endo.rₜ
    t_e             = sum(ωᵃ .* endo.𝐨)
    model_exit      = endo.S / max(t_e, 1e-8)
    ext_finance     = sum(ωᵃ .* max.(endo.𝐤 .- a⃗', 0.0) .* endo.𝐨)
    model_ext_fin   = ext_finance / max(endo.Kᵈ, 1e-8)

    emp_levels      = endo.𝐥[endo.𝐨]
    emp_masses      = ωᵃ[endo.𝐨]
    sort_idx_emp    = sortperm(emp_levels, rev=true)
    sorted_emp      = emp_levels[sort_idx_emp]
    sorted_mass_emp = emp_masses[sort_idx_emp]
    target_ent_mass = 0.10 * sum(sorted_mass_emp)
    total_emp       = sum(sorted_emp .* sorted_mass_emp)
    top10_emp_sum   = 0.0; cum_mass_emp = 0.0
    for i in eachindex(sorted_mass_emp)
        if cum_mass_emp + sorted_mass_emp[i] <= target_ent_mass
            top10_emp_sum += sorted_emp[i] * sorted_mass_emp[i]
            cum_mass_emp  += sorted_mass_emp[i]
        else
            top10_emp_sum += sorted_emp[i] * (target_ent_mass - cum_mass_emp)
            break
        end
    end
    model_top10_emp = total_emp > 0.0 ? top10_emp_sum / total_emp : 0.0

    earnings_levels = zeros(Nᶻ * Nᵃ); earnings_masses = zeros(Nᶻ * Nᵃ)
    idx = 1
    for ia in 1:Nᵃ, iz in 1:Nᶻ
        earnings_levels[idx] = endo.𝐨[iz, ia] ? endo.Π[iz, ia] : endo.wₜ
        earnings_masses[idx] = ωᵃ[iz, ia]; idx += 1
    end
    sort_idx_earn    = sortperm(earnings_levels, rev=true)
    sorted_earn      = earnings_levels[sort_idx_earn]
    sorted_mass_earn = earnings_masses[sort_idx_earn]
    target_pop_mass  = 0.05 * sum(ωᵃ)
    total_earn       = sum(sorted_earn .* sorted_mass_earn)
    top5_earn_sum    = 0.0; cum_mass_earn = 0.0
    for i in eachindex(sorted_mass_earn)
        if cum_mass_earn + sorted_mass_earn[i] <= target_pop_mass
            top5_earn_sum += sorted_earn[i] * sorted_mass_earn[i]
            cum_mass_earn += sorted_mass_earn[i]
        else
            top5_earn_sum += sorted_earn[i] * (target_pop_mass - cum_mass_earn)
            break
        end
    end
    model_top5_earn = total_earn > 0.0 ? top5_earn_sum / total_earn : 0.0

    # B. Row labels, US data, original model, and the new column values
    rows = [
        ("Top 10\\% employment",                              "%.2f", 0.69, 0.69, model_top10_emp),
        ("Top 5\\% earnings share",                           "%.2f", 0.30, 0.30, model_top5_earn),
        ("Establishment exit rate (annual)",                  "%.2f", 0.10, 0.10, model_exit),
        ("Real interest rate (annual)",                       "%.4f", 0.02, 0.02, model_r),
        ("Credit market instruments to non-financial assets", "%.4f", 0.70, 0.70, model_ext_fin),
    ]

    # C. Read existing file (to preserve other column), or initialise with placeholders
    mkpath(dirname(path))
    existing_vals = Dict{String, Tuple{String, String}}()  # label -> (grid_cell, interp_cell)
    if isfile(path)
        for line in eachline(path)
            m = match(r"^(.*?)\s*&\s*([^&]*?)\s*&\s*([^&]*?)\s*&\s*([^&\\]*?)\s*\\\\", line)
            if m !== nothing
                lbl = strip(m.captures[1])
                existing_vals[lbl] = (strip(m.captures[3]), strip(m.captures[4]))
            end
        end
    end

    # D. Assemble body rows
    body = IOBuffer()
    for (label, fmt, us, orig, newval) in rows
        new_cell = Printf.format(Printf.Format(fmt), newval)
        grid_cell, interp_cell = get(existing_vals, label, ("", ""))
        if variant == :grid
            grid_cell = new_cell
        else
            interp_cell = new_cell
        end
        us_cell   = Printf.format(Printf.Format(fmt), us)
        orig_cell = Printf.format(Printf.Format(fmt), orig)
        println(body, "$label & $us_cell & $orig_cell & $grid_cell & $interp_cell \\\\")
    end

    # E. Write the file
    open(path, "w") do io
        println(io, "\\begin{tabular}{l c c c c}")
        println(io, "\\hline")
        println(io, "Moment & US data & Original model & On the grid & Interpolated \\\\")
        println(io, "\\hline")
        print(io, String(take!(body)))
        println(io, "\\hline")
        println(io, "\\end{tabular}")
    end
    println("Saved calibration table ($variant column) to $path")
end