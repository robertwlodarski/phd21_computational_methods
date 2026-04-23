function fnSteadyStateTable(endo, params; path = "tables/ss_table.tex")
    @unpack L, c = params

    # A. Compute SS values from Endo
    f̂   = endo.f
    ŝ   = endo.S / L
    Û   = (L - endo.N) / L
    V̂   = endo.M / endo.q
    θ̂   = V̂ / (L - endo.N)

    # B. Data and original model targets
    rows = [
        ("Job-finding rate, \$f\$",         0.1125, 0.1125, f̂),
        ("Inflow rate, \$s\$",              0.0078, 0.0078, ŝ),
        ("Unemployment rate, \$U/L\$",      0.065,  0.065,  Û),
        ("Vacancy-unemp. ratio, \$\\theta\$", 0.72, 0.72,  θ̂),
        ("Vacancy rate, \$V/L\$",           "—",    "—",    V̂/L),
    ]

    # C. Terminal output
    println("="^65)
    @printf "%-30s %10s %10s %10s\n" "Variable" "Data" "EM (2013)" "Replication"
    println("-"^65)
    for (name, d, o, r) in rows
        ds = d isa String ? d : @sprintf("%.4f", d)
        os = o isa String ? o : @sprintf("%.4f", o)
        @printf "%-30s %10s %10s %10.4f\n" name ds os r
    end
    println("="^65)

    # D. LaTeX output
    open(path, "w") do io
        println(io, "\\begin{tabular}{lccc}")
        println(io, "\\toprule")
        println(io, " & Data & EM (2013) & Replication \\\\")
        println(io, "\\midrule")
        for (name, d, o, r) in rows
            ds = d isa String ? d : @sprintf("%.4f", d)
            os = o isa String ? o : @sprintf("%.4f", o)
            @printf(io, "%s & %s & %s & %.4f \\\\\n", name, ds, os, r)
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
    println("LaTeX table saved to $path")
end

function fnCyclicalElasticityTable(simu, lee, params; λ_hp = 100000, path = "tables/elasticity_table.tex")
    @unpack L = params
    T = length(simu.p⃗̂)

    # A. Construct weekly series
    N⃗  = lee.N⃗
    f⃗  = simu.f⃗
    S⃗  = simu.S⃗
    M⃗  = simu.M⃗
    q⃗  = simu.q⃗
    Y⃗  = simu.Y⃗
    U⃗  = L .- N⃗
    V⃗  = M⃗ ./ q⃗
    θ⃗  = V⃗ ./ U⃗
    s⃗  = S⃗ ./ L
    yₙ = Y⃗ ./ N⃗

    # B. Quarterly averaging (13 weeks)
    Tq = T ÷ 13
    avg(x) = [mean(x[(t-1)*13+1:t*13]) for t in 1:Tq]
    f⃗q  = avg(f⃗)
    s⃗q  = avg(s⃗)
    V⃗q  = avg(V⃗)
    θ⃗q  = avg(θ⃗)
    yₙq = avg(yₙ)

    # C. HP filter log series
    hp(x) = hp_filter(log.(x), λ_hp)[1]
    ỹ  = hp(yₙq)
    f̃  = hp(f⃗q)
    s̃  = hp(s⃗q)
    Ṽ  = hp(V⃗q)
    θ̃  = hp(θ⃗q)

    # D. Elasticities (OLS slope of x̃ on ỹ)
    β(x) = dot(x, ỹ) / dot(ỹ, ỹ)
    βf = β(f̃)
    βs = β(s̃)
    βV = β(Ṽ)
    βθ = β(θ̃)

    # E. Data and original model values
    rows = [
        ("Job-finding rate, \$f\$",             2.65,  2.55,  βf),
        ("Inflow rate, \$s\$",                  -1.89, -1.64, βs),
        ("Vacancies, \$V\$",                    2.91,  2.47,  βV),
        ("Tightness, \$\\theta = V/U\$",        6.44,  6.37,  βθ),
    ]

    # F. Terminal output
    println("="^70)
    @printf "%-30s %10s %10s %12s\n" "Variable" "Data" "EM (2013)" "Replication"
    println("-"^70)
    for (name, d, o, r) in rows
        @printf "%-30s %10.2f %10.2f %12.2f\n" name d o r
    end
    println("="^70)

    # G. LaTeX output
    open(path, "w") do io
        println(io, "\\begin{tabular}{lccc}")
        println(io, "\\toprule")
        println(io, " & Data & EM (2013) & Replication \\\\")
        println(io, "\\midrule")
        for (name, d, o, r) in rows
            @printf(io, "%s & %.2f & %.2f & %.2f \\\\\n", name, d, o, r)
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
    println("LaTeX table saved to $path")
end