using JSON
using ArbNumerics
using Plots
using FFMPEG

setprecision(1024)

# Установим Plots backend
default(fmt = :png)

# Парсинг строк вида "val +/- err" и "+- err", assume that the value is 0 in secode case
function parse_arbreal(s::String)
    s = strip(s)

    # Match full form: [value +/- radius]
    full_regex = r"^([-\d.eE+]+)\s*\+/-\s*([\d.eE+-]+)$"
    m = match(full_regex, s)
    if m !== nothing
        val = tryparse(BigFloat, m.captures[1])  # Try parsing as BigFloat
        if val === nothing
            val = parse(BigFloat, m.captures[1])  # If it fails, fall back to standard parsing
        end
        return ArbReal(val)
    end

    # Match radius-only form: [+/- radius]
    radius_only_regex = r"^\+/-\s*([\d.eE+-]+)$"
    m = match(radius_only_regex, s)
    if m !== nothing
        rad = tryparse(BigFloat, m.captures[1])  # Try parsing as BigFloat
        if rad === nothing
            rad = parse(BigFloat, m.captures[1])  # If it fails, fall back to standard parsing
        end
        return ArbReal(BigFloat(0))  # Set the midpoint to 0 for this case
    end

    # Try parsing directly if no matches found
    val = tryparse(BigFloat, s)
    if val !== nothing
        return ArbReal(val)
    end

    error("Failed to parse ArbReal from: '$s'")
end



# Парсинг строк вида "[a +/- b] + [c +/- d]j"
function parse_arbcomplex(s::String)
    regex = r"\[([^\[\]]+?)\]\s*\+\s*\[([^\[\]]+?)\]j"
    m = match(regex, s)
    if m === nothing
        error("Failed to parse ArbComplex from: $s")
    end
    real_val = parse_arbreal(String(m.captures[1]))
    imag_val = parse_arbreal(String(m.captures[2]))
    return ArbComplex(real_val, imag_val)
end

# Парсинг строки-массива ArbComplex-ов
function parse_arbcomplex_array(str::String)
    inner = strip(str, ['[', ']'])

    if isempty(inner)
        return ArbComplex[]
    end

    # Match entries like: [real ± err] + [imag ± err]j
    pattern = r"\[[^\[\]]+\]\s*[\+\-]\s*\[[^\[\]]+\]j"
    matches = collect(eachmatch(pattern, inner))

    return [parse_arbcomplex(String(m.match)) for m in matches]
end



# Основной загрузчик JSON
function load_coeffs(json_file::String)
    raw_data = JSON.parsefile(json_file)
    parsed = Dict{Int, Vector{ArbComplex}}()

    for (k, v_str) in raw_data
        idx = parse(Int, k)

        if occursin(r"\[.*\].*\[.*\]j", v_str)
            parsed[idx] = parse_arbcomplex_array(v_str)
        else
            parsed[idx] = [ArbComplex(parse_arbreal(String(strip(v_str, ['[', ']']))), 0)]
        end
    end

    return parsed
end

# RR функция
function RR(s::ArbComplex, a_coeffs::Vector{ArbComplex})
    total = ArbComplex(0, 0)
    for k in 1:length(a_coeffs)
        total += a_coeffs[k] * (k ^ -s)
    end
    return total
end

# Генерация видео
function generate_RR_evolution_gif(coeff_dict::Dict{Int, Vector{ArbComplex}}, output_path::String)
    x_vals = range(-6, stop=0, length=100)
    anim = @animate for t in 1:length(coeff_dict)
        a_coeffs = coeff_dict[t]
        R_real_vals = [RR(ArbComplex(x, 0), a_coeffs) for x in x_vals]
        
        plot(x_vals, R_real_vals,
             label = "Re(R(s))",
             xlabel = "x (Re(s))",
             ylabel = "Re(R(s))",
             title = "Odds amount $t",
             legend = :topright,
             ylim = (-10, 10),
             lw = 2)
    end

    gif(anim, output_path, fps = 5)
end

function generate_RR_evolution_video(coeff_dict::Dict{Int, Vector{ArbComplex}}, output_path::String)
    x_vals = range(-6, stop=0, length=100)
    anim = @animate for t in sort(collect(keys(coeff_dict)))
        a_coeffs = coeff_dict[t]
        R_real_vals = [RR(ArbComplex(x, 0), a_coeffs) for x in x_vals]
        
        plot(x_vals, R_real_vals,
             label = "Re(R(s))",
             xlabel = "x (Re(s))",
             ylabel = "Re(R(s))",
             title = "Odds amount $t",
             legend = :topright,
             ylim = (-10, 10),
             lw = 2)
    end

    # Save as MP4 using ffmpeg
    mp4_path = output_path * ".mp4"
    mp4(anim, mp4_path, fps = 5)
end

function generate_coeffs_evolution_gif(coeff_dict::Dict{Int, Vector{ArbComplex}}, output_path::String)
    anim = @animate for t in 1:length(coeff_dict)
        a_coeffs = coeff_dict[t]
        t_values = 1:length(coeff_dict[t])
        R_real_vals = [real(a_coeffs[t]) for t in t_values]
        
        scatter(t_values, R_real_vals,
             label = "Re(X)",
             xlabel = "T",
             ylabel = "Re(X)",
             title = "Odds $t",
             legend = :topright,
             ylim = (-2, 2),
             lw = 2)
    end

    gif(anim, output_path, fps = 5)

end

function generate_coeffs_evolution_video(coeff_dict::Dict{Int, Vector{ArbComplex}}, output_path::String)
    anim = @animate for t in 1:length(coeff_dict)
        a_coeffs = coeff_dict[t]
        t_values = 1:length(coeff_dict[t])
        R_real_vals = [real(a_coeffs[t]) for t in t_values]
        
        scatter(t_values, R_real_vals,
             label = "Re(X)",
             xlabel = "T",
             ylabel = "Re(X)",
             title = "Odds $t",
             legend = :topright,
             ylim = (-2, 2),
             lw = 2)
    end

    mp4_path = output_path * ".mp4"
    mp4(anim, mp4_path, fps = 5)

end

# Запуск
zeta_series_filename = joinpath(@__DIR__, "assets", "zeta_series.json")
coeffs = load_coeffs(zeta_series_filename)

generate_coeffs_evolution_gif(coeffs, "output/coeffs_evolution.gif")
generate_coeffs_evolution_video(coeffs, "output/coeffs_evolution")


# generate_gif(coeffs, "output/R_function_evolution.gif")
# generate_video(coeffs, "output/R_function_evolution")
