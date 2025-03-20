using Plots
using LinearAlgebra
using Printf
using AbstractAlgebra

# Определение класса DirichletSeries
struct DirichletSeries
    coefs::Vector{Complex{BigFloat}}

    function DirichletSeries(coefs::Vector{Complex{BigFloat}})
        new(coefs)
    end
end

function (R::DirichletSeries)(s::Complex{BigFloat})
    value = Complex{BigFloat}(0, 0)
    for i in 1:size(R.coefs, 1)
        value += R.coefs[i, 1] * (i ^ -s)
    end
    return value
end

# Читаем строки из файла и преобразуем их в массив комплексных чисел
function read_zeta_zeros(filename)
    zeros = Complex{BigFloat}[]
    open(filename, "r") do file
        for line in eachline(file)
            imag_part = parse(BigFloat, strip(line))  # Парсим строку в число   

            push!(zeros, Complex{BigFloat}(0.5, imag_part))  # Создаём комплексное число 0.5 + i*imag_part
            push!(zeros, Complex{BigFloat}(0.5, -imag_part))  # Создаём комплексное число 0.5 - i*imag_part
        end
    end
    return zeros
end

# Загружаем данные из файла
zeros_filename = joinpath(@__DIR__, "assets", "zeta_zeros.txt")
nont_zeros = read_zeta_zeros(zeros_filename)

if isempty(nont_zeros)
    error("Массив `nont_zeros` пуст! Возможно, не загружены нули дзета-функции.")
end

n = length(nont_zeros) + 1
println("n = $n")

A = Matrix{Complex{BigFloat}}(undef, n, n)


# Добавляем строку для a₁ = 1
A[1, :] .= Complex{BigFloat}(0, 0)
A[1, 1] = Complex{BigFloat}(1, 0)

# Заполнение матрицы коэффициентов
for (i, nont_zero) in enumerate(nont_zeros)
    println("nont_zero is: ", nont_zero)

    for j in 1:n
        A[i + 1, j] = exp(-nont_zero * log(BigFloat(j)))
    end
end

# for i in 1:(n - 1)
#     row_str = join(["(" * @sprintf("%.4f", real(x)) * " + " * @sprintf("%.4f", imag(x)) * "im)" for x in A[i, :]], ", ")
#     println("before adding 1 to A[$i] is: [", row_str, "]")
# end


# Создание вектора B с комплексными числами
B = Complex{BigFloat}[0 for _ in 1:n]  # Заполняем вектор нулями
B[1] = Complex{BigFloat}(1, 0)  # Устанавливаем последний элемент в 1

for i in 1:n
    println("B[$i] is: ", B[i])
end

# for i in 1:n
#     row_str = join(["(" * @sprintf("%.4f", real(x)) * " + " * @sprintf("%.4f", imag(x)) * "im)" for x in A[i, :]], ", ")
#     println("A[$i] is: [", row_str, "]")
# end

λ = 1e-6  # Параметр регуляризации
X = inv(A' * A + λ * I) * A' * B  # Регуляризация Тихонова

for i in 1:n
    println("X[$i] is: [", X[i, :], "]")
end

R = DirichletSeries(X)

# Вывод значений
for x in [-10, -8, -6, -5, -4, -3, -2, -1, 0]
    println("R($x) = $(R(Complex{BigFloat}(x, 0)))")
end
println("R(0.5 + 14.1347i) = $(R(Complex{BigFloat}(0.5, 14.13472514173469379045)))")
println("R(1 + 2πi/ln2) = $(R(Complex{BigFloat}(1, 2 * pi / log(2))))")

# График
x_vals = range(-3, 0, length=100)
y_vals = [real(R(Complex{BigFloat}(x, 0))) for x in x_vals]

display(plot(x_vals, y_vals, label="R(x)"))
hline!([0], linestyle=:dash, label="y=0")

savefig("output/output_plot.png")


# Построение графика на комплексной плоскости
complex_vals = [R(nont_zero) for nont_zero in nont_zeros]  # Замените на вашу логику получения комплексных значений

# Разделяем на действительную и мнимую части для графика
real_parts = real.(complex_vals)
imag_parts = imag.(complex_vals)

# Построение графика на комплексной плоскости
scatter(real_parts, imag_parts, label="Zeros on Complex Plane", xlabel="Re", ylabel="Im")

# Сохранение второго графика
savefig("output/complex_plane_plot.png")


# Извлекаем действительную и мнимую части
real_parts = [real(X[i][1]) for i in 1:length(X)]
imag_parts = [imag(X[i][1]) for i in 1:length(X)]
t_values = 1:length(X)  # Time indices

# Создаём scatter plot (точки на плоскости)
scatter(t_values, real_parts, label="X points", marker=:circle, markersize=4)

xlabel!("T")
ylabel!("Re(X)")
title!("Scatter Plot of X")

savefig("output/coefs_scatter.png")