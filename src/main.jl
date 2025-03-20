using ArbNumerics
using Printf
using Plots
using LinearAlgebra

# Читаем все строки из файла и преобразуем их в массив комплексных чисел
function read_zeta_zeros(filename)
    zeros = Vector{ArbComplex}([])

    open(filename, "r") do file
        for line in eachline(file)
            imag_part = ArbReal(strip(line))  # Парсим строку в число

            push!(zeros, ArbComplex(0.5, imag_part))  # Создаём комплексное число 0.5 + i*imag_part
            push!(zeros, ArbComplex(0.5, -imag_part))  # Создаём комплексное число 0.5 + i*imag_part
        end

        # for line in eachline(file)
        #     imag_part = ArbReal(strip(line))  # Парсим строку в число

        #     push!(zeros, ArbComplex(0.5, -imag_part))  # Создаём комплексное число 0.5 + i*imag_part
        # end

    end
    return zeros
end

# Function for matrix-vector multiplication
function matvecmul(A::Matrix{ArbComplex}, B::Vector{ArbComplex})
    result = Vector{ArbComplex}(undef, size(A, 1))  # Create a result vector
    for i in 1:size(A, 1)
        sum = ArbComplex(0, 0)  # Start with zero for complex sum
        for j in 1:size(A, 2)
            sum += A[i, j] * B[j]
        end
        result[i] = sum  # Store the result
    end
    return result
end


# Custom Gaussian Elimination for ArbComplex types with partial pivoting
function gaussian_elimination(A::Matrix{ArbComplex{128}}, B::Vector{ArbComplex{128}})
    n = size(A, 1)
    
    # Create augmented matrix [A | B]
    Augmented = hcat(A, B)
    
    # Forward elimination with partial pivoting
    for i in 1:n
        # Find the pivot row (maximum element in the current column)
        pivot_row = argmax(abs.(Augmented[i:n, i]))[1] + i - 1
        if iszero(Augmented[pivot_row, i])
            error("Matrix is singular and cannot be solved.")
        end
        
        # Swap rows if necessary
        if pivot_row != i
            Augmented[i, :], Augmented[pivot_row, :] = Augmented[pivot_row, :], Augmented[i, :]
        end
        
        # Normalize the pivot row
        pivot = Augmented[i, i]
        for k in 1:n+1
            Augmented[i, k] /= pivot
        end
        
        # Eliminate other rows
        for j in i+1:n
            factor = Augmented[j, i]
            for k in 1:n+1
                Augmented[j, k] -= factor * Augmented[i, k]
            end
        end
    end

    # Back substitution
    X = Vector{ArbComplex{128}}(undef, n)
    for i in n:-1:1
        # Start with the last element in the augmented matrix
        X[i] = Augmented[i, n+1]
        
        # Subtract the known values of X from the right-hand side
        for j in i+1:n
            X[i] -= Augmented[i, j] * X[j]
        end
    end
    
    return X
end



# Загружаем данные из файла
zeros_filename = joinpath(@__DIR__, "assets", "zeta_zeros.txt")
nont_zeros = read_zeta_zeros(zeros_filename)

if isempty(nont_zeros)
    error("Массив `nont_zeros` пуст! Возможно, не загружены нули дзета-функции.")
end

# n будет длиной массива nont_zeros + 1
n = length(nont_zeros) + 1
println("n = $n")

# Создаем матрицу A размером n x n
A = Matrix{ArbComplex{128}}(undef, n, n)

# Заполнение матрицы A с учетом нулей дзета-функции
for (i, nont_zero) in enumerate(nont_zeros)
    s_i = nont_zero

    for k in 1:n
    
        println("s_i is: ", s_i, "  k is: ", k, "  k ^ s_i is: ", k ^ s_i)
        A[i, k] = k ^ -s_i
    end
end

# Добавление строки для a₁ = 1
A[n, :] .= ArbComplex{128}(0, 0)
A[n, 1] = ArbComplex{128}(1, 0)

# Создание вектора B с комплексными числами
B::Vector{ArbComplex{128}} = Vector{ArbComplex{128}}(undef, n)
B[:] .= ArbComplex{128}(0, 0)
B[end] = ArbComplex{128}(1, 0)

for i in 1:(n - 1)
    row_str = join(["(" * @sprintf("%.4f", real(x)) * " + " * @sprintf("%.4f", imag(x)) * "im)" for x in A[i, :]], ", ")
    println("before adding 1 to A[$i] is: [", row_str, "]")
end

# Параметр регуляризации
# λ = ArbReal(1e-6)

# Создаём идентичную матрицу I для регуляризации (с типом ArbComplex)
# II = Matrix{ArbComplex}(I, n, n)

# A_transposed = Matrix(A')  # Convert Adjoint to a regular matrix

# Apply regularization
# A = A_transposed * A + λ * II
# B = matvecmul(A_transposed, B)

# Q, R = qr(A)

# Solve using QR decomposition (R is upper triangular)
# X = R \ (Q' * B)

# Solve the system using Gaussian elimination
X = gaussian_elimination(A, B)


# Получаем коэффициенты
a_coeffs = [X[i] for i in 1:n]
for k in 1:n
    println("a_$k = $(a_coeffs[k])")
end

# Функция для вычисления R(s)
function RR(s, a_coeffs)
    total = ArbComplex(0, 0)  # Инициализируем комплексное число
    for k in 1:length(a_coeffs)
        total += a_coeffs[k] * (k ^ -s)
    end
    return total
end

# Генерируем точки x от -6 до 0
x_vals = range(-6, stop=0, length=100)
R_real_vals = [real(RR(ArbComplex(x, 0), a_coeffs)) for x in x_vals]

# Строим график
plot(x_vals, R_real_vals, label="Re(R(s))")
xlabel!("x")
ylabel!("Re(R(s))")
title!("График Re(R(s)) для s = x + 0i, x ∈ [-6, 0]")
# grid(true)

# Сохранение графика
savefig("output/output_plot.png")

# Построение графика на комплексной плоскости (если нужно)
complex_vals = [RR(nont_zero, a_coeffs) for nont_zero in nont_zeros]

# Разделяем на действительную и мнимую части для графика
real_parts = real.(complex_vals)
imag_parts = imag.(complex_vals)

# Построение графика на комплексной плоскости
scatter(real_parts, imag_parts, label="Zeros on Complex Plane", xlabel="Re", ylabel="Im")

# Сохранение второго графика
savefig("output/complex_plane_plot.png")

real_parts = [real(X[i][1]) for i in 1:length(X)]
imag_parts = [imag(X[i][1]) for i in 1:length(X)]
t_values = 1:length(X)  # Time indices

# Создаём scatter plot (точки на плоскости)
scatter(t_values, real_parts, label="X points", marker=:circle, markersize=4)

xlabel!("T")
ylabel!("Re(X)")
title!("Scatter Plot of X")

savefig("output/coefs_scatter.png")