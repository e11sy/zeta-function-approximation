using ArbNumerics
using Printf
using Plots
using LinearAlgebra

setprecision(1024)

# Читаем все строки из файла и преобразуем их в массив комплексных чисел
function read_zeta_zeros(filename)
    zeros = Vector{ArbComplex}([])

    open(filename, "r") do file
        for line in eachline(file)
            imag_part = ArbReal(strip(line))  # Парсим строку в число

            push!(zeros, ArbComplex(0.5, imag_part))  # Создаём комплексное число 0.5 + i*imag_part
            push!(zeros, ArbComplex(0.5, -imag_part))  # Создаём комплексное число 0.5 + i*imag_part
        end
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
# function gaussian_elimination(A::Matrix{ArbComplex{1024}}, B::Vector{ArbComplex{1024}})
#     n = size(A, 1)
    
#     # Create augmented matrix [A | B]
#     println("Type of A: ", typeof(A))
#     println("Type of B: ", typeof(B))
#     Augmented = hcat(copy(A), copy(B))
#     println("Type of Augmented: ", typeof(Augmented))

    
#     # Forward elimination with partial pivoting
#     for i in 1:n
#         # Find the pivot row (maximum element in the current column)
#         _, rel_index = findmax(abs.(Augmented[i:n, i]))
#         pivot_row = rel_index + i - 1

#         # Swap rows if necessary
#         if pivot_row != i
#             Augmented[i, :], Augmented[pivot_row, :] = Augmented[pivot_row, :], Augmented[i, :]
#         end
        
#         # Normalize the pivot row
#         pivot = Augmented[i, i]

#         if iszero(Augmented[pivot_row, i])
#             error("Matrix is singular and cannot be solved.")
#         end
        
#         if abs(Augmented[i, i]) < 1e-100
#             println("WARNING: pivot at $i is nearly zero: ", Augmented[i, i])
#         end

#         for k in 1:n+1
#             Augmented[i, k] /= pivot

#             if isnan(real(Augmented[i, k])) || isnan(imag(Augmented[i, k]))
#                 error("NaN detected after normalization at row $i, column $i Value: ", Augmented[i, k])
#             end

#         end
        
#         # Eliminate other rows
#         for j in i+1:n
#             factor = Augmented[j, i]
#             for k in 1:n+1
#                 Augmented[j, k] -= factor * Augmented[i, k]
#             end
#         end
#     end

#     # Back substitution
#     X = Vector{ArbComplex{1024}}(undef, n)
#     for i in n:-1:1
#         # Start with the last element in the augmented matrix
#         X[i] = Augmented[i, n+1]
        
#         # Subtract the known values of X from the right-hand side
#         for j in i+1:n
#             X[i] -= Augmented[i, j] * X[j]
#         end
#     end
    
#     return X
# end

function gaussian_elimination(A::Matrix{ArbComplex{P}}, B::Vector{ArbComplex{P}}; ε=ArbReal("1e-40")) where P
    n = size(A, 1)
    Aug = hcat(A, B)  # Augmented matrix
    row_perm = collect(1:n)
    col_perm = collect(1:n)

    for i in 1:n
        # Full pivoting: search maximum absolute value in submatrix
        pivot_val = zero(ArbReal)
        pivot_row, pivot_col = i, i

        for r in i:n, c in i:n
            a_rc = Aug[r, c]
            mag = abs(a_rc)
            if mag > pivot_val
                pivot_val = mag
                pivot_row, pivot_col = r, c
            end
        end

        @info "pivot[$i] = $(Aug[pivot_row, pivot_col])"

        # Threshold check
        if pivot_val < ε
            error("Pivot too small at step $i (|pivot| = $pivot_val), system may be singular or ill-conditioned.")
        end

        # Swap rows
        if pivot_row != i
            tmp_row = copy(Aug[i, :])
            Aug[i, :] .= Aug[pivot_row, :]
            Aug[pivot_row, :] .= tmp_row
            row_perm[i], row_perm[pivot_row] = row_perm[pivot_row], row_perm[i]
        end

        # Swap columns (incl. col_perm tracking)
        if pivot_col != i
            tmp_col = copy(Aug[:, i])
            Aug[:, i] .= Aug[:, pivot_col]
            Aug[:, pivot_col] .= tmp_col
            col_perm[i], col_perm[pivot_col] = col_perm[pivot_col], col_perm[i]
        end

        # Elimination
        for j in i+1:n
            f = Aug[j, i] / Aug[i, i]
            for k in i:n+1
                Aug[j, k] -= f * Aug[i, k]
            end
        end
    end

    # Back substitution
    x = zeros(ArbComplex{P}, n)
    for i in n:-1:1
        x[i] = Aug[i, end]
        for j in i+1:n
            x[i] -= Aug[i, j] * x[j]
        end
        x[i] /= Aug[i, i]
    end

    # Undo column swaps
    x_corrected = similar(x)
    for i in 1:n
        x_corrected[col_perm[i]] = x[i]
    end

    return x_corrected
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
A = Matrix{ArbComplex{1024}}(undef, n, n)

# Заполнение матрицы A с учетом нулей дзета-функции
for (i, nont_zero) in enumerate(nont_zeros)
    s_i = nont_zero

    for k in 1:n
    
        # println("s_i is: ", s_i, "  k is: ", k, "  k ^ s_i is: ", k ^ s_i)
        A[i, k] = k ^ -s_i
        # println(A[i, k])
    end
end

# Добавление строки для a₁ = 1
A[n, :] .= ArbComplex{1024}(0, 0)
A[n, 1] = ArbComplex{1024}(1, 0)

println(typeof(A))

# Создание вектора B с комплексными числами
B::Vector{ArbComplex{1024}} = Vector{ArbComplex{1024}}(undef, n)
B[:] .= ArbComplex{1024}(0, 0)
B[end] = ArbComplex{1024}(1, 0)

for i in 1:(n - 1)
    row_str = join(["(" * @sprintf("%.4f", real(x)) * " + " * @sprintf("%.4f", imag(x)) * "im)" for x in A[i, :]], ", ")
    # println("before adding 1 to A[$i] is: [", row_str, "]")
end

if any(x -> isnan(real(x)) || isnan(imag(x)), A)
    error("Matrix A contains NaN values!")
end

# Solve the system using Gaussian elimination
X = gaussian_elimination(A, B)


# Получаем коэффициенты
a_coeffs = [X[i] for i in 1:n]
println("a_1 = $(a_coeffs[1])")
# for k in 1:n
# end

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