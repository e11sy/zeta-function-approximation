module ArbGMRES

using LinearAlgebra
using ArbNumerics

export gmres

function gmres(
    A::AbstractMatrix{ArbComplex{P}},
    b::AbstractVector{ArbComplex{P}};
    restart::Int = 30,
    tol::ArbReal{P} = ArbReal{P}("1e-40"),
    maxiter::Int = 200,
    verbose::Bool = true,
) where {P}

    n = length(b)
    x = zeros(ArbComplex{P}, n)
    r = b - A * x
    β = norm(r)

    if β < tol
        return x
    end

    for outer = 1:div(maxiter, restart)
        V = Matrix{ArbComplex{P}}(undef, n, restart + 1)
        H = zeros(ArbComplex{P}, restart + 1, restart)
        cs = zeros(ArbReal{P}, restart)
        sn = zeros(ArbReal{P}, restart)
        e₁ = zeros(ArbComplex{P}, restart + 1)
        e₁[1] = β

        V[:, 1] = r / β

        for j in 1:restart
            w = A * V[:, j]

            # Arnoldi: Modified Gram-Schmidt
            for i in 1:j
                H[i, j] = dot(V[:, i], w)
                w -= H[i, j] * V[:, i]
            end

            # Optional reorthogonalization
            for i in 1:j
                hij2 = dot(V[:, i], w)
                H[i, j] += hij2
                w -= hij2 * V[:, i]
            end

            H[j+1, j] = norm(w)

            if isnan(real(H[j+1, j])) || isnan(imag(H[j+1, j])) || abs(H[j+1, j]) < tol
                if verbose
                    @warn "Arnoldi breakdown at step $j, H = $(H[j+1, j])"
                end
                break
            end

            V[:, j+1] = w / H[j+1, j]

            # Apply Givens rotations to the new column
            for i in 1:j-1
                temp = cs[i]*H[i,j] + sn[i]*H[i+1,j]
                H[i+1,j] = -sn[i]*H[i,j] + cs[i]*H[i+1,j]
                H[i,j] = temp
            end

            # New Givens rotation
            rj = sqrt(abs2(real(H[j,j])) + abs2(real(H[j+1,j])))
            cs[j] = real(H[j,j])/rj
            sn[j] = real(H[j+1,j])/rj

            H[j,j] = ArbComplex{P}(cs[j]) * H[j,j] + ArbComplex{P}(sn[j]) * H[j+1,j]
            H[j+1,j] = 0

            e₁[j+1] = -sn[j] * e₁[j]
            e₁[j] = cs[j] * e₁[j]

            resid = abs(e₁[j+1])
            if verbose
                # @info "GMRES iter $(outer*restart + j - restart): residual norm = $resid"
            end

            if isnan(resid)
                error("NaN encountered in residual. Matrix may be too ill-conditioned.")
            end

            if resid < tol
                y = H[1:j, 1:j] \ e₁[1:j]
                x += V[:, 1:j] * y
                return x
            end
        end

        # Restart
        y = H[1:restart, 1:restart] \ e₁[1:restart]
        x += V[:, 1:restart] * y
        r = b - A * x
        β = norm(r)

        if β < tol
            return x
        end
    end

    return x
    # error("GMRES did not converge after $maxiter iterations.")
end

end # module
