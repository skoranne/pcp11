
function kaporin_grad(M::AbstractMatrix{T}, ε::Real) where {T}
    r   = size(M, 1)
    G   = Symmetric(M * M' + ε * I)
    invG = inv(G)                       # r×r dense inverse – cheap because r is tiny
    trG = tr(G)

    # 2/r * ( M / tr(G)  -  invG * M )
    return (2.0 / r) * (M / trG .- invG * M)
end
using Zygote

# --------------------------------------------------------------
# 2. Kaporin proxy (exact, cheap for r ≤ 64) ------------------
# --------------------------------------------------------------
"""
    kaporin_proxy(M::AbstractMatrix{T}, ε::Real) where {T<:Real}

Continuous Kaporin rank‑diversity proxy

    K(M) = log( tr(G)/r ) - (1/r)·logdet(G)

with G = M·M' + ε·I   (row‑wise Gram for A, column‑wise Gram for B).

The function returns a *scalar* (Float64).  It is differentiable via
Zygote because it only uses `cholesky` and elementary ops.
"""

# ---- original (forward) definition ---------------------------------
function kaporin_proxy(M::AbstractMatrix{T}, ε::Real) where {T<:Real}
    r = size(M, 1)
    G = Symmetric(M * M' + ε * I)          # SPD
    L = cholesky(G).L                     # ← the part Zygote cannot differentiate
    trG   = tr(G)
    logdetG = 2.0 * sum(log.(diag(L)))    # log‑det from Cholesky
    return log(trG / r) - logdetG / r
end

# ---- custom backward (adjoint) ------------------------------------
Zygote.@adjoint function kaporin_proxy(M::AbstractMatrix{T}, ε::Real) where {T}
    y = kaporin_proxy(M, ε)                # forward value
    function back(Δ)                       # Δ is the upstream scalar gradient (normally 1.0)
        # ∂K/∂M  (matrix)  multiplied by Δ (scalar)
        gradM = Δ * kaporin_grad(M, ε)
        # ε is a scalar constant – no gradient needed
        return (gradM, zero(ε))
    end
    return y, back
end

using LinearAlgebra, Statistics, Random, Zygote, Plots

# --------------------------------------------------------------
# 1. Kaporin proxy + custom adjoint (copied from the block above)
# --------------------------------------------------------------
# (the two functions `kaporin_proxy` and `kaporin_grad` + the @adjoint)

# --------------------------------------------------------------
# 2. Helper: effective rank (entropy‑based)
# --------------------------------------------------------------
function effective_rank(M::AbstractMatrix)
    s = svdvals(M) .+ eps()
    p = s / sum(s)
    return exp(-sum(p .* log.(p)))
end

# --------------------------------------------------------------
# 3. Core simulation (identical to the previous script, except
#    the loss closure now uses the *differentiable* Kaporin term)
# --------------------------------------------------------------
function simulate_rka(; d_out=128, d_in=64, r=8, batch=32,
                     noise_std=0.01, opt=:adam, η=1e-2,
                     β1=0.9, β2=0.999, λ=0.0, ε=1e-4,
                     steps=500, seed=42, record_every=5)

    rng = MersenneTwister(seed)

    # base weight (fixed)
    W0 = randn(rng, d_out, d_in) / sqrt(d_in)

    # one mini‑batch that we reuse (synthetic regression problem)
    X = randn(rng, batch, d_in)
    Y = X * W0' .+ noise_std * randn(rng, batch, d_out)

    # low‑rank factors (small random init)
    A = randn(rng, r, d_in) * 0.01
    B = randn(rng, d_out, r) * 0.01

    # Adam state
    mA = zeros(size(A)); vA = zeros(size(A))
    mB = zeros(size(B)); vB = zeros(size(B))

    # containers for diagnostics
    nrec = div(steps, record_every) + 1
    loss_hist = Vector{Float64}(undef, nrec)
    kapA_hist = Vector{Float64}(undef, nrec)
    kapB_hist = Vector{Float64}(undef, nrec)
    erA_hist  = Vector{Float64}(undef, nrec)
    erB_hist  = Vector{Float64}(undef, nrec)

    # ----- loss + metrics (closure for Zygote) -----------------
    function loss_and_metrics(A_, B_)
        ΔW = B_ * A_                       # low‑rank update
        Ŷ  = X * (W0 .+ ΔW)'               # predictions
        mse = sum(abs2, Ŷ .- Y) / batch
        kapA = kaporin_proxy(A_, ε)        # row‑wise Gram
        kapB = kaporin_proxy(B_', ε)       # column‑wise Gram (B' is r×d_out)
        total = mse + λ * (kapA + kapB)
        return total, (mse, kapA, kapB)
    end

    # ----- optimisation loop ------------------------------------
    rec = 1
    for t in 1:steps
        # forward + pull‑back
        (ℓ, (mse, kapA, kapB)), back = Zygote.pullback(() -> loss_and_metrics(A, B), A, B)
        ∂A, ∂B = back((one(ℓ), nothing, nothing))   # gradients w.r.t. A and B

        # Adam update
        if opt == :adam
            mA .= β1 .* mA .+ (1-β1) .* ∂A
            vA .= β2 .* vA .+ (1-β2) .* (∂A .^ 2)
            m̂A = mA ./ (1 - β1^t)
            v̂A = vA ./ (1 - β2^t)
            A .-= η .* m̂A ./ (sqrt.(v̂A) .+ 1e-8)

            mB .= β1 .* mB .+ (1-β1) .* ∂B
            vB .= β2 .* vB .+ (1-β2) .* (∂B .^ 2)
            m̂B = mB ./ (1 - β1^t)
            v̂B = vB ./ (1 - β2^t)
            B .-= η .* m̂B ./ (sqrt.(v̂B) .+ 1e-8)
        else            # simple SGD fallback
            A .-= η .* ∂A
            B .-= η .* ∂B
        end

        # record diagnostics
        if t % record_every == 0 || t == 1
            loss_hist[rec] = ℓ
            kapA_hist[rec] = kapA
            kapB_hist[rec] = kapB
            erA_hist[rec]  = effective_rank(A)
            erB_hist[rec]  = effective_rank(B')
            rec += 1
        end
    end

    return Dict(
        :loss   => loss_hist,
        :kapA   => kapA_hist,
        :kapB   => kapB_hist,
        :erA    => erA_hist,
        :erB    => erB_hist,
        :A_final=> copy(A),
        :B_final=> copy(B),
        :settings=> (d_out=d_out, d_in=d_in, r=r,
                    λ=λ, ε=ε, η=η, opt=opt, steps=steps)
    )
end

# --------------------------------------------------------------
# 4. Demo: compare “no regulariser” vs “with Kaporin”
# --------------------------------------------------------------
res_no  = simulate_rka(λ=0.0, opt=:adam, steps=800, record_every=5)
res_reg = simulate_rka(λ=0.1, opt=:adam, steps=800, record_every=5)

# --------------------------------------------------------------
# 5. Quick visualisation (requires Plots.jl)
# --------------------------------------------------------------
function plot_trajectories(no, yes; title_suffix="")
    p1 = plot(no[:loss], label="no‑reg", lw=2, xlabel="record step",
              ylabel="total loss", title="Training loss $title_suffix")
    plot!(p1, yes[:loss], label="Kaporin λ=0.1", lw=2, ls=:dash)

    p2 = plot(no[:kapA], label="kapA (no‑reg)", lw=2,
              xlabel="record step", ylabel="Kaporin proxy",
              title="Kaporin A‑side $title_suffix")
    plot!(p2, yes[:kapA], label="kapA (λ=0.1)", lw=2, ls=:dash)

    p3 = plot(no[:kapB], label="kapB (no‑reg)", lw=2,
              xlabel="record step", ylabel="Kaporin proxy",
              title="Kaporin B‑side $title_suffix")
    plot!(p3, yes[:kapB], label="kapB (λ=0.1)", lw=2, ls=:dash)

    p4 = plot(no[:erA], label="eff‑rank A (no‑reg)", lw=2,
              xlabel="record step", ylabel="effective rank",
              title="Effective rank of A $title_suffix")
    plot!(p4, yes[:erA], label="eff‑rank A (λ=0.1)", lw=2, ls=:dash)

    p5 = plot(no[:erB], label="eff‑rank B (no‑reg)", lw=2,
              xlabel="record step", ylabel="effective rank",
              title="Effective rank of B $title_suffix")
    plot!(p5, yes[:erB], label="eff‑rank B (λ=0.1)", lw=2, ls=:dash)

    plot(p1, p2, p3, p4, p5; layout=(3,2), size=(960,800), legend=:bottomright)
end

plot_trajectories(res_no, res_reg, title_suffix="(r=8, Adam)")
