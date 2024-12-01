using DataIntegration

using Test: @test

# ----------------------------------------------------------------------------------------------- #

using DataFrames

using Omics

# ---- #

ts_ = Tuple(
    pkgdir(DataIntegration, "data", ts) for
    ts in ("GSE18520.tsv", "GSE66957.tsv", "GSE69428.tsv")
)

# ---- #

ud = lastindex(ts_)

fe___ = Vector{Vector{String}}(undef, ud)

us_ = Vector{Int}(undef, ud)

ma_ = Vector{Matrix{Float64}}(undef, ud)

for id in 1:ud

    da = Omics.Table.rea(ts_[id])

    fe___[id] = da[!, 1]

    us_[id] = size(da, 2) - 1

    ma_[id] = Matrix(da[!, 2:end])

end

us = sum(us_)

# ---- #

fe_ = intersect(fe___...)

uf = lastindex(fe_)

# ---- #

for id in 1:ud

    ma_[id] = ma_[id][indexin(fe_, fe___[id]), :]

end

# ---- #

me_ = Vector{Float64}(undef, uf)

va_ = Vector{Float64}(undef, uf)

for ie in 1:uf

    for id in 1:ud

        nu_ = ma_[id][ie, :]

        ua = us_[id]

        me = sum(nu_) / ua

        me_[ie] += me * ua / us

        va_[ie] += sum(nu -> (nu - me)^2, nu_) / us

    end

end

@test isapprox(me_[1:3], [4.73804444, 5.95952282, 3.85754567])

@test isapprox(va_[1:3], [0.16001282, 0.20522794, 0.74683924])

# ---- #

me_
