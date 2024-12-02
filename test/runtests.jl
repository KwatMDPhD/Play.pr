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

st_ = deepcopy(ma_)

for ie in 1:uf

    me = 0.0

    va = 0.0

    for id in 1:ud

        nu_ = ma_[id][ie, :]

        u1 = us_[id]

        m1 = sum(nu_) / u1

        v1 = sum(nu -> (nu - m1)^2, nu_)

        me += m1 * u1 / us

        va += v1 / us

    end

    me_[ie] = me

    va_[ie] = va

    for id in 1:ud

        st = st_[id]

        for is in 1:us_[id]

            st[ie, is] = (st[ie, is] - me) / sqrt(va)

        end

    end

end

@test isapprox(
    me_[vcat(1:3, (end - 2):end)],
    [4.73804444, 5.95952282, 3.85754567, 8.1264884, 7.05907568, 8.35037346],
)

@test isapprox(
    va_[vcat(1:3, (end - 2):end)],
    [0.16001282, 0.20522794, 0.74683924, 0.44558136, 0.44762729, 0.67700222],
)

@test isapprox(st_[1][1, 1:3], [-1.49485243, -0.37141835, -0.52958901])

@test isapprox(st_[end][end, (end - 2):end], [0.7445448, 0.25915602, 0.32268557])
