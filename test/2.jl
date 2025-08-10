using StatsBase: mean, var

using Nucleus

using Play

# ---- #

const BA_ = "GSE18520.tsv", "GSE66957.tsv", "GSE69428.tsv"

const U3 = lastindex(BA_)

# ---- #

const ST__ = Vector{Vector{String}}(undef, lastindex(BA_))

const UM_ = similar(ST__, Int)

const N_ = similar(UM_, Matrix{Float64})

for nd in 1:U3

    _, ST__[nd], _, N = Nucleus.Table.ge(Nucleus.Table.rea(joinpath(Play.IN, BA_[nd])))

    UM_[nd] = size(N, 2)

    N_[nd] = N

end

# ---- #

const ST_ = reduce(intersect, ST__)

const U1 = lastindex(ST_)

@assert U1 === 17126

for nd in 1:U3

    N_[nd] = N_[nd][indexin(ST_, ST__[nd]), :]

end

@assert UM_ == [63, 69, 29]

# ---- #

# Global mean
const ME_ = Vector{Float64}(undef, U1)

# Global variance
const VA_ = similar(ME_)

const Z_ = map(similar, N_)

const PR = inv(sum(UM_))

@assert PR === 0.006211180124223602

const M = Matrix{Float64}(undef, U1, U3)

const V = similar(M)

for i1 in 1:U1

    su = va = 0

    for i3 in 1:U3

        nu_ = N_[i3][i1, :]

        um = UM_[i3]

        me = sum(nu_) / um

        su += me * um

        va += sum(nu -> (nu - me)^2, nu_)

    end

    ME_[i1] = su *= PR

    VA_[i1] = va *= PR

    for i3 in 1:U3

        N = N_[i3]

        Z = Z_[i3]

        um = UM_[i3]

        for i2 in 1:um

            Z[i1, i2] = (N[i1, i2] - su) / sqrt(va)

        end

        zs_ = Z[i1, :]

        M[i1, i3] = sum(zs_) / um

        V[i1, i3] = var(zs_; corrected = false)

    end

end

# ---- #

const IN_ = vcat(1:3, (U1 - 2):U1)

@assert isapprox(
    ME_[IN_],
    [4.73804444, 5.95952282, 3.85754567, 8.1264884, 7.05907568, 8.35037346],
)

@assert isapprox(
    VA_[IN_],
    [0.16001282, 0.20522794, 0.74683924, 0.44558136, 0.44762729, 0.67700222],
)

@assert isapprox(Z_[1][1, 1:3], [-1.49485243, -0.37141835, -0.52958901])

@assert isapprox(Z_[end][end, (end - 2):end], [0.7445448, 0.25915602, 0.32268557])

@assert isapprox(
    M[IN_, :],
    [
        -0.56819588 -0.09791935 1.24568143 -0.93382787 -1.3277284 -0.05480174
        1.34767285 1.75551025 -1.03638821 2.54451507 1.6699258 0.21199566
        -1.97217538 -3.96418236 -0.24024632 -4.02553049 -1.08889624 -0.38535142
    ]',
)

@assert isapprox(
    V[IN_, :],
    [
        0.52503647 0.32167221 0.14062541 0.80817126 0.07214159 0.5086251
        1.52436832 1.80318174 2.0726262 1.11581255 2.03429736 1.3820283
        0.78418236 0.56259004 0.31480659 1.1411774 0.554778 1.15850571
    ]',
)

# ---- #

m1_ = map(mean, eachcol(M))

v1_ = map(var, eachcol(M))

@assert isapprox(m1_, [-0.2632684539878588, 0.8935639505524893, -1.5541379305822973])

@assert isapprox(
    v1_,
    [2.6366798372981464, 2.1904522505055417, 1.8341420904760084];
    atol = 1e-3,
)

# ---- #

m2_ = map(mean, eachcol(V))

v2_ = map(var, eachcol(V))

p1_ = map((m2, v2) -> (m2^2 + 2 * v2) / v2, m2_, v2_)

p2_ = map((m2, v2) -> (m2^3 + m2 * v2) / v2, m2_, v2_)

@assert isapprox(
    p1_,
    [3.649330176189566, 11.065529845342768, 4.090202496497435];
    atol = 1e-3,
)

@assert isapprox(
    p2_,
    [1.2237704407501055, 14.305946682918783, 3.604970820054511];
    atol = 1e-3,
)
