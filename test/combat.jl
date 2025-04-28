using Nucleus

using Play

# ---- #

ts_ =
    Tuple(joinpath(Play.IN, ba) for ba in ("GSE18520.tsv", "GSE66957.tsv", "GSE69428.tsv"))

# ---- #

u1 = lastindex(ts_)

st__ = Vector{Vector{String}}(undef, u1)

um_ = Vector{Int}(undef, u1)

F_ = Vector{Matrix{Float64}}(undef, u1)

for id in 1:u1

    A = Nucleus.Table.rea(ts_[id])

    st__[id] = A[!, 1]

    um_[id] = size(A, 2) - 1

    F_[id] = Matrix(A[!, 2:end])

end

u3 = sum(um_)

# ---- #

fe_ = reduce(intersect, st__)

u2 = lastindex(fe_)

# ---- #

for id in 1:u1

    F_[id] = F_[id][indexin(fe_, st__[id]), :]

end

# ---- #

me_ = Vector{Float64}(undef, u2)

va_ = similar(me_)

st_ = map(similar, F_)

ga = Matrix{Float64}(undef, u2, u1)

de = similar(ga)

for ie in 1:u2

    me = 0.0

    va = 0.0

    for id in 1:u1

        nu_ = F_[id][ie, :]

        u1 = um_[id]

        m1 = sum(nu_) / u1

        v1 = sum(nu -> (nu - m1)^2, nu_)

        me += m1 * u1 / u3

        va += v1 / u3

    end

    me_[ie] = me

    va_[ie] = va

    for id in 1:u1

        ma = F_[id]

        st = st_[id]

        ua = um_[id]

        for is in 1:ua

            st[ie, is] = (ma[ie, is] - me) / sqrt(va)

        end

        nu_ = st[ie, :]

        ga[ie, id] = sum(nu_) / ua

        de[ie, id] = var(nu_; corrected = false)

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

# ---- #

@test isapprox(
    ga[vcat(1:3, (end - 2):end), :],
    [
        -0.56819588 -0.09791935 1.24568143 -0.93382787 -1.3277284 -0.05480174
        1.34767285 1.75551025 -1.03638821 2.54451507 1.6699258 0.21199566
        -1.97217538 -3.96418236 -0.24024632 -4.02553049 -1.08889624 -0.38535142
    ]',
)

@test isapprox(
    de[vcat(1:3, (end - 2):end), :],
    [
        0.52503647 0.32167221 0.14062541 0.80817126 0.07214159 0.5086251
        1.52436832 1.80318174 2.0726262 1.11581255 2.03429736 1.3820283
        0.78418236 0.56259004 0.31480659 1.1411774 0.554778 1.15850571
    ]',
)

# ---- #

gm_ = map(mean, eachcol(ga))

gv_ = map(var, eachcol(ga))

@test isapprox(gm_, [-0.2632684539878588, 0.8935639505524893, -1.5541379305822973])

@test isapprox(
    gv_,
    [2.6366798372981464, 2.1904522505055417, 1.8341420904760084];
    rtol = 1e-4,
)

dm_ = map(mean, eachcol(de))

dv_ = map(var, eachcol(de))

p1_ = map((dm, dv) -> (dm^2 + 2 * dv) / dv, dm_, dv_)

p2_ = map((dm, dv) -> (dm^3 + dm * dv) / dv, dm_, dv_)

@test isapprox(p1_, [3.649330176189566, 11.065529845342768, 4.090202496497435]; rtol = 1e-4)

@test isapprox(
    p2_,
    [1.2237704407501055, 14.305946682918783, 3.604970820054511];
    rtol = 1e-4,
)

# ---- #

for id in 1:u1

    sa_ = st_[id]

    ga_ = ga[:, id]

    de_ = de[:, id]

    gs_ = Float64[]

    ds_ = Float64[]

    for ie in 1:u2

    end

end
