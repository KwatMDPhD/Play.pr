using Distributions: Binomial

using GLM: StatisticalModel, coef

using Nucleus

using Omics

using Play

# ---- #

# TODO: Avoid DataFrame
const A =
    Nucleus.Table.rea(pkgdir(Nucleus, "data", "Table", "_.tsv"); missingstring = "NA")[
        !,
        ["name", "survived", "sex", "age", "fare"],
    ]

filter!(an_ -> all(!ismissing, an_), A)

# ---- #

const S1 = "Passenger"

const S1_ = A[!, 1]

const U1 = lastindex(S1_)

for st in ("Connolly, Miss. Kate", "Kelly, Mr. James")

    S1_[findlast(==(st), S1_)] *= " (2)"

end

@assert allunique(S1_)

# ---- #

const S2 = "Survival"

const NU_ = A[!, 2]

# ---- #

const S2_ = "Sex", "Age", "Fare"

const NU__ = convert(Vector{Int}, map(st -> if st == "female"
    0
elseif st == "male"
    1
end, A[!, 3])),
convert(Vector{Float64}, A[!, 4]),
convert(Vector{Float64}, A[!, 5])

const U2 = lastindex(S2_)

# ---- #
# TODO: Use Binary

const GL_ = Vector{StatisticalModel}(undef, U2)

for i1 in 1:U2

    st = S2_[i1]

    nu_ = NU__[i1]

    GL_[i1] = gl = Omics.LinearModel.make(nu_, Binomial(), NU_)

    gr_ = Omics.Grid.make(nu_, U1)

    e1_, e2_, e3_ = Omics.LinearModel.make(gl, gr_)

    Omics.LinearModelPlot.writ(
        joinpath(Play.OU, "titanic.$st.html"),
        nu_,
        NU_,
        gr_,
        e1_,
        e2_,
        e3_,
        Dict(
            "title" => Dict("text" => coef(gl)[2]),
            "yaxis" => Dict("title" => Dict("text" => S2)),
            "xaxis" => Dict("title" => Dict("text" => st)),
        ),
    )

end

# ---- #

const PR = sum(isone, NU_) / U1

const P = [Omics.LinearModel.make(GL_[i1], NU__[i1][i2])[1] for i1 in 1:U2, i2 in 1:U1]

for nd in Nucleus.Extreme.index(map(po_ -> Omics.Evidence.make(PR, po_), eachcol(P)), 4)

    po_ = P[:, nd]

    for um in 0:U2

        Omics.EvidencePlot.writ(
            joinpath(Play.OU, "titanic.$nd.$um.html"),
            PR,
            S2_[1:um],
            po_[1:um],
        )

    end

end
