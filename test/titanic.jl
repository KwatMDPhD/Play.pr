using DataFrames: completecases, disallowmissing!

using GLM: StatisticalModel

using Omics

using Play

# ---- #

const TB = Omics.Table.rea(
    pkgdir(Omics, "data", "Table", "titanic.tsv");
    select = ["name", "survived", "sex", "age", "fare"],
    missingstring = "NA",
)

TB = disallowmissing!(TB[completecases(TB), :])

# ---- #

const NS = "Passenger"

const SA_ = TB[!, "name"]

# ---- #

for na in ("Connolly, Miss. Kate", "Kelly, Mr. James")

    SA_[findlast(==(na), SA_)] *= " (2)"

end

@assert allunique(SA_)

# ---- #

const NT = "Survival"

const TA_ = TB[!, "survived"]

# ---- #

const NF_ = ("Sex", "Age", "Fare")

const DA = stack(
    (map(se -> se == "female" ? 0 : 1, TB[!, "sex"]), TB[!, "age"], TB[!, "fare"]);
    dims = 1,
)

# ---- #

const UF, US = size(DA)

# ---- #

const GL_ = Vector{StatisticalModel}(undef, UF)

for ie in 1:UF

    nf = NF_[ie]

    fe_ = DA[ie, :]

    GL_[ie] = gl = Omics.GeneralizedLinearModel.fit(TA_, fe_)

    io_ = sortperm(fe_)

    fe_ = fe_[io_]

    Omics.GeneralizedLinearModel.plot(
        joinpath(p.OU, "titanic.$nf.html"),
        NS,
        SA_[io_],
        NT,
        TA_[io_],
        nf,
        fe_,
        Omics.GeneralizedLinearModel.predic(gl, range(fe_[1], fe_[end], US))...,
    )

end

# ---- #

const PR = sum(isone, TA_) / lastindex(TA_)

# ---- #

const PO = stack((
    map(ie -> Omics.GeneralizedLinearModel.predic(GL_[ie], sa_[ie])[1], 1:UF) for
    sa_ in eachcol(DA)
))

# ---- #

for is in Omics.Extreme.get(map(po_ -> Omics.Evidence.ge(PR, po_), eachcol(PO)), 4)

    po_ = PO[:, is]

    for uf in 0:UF

        Omics.Evidence.plot(
            joinpath(Play.OU, "titanic.$is.$uf.html"),
            SA_[is],
            PR,
            map(ie -> "$(NF_[ie]) = $(DA[ie, is])", 1:uf),
            po_[1:uf],
        )

    end

end
