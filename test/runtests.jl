using DataIntegration

using Test: @test

# ----------------------------------------------------------------------------------------------- #

using DataFrames

using Omics

# ---- #

da = pkgdir(DataIntegration, "data")

e1 = Omics.Table.rea(joinpath(da, "GSE18520.tsv"), )

e2 = Omics.Table.rea(joinpath(da, "GSE66957.tsv"), )

e3 = Omics.Table.rea(joinpath(da, "GSE69428.tsv"), )

# ---- #

ex = innerjoin(e1, e2, e3; on = "gene_symbol")
