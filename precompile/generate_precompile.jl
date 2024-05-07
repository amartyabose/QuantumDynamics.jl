using SnoopCompile

inf_timing = @snoopi include("snoop/quapi.jl")
pcs = SnoopCompile.parcel(inf_timing)
SnoopCompile.write("precompile_sources", pcs)