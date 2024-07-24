using SnoopCompile

inf_timing = @snoopi_deep include("snoop/quapi.jl")
ts, pcs = SnoopCompile.parcel(inf_timing)
SnoopCompile.write("precompile_sources", pcs)