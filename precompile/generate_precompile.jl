using SnoopCompile

inf_timing = @snoopi tmin = 0.001 include("snoop/snoop.jl")
pcs = SnoopCompile.parcel(inf_timing)
SnoopCompile.write("precompile_sources", pcs)
cp("precompile_sources/precompile_QuantumDynamics.jl", "precompile.jl"; force=true)