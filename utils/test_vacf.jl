include("../molecule_tools/read_xyz.jl")
include("correlation_functions.jl")
using .Correlation_Functions

header, labels, geoms = read_xyz("data/w6_dynamics_velocity_production_full.xyz")

series = form_atomic_time_series(geoms)

signal = Correlation_Functions.autocorrelation_function(series, 5000)
Correlation_Functions.plot_correlation_function(signal, 2.5)
freqs, vdos = Correlation_Functions.vdos(signal, 2.5)
Correlation_Functions.plot_vdos(freqs, vdos)
