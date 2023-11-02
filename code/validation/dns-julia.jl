using BoundaryLayerDynamics, MPI

function closed_channel_dns(; measure_only = false)

    outdir = abspath("$(dirname(PROGRAM_FILE))/../../data/dns-julia")

    ν = 3.5e-4
    η = 0.97 # grid stretching factor
    grid_size = (256,192,96)
    p = 2 # order of polynomial for IC, should be even
    u0(x,y,z) = (1 - (z-1)^p)*(p+1)/p

    domain = Domain((4π, 2π, 2), SmoothWall(), SmoothWall(),
                    SinusoidalMapping(η, :symmetric))
    processes = incompressible_flow(ν, constant_flux = 1)
    abl = Model(grid_size, domain, processes)
    initialize!(abl, vel1 = u0)

    # manually add noise with appropriate intensity
    BoundaryLayerDynamics.State.add_noise!(abl.state.vel1, 1e-2)

    turnovers = 20
    dt = 4.48e-3 # time step estimate
    steps_per_turnover = 3_500 # time ~15.7 / dt, based on uτ of LM2015
    nt = measure_only ? 100 : turnovers * steps_per_turnover
    sfreq = div(steps_per_turnover, 10) * dt
    snapshots = Snapshots(path = joinpath(outdir, "snapshots"),
                          frequency = sfreq, precision = Float32)

    evolve!(abl, dt * nt, dt = dt, output = snapshots, verbose = true)

end

MPI.Init()
closed_channel_dns(measure_only = haskey(ENV, "MEASURE"))
