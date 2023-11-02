using BoundaryLayerDynamics, MPI

function open_channel_les(; measure_only = false)

    outdir = abspath("$(dirname(PROGRAM_FILE))/../../data/les-julia")

    Re = 1e8
    z0 = 1e-4
    grid_size = (64, 64, 64)

    domain = Domain((2π, 4/3*π, 1), RoughWall(z0), FreeSlipBoundary())
    processes = incompressible_flow(1/Re, sgs_model =
            StaticSmagorinskyModel(Cs=0.1, wall_damping=true))
    abl = Model(grid_size, domain, processes)
    initialize!(abl, vel1 = (x,y,z) -> 1/0.4 * log(z/z0))

    # manually add noise with limited intensity
    BoundaryLayerDynamics.State.add_noise!(abl.state.vel1, 1e-3)

    turnovers = 100
    steps_per_turnover = 10_000
    dt = 1/steps_per_turnover
    nt = measure_only ? 100 : turnovers * steps_per_turnover
    sfreq = div(steps_per_turnover, 2) * dt
    snapshots = Snapshots(path = joinpath(outdir, "snapshots"),
                          frequency = sfreq, precision = Float32)
    evolve!(abl, dt * nt, dt = dt, method = AB2(),
            output = snapshots, verbose = true)

end

MPI.Init()
open_channel_les(measure_only = haskey(ENV, "MEASURE"))
