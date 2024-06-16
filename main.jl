using FFTW
using Distributions
using GLMakie

function initialise_sources(N)
    l = rand(Uniform(-1.0, 1.0), N)
    m = rand(Uniform(-1.0, 1.0), N)
    brightness = rand(Uniform(0.0, 10.0), N)

    sources = stack([l, m, brightness], dims=1)
    return sources
end

function initialise_baselines(N)
    u = rand(-200:200, N)
    v = rand(-200:200, N)

    baselines = stack([u, v], dims=1)
    return baselines
end

function calculate_vis(baselines, sources)

    vis = Vector{ComplexF64}(undef, length(baselines))
    # Calculate visibility for each baseline
    for i in range(1, length(baselines))
        vis[i] = sum(sources[3, :] .* exp.(2 * pi * im * (baselines[1, :] .* sources[1, :] .+ baselines[2, :] .* sources[2, :])))
    end
    println(display(vis))
end

function gridding(baselines, vis)
    # Create the grid 
    u_max = 50
    du = 0.5

end

# REMEMBER COLUMN MAJOR!
function main()
    baselines = initialise_baselines(10)
    sources = initialise_sources(10)

    println("BASELINES")
    println(display(baselines))

    println("SOURCES")
    println(display(sources))
    calculate_vis(baselines, sources)

    f = Figure()
    ax = Axis(f[1, 1], limits=(-1, 1, -1, 1), aspect=AxisAspect(1), xlabel="l", ylabel="m", title="sources")
    scatter!(ax, sources[1, :], sources[2, :]; color=sources[3, :])
    Colorbar(f[1, 2], limits=(minimum(sources[3, :]), maximum(sources[3, :])))
    wait(display(f))


    f2 = Figure()
    ax2 = Axis(f2[1, 1], limits=(-200, 200, -200, 200), aspect=AxisAspect(1), xlabel="u", ylabel="v", title="baselines")
    scatter!(ax2, baselines[1, :], baselines[2, :])
    wait(display(f2))
end

main()
