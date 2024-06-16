using Base: _linspace
using FFTW
using Distributions
using GLMakie

function initialise_sources(N)
    l = rand(Uniform(-1.0, 1.0), N)
    m = rand(Uniform(-1.0, 1.0), N)
    brightness = rand(Uniform(10.0, 50.0), N)

    sources = stack([l, m, brightness], dims=1)
    return sources
end

function initialise_baselines(N)
    u = rand(-100:100, N)
    v = rand(-100:100, N)

    baselines = stack([u, v], dims=1)
    return baselines
end

function calculate_vis(baselines, sources)
    vis = zeros(ComplexF64, size(baselines)[2])
    # Calculate visibility for each baseline
    for i in range(1, size(sources)[2])
        vis .+= sources[3, i] .* exp.(-2 .* pi .* im .* (baselines[1, :] .* sources[1, i] .+ baselines[2, :] .* sources[2, i]))
    end

    return vis
end

function find_uv_index(u, v, u_vec, v_vec)
    return argmin(abs.(u_vec .- u)), argmin(abs.(v_vec .- v))
end

function gridding(baselines, vis)
    # Create the grid 
    u_max = 100
    du = 0.5

    n_grid = floor(Int, u_max / du)
    u_vec = range(-u_max, u_max, n_grid + 1)
    v_vec = range(-u_max, u_max, n_grid + 1)

    gridded_vis = zeros(ComplexF64, (n_grid + 1, n_grid + 1))
    weights = zeros((n_grid + 1, n_grid + 1))

    for i in range(1, size(baselines)[2])
        u_ind, v_ind = find_uv_index(baselines[1, i], baselines[2, i], u_vec, v_vec)
        gridded_vis[u_ind, v_ind] += vis[i]
        weights[i] += 1
    end

    gridded_vis[findall(weights .!= 0)] ./= weights[findall(weights .!= 0)]
    return (gridded_vis)
end

function imaging(vis)
    image = circshift(fftshift(ifft(circshift(ifftshift(vis, 1), (1, 1)))), (-1, 0))

    return image
end

# REMEMBER COLUMN MAJOR!
# BUG: I'm pretty sure something is still wrong here, I think in the calculate_vis function.
function main()
    baselines = initialise_baselines(10000)
    sources = initialise_sources(10000)

    println("BASELINES")
    println(display(baselines))

    println("SOURCES")
    println(display(sources))

    vis = calculate_vis(baselines, sources)
    println("VIS")
    println(display(vis))

    gridded_vis = gridding(baselines, vis)

    f = Figure()
    ax = Axis(f[1, 1], limits=(-1, 1, -1, 1), aspect=AxisAspect(1), xlabel="l", ylabel="m", title="sources")
    scatter!(ax, sources[1, :], sources[2, :]; color=sources[3, :])
    # Colorbar(f[1, 1][1, 2], limits=(minimum(sources[3, :]), maximum(sources[3, :])))

    ax2 = Axis(f[1, 2], limits=(-200, 200, -200, 200), aspect=AxisAspect(1), xlabel="u", ylabel="v", title="baselines")
    scatter!(ax2, baselines[1, :], baselines[2, :])

    ax3 = Axis(f[2, 1], aspect=AxisAspect(1), xlabel="u", ylabel="v", title="visibilities")
    heatmap!(ax3, abs.(gridded_vis))

    image = imaging(gridded_vis)

    print(display(image))
    ax4 = Axis(f[2, 2], aspect=AxisAspect(1), xlabel="u", ylabel="v", title="image")
    heatmap!(ax4, abs.(image))

    wait(display(f))
    # wait(display(f3))
end

main()
