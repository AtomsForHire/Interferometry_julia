using FFTW
using Distributions
using GLMakie

function initialise_sources(N)
    if N == 1
        l = 0
        m = 0
        brightness = 10

        sources = [l m brightness]'
    else
        l = rand(Uniform(-1.0, 1.0), N)
        m = rand(Uniform(-1.0, 1.0), N)
        # brightness = rand(Uniform(10.0, 50.0), N)
        brightness = ones(N)

        sources = stack([l, m, brightness], dims=1)
    end
    return sources
end

function initialise_baselines(N)
    u = rand(-50:50, N)
    v = rand(-50:50, N)

    baselines = stack([u, v], dims=1)
    return baselines
end

function calculate_vis_analytic(sources)
    u_max = 100
    du = 0.5
    n_grid = floor(Int, u_max / du)
    u_vec = range(-u_max, u_max, n_grid + 1)
    v_vec = range(-u_max, u_max, n_grid + 1)

    u_mat = u_vec' .* ones(n_grid + 1)
    v_mat = v_vec .* ones(n_grid + 1)'

    vis = zeros(ComplexF64, (n_grid + 1, n_grid + 1))

    for i in range(1, size(sources)[2])
        result = exp.(-2 * pi * im * (u_mat .* sources[1, i] .+ v_mat .* sources[2, i]))
        vis += sources[3, i] .* result
    end

    return vis
end

function calculate_vis(baselines, sources)
    vis = zeros(ComplexF64, size(baselines)[2])
    # Calculate visibility for each baseline
    for i in range(1, size(sources)[2])
        vis += sources[3, i] .* exp.(-2 .* pi .* im .* (baselines[1, :] .* sources[1, i] .+ baselines[2, :] .* sources[2, i]))
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
    baselines = initialise_baselines(1000)
    sources = initialise_sources(1)

    # println("BASELINES")
    # println(display(baselines))

    # println("SOURCES")
    # println(display(sources))

    vis = calculate_vis(baselines, sources)
    # println("VIS")
    # println(display(vis))

    gridded_vis = gridding(baselines, vis)
    # println("GRIDDED VIS")
    # println(display(gridded_vis))
    analytic_vis = calculate_vis_analytic(sources)
    # println("ANALYTIC VIS")
    # println(display(analytic_vis))

    f = Figure()
    ax = Axis(f[1, 1], limits=(-1, 1, -1, 1), aspect=AxisAspect(1), xlabel="l", ylabel="m", title="sources")
    scatter!(ax, sources[1, :], sources[2, :]; color=sources[3, :])
    # Colorbar(f[1, 1][1, 2], limits=(minimum(sources[3, :]), maximum(sources[3, :])))

    ax2 = Axis(f[2, 1], limits=(-50, 50, -50, 50), aspect=AxisAspect(1), xlabel="u", ylabel="v", title="baselines")
    scatter!(ax2, baselines[1, :], baselines[2, :])

    ax3 = Axis(f[1, 2], aspect=AxisAspect(1), xlabel="u", ylabel="v", title="analytic visibilities")
    heatmap!(ax3, abs.(analytic_vis))

    ax4 = Axis(f[2, 2], aspect=AxisAspect(1), xlabel="u", ylabel="v", title="visibilities")
    heatmap!(ax4, abs.(gridded_vis))

    image = imaging(gridded_vis)
    analytic_image = imaging(analytic_vis)

    ax5 = Axis(f[1, 3], aspect=AxisAspect(1), xlabel="u", ylabel="v", title="analytic image")
    heatmap!(ax5, abs.(analytic_image))

    ax6 = Axis(f[2, 3], aspect=AxisAspect(1), xlabel="u", ylabel="v", title="real image")
    heatmap!(ax6, abs.(image))

    save("output.png", f)
end

@time main()
