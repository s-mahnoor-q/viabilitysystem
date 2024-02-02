# parameter 1 is the frequency, parameter 2 is the peak
# 2-node display compatibility
using GLMakie
using Graphs
using GraphMakie
using Makie.Colors
using JLD2
using DynamicalSystems
using Attractors
using BifurcationKit

recordFromSolution(x, p) = (u1 = x[1], u2 = x[2])

struct ODENetworkModel
    num_nodes::Int64 
    adj_matrix::Array{Int64}
    par_adj_matrix::Array{Int64}
    parameters::Array{Float64}
    p_index::Int64
    p_start::Float64
    p_end::Float64
    u₀::Array{Float64}
end

function input(t, freq, peak = 2.0)
    return peak .* (0.5 .* sin.(freq .* t) .+ 0.5)
end

# create equations of motion
function create_viability_model(odenetworkmodel::ODENetworkModel, u₀ = odenetworkmodel.u₀; return_eqs = false, input_f = false)
    #should be in place
    function model!(du, u, p, t)
        Φ(x, ϵ = 0.05, θ = 0.5) = @. 1.0 / (1.0 + exp(-(x - θ)/ϵ))
        
        freq = p[1]; peak = p[2]
        
        F = input(t, freq, peak)

        for i in range(1, odenetworkmodel.num_nodes)
            du[i] = -u[i]
            input = 0.0

            for j in range(1, odenetworkmodel.num_nodes)
                if odenetworkmodel.adj_matrix[i,j] == 2
                    input += u[j]
                end
            end

            # if there were no OR relationships, reset input ready for multiplication
            if input == 0.0 && odenetworkmodel.adj_matrix[i,:] != zeros(Int64, odenetworkmodel.num_nodes)
                input = 1.0
            end

            for j in range(1, odenetworkmodel.num_nodes)
                if odenetworkmodel.adj_matrix[i,j] == 1
                    input *= u[j]
                end
            end
            
            # println("here")
            if odenetworkmodel.par_adj_matrix[1,i] != 0
                if input_f
                    input *= F
                end
            end

            du[i] += Φ(input)
        end

        return du
    end

    if return_eqs
        return model!
    end

    ds = ContinuousDynamicalSystem(model!, u₀, odenetworkmodel.parameters)

    return ds
end

function parameter_labels(num_nodes, p_matrix)
    labels = ["" for i in 1:num_nodes]

    for i in 1:size(p_matrix,1)
        for node in 1:num_nodes
            if p_matrix[i,node] != 0
                labels[node] = labels[node] * "p$i"
            end
        end
    end
    return labels
end

function create_UI(ODENetworkModel)
    resolution = 0.1; # has been removed to allow the trajectory to decide on step sizes.
    simulation_time = 20.0
    g = SimpleDiGraph(ODENetworkModel.num_nodes)
    fig = Figure(size = (1920, 1080))

    ds = create_viability_model(ODENetworkModel; input_f = true)

    X, t = DynamicalSystems.trajectory(ds, simulation_time
        # ; Δt = resolution
    )

    network_and_init = fig[1,1] = GridLayout()

    # create graph visualisation
    ax1 = Axis(network_and_init[1,1], title = "$(ODENetworkModel.num_nodes) Node System", width = 250, height = 250, xrectzoom = false, yrectzoom = false)
    println("Created network...")
    edge_colours = []
    for i in 1:ODENetworkModel.num_nodes
        for j in 1:ODENetworkModel.num_nodes
            if ODENetworkModel.adj_matrix[i,j] == 1
                add_edge!(g, i, j)
                push!(edge_colours,:black)
            elseif ODENetworkModel.adj_matrix[i,j] == 2
                add_edge!(g, i, j)
                push!(edge_colours,:red)
            end
        end
    end
    
    p = graphplot!(ax1, g; edge_color = edge_colours,
        ilabels = [i for i in 1:ODENetworkModel.num_nodes], 
        nlabels = parameter_labels(ODENetworkModel.num_nodes, ODENetworkModel.par_adj_matrix),
        curve_distance=0.5, curve_distance_usage=true)
    
    offsets = (p[:node_pos][]) * 0.2
    p.nlabels_offset[] = offsets

    xmin = p[:node_pos][][1][1]; xmax = p[:node_pos][][1][1]
    ymin = p[:node_pos][][1][2]; ymax = p[:node_pos][][1][2]

    for i in 1:ODENetworkModel.num_nodes
        point = p[:node_pos][][i]

        point[1] < xmin ? xmin = point[1] : nothing
        point[1] > xmax ? xmax = point[1] : nothing

        point[2] < ymin ? ymin = point[2] : nothing
        point[2] > ymax ? ymax = point[2] : nothing
    end
    fill = 0.75
    ax1.limits[] = ((xmin - fill, xmax + fill), (ymin - fill, ymax + fill))
    # hidespines!(ax1)
    hidedecorations!(ax1)

    # create initial condition text boxes and rerun button
    textboxes_frame = network_and_init[1,2] = GridLayout()
    initial_conditions = []
    Label(textboxes_frame[1,1], "Initial Conditions (Float64)", width = 150)
    for i in 1:ODENetworkModel.num_nodes
        tb = Textbox(textboxes_frame[i+1,1], placeholder = "$(ODENetworkModel.u₀[i])", width = 100, validator = Float64)
        push!(initial_conditions, tb)
        tb.stored_string = ODENetworkModel.u₀[i]
    end
    rerun_button = Button(textboxes_frame[ODENetworkModel.num_nodes + 2,1], label = "Rerun", width = 100)

    rand_button = Button(textboxes_frame[ODENetworkModel.num_nodes + 3,1], label = "Random", width = 100)

    timeseries_box = fig[1,2] = GridLayout()
    timeseries = Axis(timeseries_box[1,1], title = "Time Series", width = 600, height = 250, xrectzoom = false, yrectzoom = false, limits = (nothing, (-0.1,ODENetworkModel.parameters[2])))
    colors = [RGBf(1, 0, 0), RGBf(1, 0.5, 0), RGBf(1, 1, 0), RGBf(0.5, 1, 0),
    RGBf(0, 1, 0), RGBf(0, 1, 0.5), RGBf(0, 1, 1), RGBf(0, 0.5, 1),
    RGBf(0, 0, 1), RGBf(0.5, 0, 1), RGBf(1, 0, 1), RGBf(1, 0, 0.5),
    RGBf(1, 0, 0), RGBf(1, 0.5, 0), RGBf(1, 1, 0)]

    trajectory_box = fig[2,2] = GridLayout()
    axes = []
    if ODENetworkModel.num_nodes == 2
        trajectory = Axis(trajectory_box[1,1], title = "Trajectory", width = 250, height = 250, limits = ((-0.1, 1.1), (-0.1, 1.1)), xrectzoom = false, yrectzoom = false)
        lines!(trajectory, X[:,1], X[:,2], color = :red)
        scatter!(trajectory, X[1,:], color = :red, markersize = 10, marker = :star6)
    else    
        trajectory = Axis3(trajectory_box[1,1], title = "Trajectory", width = 250, height = 250, limits = ((-0.1, 1.1), (-0.1, 1.1), (-0.1, 1.1)))
        axis_opts = trajectory_box[1,2] = GridLayout()
        variables = [i for i in 1:ODENetworkModel.num_nodes]
        axis3d_1 = Menu(axis_opts[1,1], options = variables, default = 1, width = 100); Label(axis_opts[1,2], "x-axis", width = 75)
        axis3d_2 = Menu(axis_opts[2,1], options = variables, default = 2, width = 100); Label(axis_opts[2,2], "y-axis", width = 75)
        axis3d_3 = Menu(axis_opts[3,1], options = variables, default = 3, width = 100); Label(axis_opts[3,2], "z-axis", width = 75)
        GLMakie.scatter!(trajectory, X[1,1], X[1,2], X[1,3], color = :red, markersize = 10, marker = :star6)
        axes = [axis3d_1, axis3d_2, axis3d_3]
    end

    function plot_2d_trajectory(X, t)
        empty!(trajectory)
        lines!(trajectory, X[:,1], X[:,2], color = :red)    
        GLMakie.scatter!(trajectory, X[1,:], color = :red, markersize = 10, marker = :star6)
    end

    function plot_3d_trajectory(X, t)
        empty!(trajectory)
        d1 = axes[1].selection[]; d2 = axes[2].selection[]; d3 = axes[3].selection[]
        scatter!(trajectory, X[1,d1], X[1,d2], X[1,d3]; color = :red, markersize = 10, marker = :star6)
        lines!(trajectory, X[:,d1], X[:,d2], X[:,d3], color = :red)
    end

    function plot_trajectory(X, t)
        if ODENetworkModel.num_nodes == 2
            plot_2d_trajectory(X, t)
        else
            plot_3d_trajectory(X, t)
        end
    end

    function plot_timeseries(X, t, freq = freq_slider.value[], peak = peak_slider.value[])
        empty!(timeseries)
        
        timeseries.limits = (nothing, (-0.1, max(1.2, peak + 0.1)))
        lines!(timeseries, t, input(t, freq, peak), linestyle = :dash, color = :gray, alpha = 0.5)

        lins = [lines!(timeseries, t, X[:,i]; color = colors[i % 15 + 1], label = "var $i") for i in 1:ODENetworkModel.num_nodes]
        Legend(timeseries_box[1,2], lins, ["Process $i" for i in 1:ODENetworkModel.num_nodes]; margin = (10, 10, 10, 10), halign = :right, valign = :top)
    end

    # rerun timeseries, trajectory
    on(rerun_button.clicks) do x
        u_new = [tb.stored_string[] isa Float64 ? tb.stored_string[] : parse(Float64, tb.stored_string[]) for tb in initial_conditions]
        
        ds = create_viability_model(ODENetworkModel, u_new; input_f = true)
        
        X, t = DynamicalSystems.trajectory(ds, simulation_time; 
            # Δt = resolution
        )
        plot_trajectory(X, t)

        plot_timeseries(X, t, freq_slider.value[], peak_slider.value[])
    end

    # random trajectory generator
    on(rand_button.clicks) do x
        u_new = [rand() for i in 1:ODENetworkModel.num_nodes]

        for i in 1:ODENetworkModel.num_nodes
            initial_conditions[i].displayed_string[] = string(round(u_new[i]; sigdigits = 3))
            println(initial_conditions[i].displayed_string)
        end

        ds = create_viability_model(ODENetworkModel, u_new; input_f = true)
        
        X, t = DynamicalSystems.trajectory(ds, simulation_time; 
            # Δt = resolution
        )
        plot_trajectory(X, t)

        plot_timeseries(X, t, freq_slider.value[], peak_slider.value[])
    end

    # create parameter sliders...
    sliders = []
    sliders_lbl = []
    slider_box = fig[2,1] = GridLayout()

    # creates slider which changes the freq of input parameter
    freq_slider = Slider(slider_box[1,2], range = 0.0:0.05:10.0, startvalue = 0.5, width = 300)
    lbl = Label(slider_box[1,3], "Frequency", width = 75)
    value = Label(slider_box[1,1], "$(freq_slider.value[])", width = 75)
    push!(sliders, freq_slider.value)
    push!(sliders_lbl, value)

    peak_slider = Slider(slider_box[2,2], range = 0.0:0.05:15.0, startvalue = 0.5, width = 300)
    lbl = Label(slider_box[2,3], "Peak", width = 75)
    value = Label(slider_box[2,1], "$(peak_slider.value[])", width = 75)
    push!(sliders, peak_slider.value)
    push!(sliders_lbl, value)
    
    println("Created slider...")

    if ODENetworkModel.num_nodes != 2
        for a in axes
            on(a.selection) do select
                plot_3d_trajectory(X, t)
            end
        end
    end

    # change frequency and peak and recalculate trajectory
    @lift begin
        sliders_lbl[1].text[] = "$(sliders[1].val)"
        sliders_lbl[2].text[] = "$(sliders[2].val)"
        set_parameter!(ds, 1, $(sliders[1]))
        set_parameter!(ds, 2, $(sliders[2]))

        X, t = DynamicalSystems.trajectory(ds, simulation_time; 
            # Δt = resolution
        )
        
        plot_timeseries(X, t, freq_slider.value[], peak_slider.value[])
        
        plot_trajectory(X, t)
    end
    
    plot_timeseries(X, t, freq_slider.value[], peak_slider.value[])
    println("Created timeseries...")

    plot_trajectory(X, t)
    println("Created trajectory...")
    
    display(fig)
    println("UI ready!")
end

function create_attractors(ODENetworkModel)
    ds = create_viability_model(ODENetworkModel)
    # set_parameter!(ds, 1, 1.33)
    # set_parameter!(ds, 2, 1.35)
   
    xg = -0.5:0.3:1.3
    grid = tuple((xg for i in 1:ODENetworkModel.num_nodes)...)
    
    mapper = AttractorsViaRecurrences(ds, grid)

    basins, attractors = basins_of_attraction(mapper, grid; show_progress = true)

    fig = Figure()

    ax = Axis3(fig[1,1]; title = "Attractors", width = 600, height = 600, limits = (0, 1, 0, 1, 0, 1))

    for key in keys(attractors)
        attr = vec(attractors[key])
        println(attr)
        a = attractors[key]
        scatter!(ax, a[:,1], a[:,2], a[:,3]; markersize = 50)
    end

    display(fig)

    return attractors
end

function create_branch(model!, ax, u0, par, opts_br, p_index)
	prob = BifurcationProblem(model!, u0, par,
		# specify the continuation parameter
		(@lens _[p_index]), record_from_solution = recordFromSolution)

    br = BifurcationKit.continuation(prob, PALC(), opts_br; verbosity = 3, bothside = true)

    println(br)

	BifurcationKit.plot!(ax, br)
end

function create_bifurcation_diagram(ODENetworkModel)
    # we know that for low F, there will be eq near [0 0 ... 0]
    # and for high F there MAY be an eq near [1 1 ... 0]
    # so we will use this for the continuation
    fig = Figure()
    ax = Axis(fig[1,1], title = "Viability Bifurcation Diagram")
    model = create_viability_model(ODENetworkModel; return_eqs = true)
    println("Created model...")

    opts_br = ContinuationPar(p_min = ODENetworkModel.p_start, p_max = ODENetworkModel.p_end, dsmax = 0.05,
		# options to detect bifurcations
		detect_bifurcation = 3, n_inversion = 8, max_bisection_steps = 25,
		# maximum number of continuation steps
		max_steps = 5000,)

        opts_br2 = ContinuationPar(p_min = ODENetworkModel.p_start, p_max = ODENetworkModel.p_end, ds = -0.01, dsmax = 0.01,
		# options to detect bifurcations
		detect_bifurcation = 3, n_inversion = 8, max_bisection_steps = 100,
		# maximum number of continuation steps
		max_steps = 5000,)

    create_branch(model, ax, zeros(ODENetworkModel.num_nodes),  ODENetworkModel.parameters, opts_br, ODENetworkModel.p_index)
	create_branch(model, ax, ones(ODENetworkModel.num_nodes),  ODENetworkModel.parameters, opts_br2, ODENetworkModel.p_index)

    display(fig)
end

function main()
    problem = nothing

    if pwd()[end-14:end] != "viabilitysystem"
        @warn "pwd() is not in the correct directory. Please change to '.../viabilitysystem' using cd() and try again."
    end
    
    println("Welcome to the Viability System Applet! Please enter the name of your problem file located in /Networks (e.g. 'two_node_params.jld2'):")
    invalid = true
    while invalid 
        try
            file = readline()
            directory = pwd() * "\\Applet_Param\\Networks\\" * file
            # println(directory)
            @load directory problem
            invalid = false
        catch y
            println("File not found. Please try again and enter a valid file name:")
            continue
        end
    end
    
    println("Program starting...")
    println(problem)
    create_UI(problem)

    # create_bifurcation_diagram(problem)
end

main()