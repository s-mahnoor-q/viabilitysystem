# implementing OR relationships
using JLD2

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

function create_new_network()
    println("Welcome to the network generator!")
    println("---------------------------------")
    println("Note that entries with 1 will be AND together, and 2 will be OR together")
    println("Please select the number of nodes in your network:")
    invalid = true
    num_nodes = nothing
    while invalid
        try
            num_nodes = parse(Int64, readline())
            invalid = false
        catch y
            @warn "Incorrect type. Please try again and enter a valid integer:"
            continue
        end
    end

    adj_matrix = zeros(Int64, num_nodes, num_nodes)


    for row in range(1, num_nodes)
        invalid = true
        while invalid
            println("Please enter the adjacency row of node $row, separated by spaces:")
            adj_row = [parse.(Int64, s) for s in split(readline())]
            
            if length(adj_row) != num_nodes
                @warn "Incorrect number of nodes. Please try again."
                continue
            else
                adj_matrix[row,:] = adj_row
                invalid = false
            end
        end
        
    end

    println("Adjacency matrix:")
    println(adj_matrix)

    
    println("Please enter the parameters of your network, separated by spaces:")
    parameters = [parse(Float64, s) for s in split(readline())]
    par_adj_matrix = zeros(Int64, length(parameters), num_nodes)

    println("")
    
    println("Entering parameter adjacency matrix: e.g. parameter 1 feeds into node 1 and 2, so has adjacency row [1 1 0 ...]")
    for i in range(1,length(parameters))
        invalid = true
        while invalid
            println("Please enter the adjacency row of parameter $i = $(parameters[i]), separated by spaces:")
            adj_row = [parse.(Int64, s) for s in split(readline())]
            
            if length(adj_row) != num_nodes
                @warn "Incorrect number of nodes. Please try again."
                continue
            else
                par_adj_matrix[i,:] = adj_row
                invalid = false
            end
        end
    end

    println("Parameters: ", parameters)

    println("Please enter the index of the parameter to bifurcate: ")
    invalid = true
    p_index = 1
    while invalid
        try
            p_index = parse(Int64, readline())
            if p_index < 1 || p_index > length(parameters)
                @warn "Incorrect index. Please try again and enter a valid integer:"
                continue
            else
                invalid = false
            end
        catch y
            @warn "Incorrect type. Please try again and enter a valid integer:"
            continue
        end
    end

    println("Index: ", p_index)

    println("Please enter the range to bifurcate from and to separated by a space: ")
    p_start, p_end = [parse(Float64, s) for s in split(readline())]

    println("*Please enter the initial conditions of your network, separated by spaces:")
    invalid = true
    u₀ = []
    while invalid
        u₀ = [parse(Float64, s) for s in split(readline())]
        if length(u₀) != num_nodes
            @warn "Incorrect number of nodes. Please try again."
            continue
        else
            invalid = false
        end
    end

    problem = ODENetworkModel(num_nodes, adj_matrix, par_adj_matrix, parameters, p_index, p_start, p_end, u₀)

    println("Would you like to save this network? (y/n)")
    response = readline()

    if response == "y"
        println("Please enter the name of the file you would like to save to (e.g. net.jld2): ")
        filename = readline()
        jldsave("C:/Users/smahn/JuliaScripts2/Pipeline/data/" * filename; problem)
        println("File saved successfully!")
    elseif response == "n"
        println("File not saved.")
    end
    
    return problem
end

create_new_network()