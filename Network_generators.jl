 
## Collection of Network Generators ## 

# cd("Dir") - change to correct directory 

# Required packages # 
using Random, Statistics , DataFrames , Inequality, CSV,  Graphs 

include("Network_properties.jl")

function ER_model(N,p)    
    graph = Graph(N)

    for u in 1:N
        for v in (u+1):N  # currently only looking at the upper triangle of the matrix (for directed might need a diff. thing)
            if rand() < p 
                add_edge!(graph, u, v)
            end
        end
    end
    return graph 
end


############################################

############################################


function WS_model(N, m, p)

    graph = Graph(N)
    # Create connections to m/2 neighbors on each side (ring lattice)
    if iseven(m)
        m_right = m/2
    else
        m_right = (m+1)/2
    end

    for i in 1:N 
                
        for j in 1:m_right

            add_edge!(graph, i,mod((i + j - 1), N) + 1) 
        end
    end
       
    edgelist = [[src(e), dst(e)] for e in collect(edges(graph))]
    
      
    for i in 1:length(edgelist)  

        if rand() < p

            (u, v) = edgelist[i]  
            
            if has_edge(graph,u,v) # to make sure we're not double-re-rewiring
                
                new_v = rand(1:N)
            
                while new_v == u || has_edge(graph,new_v,u) || has_edge(graph,u,new_v)
                        
                    new_v = rand(1:N)
                end     
                   
                rem_edge!(graph, u, v)
        
                add_edge!(graph,u, new_v) 
            end             
        end
    end
    return graph
end


######################################

######################################


function BA_model(N, m)

   
    graph = Graph(m)
    
   
    for source in 1:m
        for target in (source+1):m
            add_edge!(graph, source,target)
        end
    end
   
    degrees = degree(graph)
    weights = Weights(degrees) 

    for u in (m + 1):N
        add_vertex!(graph)

        list = sample(1:length(degrees), Weights(weights), m,replace=false)
    
        for each in list 

            add_edge!(graph, u, each)
            degrees[each] += 1
           
        end
        push!(degrees, m)
        weights = Weights(degrees)
    end
    return graph 
end


############################################

###       CONFIGURATION MODEL            ###

############################################

function Configuration_Model(degree_sequence,max_attempts,seed)

    Random.seed!(seed)
    total_degree = sum(degree_sequence)
    @assert iseven(total_degree) "Sum of degrees must be even." 
    
    if !iseven(total_degree)
        x = rand(degree_sequence)
        degree_sequence[x] = degree_sequence[x] +=1
    end

    total_degree = sum(degree_sequence)
    @assert iseven(total_degree) "Sum of degrees must be even." 

    # Dic of nodes and edges to keep track of remaining degrees

    nodes_dict = OrderedDict{Int, Int}()
    for (node, degree) in enumerate(degree_sequence)
        nodes_dict[node] = degree
    end

    graph = Graph(length(degree_sequence))

    
    remaining_nodes = collect(keys(nodes_dict))

    iter_attempts = 0

    while sum(values(nodes_dict)) > 0

        candidates = filter(x -> nodes_dict[x] > 0, remaining_nodes)

        if isempty(candidates)
            break
        end
    
        original_node = rand(candidates)       
        node_friends = neighbors(graph, original_node)
    
        
        neighbour_options = filter(x -> !(x in node_friends) && x !=original_node, candidates)
        
       
        num_to_select =  min(length(neighbour_options), nodes_dict[original_node])
        selected_nodes = sample(collect(neighbour_options), num_to_select, replace=false)

        for each_node in selected_nodes

            add_edge!(graph, original_node, each_node)
            nodes_dict[original_node] -= 1
            nodes_dict[each_node] -= 1
        
        end

        filter!(x ->nodes_dict[x] != 0, candidates)

        # to avoid infinite loop

        if nodes_dict[original_node] > 0
            iter_attempts += 1
    
        end

        if iter_attempts > max_attempts
            break
        end

    
    end

    # In case it couldn't match all the stubs 
    leftover_candidates = filter(x -> nodes_dict[x] > 0, remaining_nodes)

    odd_nodes = filter(node -> nodes_dict[node] % 2 == 1 , leftover_candidates)
    
    if !isempty(odd_nodes) 

        remaining_edgelist = [[src(e), dst(e)] for e in collect(edges(graph))]

        odd_nodes_neighbours = unique(reduce(vcat, [neighbors(graph, node) for node in odd_nodes]))
        
        
        odd_filtered_edgelist = [(src, dst) for (src, dst) in remaining_edgelist if !(src in odd_nodes || dst in odd_nodes || src in odd_nodes_neighbours || dst in odd_nodes_neighbours)]

        while !isempty(odd_nodes)

            last_nodes = sample(odd_nodes,2,replace=false)

            odd_edges_options = sample(odd_filtered_edgelist,1,replace=false)
        
            last_candidates = last_nodes[1],last_nodes[2], odd_edges_options[1][1], odd_edges_options[1][2]
        
            while length(last_candidates) != length(unique(last_candidates))
        
                odd_edges_options = sample(odd_filtered_edgelist,1,replace=false)
        
                last_candidates = last_nodes[1],last_nodes[2], odd_edges_options[1][1], odd_edges_options[1][2]
            end
        
            add_edge!(graph, last_nodes[1], odd_edges_options[1][1])
            add_edge!(graph, last_nodes[2], odd_edges_options[1][2])
        
            rem_edge!(graph,odd_edges_options[1][1],odd_edges_options[1][2])

            nodes_dict[last_nodes[1]] -= 1 
            nodes_dict[last_nodes[2]] -= 1
            
            odd_nodes = filter(node -> nodes_dict[node] % 2 == 1 , odd_nodes)
        end

    end
   
    even_nodes = filter(x -> nodes_dict[x] > 0, remaining_nodes)
    
    while !isempty(even_nodes)

        each_even_node = sample(even_nodes,1,replace=false)[1]
    
        
        num_edges = div(nodes_dict[each_even_node], 2)
        
        while num_edges > 0
    
            last_edgelist = [[src(e), dst(e)] for e in collect(edges(graph))]
            
            even_nodes_neighbours = neighbors(graph, each_even_node)
    
            even_filtered_edgelist = [(src, dst) for (src, dst) in last_edgelist if !(src in even_nodes || dst in even_nodes || src in even_nodes_neighbours || dst in even_nodes_neighbours)]
    
            even_edges_options = sample(even_filtered_edgelist,1,replace=false)
            
            add_edge!(graph, each_even_node, even_edges_options[1][1])
            add_edge!(graph, each_even_node, even_edges_options[1][2])
    
            rem_edge!(graph,even_edges_options[1][1], even_edges_options[1][2])
            nodes_dict[each_even_node] -= 2  
    
            num_edges = div(nodes_dict[each_even_node], 2)
    
        end
    
        even_nodes = filter(x -> nodes_dict[x] != 0 , even_nodes)
      
    end
    
    return graph
end

