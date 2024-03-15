
## NETWORK PROPERTY CALCULATIONS ##

# Dataframe structure to store all calculations

df_properties =  DataFrame(
    n_nodes = Int64[],
    n_edges = Int64[],
    giant_component = Int64[],
    init_avg_deg = Any[],
    p_value = Any[],
    Degree = String[],
    Gini_degree = Float64[],
    Closeness = String[],
    Gini_close = Float64[],
    Betweenness = Float64[],
    Betweenness_per = String[], 
    Gini_between = Float64[],
    Eigenvector = String[],
    Gini_eigen = Float64[],
    Density = Float64[],
    Assortativity = Float64[],
    Clustering_coeff = Float64[],
    Global_transit = Float64[],
    Geodesic_mean = Float64[],
    Geodesic_gini = Float64[],
    Diameter = Float64[]
)

## Find the largest connected component before taking any measures ##
function main_component(g)
    c = connected_components(g)
    _, i = findmax(length.(c))
    g[c[i]]
end

# obtain Q1 - Q3 and round them up
function calculate_percentiles(var_input) 
    IQ1 = quantile(var_input, 25/100)
    IQ3 = quantile(var_input, 75/100)
    rounded_IQ1 = round(IQ1, digits = 3)
    rounded_IQ3 = round(IQ3, digits=3)

    return string.(rounded_IQ1, " - ",  rounded_IQ3)
end

# Calculate the geodesic
function calculate_geodesic(var_input)
    paths = Float64[] # have to specify float as if you leave it as Any Gini fn doesn't work 
    nodes = collect(1:nv(var_input))
    for each in nodes
       a = dijkstra_shortest_paths(var_input,each)
       append!(paths,a.dists) # use append so we have a 1 row vector and not n-length one
    end 
    return paths   
end

## Calculate centrality measures & append them ## 
function centrality_calc(network)

    largest_c = main_component(network)

    n_nodes = nv(network)
    n_edges = ne(network)
    giant_component = nv(largest_c)
    init_avg_deg = missing 
    p_value = missing # replace where appropriate
    degree_cent = calculate_percentiles(degree(largest_c))
    gini_deg = round(gini(degree(largest_c)),digits =3)
    close_cent = calculate_percentiles(closeness_centrality(largest_c)) 
    gini_close = round(gini(closeness_centrality(largest_c)),digits =3)
    betweenness_cent = round(mean(betweenness_centrality(largest_c)), digits=3)
    betweenness_cent_per = calculate_percentiles(betweenness_centrality(largest_c)) #  essential to get the right estimate
    gini_between = round(gini(betweenness_centrality(largest_c)),digits =3)
    eigen_cent = calculate_percentiles(eigenvector_centrality(largest_c))
    gini_eigen = round(gini(eigenvector_centrality(largest_c)),digits =3)
    density_value = round(density(largest_c),digits =3)
    assort = round(assortativity(largest_c),digits =3)
    clustering = round(mean(local_clustering_coefficient(largest_c)),digits =3)
    global_transit = round(global_clustering_coefficient(largest_c),digits =3)
    # modularity_measure = round(modularity(largest_c), digits = 3)
    geodesic_mean = round(mean(calculate_geodesic(largest_c)),digits =3)
    geodesic_gini =	round(gini(calculate_geodesic(largest_c)),digits =3)
    diameter_v = round(diameter(largest_c),digits =3)
    
    # append to dataframe
    push!(df_properties, (n_nodes,n_edges,giant_component, init_avg_deg, p_value, degree_cent, gini_deg, close_cent, gini_close, betweenness_cent,betweenness_cent_per, gini_between, 
    eigen_cent, gini_eigen, density_value, assort, clustering, global_transit, geodesic_mean,geodesic_gini, diameter_v)) 

    return df_properties
end


