library(sf)
library(FNN)
library(dbscan)
library(igraph)
library(data.table)

##### 1. Use 2D points to create a sparse, weighted, undirected graph #############################################################

# called by the dbscan_local_init function to create a graph
# called by split_area to split study area into more point density homogeneous parts
# requires as input:
# a matrix of coordinates (projected),
# a vector of IDs (length equal to nrow of coords),
# the maximum distance that a point should have at least one neighbor,
# the minimum number of points in the graph,
# the number of nearest neighbors to produce a sparse matrix
# and whether to decompose the graph into self-connected subgraphs

points2graph <- function(coords, id, max_dist, min_pnts, k, graph_decompose = T){
  
  if (nrow(coords) != length(id)){
    stop("the number of rows for coords and the length of id should be equal")
  }
  
  if (nrow(coords) <= k){
    k <- nrow(coords) - 1
  }
  
  
  # k-nearest neighbors
  knn <- get.knn(coords, k = k)
  
  # matrix with three columns,
  # col1 -> point id, col2 -> neighboring point id, col3 -> euclidean distance
  vertices <- do.call(rbind, lapply(1:nrow(coords),
                                    function(x) cbind(id[x],
                                                      id[knn$nn.index[x, ]],
                                                      knn$nn.dist[x, ])))
  
  if (any(is.na(vertices))){
    stop("NAs were created for vertices")
  }
  
  # Create undirected graph
  G <- graph.data.frame(vertices[,1:2], directed = F)
  
  # weight by distance
  E(G)$weight <- as.numeric(vertices[,3])
  
  # remove loops
  G <- simplify(G, edge.attr.comb="min")
  
  # remove edges with distance > max_dist
  G <- delete.edges(G, which(E(G)$weight >= max_dist))
  
  # delete points without neighbors
  idx_rm <- attributes(which(degree(G) < 1))$names
  if (length(idx_rm) > 0){
    G <- delete.vertices(G, idx_rm)
  }
  
  # if we want to split the graph into self-connected subgraphs
  if (graph_decompose){
    # return a list of those subgraphs that have at least as many points as min_pnts
    return(decompose.graph(G, min.vertices = min_pnts))
    }
  
  # otherwise return a list containing the graph
  list(G)
  
}

##### 2. Split by proximity ####################################################################################################

# example (assuming point data shapefile my_data.shp):
# shp <- st_read("my_data.shp", stringAsFactors = F)
# coords <- do.call(rbind, unclass(shp$geometry))
# shp$graph_id <- split_area(coords)

split_area <- function(coords, max_dist = 300, min_pnts = 10, k = 15){
  
  DT <- data.table(id = 1:nrow(coords))
  setkey(DT, id)
  
  graph_list <- points2graph(coords, DT$id, max_dist, min_pnts, k, graph_decompose = T)
  
  lookup <- as.data.frame(do.call(rbind, lapply(seq_along(graph_list), 
                                                function(graph_id) cbind(graph_id, as.integer(V(graph_list[[graph_id]])$name)))))
  
  setDT(lookup)
  setnames(lookup, names(lookup), c("graph_id", "id"))
  setkey(lookup, id)
  
  DT <- lookup[DT, on = "id"]
  DT[is.na(graph_id), graph_id := 0]
  
  DT$graph_id
  
}


##### 3. Get eps value as x% of k-NN distance ####################################################################################

calc_eps_clust <- function(coords_mat, k = 4, pct = 0.95){
  
  eps <- as.integer(round(quantile(knn.dist(coords_mat, k), pct)))
  
  eps
}

##### 4. Apply DBSCAN ##############################################################################################################

get_clustering <- function(coords_mat, eps_max = Inf, eps_min = -Inf, min_pnts = 10, border_points = T){
  
  if (nrow(coords_mat) < min_pnts){
    out <- list(eps_temp = -99, eps = -99, clustering = 0)
    return(out)
  }
  
  # get eps as 95% of 4-NN distance
  eps <- calc_eps_clust(coords_mat = coords_mat)

  if (eps > eps_max){
    eps <- eps_max
  } else if (eps < eps_min){
    eps <- eps_min
  }
  
  clustering_temp <- dbscan::dbscan(coords_mat, eps = eps, minPts = min_pnts, borderPoints = border_points)$cluster

  clusters_nr <- table(clustering_temp)
  
  # identify clusters with less than min_pnts members
  clusters_below_threshold <- as.integer(attributes(which(clusters_nr < min_pnts))$name)
  
  if (length(clusters_below_threshold) > 0){
    clustering_temp[clustering_temp %in% clusters_below_threshold] <- 0
  }

  list(eps = eps, clustering = clustering_temp)

}

##### 5. Get cluster id with density closer to overall density ###################################################################

check_clusters_density_nearest <- function(coords_mat, clustering_cl, eps_overall, k){
  
  # identify outliers
  idx_rm <- which(clustering_cl == 0)
  
  # if there are outliers remove them
  if (length(idx_rm) > 0){
    clustering_cl <- clustering_cl[-idx_rm]
    coords_mat <- coords_mat[-idx_rm, ]
  }
  
  # cluster names
  unique_cl <- unique(clustering_cl)
  
  # crosstab clusters
  clusters_nr <- table(clustering_cl)
  
  # eps values for each cluster
  all_eps <- sapply(seq_along(unique_cl), 
                    function(x) calc_eps_clust(coords_mat[clustering_cl == names(clusters_nr)[x], ], k))
  
  # id of the cluster with density (eps) closer to the overall eps
  idx_sele <- which.min(sapply(all_eps, function(x) abs(x - eps_overall)))
  
  # output is a list with the name of the cluster, the eps of the cluster,
  # the standard deviation of all eps, and the names of the clusters
  list(cl_sele = as.integer(names(clusters_nr)[idx_sele]),
       eps_cl_sele = all_eps[idx_sele],
       eps_sd = sd(all_eps),
       all_eps = all_eps,
       all_cl = unique_cl)
}

##### 6. Get neighboring points from a graph of points ###########################################################################

get_neighbors <- function(the_graph, point_ids){
  
  # get all point ids that share an edge with the input point_ids
  all_neighbors <- unique(unlist(lapply(as.character(point_ids), 
                                        function(pid) attributes(neighbors(the_graph, pid))$names)))
  
  all_neighbors
}

##### 7. Apply local DBSCAN ######################################################################################################

# this is a closure (a function that returns a function),
# when run from a data table using the by argument to loop over subgraphs
# we need to keep track of cluster names, otherwise the cluster name will start from 1 for each subgraph
# example:
# call dbscan_local_init with a starting value for the clustering (e.g. 0)
# dbscan_local <- dbscan_local_init(0)
# this returns a function with one parameter that is required,
# an object of class sfc (usually a column named geometry, see sf library)
# you can then use the dbscan_local function as follows (assuming that there is a column of subgraph ids in the data table)
# dt[, dbscan_local(geometry), by = graph_id]
# which will add four columns to the data table by reference,
# i.e. the clustering, the global eps, the local eps and the iteration id

dbscan_local_init <- function(cluster_id){
  
  max_cl <- cluster_id
  
  function(sfc_obj, eps_max = Inf, eps_min = -Inf, min_pnts = 10, border_points = T, std_thresh = 1, max_dist = 300, k = 15){
  
    # container for the output, of class data.table
    out_dt <- data.table(do.call(rbind, unclass(sfc_obj)), seq_along(sfc_obj), 0, 0, 0, 0)
    setnames(out_dt, 1:7, c("easting", "northing", "id", "clustering", "eps_overall", "eps_local", "iteration"))
    
    # build a sparse weighted (by euclidean distance) undirected graph from the input point data,
    # sparse edgelist is created using k-nearest neighbors
    # min_pnts is actually not required here but it is used when splitting a study area
    # remove edges with euclidean distance > max_dist, the point ids are added as attribute of the vertices
    graph_of_points <- points2graph(coords = as.matrix(out_dt[, 1:2]), id = out_dt$id, max_dist = max_dist,
                                    min_pnts = min_pnts, k = k, graph_decompose = F)
    
    setkey(out_dt, id)
    
    # count iterations of the algorithm
    count_iter <- 0
    
    # these are the points to cluster
    ids_to_visit <- out_dt$id
    
    # for as long as the number of points to cluster > min_pnts
    while (length(ids_to_visit) > min_pnts){
      
      count_iter <- count_iter + 1
      
      coords_mat <- as.matrix(out_dt[ids_to_visit, 1:2])
      
      # initial (global) dbscan
      dbscan_temp <- get_clustering(coords_mat = coords_mat, eps_max = eps_max, eps_min = eps_min,
                                    min_pnts = min_pnts, border_points = border_points)
      # crosstab clusters
      clusters_nr <- table(dbscan_temp$clustering)
      
      # remove outliers (cluster id == 0)
      clusters_nr <- clusters_nr[names(clusters_nr) != 0]
      
      # if there are two or more clusters, check how similar their densities are
      if (length(clusters_nr) > 1){
        
        # first get information about clusters' density (eps) and identify cluster with eps closer to overall eps
        cl_density <- check_clusters_density_nearest(coords_mat, dbscan_temp$clustering, dbscan_temp$eps, k=4)
        
        # get the point ids of the selected cluster
        cl_sele_id <- ids_to_visit[dbscan_temp$clustering == cl_density$cl_sele]
        
        # get all neighboring point ids (i.e. those that share an edge in the graph with the selected cluster points)
        neigh_pnts_id <- as.integer(get_neighbors(graph_of_points[[1]], cl_sele_id))
        
        # using the point ids get the neighboring cluster names
        neigh_cl <- unique(dbscan_temp$clustering[ids_to_visit %in% neigh_pnts_id])
        
        # remove from those the selected cluster name and the outliers name
        neigh_cl <- neigh_cl[! neigh_cl %in% c(cl_density$cl_sele, 0)]
        
        # this is the threshold to identify clusters with similar density, i.e.
        # standard deviation threshold (default is 1) * standard deviation of the eps values in the clustering
        eps_thresh <- std_thresh * cl_density$eps_sd
        
        # select clusters with similar density, i.e. eps in the range:
        # eps_cl <= eps_selected_cl + std_thresh * eps standard deviation & eps_cl <= eps_selected_cl - std_thresh * eps standard deviation
        neigh_cl <- neigh_cl[neigh_cl %in% cl_density$all_cl[which((cl_density$all_eps <= cl_density$eps_cl_sele + eps_thresh) & (cl_density$all_eps >= cl_density$eps_cl_sele - eps_thresh))]]
        
        # if there are clusters with similar density
        if (length(neigh_cl) > 0){
          
          # get their point ids
          neigh_cl_id <- ids_to_visit[dbscan_temp$clustering %in% neigh_cl]
          
          # and then include all the point ids that neighbor those clusters
          neigh_pnts_id <- unique(c(neigh_pnts_id, as.integer(get_neighbors(graph_of_points[[1]], neigh_cl_id))))
          
          # subset those points to be either members of the selected clusters or neighboring outliers
          sele_ids <- ids_to_visit[ids_to_visit %in% neigh_pnts_id & dbscan_temp$clustering %in% c(cl_density$cl_sele, neigh_cl, 0)]
          
          area_of_study <- out_dt[id %in% sele_ids]
          
          # if there aren't any clusters with similar density
        } else {
          
          # get the points that are members of the selected cluster or neighboring outliers
          sele_ids <- ids_to_visit[ids_to_visit %in% neigh_pnts_id & dbscan_temp$clustering %in% c(cl_density$cl_sele, 0)]
          
          area_of_study <- out_dt[id %in% sele_ids]
          
        }
        
        # apply dbscan in the area of homogeneous point density
        dbscan_fin <- get_clustering(coords_mat = area_of_study[,1:2], eps_max = eps_max, eps_min = eps_min,
                                     min_pnts = min_pnts, border_points = border_points)
        
        # reclassify cluster names
        dbscan_fin$clustering <- ifelse(dbscan_fin$clustering != 0, dbscan_fin$clustering + max_cl, 0)
        
        area_of_study[, c("clustering", "eps_overall", "eps_local", "iteration") := list(dbscan_fin$clustering, dbscan_temp$eps, dbscan_fin$eps, count_iter)]
        
        out_dt[area_of_study, c("clustering", "eps_overall", "eps_local", "iteration") := list(i.clustering, i.eps_overall, i.eps_local, i.iteration), on = "id"]
        
        ids_to_visit <- ids_to_visit[! ids_to_visit %in% area_of_study$id]
        
        # if there are less than 2 clusters, just use dbscan_temp
      } else {
        
        dbscan_temp$clustering <- ifelse(dbscan_temp$clustering != 0, dbscan_temp$clustering + max_cl, 0)
        
        out_dt[id %in% ids_to_visit, c("clustering", "eps_overall", "eps_local", "iteration") := list(dbscan_temp$clustering, dbscan_temp$eps, dbscan_temp$eps, count_iter)]
        
        ids_to_visit <- 0
      }
      
      max_cl <<- max(out_dt$clustering)
    }
    
    list(out_dt$clustering, out_dt$eps_overall, out_dt$eps_local, out_dt$iteration)
  
  }
  
}
