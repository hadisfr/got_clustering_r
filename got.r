#!/usr/bin/env Rscript

# Algorithmic Graph Theory ~ F 97 :: CA2 ~ Q2
# by Hadi Safari (hadi@hadisafari.ir)

library(igraph)

# We could use `decompose` to make connected components seperate

################

find_clusters <- function(g) {
    E(g)$width <- E(g)$weight / max(E(g)$weight) * 16
    coords <- layout_with_fr(g)
    modularities <- data.frame()
    make_plot <- function(cluster, name) {
        print(paste("clustering using", name))
        png(paste(name, '.png', sep = ""), units="in", width=20, height=20, res=300)
        if (is.null(cluster))
            plot(g, layout = coords, vertex.size = 10)
        else {
            modularities <<- rbind(modularities, modularity(cluster))
            rownames(modularities) <<- c(head(rownames(modularities), -1), name)
            plot(cluster, g, layout = coords, vertex.size = 10)
        }
        dev.off()
    }
    make_plot(NULL, "pure_graph")
    make_plot(cluster_fast_greedy(g), "cluster_fast_greedy")
    make_plot(cluster_edge_betweenness(g), "cluster_edge_betweenness")
    make_plot(cluster_infomap(g), "cluster_infomap")
    make_plot(cluster_leading_eigen(g, options = setNames(c(1250), c("maxiter"))), "cluster_leading_eigen")
        # See https://github.com/igraph/igraph/issues/512#issuecomment-63691590 or increase maxiter
        # for "ARPACK error, Maximum number of iterations reached" in cluster_leading_eigen
    make_plot(cluster_louvain(g), "cluster_louvain")
    make_plot(cluster_optimal(g), "cluster_optimal")
    make_plot(cluster_spinglass(g), "cluster_spinglass")
    make_plot(cluster_walktrap(g), "cluster_walktrap")
    colnames(modularities) <- "modularity"
    return(modularities)
}

########

get_most_similar_vertices <- function(g, v_name) {
    v <- V(g)[v_name]
    v_index <- as.numeric(v)
    methods <- c("jaccard", "dice", "invlogweighted")
    df <- do.call(rbind.data.frame, lapply(methods, function(method) {
        inx <- which.max(similarity(g, method = method)[, v_index][-v_index])
        return(V(g)[inx + ifelse(inx < v_index, 0, 1)]$name)
    }))
    rownames(df) <- methods
    colnames(df) <- "Most Similar Vertex"
    return(df)
}

got <- read_graph("got.graphml", format = "graphml")
find_clusters(got)
get_most_similar_vertices(got, "Varys")
