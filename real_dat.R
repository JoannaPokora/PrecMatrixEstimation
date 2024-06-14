library(igraph)
library(dplyr)
library(tidyr)
library(glasso)
library(MetaboCrates)
library(cvCovEst)

real.dat <- readxl::read_excel("2023-10-12_Conc.xlsx", skip = 1)
metabo.dat <- read_data("2023-10-12_Conc.xlsx")
metabo.dat <- complete_data(metabo.dat, "halflimit", "limit", "limit")
dat <- attr(metabo.dat, "completed")

pathway <- read.csv("pathway.csv")
metabo.code <- read.csv("kegg_hmdb.csv") %>%
  select(c(1,3)) %>%
  unique()

first.lod.row <- min(which(grepl("LOD", real.dat[["Measurement time"]])))

metabolites <- attr(metabo.dat, "metabolites")

metabo.id <- real.dat[1:(first.lod.row-1),] %>%
  select(all_of(metabolites)) %>%
  pivot_longer(everything(),
               names_to = "metabo.name",
               values_to = "hmdb.code") %>%
  inner_join(metabo.code, by = join_by(hmdb.code == hmdb_id)) %>%
  select(metabo.name, kegg_cpd_id) %>%
  unique()

metabo.vals <- dat %>%
  select(all_of(metabolites))

lambdas <- c(0.005, 0.01, 0.5, 1)
glasso.prec <- lapply(lambdas, function(lambda){
  glasso(cov(metabo.vals), lambda)
  })

glasso.disc <- setNames(sapply(glasso.prec, function(glasso.prec.lam){
  prec.matr <- glasso.prec.lam[["wi"]]
  sum(prec.matr[lower.tri(prec.matr)] != 0)
  }), lambdas)
                      
shrink.prec <- solve(linearShrinkLWEst(metabo.vals))

shrink.disc <- sum(shrink.prec[lower.tri(shrink.prec)] != 0)

get_edges <- function(prec.matr, metabolites){
  as.data.frame(which(prec.matr != 0, arr.ind = TRUE)) %>%
    filter(row > col) %>%
    mutate(row = metabolites[row], col = metabolites[col]) %>%
    unique() %>%
    as.matrix()
}

shrink.edges <- get_edges(shrink.prec, metabolites)

glasso.edges <- lapply(glasso.prec,
                       function(res) get_edges(res[["wi"]], metabolites))

get_est_edges_in_pathway <- function(est.edges, pathway, metabo.id){
  edg <- apply(est.edges, 1, function(est.edge){
    kegg.id.x <- metabo.id[which(metabo.id[,1] == est.edge[1]), 2]
    kegg.id.y <- metabo.id[which(metabo.id[,1] == est.edge[2]), 2]
    kegg.id <- expand_grid(x = kegg.id.x, y = kegg.id.y)
    apply(kegg.id, 1, function(est.edge){
      apply(pathway, 1, function(path.edge){
        if(all(path.edge %in% est.edge)) return(path.edge)
        else NULL
      })
    })
  })
  
  edg.unl <- unlist(edg)
  data.frame(x = edg.unl[seq(1, length(edg.unl)/2-1, 2)],
             y = edg.unl[seq(2, length(edg.unl)/2, 2)]) %>%
    unique()
}

shrink.path.edges <- get_est_edges_in_pathway(shrink.edges, pathway, metabo.id)

glasso.path.edges <- lapply(glasso.edges, function(glasso.edg){
  get_est_edges_in_pathway(glasso.edg, pathway, metabo.id)
})

pathway.graph <- pathway %>%
  as.matrix() %>%
  graph_from_edgelist()

pathway.layout <- layout_with_kk(pathway.graph)

pathway %>%
  as.matrix() %>%
  graph_from_edgelist() %>%
  plot(layout = pathway.layout, edge.arrow.size = 0.4, edge.color = "#0b2206",
       vertex.color = "#30CB14", vertex.size = 8,
       vertex.frame.color = "#0b2206", vertex.label = NA)

plot_found_pathway <- function(est.edges, pathway.graph, pathway.layout){
  E(pathway.graph)$found <- do.call(paste0, pathway) %in%
    do.call(paste0, est.edges)
  
  V(pathway.graph)$found <- names(V(pathway.graph)) %in%
    unlist(est.edges)
  
  pathway.graph %>%
    plot(layout = pathway.layout, edge.arrow.size = 0,
         edge.color = c("#A5A5A5", "#0b2206")[1+E(pathway.graph)$found],
         vertex.color = c("#A5A5A5", "#30CB14")[1+V(pathway.graph)$found],
         vertex.frame.color = c("#A5A5A5", "#0b2206")[1+V(pathway.graph)$found],
         vertex.size = 8, vertex.label = NA) 
}

plot_found_pathway(shrink.path.edges, pathway.graph, pathway.layout)

lapply(glasso.path.edges, function(glasso.edg){
  plot_found_pathway(glasso.edg, pathway.graph, pathway.layout)
})
