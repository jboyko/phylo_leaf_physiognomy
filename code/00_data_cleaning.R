setwd("/Users/jboyko/University\ of\ Michigan\ Dropbox/James\ Boyko/James\ Boyko’s\ files/Home/phylo_leaf_physiognomy")

# imports
library(ape)

# load and clean data
dat <- read.csv("data/RoyerLeafShapeClimateDataFixedNames_June2012.csv")
dat$genusSpecies[dat$genusSpecies == " Dialyanthera sp."] <- "Dialyanthera sp."
dat$genusSpecies <- gsub(" ", "_", dat$genusSpecies)
genera_names <- unlist(lapply(strsplit(dat$genusSpecies, "_"), function(x) x[1]))

# tree organization
tree <- ladderize(read.tree("data/best_wcvp.tre_dated"))
split_names <- sapply(tree$tip.label, function(x) strsplit(x, split = "_"))
for(i in 1:length(split_names)){
  if(length(split_names[[i]]) > 4){
    split_names[[i]] <- c(split_names[[i]][1:3], 
      paste(split_names[[i]][4:length(split_names[[i]])], collapse = "_"))
  }
}
name_table <- as.data.frame(do.call(rbind, split_names))
colnames(name_table) <-  c("order", "family", "genus", "species")
tree_sp <- tree
tree_sp$tip.label <- paste(name_table[,3], name_table[,4], sep = "_")
to_keep <- rownames(name_table)[!duplicated(name_table$genus)]
tree_gen <- keep.tip(tree, to_keep)
tree_gen$tip.label <- sapply(tree_gen$tip.label, 
  function(x) gsub("^[^_]+_[^_]+_([^_]+)_.*$", "\\1", x))

# site, species, and genera averages
dat_sp <- aggregate(dat[,c(9:ncol(dat))], 
  by = list(dat$genusSpecies, genera_names, dat$Family, dat$Order), 
  FUN = mean, na.rm=TRUE)
dat_sp <- dat_sp[!dat_sp$Group.4 == "unknown",]
dat_site <- aggregate(dat[,c(9:ncol(dat))], by = list(dat$Site), FUN = mean, na.rm=TRUE)
# write.csv(dat_site, "data/dat_site.csv", row.names = FALSE)

# prune trees and add tips to match data
for(i in 1:nrow(dat_sp)){
  cat("\r", round(i/nrow(dat_sp)*100, 1), "%")
  current_sp <- dat_sp$Group.1[i]
  if(current_sp %in% tree_sp$tip.label) next
  target_node <- NULL
  genus_tips <- tree_sp$tip.label[name_table$genus == dat_sp$Group.2[i]]
  if(length(genus_tips) >= 2) {
    target_node <- getMRCA(tree_sp, genus_tips)
  } else if(length(genus_tips) == 1) {
    target_node <- which(tree_sp$tip.label == genus_tips)
  }
  if(is.null(target_node)) {
    fam_tips <- tree_sp$tip.label[name_table$family == dat_sp$Group.3[i]]
    if(length(fam_tips) >= 2) target_node <- getMRCA(tree_sp, fam_tips)
  }
  if(is.null(target_node)) {
    ord_tips <- tree_sp$tip.label[name_table$order == dat_sp$Group.4[i]]
    if(length(ord_tips) >= 2) target_node <- getMRCA(tree_sp, ord_tips)
  }
  if(!is.null(target_node)) {
    h <- max(nodeHeights(tree_sp))
    node_dist <- nodeheight(tree_sp, target_node)
    edge_len <- h - node_dist
    tree_sp <- bind.tip(tree_sp, tip.label = current_sp, 
      where = target_node, edge.length = edge_len)
  }
}

# write.tree(tree_sp, file = "data/tre_add.tre")
tree_pruned <- keep.tip(tree_sp, dat_sp$Group.1)
# write.tree(tree_pruned, file = "data/tre_pruned.tre")

dat_pruned <- dat_sp[match(tree_pruned$tip.label, dat_sp$Group.1),]
# write.csv(dat_species, file = "data/data_cleaned.csv", row.names = FALSE)

write.csv(dat_pruned, file = "data/data_species.csv", row.names = FALSE)


