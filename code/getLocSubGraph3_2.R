###dependent function 1
##identify subpathways with s>5 from all pathways
getLocSubGraph3_2<-function(cor_pvalue,dif_pvalue,value,moleculeList,graphList,type = "gene",s = 5, method = "shortestPaths", ignoreAmbiguousEnzyme = TRUE) 
{
    #cor_pvalue--Pvalue of correlation coefficient
	#hyper_pvalue--Pvalue of Hypergeometric  
	#value--the threshold m of formula scores
	#moleculeList--DE genes and ceRNA nodes
	#graphList--the reconstructed pathway graphs
    if (typeof(moleculeList) != "character") {
        print("warning: your moleculeList must be 'character' vector. Because the type of your current moleculeList is not correct, it has been conveted arbitrarily using the function as.character().")
        moleculeList <- as.character(moleculeList)
    }
    if (method != "shortestPaths" && method != "allShortestPaths") 
        stop("the argument should be shortestPaths or allShortestPaths.")
    if (!exists("k2ri")) 
        initializeK2ri()
    subgraphList <- list()
    kk <- 0
    graphGeneNodeList <- character()
    graphCompoundNodeList <- character()
    graphGeneCompoundNodeList <- character()
    beforeOrg <- graphList[[1]]$org
    if (type == "gene" || type == "gene_compound") {
        if (beforeOrg == "ko") {
            graphGeneNodeList <- getKOFromGene(moleculeList)
        }
        else if (beforeOrg == "ec") {
            graphGeneNodeList <- getEnzymeFromGene(moleculeList, 
            ignoreAmbiguousEnzyme = ignoreAmbiguousEnzyme)
        }
        else if (unlist(strsplit(beforeOrg, "_"))[1] == getOrgAndIdType()[1]) {
            if (is.na(unlist(strsplit(beforeOrg, "_"))[2])) {
                graphGeneNodeList <- getKGeneFromGene(moleculeList)
            }
            else {
                if (unlist(strsplit(beforeOrg, "_"))[2] != getOrgAndIdType()[2]) {
                  stop("organism error")                                                                                                                                                   
                }
                else {
                  graphGeneNodeList <- moleculeList
                }
            }
        }
        else {
            stop(paste("graph ", i, "  error: it is not ec, ko, or org graph.",sep = ""))
        }
    }
    if (type == "compound" || type == "gene_compound") {
        graphCompoundNodeList <- paste("cpd:", moleculeList,sep = "")
    }
    if (type == "gene_compound") {
        graphMoleculeNodeList <- c(graphGeneNodeList, graphCompoundNodeList)
    }
    else if (type == "gene") {
        graphMoleculeNodeList <- graphGeneNodeList
    }
    else if (type == "compound") {
        graphMoleculeNodeList <- graphCompoundNodeList
    }
    if (length(graphList) > 0) {
        for (i in 1:length(graphList)) 
            {
              print(i)
              currentOrg <- graphList[[i]]$org
              if (currentOrg != beforeOrg) 
			  {
                if (type == "gene" || type == "gene_compound") {
                  if (currentOrg == "ko") {
                    graphGeneNodeList <- getKOFromGene(moleculeList)
                  }
                  else if (currentOrg == "ec") {
                    graphGeneNodeList <- getEnzymeFromGene(moleculeList,ignoreAmbiguousEnzyme = ignoreAmbiguousEnzyme)
                  }
                  else if (unlist(strsplit(currentOrg, "_"))[1] == getOrgAndIdType()[1]) {
                    if (is.na(unlist(strsplit(beforeOrg, "_"))[2])) {
                      graphGeneNodeList <- getKGeneFromGene(moleculeList)
                    }
                    else { 
                      if (unlist(strsplit(beforeOrg, "_"))[2] != getOrgAndIdType()[2]) {
                        stop("organism error")
                      }
                      else {
                        graphGeneNodeList <- moleculeList
                      }
                    }
                  }
                  else {
                    stop(paste("graph ", i, "  error: it is not ec, ko, or org graph.",sep = ""))
                  }
                }
                if (type == "compound" || type == "gene_compound") {
                  graphCompoundNodeList <- paste("cpd:", moleculeList,sep = "")
                }
                if (type == "gene_compound") {
                  graphMoleculeNodeList <- c(graphGeneNodeList,graphCompoundNodeList)
                }
                else if (type == "gene") {
                  graphMoleculeNodeList <- graphGeneNodeList
                }
                else if (type == "compound") {
                  graphMoleculeNodeList <- graphCompoundNodeList
                }
                beforeOrg <- currentOrg
            }
            nodes <- character()
            directed <- is.directed(graphList[[i]])
            hit <- sapply(V(graphList[[i]])$names, function(x) ifelse(any(unlist(strsplit(x,"[ ;]")) %in% graphMoleculeNodeList), TRUE, FALSE))     
            if (any(hit) == TRUE) {
                nodes <- as.character(V(graphList[[i]])[hit])
            }          #if moleculeList (DE genes and ceRNA nodes) appear in graphList (pathway), then the moleculeList map to nodes, and continue
            used_nodes <- c()
            subpathways_list <- list()
            subpathways_number <- 0
            all_shortest_paths_length <- shortest.paths(graphList[[i]], mode = "out")
            while (length(used_nodes) < length(nodes)) {
                shortest_path_set_in_subpathways <- list()
                unused_nodes <- setdiff(nodes, used_nodes)
                available_nodes <- unused_nodes[1]
                subpathway_nodes <- unused_nodes[1]
                while (length(available_nodes) > 0)
                    {
                        current_node <- available_nodes[1]
                        unused_nodes <- setdiff(nodes, used_nodes)
                        other_nodes <- setdiff(unused_nodes, current_node)
                        if (length(other_nodes) < 1 || length(current_node) <1) 
				         { }
                        else 
				            {    
                                #traverse each node of interest and every other node, compute and output all edge weights for all pathways					
                                shortest_path_set <- getOneNodePath3_2(cor_pvalue,dif_pvalue,value,current_node,other_nodes, graphList[[i]], all_shortest_paths_length, method = method) 
                                if(length(shortest_path_set)<1) 
                                {}
                                else
                                    {
                                        new_hit_nodes <- setdiff(intersect(unlist(shortest_path_set),nodes), used_nodes) 
                                        subpathway_nodes <- union(subpathway_nodes,new_hit_nodes) 
                                        available_nodes <- union(available_nodes,new_hit_nodes)
                                        shortest_path_set_in_subpathways <- c(shortest_path_set_in_subpathways,shortest_path_set)
                                    }
                            }
                        used_nodes <- union(used_nodes, current_node)
                        available_nodes <- setdiff(available_nodes,current_node)
                    }
                subpathways_number <- subpathways_number + 1
                subpathways <- list()
                subpathways[[1]] <- subpathway_nodes
                subpathways[[2]] <- shortest_path_set_in_subpathways
                subpathways_list[[subpathways_number]] <- subpathways
            }
            for (k in seq(subpathways_list)) {
                entry1 <- c()
                entry2 <- c()
                if (length(subpathways_list[[k]][[2]]) > 0) {
                  V <- as.integer(unique(unlist(subpathways_list[[k]][[2]])))
                  if (length(V) >= s) {
                    kk <- kk + 1
                    subgraphList[[kk]] <- induced.subgraph(graphList[[i]], V) 
                     
                    subgraphList[[kk]]$number <- paste(subgraphList[[kk]]$number,k, sep = "_") 
                      
                    names(subgraphList)[kk] <- subgraphList[[kk]]$number
                    }
                }
                else {
                  V <- as.integer(unique(subpathways_list[[k]][[1]]))
                  if (length(V) >= s) 
                    {
                       kk <- kk + 1
                       subgraphList[[kk]] <- induced.subgraph(graphList[[i]], V) 
                     
                       subgraphList[[kk]]$number <- paste(subgraphList[[kk]]$number,k, sep = "_") 
                      
                       names(subgraphList)[kk] <- subgraphList[[kk]]$number
                    }
                }
            }
        }
    }
    return(subgraphList)
}

