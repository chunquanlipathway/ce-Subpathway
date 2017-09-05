###dependent function 2
##traverse each node of interest and every other node, compute and output all edge weights for all pathways
getOneNodePath3_2<-function(cor_pvalue,dif_pvalue,value,current_node,other_nodes, pathway,all_shortest_paths_length,method = "shortestPaths") 
{
                #cor_pvalue--Pvalue of correlation coefficient
	            #dif_pvalue--Pvalue of DEgenes 
				#value--the threshold m of formula scores
				#current_node--a node of interest
				#other_nodes--the other nodes of interest related to the current node in the pathway
				#pathway--a reconstructed pathway graph
				#all_shortest_paths_length--the shortest path of the reconstructed pathway graph
                current_node_length <- length(current_node)
  		        other_nodes_length <- length(other_nodes)
   		        if (other_nodes_length < 1) 
                        warning("should have one node in other_nodes")
                if (current_node_length < 1) 
                        warning("should have one current_node")
                shortest_path_set <- list()
                current_node <- as.numeric(current_node)
                other_nodes <- as.numeric(other_nodes)                   
                if (other_nodes_length > 0)
				{   
                    for (i in 1:other_nodes_length) 
				    {
						#nested function, exist_path
                        withpath_node<-exist_path(current_node, other_nodes[i],all_shortest_paths_length) 
                        if (method == "shortestPaths") 
                        { 
                            if(withpath_node!="Inf"&&as.numeric(withpath_node)!=0)
                            {
                                shortest_path <- get.shortest.paths(pathway, current_node, other_nodes[i], mode = "out") 
                            }
                            else
                            {
                                shortest_path <- list()
                            }       
                        }
                            else if(withpath_node!="Inf"&&as.numeric(withpath_node)!=0)
				            {                                
                                shortest_path <- get.all.shortest.paths(pathway,current_node, other_nodes[i], mode = "out")   
                            }
                            else
                            {
                                shortest_path <- list()
                            }						 
                        if (length(shortest_path)>0)
					    {
                            shortest_path_length <- length(unlist(shortest_path[[1]]))-1  
							#nested function, nodep
                            pvalue<-nodep(current_node,other_nodes[i],cor_pvalue,dif_pvalue,pathway)                             
                            #nested function, merged_formula
                            if (as.numeric(pvalue)<0.05)
                            { 	
							    thp<-as.numeric(pvalue)								
                                score<-merged_formula(shortest_path_length,thp)                                                             
                                if (score>value)####threshold
                                {
                                    output_shortest_path_length<-as.numeric(shortest_path_length)
									shortest_path_set<- c(shortest_path_set,shortest_path)
                                    edge<-cbind(V(pathway)[current_node]$names,V(pathway)[other_nodes[i]]$names,thp,output_shortest_path_length,score,pathway$title)								
                                }
                            }
                        }
                    }              
                }  
                return(shortest_path_set)  
}
#if the shortest path between two nodes exist, then continue
exist_path<-function(current_node, other_nodes,all_shortest_paths_length)
{              
               result<-all_shortest_paths_length[current_node,other_nodes]
               return(as.character(result))
}
#convert the KO nodes to KOKO pairs and get the corresponding minimum pvalue occurred in the Pvalue table
nodep<-function(current_node,other_nodes,cor_pvalue,dif_pvalue,pathway)   
{ 
    current_node<-V(pathway)[current_node]$names   #one node containing one or more ko 
    other_nodes<-V(pathway)[other_nodes]$names     #one node containing one or more ko
    current_node_split<-unlist(strsplit(current_node,"[ ;]"))   #one or more ko splitted by ; at one node
	other_nodes_split<-unlist(strsplit(other_nodes,"[ ;]"))   #one or more ko splitted by ; at one node
    if(length(current_node_split)>0)
    {  
        pp<-c()
        for(i in 1:length(current_node_split))
        {                       
            if(length(other_nodes_split)>0)
            { 
                p<-psig(current_node_split[i],other_nodes_split,cor_pvalue,dif_pvalue)  #input some one ko at current node and one or more ko at other nodes to psig function
                pp<-c(pp,p)  
            }
            if(length(other_nodes_split)==0)
            { 
                pp<-c(pp,1)
            }
        }      
        p_value<-min(as.double(pp))   #get the minimum pvalue of some one ko pair as the pvalue of current_node and other_nodes		 
    } 
    if(length(current_node_split)==0)  
        {p_value<-1} 
    return(p_value)       #return "igraphID igraphID pvalue"
}
#get the minimum pvalue of some one ko pair as the final pvalue 
#only one ko at current node, one or more ko at other nodes
psig<-function(ko1,ko2,cor_pvalue,dif_pvalue)
{   
    finalp<-c()
    if(length(ko2)==1)    #if only one ko at other node (ko2)
    { 
        #for correlation pvalue
		loca1<-which(cor_pvalue[,1] %in% as.character(ko1))
        loca2<-which(cor_pvalue[,2] %in% as.character(ko2))
        inter1<-intersect(loca1,loca2)
        loca3<-which(cor_pvalue[,1] %in% as.character(ko2))
        loca4<-which(cor_pvalue[,2] %in% as.character(ko1))
        inter2<-intersect(loca3,loca4)
        if(length(inter1)>0 || length(inter2)>0)  #if they have the same location, then they are corrlation pairs in cor_pvalue table
	    {
		    inter<-c(inter1,inter2)
	        p_exp<-min(as.double(cor_pvalue[inter,3]))  #if both inter1 and inter2 exist, get the correlation pairs with the minimum pvalue 
	    }
        if(length(inter1)==0 && length(inter2)==0)  #when both inter1 and inter2 don't exist, give pvalue=1
           {p_exp<-1}
		   
	 	#for dif pvalue  
        loca1<-which(dif_pvalue[,1] %in% as.character(ko1))
        loca2<-which(dif_pvalue[,1] %in% as.character(ko2))      
        if(length(loca1)>0 && length(loca2)>0)
	    {
            loca<-c(loca1,loca2)
		    p_dif<-min(as.double(dif_pvalue[loca,2]))
	    }
        if(length(loca1)==0 || length(loca2)==0)
           {p_dif<-1}
		   
	    if(p_exp < 1 )                                
           {finalp<-p_exp*p_dif}
		if(p_exp==1  )
		   {finalp<-1 }
    }
    if(length(ko2)>1)   #if more than one ko at other node (ko2)
    { 
      store<-c()
      for (j in 1:length(ko2))
      { 
	    #for correlation pvalue
        loca1<-which(cor_pvalue[,1] %in% as.character(ko1))
        loca2<-which(cor_pvalue[,2] %in% as.character(ko2[j]))
        inter1<-intersect(loca1,loca2)
        loca3<-which(cor_pvalue[,1] %in% as.character(ko2[j]))
        loca4<-which(cor_pvalue[,2] %in% as.character(ko1))
        inter2<-intersect(loca3,loca4)
	    if(length(inter1)>0 || length(inter2)>0)
	    { 
		    inter<-c(inter1,inter2)
	        p_exp<-min(as.double(cor_pvalue[inter,3]))
	  	}
      	if(length(inter1)==0 && length(inter2)==0)
       	   {p_exp<-1}
		   
	    #for dif pvalue 
        loca1<-which(dif_pvalue[,1] %in% as.character(ko1))
        loca2<-which(dif_pvalue[,2] %in% as.character(ko2[j]))        
	    if(length(loca1)>0 && length(loca2)>0 )
	  	{ 
            loca<-c(loca1,loca2)
	  	    p_dif<-min(as.double(dif_pvalue[loca,2]))
	   	}
        if(length(loca1)==0 || length(loca2)==0)
           {p_dif<-1}
		   
		if(p_exp < 1 )   
           {p<-p_exp*p_dif}
		if(p_exp==1 )
		   {p<-1 }
        store<-c(store,p)  #store current node and several other nodes and the corresponding pvalues in several lines
      }
      finalp<-min(as.double(store))  #get the minimum pvalue of some one ko pair as the final pvalue 
    }    
    return(finalp)
}
#merge the shortest path length and the final pvalue to form formula  
merged_formula<-function(shortest_path_length,finalp)
{                    
    shortest_path_length<-as.numeric(shortest_path_length)
    finalp<-as.numeric(finalp)
	final_zscore<-qnorm(finalp,mean=0,sd=1,lower.tail = FALSE)  #convert p to z-score statistics
	merged_formula_value<-exp(-shortest_path_length/final_zscore)   #compute the formlula value
	return(as.numeric(merged_formula_value))
}
