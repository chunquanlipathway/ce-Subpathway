###Depends：R2.15.2；igraph package；iSubpathwayMiner package
##compute all edge weight for all pathways 
library(iSubpathwayMiner);
library(igraph);

#dependent functions
source("getLocSubGraph3_2.R")
source("getOneNodePath3_2.R")

#get all metabolic and non-metabolic pathways in KEGG 
#and convert to simple graphs with KOs as nodes 
#KO nodes could keep pathway structure
g1<-getMetabolicKOKOUGraph();
g2<-getNonMetabolicKOKOUGraph();
graphList<-c(g1,g2);

#read DE genes and ceRNA nodes from files
input1<-read.csv("diffgenes_qvalue0.01.csv",header=TRUE)  #input, DE genes
symbolList<-getGeneFromSymbol(input1[,1])
input2<-read.csv("ceRNA nodes.csv",header=TRUE)  #input, ceRNA nodes
ceRNAList<-getGeneFromSymbol(input2[,1])
moleculeList<-unique(c(as.character(symbolList),as.character(ceRNAList)))  #union of DE genes and ceRNA nodes  

#get Pvalue1 of PCC for gene expression profile
Pvalue1<-read.csv("ceRNA pairs.csv",header=TRUE)            #input, gene gene Pvalue
#get ko number of mRNAs and the corresponding Pvalue from Pvalue1 result
cor_pvalue<-c() 
for(i in 1:dim(Pvalue1)[1]) 
{		        
    mRNA11<-getKOFromSymbol(Pvalue1[i,1])	
    mRNA12<-getKOFromSymbol(Pvalue1[i,2])		
	if(length(mRNA11)>0 && length(mRNA12)>0)   
	{     	 
	    for(j in 1:length(mRNA11))
	        {   
			    new_Pvalue<-cbind(mRNA11[j],mRNA12,Pvalue1[i,3],Pvalue1[i,4])	
		        cor_pvalue<-rbind(cor_pvalue,new_Pvalue) #ko ko Pvalue 
	        }             
    }	                   			
}

#get Pvalue2 of DE genes
Pvalue2<-input1          #gene gene Pvalue
#get ko number of mRNAs and the corresponding Pvalue from Pvalue2 result
dif_pvalue<-c()
for(i in 1:dim(Pvalue2)[1]) 
{		        
    mRNA<-getKOFromSymbol(Pvalue2[i,1])		
	if(length(mRNA)>0)   
	{     
	    new_Pvalue<-cbind(as.character(mRNA),Pvalue2[i,2])	
        dif_pvalue<-rbind(dif_pvalue,new_Pvalue) #ko ko Pvalue	            
    }	                   			
}

#delete the rows that the first column gene and the second column gene are mapped to the same ko number
#for Pvalue1 of PCC 
delete_location_cor<-c()        
for(i in 1:dim(cor_pvalue)[1])
{ 
   if(length(intersect(cor_pvalue[i,1],cor_pvalue[i,2]))>0)
    { 
     location<-i
     delete_location_cor<-c(delete_location_cor,location)
    }
}
cor_pvalue<-cor_pvalue[-delete_location_cor,] 

#get threshold m and identify significant subpathways(P<0.05) 
formulascores<-read.table("ce-scores.txt",header=F,sep="\t")     #input, all CorSP scores for all metabolic pathways from the output of Subpathway-CorSP code1
sortedformulascores<-sort(formulascores[,1])
m_loca<-floor(0.75*length(sortedformulascores))
m<-sortedformulascores[m_loca]
subGraphList<-getLocSubGraph3_2(cor_pvalue,dif_pvalue,value=m,moleculeList, graphList, type = "gene",s = 5, method = "shortestPaths", ignoreAmbiguousEnzyme = TRUE) 
ann<-identifyGraph(moleculeList,subGraphList,type="gene",background=getPrioBackground(type="gene"))
result<-printGraph(ann,detail=TRUE)
#write.table(result,file="all subpathways.txt",row.names=FALSE,sep="\t")
result_risk<-result[which(result[,"pvalue"]<0.05),]
write.csv(result_risk,file="risk subpathways.csv",row.names=FALSE)     #output, the significant subpathways

##可视化一个子通路
#plotAnnGraph("path:04510_12",subGraphList,ann,gotoKEGG=TRUE)
#plotAnnGraph("path:04510_12",subGraphList,ann,vertex.label=V(subGraphList["04910_2"][[1]])$names,gotoKEGG=TRUE)



