#########################################
##     Correcting Under-Reporting      ##
##              Functions              ##
#########################################
library(nimble)
# Function to compute kronecker production
kronecker_prod <- function(A, B) {
  C = matrix(data=0, nrow=nrow(A)*nrow(B), ncol = ncol(A)*ncol(B))
  m = nrow(A)
  n = ncol(A)
  p = nrow(B)
  q = ncol(B)
  for (r in 1:m) {
    for (s in 1:n) {
      for (v in 1:p) {
        for (w in 1:q){
          C[p*(r-1)+v,q*(s-1)+w] = A[r,s]*B[v,w];
        }
      }
      
    }
  }
  C
}

#A1 = matrix(data=c(1,3,2,4),nrow=2)
#B1 = matrix(data=c(0,6,5,7),nrow=2) 

#kronecker_prod(A1,B1)

#A2 = matrix(data=c(1,-2,-4,3,7,3),nrow=2)
#B2 = matrix(data=c(8,1,2,1,
#                   -9,-3,8,2,
#                   -6,-4,-8,-5,
#                   5,7,-3,-1),nrow=4)
#kronecker_prod(A2,B2)


# This function was found online at:
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

tran_var=function(x){(x-mean(x))/sd(x)}

#######################mungeCardata4stan############################
mungeCARdata4stan = function(adjBUGS,numBUGS) {
  N = length(numBUGS);
  nn = numBUGS;
  N_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:N) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
}
