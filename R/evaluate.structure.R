#' Produce a multidimensional dimension matrix
#'
#' @param X a data matrix
#' @param PrintGraph an ElPiGraph object
#' @param Start the Starting node. If NULL (the default), a random end point will be selected
#'
#' @return a list with two elements
#' \itemize{
#'  \item Branches, a list containing the (ordered) nodes composing the different branches
#'  \item DimMatrix, a matrix contianing the normalized position of points acoross branches. The columns
#'  correpsond to the different branches definezd by Branches. If a points is not assigned to a branch NA is reported.
#' }
#' @export
#'
#' @examples
BranchingDimension <- function(X, PrintGraph, Start = NULL) {

  # Define supporting structures
  Net <- ElPiGraph.R::ConstructGraph(PrintGraph = PrintGraph)
  ProjStruct <- project_point_onto_graph(X = X, NodePositions = PrintGraph$NodePositions, Edges = PrintGraph$Edges$Edges)

  # Fix starting point if needed
  if(is.null(Start)){
    Start <- sample(which(igraph::degree(Net)==1), 1)
  }

  # Get the branches
  AllBr <- GetSubGraph(Net = Net, Structure = 'branches')
  FixedBr <- list()

  RetStruct <- matrix(NA, nrow = nrow(X), ncol = length(AllBr))
  colnames(RetStruct) <- names(AllBr)

  for(i in 1:length(AllBr)){

    Paths <- igraph::get.shortest.paths(graph = Net, from = Start,
                                        to = names(AllBr[[i]][c(1, length(AllBr[[i]]))]))$vpath

    if( length(Paths[[1]]) < length(Paths[[2]]) ){
      FixedBr[[i]] <- AllBr[[i]]
    } else {
      FixedBr[[i]] <- rev(AllBr[[i]])
    }
    Pt <- getPseudotime(ProjStruct = ProjStruct, NodeSeq = FixedBr[[i]])

    RetStruct[which(!is.na(Pt$Pt)), i] <- Pt$Pt[!is.na(Pt$Pt)]/Pt$PathLen

  }

  return(list(Branches = FixedBr, DimMatrix = RetStruct))

}



