#' Measure smoothness of the differnt feature of the data
#'
#' @param X
#' @param Paths
#' @param TargetPG
#' @param SmoothMode
#' @param CollMode
#' @param Partition
#' @param PrjStr
#' @param Randomize
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
MeasureSmoothness <- function(X,
                              Paths,
                              TargetPG,
                              SmoothMode = "lowess",
                              CollMode = "var",
                              Partition = NULL,
                              PrjStr = NULL,
                              Randomize = 0,
                              ...) {

  # X = tree_data
  # Paths = lapply(SelPaths[1:4], function(x){names(x)})
  # TargetPG = TreeEPG[[1]]
  # Partition = PartStruct$Partition
  # PrjStr = ProjStruct
  # Mode = "lowess"
  # Randomize = 0

  if(is.null(Partition)){
    Partition <- PartitionData(X = X, NodePositions = TargetPG$NodePositions)$Partition
  }

  if(is.null(PrjStr)){
    PrjStr <- project_point_onto_graph(X = X,
                                           NodePositions = TargetPG$NodePositions,
                                           Edges = TargetPG$Edges$Edges,
                                           Partition = Partition)
  }

  AllPt <- lapply(Paths, function(x){
    getPseudotime(ProjStruct = PrjStr, NodeSeq = as.integer(x))
  })


  AllData <- lapply(AllPt, function(x){

    Pt_rest <- x$Pt
    X_rest <- X[!is.na(Pt_rest),]

    Pt_rest <- Pt_rest[!is.na(Pt_rest)]
    Ordered_Pt <- order(Pt_rest)

    X_rest <- X_rest[Ordered_Pt,]
    Pt_rest <- Pt_rest[Ordered_Pt]


    if(SmoothMode == "lowess"){

      Diff <- apply(X_rest, 2, function(y){
        lowess(Pt_rest, y, ...)$y
      })

      if(CollMode == "sum"){
        return(apply(abs(Diff), 2, sum))
      }

      if(CollMode == "var"){
        return(apply(Diff, 2, var))
      }

      if(CollMode == "cv"){
        return(apply(Diff, 2, var)/abs(apply(Diff, 2, mean)))
      }

    }



    if(SmoothMode == "akima"){

      Diff <- apply(X_rest, 2, function(y){
        akima::aspline(Pt_rest, y, ...)$y
      })

      if(CollMode == "sum"){
        return(apply(abs(Diff), 2, sum))
      }

      if(CollMode == "var"){
        return(apply(Diff, 2, var))
      }

      if(CollMode == "cv"){
        return(apply(Diff, 2, var)/abs(apply(Diff, 2, mean)))
      }

    }


  })

  return(do.call(rbind, AllData))

}
