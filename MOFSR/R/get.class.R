#' @title Get Cluster Assignments
#' @description Extract cluster assignments from Consensus Clustering results for a specific number of clusters.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @param cc.res Consensus clustering results from ConsensusClusterPlus.
#' @param k Integer. The number of clusters to extract.
#' @return A data frame with the following columns:
#'   - ID: The sample identifier.
#'   - Cluster: The assigned cluster label, prefixed by 'C'.
#' @examples
#' data <- mtcars
#' cc_res <- RunCC(data)
#' clu <- get.class(cc_res, 2)
#' clu
#'
#' @export
get.class <- function(cc.res, k) {
  clu <- cc.res[[k]][["consensusClass"]]
  return(data.frame(
    ID = names(clu),
    Cluster = paste0("C", clu)
  ))
}
