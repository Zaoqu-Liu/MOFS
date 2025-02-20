#' @title Get Binary Clusters from Clustering Results
#' @description This function extracts binary cluster assignments from multiple clustering results.
#' @param res A list containing clustering results
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#' @return A data frame where each row represents a binary encoding of cluster assignments across different methods.
#' @export
get.binary.clusters <- function(res) {
  resm <- lapply(res, function(d) {
    d2 <- model.matrix(~ 0 + d$Cluster2) %>% as.data.frame()
    rownames(d2) <- rownames(d)
    colnames(d2) <- substr(colnames(d2), 11, stop = 100)
    d2 <- as.data.frame(d2)
    return(d2)
  }) %>%
    BiocGenerics::Reduce(cbind, .) %>%
    t() %>%
    as.data.frame()
  return(resm)
}
