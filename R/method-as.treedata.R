#' @title convert dataframe contained hierarchical relationship or other classes to treedata class
#' @param data data.frame, such like the tax_table of phyloseq.
#' @param type character, the type of datasets, default is "species", if the dataset is not about species,                                                                                                         #' such as dataset of kegg function, you should set it to "others".
#' @param include.rownames logical, whether to set the row names as the tip labels, default is FALSE.
#' @param ..., additional parameters.
#' @return treedata class.
#' @author Shuangbin Xu
#' @importFrom tibble as_tibble
#' @importFrom tidytree treedata
#' @export
#' @examples
#' \dontrun{
#'   data(hmp_aerobiosis_small)
#'   head(taxda)
#'   treedat <- convert_to_treedata(taxda, include.rownames = FALSE)
#' }
convert_to_treedata <- function(data, type="species", include.rownames=FALSE, ...){
    if (!"fillNAtax" %in% names(attributes(data))){
        data <- fillNAtax(data, type=type)
    }
    if (include.rownames){
        data$OTU <- rownames(data)
    }
    data <- data.frame(Root=rep("r__root", nrow(data)), data)
    datalist <- list()
    clnm <- colnames(data)
    for (i in seq_len(ncol(data)-1)){
    	tmpdat <- data[,c(i, i+1)]
    	colnames(tmpdat) <- c("parent", "child")
        tmpdat %<>% dplyr::mutate(nodeClass=clnm[i+1], nodeDepth=i) %>%
                    dplyr::distinct()
    	datalist[[i]] <- tmpdat
    }
    datalist <- do.call("rbind", datalist)
    datalist <- datalist[!duplicated(datalist),]
    isTip <- !as.vector(datalist$child) %in% as.vector(datalist$parent)
    index <- rep(NA, length(isTip))
    index[isTip] <- seq(1,sum(isTip))
    index[!isTip] <- seq(sum(isTip)+2,length(isTip)+1)
    mapping <- data.frame(node=index, labelnames=as.vector(datalist$child), isTip)
    indxx <- match(mapping$labelnames, datalist$child)
    mapping$nodeClass <- datalist[indxx, "nodeClass"]
    mapping$nodeDepth <- datalist[indxx, "nodeDepth"]
    parentnode <- mapping[match(as.vector(datalist$parent), as.vector(mapping$labelnames)),]$node 
    childnode <- mapping[match(as.vector(datalist$child), as.vector(mapping$labelnames)),]$node
    edges <- cbind(parentnode, childnode) 
    colnames(edges) <- NULL
    edges[is.na(edges)] <- sum(isTip) + 1
    root <- data.frame(node=sum(isTip)+1,labelnames="r__root",
    		       isTip=FALSE, nodeClass="Root", nodeDepth=0)
    mapping <- rbind(root, mapping)
    mapping <- mapping[order(mapping$node),]
    node.label <- as.vector(mapping$labelnames)[!mapping$isTip]
    tip.label <- as.vector(mapping$labelnames)[mapping$isTip]
    mapping <- mapping[,colnames(mapping) %in% c("node", "nodeClass", "nodeDepth")]
    taxphylo <- structure(list(edge=edges, node.label=node.label,
                               tip.label=tip.label,
                               Nnode = length(node.label)), class="phylo")
    res <- treedata(phylo=taxphylo, data=as_tibble(mapping))
    return(res)
}

#' convert taxonomyTable to treedata
#'
#' @title as.treedata
#' @param tree object, This is for taxonomyTable class, 
#' so it should be a taxonomyTable.
#' @param include.rownames logical, whether to set the rownames of 
#' taxonomyTable to tip labels, default is FALSE.
#' @param ... additional parameters. 
#' @method as.treedata taxonomyTable
#' @rdname as.treedata
#' @export
#' @examples
#' \dontrun{
#'   data(test_otu_data)
#'   test_otu_data %<>% as.phyloseq()
#'   tree <- as.treedata(phyloseq::tax_table(test_otu_data), include.rownames = TRUE)
#' }
as.treedata.taxonomyTable <- function(tree, include.rownames = FALSE,...){
    convert_to_treedata(data.frame(tree, check.names=FALSE), include.rownames = include.rownames)
}

#' @method as.treedata treedata
#' @export
as.treedata.treedata <- function(tree, ...){
    return(tree)
}

# #' @method as.treedata tbl_mpse
# #' @importFrom dplyr left_join
# #' @export
# as.treedata.tbl_mpse <- function(tree, use_taxatree=TRUE, tiplevel="OTU", ...){
#     tr <- attr(tree, "otutree")
#     taxavar <- attr(tree, "taxavar")
#     if (!is.null(tr) && !use_taxatree){
#         treeda <- tr %>% as_tibble()
#         extrada <- tree %>% nest()
#         tiplevel <- "OTU"
#     }else{
#         if (use_taxatree && !is.null(taxavar)){
#             taxavar <- c(taxavar[taxavar!="OTU"], "OTU")
#             if (!is.null(taxavar)){
#                 indx <- which(taxavar==tiplevel)
#                 n <- length(taxavar)
#                 flag <- indx < n
#                 if (flag){
#                     taxavar <- taxavar[!taxavar %in% taxavar[indx+1:n]]
#                 }
#                 treeda1 <- tree %>% select(taxavar) %>% distinct()
#                 treeda <- convert_to_treedata(data=treeda1) %>% as_tibble() %>% select(-"nodeClass")
#                 extrada <- tree %>% as_tibble() %>% pivot_longer(cols=taxavar, names_to="nodeClass", values_to=tiplevel) 
#                 othervars <- colnames(tree)[!colnames(tree) %in% taxavar]
#                 params <- lapply(othervars, function(i)i)
#                 names(params) <- othervars
#                 params[[".data"]] <- extrada %>% as_tibble()
#                 extrada <- do.call("nest", params)
#             }else{
#                 stop("The tax table is empty in the object!")
#             }
#         }else{
#             stop("The tree slot is empty, you can use the taxa tree via set use_taxatree=TRUE")
#         }
#     }
#     treeda %<>% left_join(extrada, by=c("label"=tiplevel)) %>% as.treedata()
#     return(treeda)
# }
# 
# #' @method as.treedata grouped_df_mpse
# #' @export
# as.treedata.grouped_df_mpse <- function(tree, use_taxatree=TRUE, tiplevel="OTU", ...){
#     tree <- tree %>% ungroup()
#     #mutatevar <- attr(tree, "mutatevar")
#     #tree <- tree %>% select(-mutatevar)
#     treeda <- as.treedata(tree=tree, use_taxatree=use_taxatree, tiplevel=tiplevel, ...)
#     return(treeda)
# }
