#' @author GhuangChuangYu
#' @importFrom grDevices colorRampPalette
#' @keywords internal
# this is from `ggtree`
get_cols <- function (n){
    col <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
             "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
             "#ccebc5", "#ffed6f")
    col2 <- c("#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
              "#a6cee3", "#006d2c", "#4d4d4d", "#8c510a", "#d73027",
              "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f")
    col3 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
              "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
              "#ffff99", "#b15928")
    colorRampPalette(col2)(n)
}

#' @keywords internal 
setfactorlevels <- function(data, factorlist){
    factornames <- intersect(colnames(data), names(factorlist))
    if (length(factornames)>0){
            for(i in factornames){
                    data[[match(i,colnames(data))]] <- factor(data[[match(i, colnames(data))]], 
                               levels=as.vector(factorlist[[match(i,names(factorlist))]]))
            }
    }
    return(data)
}

#' @keywords internal
get_otudata <- function(obj){
    otudata <- obj@otu_table
    otudata <- data.frame(otudata, check.names=FALSE)
    if(obj@otu_table@taxa_are_rows){
        otudata <- data.frame(t(otudata), check.names=FALSE)
    }
    return (otudata)
}

#' @keywords internal
checkotu <- function(obj){
    if (is.null(obj@otu_table)){
    	stop("The otu table is empty!")
    }else{
    	otuda <- get_otudata(obj)
    	return(otuda)
    }
}

#' @keywords internal
checksample <- function(obj){
    if (is.null(obj@sam_data)){
    	stop("The sample_data is empty")
    }else{
    	sampleda <- get_sample(obj)
    	return(sampleda)
    }
}

#' @keywords internal.
get_sample <- function(obj){
    if (is.null(obj@sam_data)){
    	sampleda <- NULL
    }else{
    	sampleda <- data.frame(obj@sam_data, check.names=FALSE)
    }
    return(sampleda)
}

# #' @keywords internal
# #taxlevelchar <- c("k", "p", "c", "o", "f", "g", "s", "st")

newtaxname <- function(x, y){
    y <- as.vector(y)
    x[y] <- paste(taxlevelchar[y], x[y], sep="__un_")
    x
}

#' @importFrom zoo na.locf
#' @keywords internal
filltaxname <- function(taxdf){#, type="species"){
   # if (type != "species"){
   #     taxlevelchar <- paste0("d", seq_len(ncol(taxdf)))
   # }else{
   #     taxlevelchar <- taxlevelchar[seq_len(ncol(taxdf))]
   # }
    tmprownames <- rownames(taxdf)
    indexmark <- apply(taxdf, 2, function(x){nchar(x, keepNA = TRUE)})<=4
    taxdf[indexmark] <- NA
    indextmp <- apply(is.na(taxdf), 1, which)
    if(length(indextmp)==0){
        taxdf <- data.frame(taxdf, check.names=FALSE)
        return(taxdf)
    }
    taxdf <- apply(taxdf, 1, na.locf)
    taxdf <- lapply(seq_len(ncol(taxdf)), function(i) taxdf[,i])
    #newtaxname <- function(x, y){
    #    y <- as.vector(y)
    #    x[y] <- paste(taxlevelchar[y], x[y], sep="__un_")
    #    x
    #}
    taxdf <- data.frame(t(mapply(newtaxname, taxdf, indextmp)), 
                        stringsAsFactors=FALSE)
    rownames(taxdf) <- tmprownames
    return(taxdf)
}

#' @keywords internal
addtaxlevel <- function(taxdf){#, type="species"){
    #if (type != "species"){
    #    taxlevelchar <- paste0("d", seq_len(ncol(taxdf)))
    #}else{
    #    taxlevelchar <- taxlevelchar[seq_len(ncol(taxdf))]
    #}
    taxlevelchar <- taxlevelchar[seq_len(length(taxdf))]
    paste(taxlevelchar, taxdf, sep="__")
}

#' @importFrom tibble column_to_rownames
#' @keywords internal
fillNAtax <- function(taxdf, type="species"){
    #taxdf <- remove_unclassfied(taxdf)
    taxdf <- remove_na_taxonomy_rank(taxdf)
    if (type!="species"){
        assign("taxlevelchar", paste0("d", seq_len(ncol(taxdf))), envir = .GlobalEnv)
    }else{
        assign("taxlevelchar", c("k", "p", "c", "o", "f", "g", "s", "st"), envir = .GlobalEnv)
    }
    if (any(is.na(taxdf[,1]))){taxdf[is.na(taxdf[,1]),1] <- "Unknown"}
    if (!(grepl("^k__", taxdf[1,1]) || grepl("^d1__", taxdf[1,1]))){
    	tmprownames <- rownames(taxdf)
    	tmpcolnames <- colnames(taxdf)
        taxdf <- t(apply(taxdf, 1, as.character))
        taxdf[is.na(taxdf)] <- ""
        taxdf <- data.frame(t(apply(taxdf, 1, addtaxlevel)),
                            stringsAsFactors=FALSE)
    	rownames(taxdf) <- tmprownames
    	colnames(taxdf) <- tmpcolnames
    }
    taxdf <- filltaxname(taxdf)
    taxdf <- repduplicatedtaxcheck(taxdf) #%>% column_to_rownames(var="rowname")
    attr(taxdf, "fillNAtax") <- TRUE 
    return(taxdf)
}

#' @keywords internal
remove_unclassfied <- function(taxdf){
    taxdf[grepl.data.frame("Unclassified|uncultured|Ambiguous|Unknown|unknown|metagenome", taxdf)] <- NA
    return(taxdf)
}

grepl.data.frame <- function(pattern, x, ...){
    y <- if (length(x)) {
             do.call("cbind", lapply(x, "grepl", pattern=pattern, ...))
         }else{
             matrix(FALSE, length(row.names(x)), 0)
         }
    if (.row_names_info(x) > 0L)
        rownames(y) <- row.names(x)
    y
}


#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @keywords internal
duplicatedtaxcheck <- function(taxdf){
    if (ncol(taxdf)==1){return(taxdf)}
    taxdf <- taxdf %>% rownames_to_column()
    for (i in ncol(taxdf):3){
    	tmp <- split(taxdf,taxdf[,i])
    	for (j in seq_len(length(tmp))){
    		flag <- length(unique(as.vector(tmp[[j]][,i-1])))
    		if (flag > 1){
    			tmp[[j]][,i] <- paste(tmp[[j]][,i],tmp[[j]][,i-1],sep="_")
    		}
    	}
    	taxdf <- do.call("rbind",c(tmp, make.row.names=FALSE)) 
    }
    return(taxdf)
}

#' @keywords internal
repduplicatedtaxcheck <- function(taxdf){
    for (i in seq_len(7)){
    	taxdf <- duplicatedtaxcheck(taxdf) %>% 
		column_to_rownames(var="rowname")
    }
    return(taxdf)
}

## #' @keywords internal
## ## reference https://rdrr.io/cran/stackoverflow/man/match.call.defaults.html
## match.call.defaults <- function(fun) {
##     if (!is.na(fun)){
##         print(args(diff_analysis.data.frame))
##         args(diff_analysis.data.frame)
##     }else{
##         call <- evalq(match.call(expand.dots=TRUE), parent.frame(1))
##         formals <- evalq(formals(), parent.frame(1))
##         for(i in setdiff(names(formals), c(names(call)))){
##             call[i] <- list(formals[[i]])
##         }
##         match.call(sys.function(sys.parent()), call)
##     }
## }

#' @keywords internal
extract_args <- function(obj, arg){
    if (!"someparams" %in% methods::slotNames(obj)){
        stop("The object don't have someparams slot!")
    }else{
        args <- obj@someparams
        argres <- args[[arg]]
        return(argres)
    }
}

#' @importFrom utils globalVariables
utils::globalVariables('taxlevelchar')

ddply <- getFromNamespace("ddply", "plyr")

#' @importFrom stats sd 
# Adapted from Rmisc
summarySE <- function (data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                       conf.interval = 0.95, .drop = TRUE){
    length2 <- function(x, na.rm = FALSE) {
        if (na.rm)
            sum(!is.na(x))
        else length(x)
    }
    datac <- ddply(data, 
                   groupvars, 
                   .drop = .drop, 
                   .fun = function(xx, col, na.rm) {
          c(N = length2(xx[, col], na.rm = na.rm), mean = mean(xx[, col], na.rm = na.rm), 
            sd = sd(xx[, col], na.rm = na.rm))
          }, measurevar, na.rm)
    datac %<>% dplyr::rename(!!measurevar:="mean")
    datac$se <- datac$sd/sqrt(datac$N)
    ciMult <- qt(conf.interval/2 + 0.5, datac$N - 1)
    datac$ci <- datac$se * ciMult
    return(datac)
}

#' @importFrom stats qt
# reference to Rmisc
CI <- function (x, ci = 0.95, na.rm=FALSE){
    a <- mean(x, na.rm = na.rm)
    s <- sd(x, na.rm = na.rm)
    if (na.rm){
        x <- x[!is.na(x)]
    }
    n <- length(x)
    error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
    return(c(upper = a + error, mean = a, lower = a - error))
}

setGeneric(
    name = "drop.tip",
    def = function( object, tip, ... )
        standardGeneric("drop.tip")
)

setMethod("drop.tip", signature(object="treedata"),
          function(object, tip, ...) {
              drop.tip.treedata(object, tip, ...)
          })

setMethod("drop.tip", signature(object="phylo"),
          function(object, tip, ...){
              ape::drop.tip(object, tip, ...)
          })

drop.tip.treedata <- function(object, tip, ...){
    params <- list(...)
    if ("interactive" %in% names(params) && params$interactive){
        message("The interactive mode is not implemented for treedata object!")
        params$interactive <- FALSE
    }
    res <- build_new_labels(tree=object)
    tree <- res$tree
    old_and_new <- res$node2old_new_lab
    if(is.character(tip)){
        tip <- old_and_new[old_and_new$old %in% tip, "new"] %>% unlist(use.names=FALSE)
    }
    params$phy <- tree
    params$tip <- tip
    new_tree <- do.call(ape::drop.tip, params)

    if (is.null(new_tree)){
        return(new_tree)
    }

    trans_node_data <- old_new_node_mapping(tree, new_tree)
    object@phylo <- build_new_tree(tree=new_tree, node2old_new_lab=old_and_new)

    update_data <- function(data, trans_node_data) {
        data <- data[match(trans_node_data$old, data$node),]
        data$node <- trans_node_data$new
        return(data)
    }

    if (nrow(object@data) > 0) {
        object@data <- update_data(object@data, trans_node_data)
    }

    if (nrow(object@extraInfo) > 0) {
        object@extraInfo <- update_data(object@extraInfo, trans_node_data)
    }
    return (object)
}

#' @importFrom treeio Nnode
build_new_labels <- function (tree){
    node2label_old <- tree %>% as_tibble() %>% dplyr::select(c("node",
        "label"))
    if (inherits(tree, "treedata")) {
        tree <- tree@phylo
    }
    tree$tip.label <- paste0("t", seq_len(Ntip(tree)))
    tree$node.label <- paste0("n", seq_len(Nnode(tree)))
    node2label_new <- tree %>% as_tibble() %>% dplyr::select(c("node",
        "label"))
    old_and_new <- node2label_old %>% dplyr::inner_join(node2label_new,
        by = "node") %>% dplyr::rename(old = "label.x", new = "label.y")
    return(list(tree = tree, node2old_new_lab = old_and_new))
}

old_new_node_mapping <- function (oldtree, newtree){
    treelab1 <- oldtree %>% as_tibble() %>% dplyr::select(c("node",
        "label"))
    treelab2 <- newtree %>% as_tibble() %>% dplyr::select(c("node",
        "label"))
    node_map <- dplyr::inner_join(treelab1, treelab2, by = "label") %>%
        dplyr::select(c("node.x", "node.y")) %>% dplyr::rename(c(old = "node.x",
        new = "node.y"))
    return(node_map)
}

build_new_tree <- function (tree, node2old_new_lab){
    treeda <- tree %>% as_tibble()
    treeda1 <- treeda %>% dplyr::filter(.data$label %in% node2old_new_lab$new)
    treeda2 <- treeda %>% dplyr::filter(!(.data$label %in% node2old_new_lab$new))
    treeda1$label <- node2old_new_lab[match(treeda1$label, node2old_new_lab$new),
        "old"] %>% unlist(use.names = FALSE)
    treeda <- rbind(treeda1, treeda2)
    tree <- treeda[order(treeda$node), ] %>% as.phylo()
    return(tree)
}
