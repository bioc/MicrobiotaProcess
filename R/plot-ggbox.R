#' @title A box or violin plot with significance test
#' @param obj object, alphasample or data.frame (row sample x column features).
#' @param sampleda data.frame, sample information if obj is data.frame, the 
#' sampleda should be provided.
#' @param factorNames character, the names of factor contained in sampleda.
#' @param indexNames character, the vector character, should be the names of 
#' features contained object.
#' @param geom character, "boxplot" or "violin", default is "boxplot".
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param compare logical, whether test the features among groups,default is TRUE.
#' @param testmethod character, the method of test, default is `wilcox.test`.
#' see also \code{\link[ggsignif]{stat_signif}}.
#' @param signifmap logical, whether the pvalue are directly written a annotaion
#' or asterisks are used instead, default is (pvalue) FALSE. see also
#' \code{\link[ggsignif]{stat_signif}}.
#' @param p_textsize numeric, the size of text of pvalue or asterisks, 
#' default is 2.
#' @param step_increase numeric, see also \code{\link[ggsignif]{stat_signif}},
#' default is 0.1.
#' @param boxwidth numeric, the width of boxplot when the geom is 'violin',
#' default is 0.2.
#' @param facetnrow integer, the nrow of facet, default is 1.
#' @param controlgroup character, the names of control group, if it was set, the other groups 
#' will compare to it, default is NULL.
#' @param comparelist list, the list of vector, default is NULL.
#' @param ... additional arguments, see also \code{\link[ggsignif]{stat_signif}}.
#' @return a 'ggplot' plot object, a box or violine plot.
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#' library(magrittr)
#' otudafile <- system.file("extdata", "otu_tax_table.txt",
#'                          package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", 
#'                     header=TRUE, row.names=1,
#'                     check.names=FALSE, skip=1, 
#'                     comment.char="")
#' samplefile <- system.file("extdata",
#'                           "sample_info.txt", 
#'                           package="MicrobiotaProcess")
#' sampleda <- read.table(samplefile,
#'                        sep="\t", header=TRUE, row.names=1)
#' otuda <- otuda[sapply(otuda, is.numeric)] %>% t() %>%
#'          data.frame(check.names=FALSE)
#' set.seed(1024)
#' alphaobj1 <- get_alphaindex(otuda, sampleda=sampleda)
#' p1 <- ggbox(alphaobj1, factorNames="group")
#' data(test_otu_data)
#' test_otu_data %<>% as.phyloseq()
#' set.seed(1024)
#' alphaobj2 <- get_alphaindex(test_otu_data)
#' class(alphaobj2)
#' head(as.data.frame(alphaobj2))
#' p2 <- ggbox(alphaobj2, factorNames="group")
#' # set factor levels.
#' p3 <- ggbox(obj=alphaobj2, factorNames="group", 
#'             factorLevels=list(group=c("M", "N", "B", "D")))
#' # set control group.
#' p4 <- ggbox(obj=alphaobj2, factorNames="group", controlgroup="B")
#'  set comparelist
#' p5 <- ggbox(obj=alphaobj2, factorNames="group", 
#'             comparelist=list(c("B", "D"), c("B", "M"), c("B", "N")))
#' }
setGeneric("ggbox", function(obj, factorNames, ...){standardGeneric("ggbox")})

#' @aliases ggbox,data.frame
#' @rdname ggbox
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin aes_string facet_wrap 
#' @importFrom ggsignif geom_signif
#' @importFrom rlang .data
#' @export
setMethod("ggbox", "data.frame", 
          function(obj, sampleda, factorNames, indexNames, geom="boxplot",
                   factorLevels=NULL, compare=TRUE, testmethod="wilcox.test",
                   signifmap=FALSE, p_textsize=2, step_increase=0.1, boxwidth=0.2,
                   facetnrow=1, controlgroup=NULL, comparelist=NULL, ...){
    if (missing(sampleda) || is.null(sampleda)){
        stop("the sampleda should be provided!")
    }
    if (missing(factorNames) || is.null(factorNames)){
        stop("the factorNames should be provided!")
    }
    sampleda <- sampleda[,match(factorNames, colnames(sampleda)), drop=FALSE]
    obj <- merge(obj, sampleda, by=0)
    rownames(obj) <- obj$Row.names
    obj$Row.names <- NULL
    obj <- obj %>% tidyr::pivot_longer(!factorNames, names_to="feature", values_to="value")
    if (!is.null(factorLevels)){
        obj <- setfactorlevels(obj, factorLevels)
    }
    if (missing(indexNames)||is.null(indexNames)){
        indexNames <- unique(as.vector(obj$feature))
    }
    obj <- obj %>% dplyr::filter(.data$feature %in% indexNames)
    if (is.null(comparelist)){
        comparelist <- get_comparelist(data=obj, classgroup=factorNames, controlgroup=controlgroup)
    }
    mapping <- aes_string(x=factorNames, y="value", fill=factorNames)
    message("The color has been set automatically, you can reset it manually by adding scale_fill_manual(values=yourcolors)")
    p <- ggplot(data=obj,mapping)
    ifelse(geom=="boxplot",p <- p + geom_boxplot(outlier.size=0.5,outlier.shape=21),
           p <- p + geom_violin(trim=FALSE)+
                geom_boxplot(outlier.size=0.5,outlier.shape=21, 
			     width=boxwidth,
                             fill="white", show.legend=FALSE))
    if (compare){
        p <- p + geom_signif(comparisons = comparelist,
			     test = testmethod,
                             map_signif_level = signifmap,
                             textsize = p_textsize, 
                             step_increase = step_increase,
                             ...)
    }
    p <- p + facet_wrap(~feature, scales="free", nrow=facetnrow) +
         theme_bw() + xlab(NULL) + ylab(NULL)
    return(p)
})

#' @aliases ggbox,alphasample
#' @rdname ggbox
#' @export
setMethod("ggbox", "alphasample", function(obj,factorNames,...){
    alphatab <- obj@alpha
    sampleda <- obj@sampleda
    p <- ggbox(obj=alphatab,
               sampleda=sampleda,
               factorNames=factorNames,
               ...)
    return(p)
})

#' @keywords intermal
set_newlevels <- function(data, newlevels, factorNames){
    oldlevels <- unique(as.vector(data[[factorNames]]))
    newlevels <- oldlevels[match(newlevels, oldlevels)]
    newlevels <- newlevels[!is.na(newlevels)]
    data[[factorNames]] <- factor(data[[factorNames]], levels=newlevels)
    return(data)
}

#' @keywords internal
get_comparelist <- function(data, classgroup, controlgroup){
    groups <- get_classlevels(sampleda=data, classgroup=classgroup)
    if (!is.null(controlgroup)){
        groups <- setdiff(groups, controlgroup)
        tmplen <- length(groups)
        comparelist <- matrix(data=c(rep(controlgroup, tmplen), groups), nrow=tmplen)
    }else{
        comparelist <- get_compareclass(classlevels=groups)
    }
    comparelist <- split(comparelist, slice.index(comparelist, 1))
    names(comparelist) <- NULL
    return(comparelist)
}

