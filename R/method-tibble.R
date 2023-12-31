#' @method as_tibble phyloseq
#' @importFrom tibble rownames_to_column 
#' @importFrom dplyr left_join
#' @importFrom tidyr pivot_longer
#' @export
as_tibble.phyloseq <- function(x, ...){
    otuda <- checkotu(x) %>% 
             tibble::as_tibble(rownames="Sample") %>%
             pivot_longer(!"Sample", names_to = "OTU", values_to = "Abundance") 
    sampleda <- get_sample(x) %>%
                avoid_conflict_names() %>%
                rownames_to_column(var="Sample")
    if (ncol(sampleda)!=0){
        otuda <- otuda %>% left_join(sampleda, by="Sample", suffix=c("", ".y"))
        samplevar <- colnames(sampleda)
    }else{
        samplevar <- NULL
    }
    
    taxada <- as.data.frame(x@tax_table)
    if (!all(dim(taxada)==0)){
        taxavar <- colnames(taxada)
        taxatree <- try_convert_taxa(taxada)
        taxada %<>% avoid_conflict_names() %>% 
                    tibble::as_tibble(rownames="OTU")
        otuda <- otuda %>% left_join(taxada, by="OTU", suffix=c("", ".y"))
        fillNAtaxflag <- TRUE
    }else{
        taxavar <- NULL
        taxatree <- NULL
        fillNAtaxflag <- NULL
    }
    if (!is.null(x@phy_tree)){
        otutree <- x@phy_tree %>% as.treedata() 
    }else{
        otutree <- NULL
    }
    attr(otuda, "samplevar") <- samplevar
    attr(otuda, "taxavar") <- taxavar
    attr(otuda, "assaysvar") <- "Abundance"
    attr(otuda, "fillNAtax") <- fillNAtaxflag
    attr(otuda, "otutree") <- otutree
    attr(otuda, "taxatree") <- taxatree
    attr(otuda, "refseq") <- x@refseq
    class(otuda) <- c("tbl_mpse", class(otuda))
    return(otuda)
}

#' @method as_tibble grouped_df_mpse
#' @export
as_tibble.grouped_df_mpse <- function(x, ...){
    res <- NextMethod()
    res <- add_attr.tbl_mpse(x1=res, x2=x)
    res <- drop_class(x=res, class=c("grouped_df_mpse", "grouped_df"))
    return(res)
}

#' @method as_tibble MPSE
#' @export
as_tibble.MPSE <- function(x, ...){
    res <- .as_tibble_MPSE(x=x, ...)
    return (res)
}

.as_tibble_MPSE <- function(x, .subset = NULL){
    .subset = rlang::enquo(.subset)
    otuda <- extract_count_data(x)
    sampleda <- x %>%
                mp_extract_sample()

    if (ncol(sampleda)>1){
        otuda <- otuda %>% left_join(sampleda, by="Sample", suffix=c("", ".y"))
        samplevar <- colnames(sampleda)
    }else{
        samplevar <- "Sample"
    }
    otumeta <-
        SummarizedExperiment::rowData(x) %>%
        avoid_conflict_names() 
    if (ncol(otumeta)>0){
        otumetavar <- colnames(otumeta)
        otumeta %<>% tibble::as_tibble(rownames="OTU") %>% modify_AsIs_list()
        otuda <- otuda %>% left_join(otumeta, by="OTU", suffix=c("", ".y"))
    }else{
        otumetavar <- NULL
    }
    
    if (!is.null(x@taxatree)){
        check_installed('purrr', "for `as_tibble()` with MPSE object.")
        taxada <- taxatree_to_tb(x@taxatree) 
        uniqnm <- setdiff(colnames(taxada), colnames(otuda))
        taxada %<>% dplyr::select(uniqnm)
        taxada <- taxada[,!vapply(taxada, rlang::is_list, logical(1)),drop=FALSE]
        taxavar <- colnames(taxada)  
        taxada %<>% tibble::as_tibble(rownames="OTU")
        otuda <- otuda %>% left_join(taxada, by="OTU", suffix=c("", ".y"))
        fillNAtaxflag <- TRUE
    }else{
        taxavar <- NULL
        fillNAtaxflag <- FALSE 
    }
    
    if (!rlang::quo_is_null(.subset)){
        otuda <- otuda %>% dplyr::select(!!.subset)
        samplevar <- samplevar[samplevar %in% colnames(otuda)]
        taxavar <- taxavar[taxavar %in% colnames(otuda)]
        otumetavar <- otumetavar[otumetavar %in% names(otuda)]
    }
    attr(otuda, "samplevar") <- samplevar
    attr(otuda, "taxavar") <- taxavar
    attr(otuda, "otumetavar") <- otumetavar
    attr(otuda, "assaysvar") <- names(x@assays)
    attr(otuda, "fillNAtax") <- fillNAtaxflag
    attr(otuda, "otutree") <- x@otutree
    attr(otuda, "taxatree") <- x@taxatree
    attr(otuda, "refseq") <- x@refseq
    attr(otuda, "internal_attr") <- x %>% attr("internal_attr")
    class(otuda) <- c("tbl_mpse", class(otuda))
    return(otuda)
}

avoid_conflict_names <- function(dat, spename=NULL){
    cls <- colnames(dat) %>% 
           gsub("^OTU$", "OTU.x", .) %>%
           gsub("^Sample$", "Sample.x", .) %>%
           gsub("^Abundance$", "Abundance.x", .)
    if(!is.null(spename)){
        cls <- gsub(paste0("^", spename, "$"), paste0(spename, ".x"), cls)
    }
    colnames(dat) <- cls
    return(dat)
}

extract_count_data <- function(SE_object){
    allda <- SummarizedExperiment::assays(SE_object) %>% as.list()
    clnm <- names(allda)

    if (clnm[1] != "Abundance" || is.null(clnm)){
        clnm[1] <- "Abundance"
    }
    if (any(nchar(clnm)==0)){
       indx <- which(nchar(clnm)==0) 
       clnm[indx] <- paste0("Abund.", seq_len(length(indx)))
    }else if (length(clnm) != length(allda)){
       clnm[seq_len(length(allda))!=1] <- paste0("Abund.", seq_len(length(allda)-length(clnm)))
    }

    da <- mapply(trans_to_longer, 
              allda, 
              clnm, 
              SIMPLIFY=FALSE) %>%
          purrr::reduce(left_join, by=c("OTU", "Sample"), suffix = c("", ".y")) %>%
          suppressMessages()
    abunclnm <- da %>% colnames() %>% magrittr::extract(3)
    da %<>% dplyr::filter(!is.na(!!as.symbol(abunclnm)))
    return (da)
}

trans_to_longer <- function(.data, name){
    if (is.null(colnames(.data))){
        .data %<>% methods::as('matrix') %>% tibble::as_tibble(rownames="OTU") %>% 
            suppressWarnings()
    }else{
        .data %<>% methods::as('matrix') %>% 
           tibble::as_tibble(rownames = "OTU", .name_repair = "minimal")
    }
    .data %<>%
        tidyr::pivot_longer(!.data$OTU, names_to = "Sample", values_to=name)
    return(.data)
}

