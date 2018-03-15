#' CERES wrapper
#'
#' @param sg_path path to matrix with sgRNA sequences as row names, samples (i.e. replicates of cell line screens) as column names
#' @param cn_path path to matrix with genomic cut sites as row names, cell lines as column names
#' @param guide_locus_path path to data.frame with column names `Guide` and `Locus`. Values must match those corresponding to row names of sg_data and cn_data
#' @param locus_gene_path path to data.frame with column names `Locus` and `Gene`
#' @param replicate_map_path path to data.frame with column names `Replicate` and `CellLine`
#' @param run_id name for the CERES run
#' @param optimize toggle optimization routine
#' @param params list of parameters; must specify or use defaults when optimize = F
#' @param fit_efficacy boolean indicating whether sgRNA offsets and activity scores are computed
#'
#' 
#' @return A list of numeric vectors and matrices containing the results of the CERES fit.
#' @examples
#' 
#' @export wrap_ceres
#'
wrap_ceres_no_path <- function(sg_path, cn_path, guide_locus_path,
                       locus_gene_path, replicate_map_path, 
                       run_id, optimize=F, params=NULL, cl_subset=NULL,
                       fit_efficacy=T){
  
  

    sg_data <- sg_path
    cn_data <- cn_path
    guide_locus <- guide_locus_path
    locus_gene <- locus_gene_path
    replicate_map <- replicate_map_path
  
  
  if(!is.null(cl_subset)){
    if(is.character(cl_subset)){
      replicate_map %<>% dplyr::filter(CellLine %in% cl_subset)
    } else if(is.numeric(cl_subset)){
      if(cl_subset > length(unique(replicate_map[["CellLine"]]))){
        stop("Number of cell lines to sample exceeds number of cell lines in data. Please supply a smaller cl_subset...")
      }
      cls_to_fit <- sample(replicate_map[["CellLine"]], cl_subset, replace = F)
      replicate_map %<>% dplyr::filter(CellLine %in% cls_to_fit)
    }
  }
  
  default_params <- list(lambda_g=0.5,
                         lambda_o=0.01,
                         lambda_s=0.001,
                         n_segments=25,
                         validation_set=0)
  
  if(!optimize & is.null(params) & interactive()){
    continue_with_defaults <- readline(
      prompt = "Continue with default hyperparameters? (y/n) "
      )
    if(tolower(continue_with_defaults) == "n"){
      stop("Please set `optimize=F` or specify hyperparameters")
    } else if(tolower(continue_with_defaults) == "y") {
      params <- default_params
    }
  }
  
  if(!optimize){
    if(is.null(params)){
      params <- list()
    }
    params[["lambda_g"]] <- ifelse(is.null(params[["lambda_g"]]),
                                   default_params[["lambda_g"]],
                                   params[["lambda_g"]])
    params[["lambda_o"]] <- ifelse(is.null(params[["lambda_o"]]),
                                   default_params[["lambda_o"]],
                                   params[["lambda_o"]])
    params[["lambda_s"]] <- ifelse(is.null(params[["lambda_s"]]),
                                   default_params[["lambda_s"]],
                                   params[["lambda_s"]])
    params[["n_segments"]] <- ifelse(is.null(params[["n_segments"]]),
                                   default_params[["n_segments"]],
                                   params[["n_segments"]])
    params[["validation_set"]] <- ifelse(is.null(params[["validation_set"]]),
                                         default_params[["validation_set"]],
                                         params[["validation_set"]])
    params[["run_name"]] <- run_id
    
    res <- run_ceres(sg_data, cn_data, guide_locus, locus_gene,
                      replicate_map, 
                      params, fit_efficacy = T)
    
  } else {
    if(is.null(params)){
      params <- list()
    }
    if(is.null(params[["n_grid"]])){
      params[["n_grid"]] <- 5
    }
    
    # Set grid of run parameters
    lambda_g <- 10^seq(-2, 0, length.out = params[["n_grid"]])
    
    # Make parameter table
    param_df <- expand.grid(lambda_g=lambda_g) %>%
      data.frame(lambda_o=default_params[["lambda_o"]], 
                 lambda_s=default_params[["lambda_s"]], 
                 n_segments=default_params[["n_segments"]], validation_set=1,
                 run_name=paste0("hyperparam", 1:nrow(.)))
    
    # Run CERES on each parameter
    params[["lambda_g"]] <- param_df[["lambda_g"]] %>% as.numeric()
    params[["lambda_o"]] <- param_df[["lambda_o"]] %>% as.numeric()
    params[["lambda_s"]] <- param_df[["lambda_s"]] %>% as.numeric()
    params[["n_segments"]] <- param_df[["n_segments"]] %>% as.integer()
    params[["validation_set"]] <- param_df[["validation_set"]] %>% as.integer()
    params[["run_name"]] <- param_df[["run_name"]] %>% as.numeric()
    
    res <- run_ceres(sg_data, cn_data, guide_locus, locus_gene,
                             replicate_map, params, fit_efficacy)
  }
  
  return(res)
  
}