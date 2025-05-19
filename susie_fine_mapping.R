args <- commandArgs(trailingOnly = TRUE)
model_number <- as.integer(args[1])

library(genio)
library(data.table)
library(susieR)
library(pgenlibr)
library(WGCNA)
library(stringr)
library(dplyr)

susie_launcher <- function (
    model.no=5, 
    file_path_R_sim = "/project/statgen/models_for_final_simulations/model_%d/model_%d_with_R_sim.rds",
    file_path_simulation = "/project/statgen/models_for_final_simulations/model_%d/model_%d_with_R_int_R_ext_and_effects.rds",
    max_iter_external_susie = 3000,
    max_iter_internal_susie = 1000,
    max_iter_simulated_susie = 1000,
    fixed_LDs = c("internal_LD", "external_LD"), 
    sim_LDs = c(10, 50, 1000, 3000, 5000),
    n_individuals = c(29957, 5000)
    ) 
  {
  
  setwd(sprintf("/project/statgen/models_for_final_simulations/model_%d", model.no))
  
  log_file <- "warnings.log"
  
  if (file.exists(log_file)) file.remove(log_file)  
  
  withCallingHandlers({
    source("/home/giulia.desanctis/final_stage/scripts/functions_for_susie.R")
    
    simulation <- readRDS(sprintf(file_path_simulation, model.no, model.no))
    R_sim.no <- ifelse(model.no %% 100 == 0, 100, model.no %% 100)
    R_simulated <- readRDS(sprintf(file_path_R_sim, R_sim.no, R_sim.no))
    
    R_internal <- simulation$R_internal
    R_external <- simulation$R_external
    sum_stats <- simulation$sum_stats
    
    # Do internal and external susie
    int_start <- Sys.time()
    tryCatch({
    do_susie(sim=simulation, R_to_use = R_internal, est_res_var=TRUE, name_of_LD="internal", model_no = model.no, num_ind = simulation$N_int, max_iter=max_iter_internal_susie)
    }, warning = function(w) {
      cat(sprintf("[%s] WARNING in model %s (internal): %s\n", Sys.time(), model.no, conditionMessage(w)), 
          file = log_file, append = TRUE)
    })
    int_end <- Sys.time()
    time_int <- int_end - int_start
    print(sprintf("Time taken internal susie: %f", time_int))
    
    ext_start <- Sys.time()
    tryCatch({
      do_susie(sim=simulation, R_to_use = R_external, est_res_var=FALSE, name_of_LD="external", model_no = model.no, num_ind = simulation$N_ext, max_iter=max_iter_external_susie)
    }, warning = function(w) {
      cat(sprintf("[%s] WARNING in model %d (external): %s\n", Sys.time(), model.no, conditionMessage(w)), 
          file = log_file, append = TRUE)
    })
    ext_end <- Sys.time()
    time_ext <- ext_end - ext_start
    print(sprintf("Time taken external susie: %f", time_ext))
    
    # Do susie with R_sim for 30k individuals
    for(i in seq_len(length(R_simulated$matrix_names_list) - 1)) {
      start_sim <- Sys.time()
      tryCatch({
        do_susie(
        sim=simulation, 
        R_to_use = R_simulated$R_sim_zscore_list[[5]], 
        est_res_var=FALSE, 
        name_of_LD="simulated", 
        model_no = model.no, 
        num_ind = simulation$N_int, 
        max_iter=max_iter_simulated_susie,
        num_phenos = (dim(R_simulated$sum_stats_list[[5]])[2]-8)/2)
      }, 
        warning = function(w) {
        cat(sprintf("[%s] WARNING in model %d (%s): %s\n", Sys.time(), model.no, R_simulated$matrix_names_list[[i]], conditionMessage(w)), 
           file = log_file, append = TRUE)
      })
  
      print(sprintf("done susie for %s", R_simulated$matrix_names_list[[i]]))
      end_sim <- Sys.time()
      time_sim <- end_sim - start_sim
      print(sprintf("Time taken %s susie: %f", R_simulated$matrix_names_list[[i]], time_sim))
    }
      
    print("Finsihed with 10 phenotypes")
    
    # Do susie with R_sim for 5k individuals
    five_k_start <- Sys.time()
    ind = length(R_simulated$matrix_names_list)
    tryCatch({
    do_susie(sim=simulation, R_to_use = R_simulated$R_sim_zscore_list[[ind]], 
             est_res_var=FALSE, name_of_LD="simulated", model_no = model.no, 
             max_iter=max_iter_internal_susie,  num_ind = (dim(R_simulated$sum_stats_list[[ind]])[2]-8)/2, 
      num_phenos=(dim(R_simulated$sum_stats_list[[ind]])[2]-8)/2)
    }, warning = function(w) {
      cat(sprintf("[%s] WARNING in model %d (%s): %s\n", Sys.time(), model.no, R_simulated$matrix_names_list[[ind]], conditionMessage(w)), 
          file = log_file, append = TRUE)
    })
    
    print(sprintf("done susie for %s", R_simulated$matrix_names_list[[ind]]))
    five_k_end <- Sys.time()
    five_k_time <- five_k_end - five_k_start
    print(sprintf("Time taken %s susie: %f", R_simulated$matrix_names_list[[ind]], five_k_time))
    
    print("Finsihed with 5000 phenotypes and 5000 individuals")
    
    # Load susie objects to do fine-mapping 
    base <- sprintf("susie_model_%d_with_", model.no)
    
    # Generate suffixes for simulated LDs
    sim_suffixes <- expand.grid(ld = sim_LDs, n = n_individuals)
    sim_suffixes <- apply(sim_suffixes, 1, function(row) {
      paste0("simulated_LD_", row["ld"], "_phenos_", row["n"], "_individuals")
    })
    
    # All model types (LD identifiers)
    ld_types <- c(fixed_LDs, sim_suffixes)
    
    # Generate full file names
    file_names <- paste0(base, ld_types, ".rds")
    
    # Filter only existing files
    existing_files <- file_names[file.exists(file_names)]
    existing_names <- gsub(paste0(base, "|\\.rds$"), "", existing_files)
    
    # Read only existing files
    models <- lapply(existing_files, readRDS)
    names(models) <- existing_names
    
    # I create a data frame to save the results 
    
    results_df <- data.frame(
      model_name    = names(models),
      coverage      = NA_real_,
      power         = NA_real_,
      coverage_num  = I(vector("list", length(models))),
      coverage_den  = I(vector("list", length(models))),
      power_num     = I(vector("list", length(models))),
      power_den     = I(vector("list", length(models))),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(models)) {
      coverage_output <- compute_coverage(effect_list = simulation$effect_list, susie_out = models[[i]])
      power_output    <- compute_power(effect_list = simulation$effect_list, susie_out = models[[i]])
      
      results_df$coverage[i] <- sum(coverage_output$coverage_num) / sum(coverage_output$coverage_den)
      results_df$power[i]    <- sum(power_output$power_num) / sum(power_output$power_den)
      
      results_df$coverage_num[[i]] <- coverage_output$coverage_num
      results_df$coverage_den[[i]] <- coverage_output$coverage_den
      results_df$power_num[[i]]    <- power_output$power_num
      results_df$power_den[[i]]    <- power_output$power_den
    }
    
    saveRDS(results_df, sprintf("susie_metrics_%d.rds", model_number))
    
  }, warning = function(w) {
    cat("WARNING: ", conditionMessage(w), "\n", file = log_file, append = TRUE)
  })
  
}


##
susie_launcher(
    model.no = model_number
  , file_path_R_sim = "/project/statgen/models_for_final_simulations/model_%d/model_%d_with_R_sim.rds"
  , file_path_simulation = "/project/statgen/models_for_final_simulations/model_%d/model_%d_with_R_int_R_ext_and_effects.rds"
  , max_iter_external_susie = 10000
  , max_iter_internal_susie = 2000
  , max_iter_simulated_susie = 5000
  , fixed_LDs = c("internal_LD", "external_LD")
  , sim_LDs = c(10, 50, 1000, 3000, 5000)
  , n_individuals = c(29957, 5000)
)














