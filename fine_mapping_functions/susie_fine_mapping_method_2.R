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
    #file_path_R_sim = "/project/statgen/models_for_final_simulations/model_%d/model_%d_with_R_sim.rds",
    file_path_R_sim = "~/test_30_05/5k_individuals/model_%d_with_R_sim.rds",
    file_path_simulation = "~/project/statgen/models_for_final_simulations/model_%d/model_%d_with_R_int_R_ext_and_effects.rds",
    max_iter_external_susie = 10000,
    max_iter_internal_susie = 2000,
    max_iter_simulated_susie = 5000,
    fixed_LDs = c("internal_LD", "external_LD"), 
    sim_LDs = c(1000, 2000, 3000, 4000, 5000),
    list_of_n_individuals = c(29957, 3000),
    n_individuals_R_sim = 3000,
    perc_phenos_with_effect = "0.1",
    effect_size = 0.01
) 
{
  
  setwd(sprintf("/scratch/giulia.desanctis/susie_results_method_2/%d_individuals/model_%d/%s_perc_phenos_with_effect/effect_%.2f", n_individuals_R_sim, model.no, perc_phenos_with_effect, effect_size))
  
  log_file <- "warnings.log"
  
  if (file.exists(log_file)) file.remove(log_file)  
  
  print(sprintf("model number: %d", model.no))
  
  withCallingHandlers({
    source("/home/giulia.desanctis/final_stage/scripts/functions_for_susie.R")
    
    simulation <- readRDS(sprintf(file_path_simulation, n_individuals_R_sim, model.no, model.no))
    #R_sim.no <- ifelse(model.no %% 100 == 0, 100, model.no %% 100)
    #R_simulated <- readRDS(sprintf(file_path_R_sim, R_sim.no, R_sim.no))
    R_simulated <- readRDS(sprintf(file_path_R_sim,  n_individuals_R_sim, model.no, perc_phenos_with_effect, effect_size, model.no))
    
    R_internal <- simulation$R_internal
    R_external <- simulation$R_external
    sum_stats <- simulation$sum_stats
    
    # Do internal susie
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
    
    # Do susie with R_sim for 30k individuals
    for(i in seq_len(length(R_simulated$matrix_names_list))) {
      start_sim <- Sys.time()
      tryCatch({
        do_susie(
          sim=simulation, 
          R_to_use = R_simulated$R_sim_zscore_list[[i]], 
          est_res_var=FALSE, 
          name_of_LD="simulated", 
          model_no = model.no, 
          num_ind = n_individuals_R_sim, 
          max_iter=max_iter_simulated_susie,
          num_phenos = (dim(R_simulated$sum_stats_list[[i]])[2]-8)/2
        )
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
    
    # Do external susie
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
    
    # Load susie objects to do fine-mapping 
    base <- sprintf("susie_model_%d_with_", model.no)
    
    # Generate suffixes for simulated LDs
    sim_suffixes <- expand.grid(ld = sim_LDs, n = n_individuals_R_sim)
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
    
    saveRDS(results_df, sprintf("susie_metrics_%d_%.2f_effect.rds", model_number, effect_size))
    
  }, warning = function(w) {
    cat("WARNING: ", conditionMessage(w), "\n", file = log_file, append = TRUE)
  })
  
}

##
susie_launcher(
    model.no = model_number
  , file_path_R_sim = "/scratch/giulia.desanctis/susie_results_method_2/%d_individuals/model_%d/%s_perc_phenos_with_effect/effect_%.2f/model_%d_with_R_sim_.rds"
  , file_path_simulation = "/scratch/giulia.desanctis/susie_results_method_2/%d_individuals/model_%d/model_%d_with_R_int_R_ext_and_effects.rds"
  , max_iter_external_susie = 20000
  , max_iter_internal_susie = 2000
  , max_iter_simulated_susie = 8000
  , fixed_LDs = c("internal_LD", "external_LD")
  , sim_LDs = c(1000, 2000, 3000, 4000, 5000)
  , n_individuals_R_sim = 1000
  , perc_phenos_with_effect = "10"
  , effect_size = 0.2
)
 













