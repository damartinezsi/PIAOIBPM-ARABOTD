```r
#==============================================================================#
# ---------------------------------------------------------------------------- #
#              Main - Simulations                  D.A.M.S.                   #
# ---------------------------------------------------------------------------- #
#==============================================================================#

# Generate clusters for parallel processing

setwd("~/Tesis/Versión final - code")
library(parallel)

# Define the number of clusters
num_cores <- detectCores() - 1 

# Create the cluster
cl <- makeCluster(num_cores)

# Export the required objects to the cluster
clusterExport(cl, varlist = c("run_simulation","generate_ellipse_image",
                              "generate_samples", "recon_trod_tensor",
                              "vecnorm", "makeDirections","unimcd",
                              "MacroPARAFAC2","Best_alpha","proj_qr"
                              ), envir = environment())

clusterEvalQ(cl, {
  library(grid)
  library(abind)
  library(rTensor)
  library(pracma)     
  library(stats)      
  library(cellWise)
  library(robustbase)
  library(rrcov)
  library(RobStatTM)
  library(StableMCD)
})

# ============================================================================ #
#                Find alpha values and explained variance                     #
# ============================================================================ #

# Set the same seed on each node for reproducibility
clusterSetRNGStream(cl, 520)
N_iter <-1000

# Run simulations and combine results
results_ba <- parLapply(cl, 1:N_iter, Best_alpha)
results_df_besta <- do.call(rbind, results_ba)

# Examine the best alpha values and retained variance
apply(results_df_besta,2,mean)

# ============================================================================ #
#        Measure decomposition time - Obtain the UCLs                         #
# ============================================================================ #

# Set the same seed on each node for reproducibility
clusterSetRNGStream(cl, 520)

results <- parLapply(cl, 1:N_iter, run_simulation)
results_df <- do.call(rbind, results)

# Stop the cluster
stopCluster(cl)

round(apply(results_df,2,mean),4)

# Alpha values for n = 60 images ---------------------------------------------

alpha_TROD <-0.053        
alpha_TROD_MCD <-0.052  

alpha_MPCA <-0.054      
alpha_MPCA_MCD <-0.05    
alpha_MPCA_MCD_80 <-0.050   

alpha_MACRO <-0.052     
alpha_MACRO_MCD <-0.052    
alpha_MACRO_MCD_80 <-0.053


# UCLs based on empirical percentiles ----------------------------------------
# T² charts ------------------------------------------------------------------

UCL_T_TROD <- quantile(results_df$max_T2_TROD, probs = 1 - (alpha_TROD / 2))
UCL_T_TROD_MCD <- quantile(results_df$max_T2_TROD_MCD, probs = 1 - (alpha_TROD_MCD / 2))

UCL_T_MPCA <- quantile(results_df$max_T2_MPCA, probs = 1 - (alpha_MPCA / 2))
UCL_T_MPCA_MCD <- quantile(results_df$max_T2_MPCA_MCD, probs = 1 - (alpha_MPCA_MCD / 2))
UCL_T_MPCA_MCD_80 <- quantile(results_df$max_T2_MPCA_MCD_80, probs = 1 - (alpha_MPCA_MCD_80 / 2))

UCL_T_MACRO <- quantile(results_df$max_T2_MACRO, probs = 1 - (alpha_MACRO / 2))
UCL_T_MACRO_MCD <- quantile(results_df$max_T2_MACRO_MCD, probs = 1 - (alpha_MACRO_MCD / 2))
UCL_T_MACRO_MCD_80 <- quantile(results_df$max_T2_MACRO_MCD_80, probs = 1 - (alpha_MACRO_MCD_80 / 2))

# Q charts -------------------------------------------------------------------

UCL_Q_TROD <- quantile(results_df$max_Q_TROD, probs = 1 - (alpha_TROD / 2))
UCL_Q_TROD_MCD <- quantile(results_df$max_Q_TROD_MCD, probs = 1 - (alpha_TROD_MCD / 2))

UCL_Q_MPCA <- quantile(results_df$max_Q_MPCA, probs = 1 - (alpha_MPCA / 2))
UCL_Q_MPCA_MCD <- quantile(results_df$max_Q_MPCA_MCD, probs = 1 - (alpha_MPCA_MCD / 2))
UCL_Q_MPCA_MCD_80 <- quantile(results_df$max_Q_MPCA_MCD_80, probs = 1 - (alpha_MPCA_MCD_80 / 2))

UCL_Q_MACRO <- quantile(results_df$max_Q_MACRO, probs = 1 - (alpha_MACRO / 2))
UCL_Q_MACRO_MCD <- quantile(results_df$max_Q_MACRO_MCD, probs = 1 - (alpha_MACRO_MCD / 2))
UCL_Q_MACRO_MCD_80 <- quantile(results_df$max_Q_MACRO_MCD_80, probs = 1 - (alpha_MACRO_MCD_80 / 2))


# ============================================================================ #
#       JOINT T²–Q CHARTS FOR ALL METHODS AND VARIANTS                        #
# ============================================================================ #

# ------------------ TROD ------------------
results_df$SIGNAL_TROD <- as.integer(
  (results_df$max_T2_TROD > UCL_T_TROD) |
    (results_df$max_Q_TROD  > UCL_Q_TROD)
)
results_df$SIGNAL_TROD_MCD <- as.integer(
  (results_df$max_T2_TROD_MCD > UCL_T_TROD_MCD) |
    (results_df$max_Q_TROD      > UCL_Q_TROD_MCD)
)

# ------------------ MPCA ------------------
results_df$SIGNAL_MPCA <- as.integer(
  (results_df$max_T2_MPCA > UCL_T_MPCA) |
    (results_df$max_Q_MPCA  > UCL_Q_MPCA)
)
results_df$SIGNAL_MPCA_MCD <- as.integer(
  (results_df$max_T2_MPCA_MCD > UCL_T_MPCA_MCD) |
    (results_df$max_Q_MPCA      > UCL_Q_MPCA_MCD)
)
results_df$SIGNAL_MPCA_MCD_80 <- as.integer(
  (results_df$max_T2_MPCA_MCD_80 > UCL_T_MPCA_MCD_80) |
    (results_df$max_Q_MPCA      > UCL_Q_MPCA_MCD_80)
)

# ------------------ MACRO-PARAFAC ------------------
results_df$SIGNAL_MACRO <- as.integer(
  (results_df$max_T2_MACRO > UCL_T_MACRO) |
    (results_df$max_Q_MACRO  > UCL_Q_MACRO)
)
results_df$SIGNAL_MACRO_MCD <- as.integer(
  (results_df$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) |
    (results_df$max_Q_MACRO      > UCL_Q_MACRO_MCD)
)
results_df$SIGNAL_MACRO_MCD_80 <- as.integer(
  (results_df$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) |
    (results_df$max_Q_MACRO      > UCL_Q_MACRO_MCD_80)
)

# ============================================================================ #
#                           Average signal rates                              #
# ============================================================================ #

cat("Average signal rate (in-control):\n")

mean_vals <- c(
  TROD        = mean(results_df$SIGNAL_TROD),
  TROD_MCD    = mean(results_df$SIGNAL_TROD_MCD),
  MPCA        = mean(results_df$SIGNAL_MPCA),
  MPCA_MCD    = mean(results_df$SIGNAL_MPCA_MCD),
  MPCA_MCD_80 = mean(results_df$SIGNAL_MPCA_MCD_80),
  MACRO       = mean(results_df$SIGNAL_MACRO),
  MACRO_MCD   = mean(results_df$SIGNAL_MACRO_MCD),
  MACRO_MCD_80= mean(results_df$SIGNAL_MACRO_MCD_80)
)

print(round(mean_vals, 3))




# ============================================================================ #
#        Create clusters to run the different simulations                    #
# ============================================================================ #

num_cores <- detectCores() - 1  

# Create the cluster
cl <- makeCluster(num_cores)

# Export the required objects to the cluster
clusterExport(cl, varlist = c("run_simulation","generate_ellipse_image",
                              "generate_samples", "recon_trod_tensor",
                              "vecnorm", "makeDirections", "unimcd",
                              "MacroPARAFAC2"), envir = environment())

clusterEvalQ(cl, {
  library(grid)
  library(abind)
  library(rTensor)
  library(pracma)     
  library(stats)      
  library(cellWise)
  library(robustbase)
  library(rrcov)
  library(RobStatTM)
})

###################################################################
library(openxlsx)

#delta_values <- seq(1, 5, by = 1)  # For delta 1
delta_values <- seq(8, 40, by = 8)  # For delta 5A

N_iter <- 1000


# ================================== 95% ======================================= 

# Set the same seed on each node for reproducibility
clusterSetRNGStream(cl, 520)

# Empty data frame to store all results
results_signal <- data.frame()

for (delta_sim in delta_values) {
    resultados <- parLapply(
    cl, 1:N_iter,
    function(i, dval) {
      run_simulation(
        i,
        delta = 0,
        n_samples = 100,
        prop_in_control = 0.95,
        region_enable = FALSE,
        noise_frac = 0.0,
        region_mode = "edge_tl",
        d4delta = dval
        )
    },
    dval = delta_sim
  )
  
  resultados_df2 <- do.call(rbind, resultados)
  
  señales_iter <- data.frame(
    delta = delta_sim,
    
    # ------------------ TROD ------------------
    señal_TROD      = as.integer(
      (resultados_df2$max_T2_TROD      > UCL_T_TROD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD)
    ),
    señal_TROD_MCD  = as.integer(
      (resultados_df2$max_T2_TROD_MCD  > UCL_T_TROD_MCD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD_MCD)
    ),
    
    # ------------------ MPCA ------------------
    señal_MPCA      = as.integer(
      (resultados_df2$max_T2_MPCA      > UCL_T_MPCA) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA)
    ),
    señal_MPCA_MCD  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD  > UCL_T_MPCA_MCD) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD)
    ),
    señal_MPCA_MCD_80  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD_80)
    ),
    
    # ------------------ MACRO-PARAFAC ------------------
    señal_MACRO     = as.integer(
      (resultados_df2$max_T2_MACRO     > UCL_T_MACRO) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO)
    ),
    señal_MACRO_MCD = as.integer(
      (resultados_df2$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD)
    ),
    señal_MACRO_MCD_80 = as.integer(
      (resultados_df2$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD_80)
    )
  )
  
  results_signal <- rbind(results_signal, señales_iter)
  cat("End of delta:", delta_sim, "\n")
}


# Calculate the detection rate for each method and delta
tasa_senal <- aggregate(. ~ delta, data = results_signal, mean)

print(tasa_senal)

write.xlsx(tasa_senal, file = "tasa_senal.xlsx", sheetName = "Resultados95", rowNames = FALSE)


# ================================== 90% ======================================= 

# Set the same seed on each node for reproducibility
clusterSetRNGStream(cl, 520)

# Empty data frame to store all results
results_signal <- data.frame()

for (delta_sim in delta_values) {
  resultados <- parLapply(
    cl, 1:N_iter,
    function(i, dval) {
      run_simulation(
        i,
        delta = 0,
        n_samples = 100,
        prop_in_control = 0.9,
        region_enable = FALSE,
        noise_frac = 0.0,
        region_mode = "edge_tl",
        d4delta = dval
      )
    },
    dval = delta_sim
  )
  
  resultados_df2 <- do.call(rbind, resultados)
  
  señales_iter <- data.frame(
    delta = delta_sim,
    
    # ------------------ TROD ------------------
    señal_TROD      = as.integer(
      (resultados_df2$max_T2_TROD      > UCL_T_TROD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD)
    ),
    señal_TROD_MCD  = as.integer(
      (resultados_df2$max_T2_TROD_MCD  > UCL_T_TROD_MCD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD_MCD)
    ),
    
    # ------------------ MPCA ------------------
    señal_MPCA      = as.integer(
      (resultados_df2$max_T2_MPCA      > UCL_T_MPCA) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA)
    ),
    señal_MPCA_MCD  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD  > UCL_T_MPCA_MCD) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD)
    ),
    señal_MPCA_MCD_80  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD_80)
    ),
    
    # ------------------ MACRO-PARAFAC ------------------
    señal_MACRO     = as.integer(
      (resultados_df2$max_T2_MACRO     > UCL_T_MACRO) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO)
    ),
    señal_MACRO_MCD = as.integer(
      (resultados_df2$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD)
    ),
    señal_MACRO_MCD_80 = as.integer(
      (resultados_df2$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD_80)
    )
  )
  
  results_signal <- rbind(results_signal, señales_iter)
  cat("End of delta:", delta_sim, "\n")
}


# Calculate the detection rate for each method and delta
tasa_senal2 <- aggregate(. ~ delta, data = results_signal, mean)

print(tasa_senal2)
write.xlsx(tasa_senal2, file = "tasa_senal2.xlsx", sheetName = "Resultados90", rowNames = FALSE)

# ================================== 80% ======================================= 

# Set the same seed on each node for reproducibility
clusterSetRNGStream(cl, 520)

# Empty data frame to store all results
results_signal <- data.frame()

for (delta_sim in delta_values) {
  resultados <- parLapply(
    cl, 1:N_iter,
    function(i, dval) {
      run_simulation(
        i,
        delta = 0,
        n_samples = 100,
        prop_in_control = 0.8,
        region_enable = FALSE,
        noise_frac = 0.0,
        region_mode = "edge_tl",
        d4delta = dval
      )
    },
    dval = delta_sim
  )
  
  resultados_df2 <- do.call(rbind, resultados)
  
  señales_iter <- data.frame(
    delta = delta_sim,
    
    # ------------------ TROD ------------------
    señal_TROD      = as.integer(
      (resultados_df2$max_T2_TROD      > UCL_T_TROD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD)
    ),
    señal_TROD_MCD  = as.integer(
      (resultados_df2$max_T2_TROD_MCD  > UCL_T_TROD_MCD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD_MCD)
    ),
    
    # ------------------ MPCA ------------------
    señal_MPCA      = as.integer(
      (resultados_df2$max_T2_MPCA      > UCL_T_MPCA) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA)
    ),
    señal_MPCA_MCD  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD  > UCL_T_MPCA_MCD) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD)
    ),
    señal_MPCA_MCD_80  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD_80)
    ),
    
    # ------------------ MACRO-PARAFAC ------------------
    señal_MACRO     = as.integer(
      (resultados_df2$max_T2_MACRO     > UCL_T_MACRO) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO)
    ),
    señal_MACRO_MCD = as.integer(
      (resultados_df2$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD)
    ),
    señal_MACRO_MCD_80 = as.integer(
      (resultados_df2$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD_80)
    )
  )
  
  results_signal <- rbind(results_signal, señales_iter)
  cat("End of delta:", delta_sim, "\n")
}

# Calculate the detection rate for each method and delta
tasa_senal3 <- aggregate(. ~ delta, data = results_signal, mean)

print(tasa_senal3)
write.xlsx(tasa_senal3, file = "tasa_senal3.xlsx", sheetName = "Resultados80", rowNames = FALSE)


# ================================== 70% ======================================= 

# Set the same seed on each node for reproducibility
clusterSetRNGStream(cl, 520)

# Empty data frame to store all results
results_signal <- data.frame()

for (delta_sim in delta_values) {
  resultados <- parLapply(
    cl, 1:N_iter,
    function(i, dval) {
      run_simulation(
        i,
        delta = 0,
        n_samples = 100,
        prop_in_control = 0.7,
        region_enable = FALSE,
        noise_frac = 0.0,
        region_mode = "edge_tl",
        d4delta = dval
      )
    },
    dval = delta_sim
  )
  
  resultados_df2 <- do.call(rbind, resultados)
  
  señales_iter <- data.frame(
    delta = delta_sim,
    
    # ------------------ TROD ------------------
    señal_TROD      = as.integer(
      (resultados_df2$max_T2_TROD      > UCL_T_TROD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD)
    ),
    señal_TROD_MCD  = as.integer(
      (resultados_df2$max_T2_TROD_MCD  > UCL_T_TROD_MCD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD_MCD)
    ),
    
    # ------------------ MPCA ------------------
    señal_MPCA      = as.integer(
      (resultados_df2$max_T2_MPCA      > UCL_T_MPCA) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA)
    ),
    señal_MPCA_MCD  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD  > UCL_T_MPCA_MCD) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD)
    ),
    señal_MPCA_MCD_80  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD_80)
    ),
    
    # ------------------ MACRO-PARAFAC ------------------
    señal_MACRO     = as.integer(
      (resultados_df2$max_T2_MACRO     > UCL_T_MACRO) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO)
    ),
    señal_MACRO_MCD = as.integer(
      (resultados_df2$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD)
    ),
    señal_MACRO_MCD_80 = as.integer(
      (resultados_df2$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD_80)
    )
  )
  
  results_signal <- rbind(results_signal, señales_iter)
  cat("End of delta:", delta_sim, "\n")
}


# Calculate the detection rate for each method and delta
tasa_senal4 <- aggregate(. ~ delta, data = results_signal, mean)

print(tasa_senal4)
write.xlsx(tasa_senal4, file = "tasa_senal4.xlsx", sheetName = "Resultados70", rowNames = FALSE)


# ================================== 60% ======================================= 

# Set the same seed on each node for reproducibility
clusterSetRNGStream(cl, 520)

# Empty data frame to store all results
results_signal <- data.frame()

for (delta_sim in delta_values) {
  resultados <- parLapply(
    cl, 1:N_iter,
    function(i, dval) {
      run_simulation(
        i,
        delta = 0,
        n_samples = 100,
        prop_in_control = 0.6,
        region_enable = FALSE,
        noise_frac = 0.0,
        region_mode = "edge_tl",
        d4delta = dval
      )
    },
    dval = delta_sim
  )
  
  resultados_df2 <- do.call(rbind, resultados)
  
  señales_iter <- data.frame(
    delta = delta_sim,
    
    # ------------------ TROD ------------------
    señal_TROD      = as.integer(
      (resultados_df2$max_T2_TROD      > UCL_T_TROD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_T_TROD)
    ),
    señal_TROD_MCD  = as.integer(
      (resultados_df2$max_T2_TROD_MCD  > UCL_T_TROD_MCD) |
        (resultados_df2$max_Q_TROD       > UCL_Q_TROD_MCD)
    ),
    
    # ------------------ MPCA ------------------
    señal_MPCA      = as.integer(
      (resultados_df2$max_T2_MPCA      > UCL_T_MPCA) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA)
    ),
    señal_MPCA_MCD  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD  > UCL_T_MPCA_MCD) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD)
    ),
    señal_MPCA_MCD_80  = as.integer(
      (resultados_df2$max_T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80) |
        (resultados_df2$max_Q_MPCA       > UCL_Q_MPCA_MCD_80)
    ),
    
    # ------------------ MACRO-PARAFAC ------------------
    señal_MACRO     = as.integer(
      (resultados_df2$max_T2_MACRO     > UCL_T_MACRO) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO)
    ),
    señal_MACRO_MCD = as.integer(
      (resultados_df2$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD)
    ),
    señal_MACRO_MCD_80 = as.integer(
      (resultados_df2$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) |
        (resultados_df2$max_Q_MACRO      > UCL_Q_MACRO_MCD_80)
    )
  )
  
  results_signal <- rbind(results_signal, señales_iter)
  cat("End of delta:", delta_sim, "\n")
}


# Calculate the detection rate for each method and delta
tasa_senal5 <- aggregate(. ~ delta, data = results_signal, mean)

print(tasa_senal5)
write.xlsx(tasa_senal5, file = "tasa_senal5.xlsx", sheetName = "Resultados60", rowNames = FALSE)










stopCluster(cl)




