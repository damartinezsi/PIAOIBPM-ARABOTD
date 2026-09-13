# ==========================================================
# Real-data analysis: Phase I image monitoring
# ==========================================================
#
# This script:
#   1. Loads the MVTec AD bottle images.
#   2. Estimates retained variance and MCD tuning parameters.
#   3. Obtains empirical UCLs using bootstrap resampling.
#   4. Evaluates image-shift and contamination scenarios.
#   5. Generates control charts for the selected methods.
#
# The original image dataset is not included in this repository.
# See data/README.md for instructions on obtaining the data.
# ==========================================================


# ==========================================================
# 1. Packages
# ==========================================================

library(magick)
library(grid)
library(abind)
library(rTensor)
library(pracma)
library(stats)
library(cellWise)
library(robustbase)
library(rrcov)
library(RobStatTM)
library(ggplot2)
library(gridExtra)
library(openxlsx)
library(parallel)


# ==========================================================
# 2. Source project functions
# ==========================================================

source("R/functions.R")


# ==========================================================
# 3. Directories
# ==========================================================

data_dir <- file.path("data", "bottle")

dir_train_good <- file.path(
  data_dir,
  "train",
  "good"
)

root_test <- file.path(
  data_dir,
  "test"
)

figures_dir <- file.path(
  "results",
  "figures"
)

tables_dir <- file.path(
  "results",
  "tables"
)

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)


# ==========================================================
# 4. Load training and test images
# ==========================================================

X_real <- load_tensor_images(dir_train_good)

X_good <- load_tensor_images(
  file.path(root_test, "good")
)

X_broken_large <- load_tensor_images(
  file.path(root_test, "broken_large")
)

X_broken_small <- load_tensor_images(
  file.path(root_test, "broken_small")
)

X_contamination <- load_tensor_images(
  file.path(root_test, "contamination")
)

X_contamination_strong <- load_tensor_images(
  file.path(root_test, "Contamination_strong")
)


# Number of Phase I images
n_images <- dim(X_real)[1]


# ==========================================================
# 5. Generate shifted versions of the in-control images
# ==========================================================

X_shifted <- array(
  0,
  dim = dim(X_real)
)

for (i in seq_len(n_images)) {
  X_shifted[i, , , ] <- shift_image(
    X_real[i, , , ],
    up = 2,
    left = 2
  )
}


X_shifted2 <- array(
  0,
  dim = dim(X_real)
)

for (i in seq_len(n_images)) {
  X_shifted2[i, , , ] <- shift_image(
    X_real[i, , , ],
    up = 1,
    left = 1
  )
}


# ==========================================================
# 6. Bootstrap estimation
# ==========================================================

N_iter <- 1000

num_cores <- max(
  1,
  detectCores() - 3
)

cl <- makeCluster(num_cores)


# Objects required by the parallel workers
clusterExport(
  cl,
  varlist = c(
    "recon_trod_tensor",
    "vecnorm",
    "makeDirections",
    "unimcd",
    "MacroPARAFAC2",
    "Best_alpha_real",
    "proj_qr",
    "Real_images",
    "X_real"
  ),
  envir = environment()
)


clusterEvalQ(
  cl,
  {
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
  }
)


# ----------------------------------------------------------
# 6.1 Estimate alpha and retained variance
# ----------------------------------------------------------

clusterSetRNGStream(
  cl,
  iseed = 520
)

results_ba <- parLapply(
  cl,
  seq_len(N_iter),
  function(i) {
    Best_alpha_real(
      i,
      X_IMAGES = X_real
    )
  }
)

results_df_besta <- do.call(
  rbind,
  results_ba
)


mean_best_alpha <- apply(
  results_df_besta,
  2,
  mean
)


print(
  round(mean_best_alpha, 4)
)


# ----------------------------------------------------------
# 6.2 Bootstrap distributions for UCL estimation
# ----------------------------------------------------------

clusterSetRNGStream(
  cl,
  iseed = 520
)

results <- parLapply(
  cl,
  seq_len(N_iter),
  function(i) {
    Real_images(
      i,
      X_IMAGES = X_real
    )
  }
)

Resultados_bootstrap_df <- do.call(
  rbind,
  results
)

stopCluster(cl)


# ==========================================================
# 7. Significance levels
# ==========================================================

alpha_TROD <- 0.051
alpha_TROD_MCD <- 0.05
alpha_TROD_MCD_80 <- 0.05

alpha_MPCA <- 0.05
alpha_MPCA_MCD <- 0.05
alpha_MPCA_MCD_80 <- 0.051

alpha_MACRO <- 0.05
alpha_MACRO_MCD <- 0.05
alpha_MACRO_MCD_80 <- 0.05


# ==========================================================
# 8. Empirical UCLs
# ==========================================================

# ----------------------------------------------------------
# T2 charts
# ----------------------------------------------------------

UCL_T_TROD <- quantile(
  Resultados_bootstrap_df$max_T2_TROD,
  probs = 1 - alpha_TROD / 2
)

UCL_T_TROD_MCD <- quantile(
  Resultados_bootstrap_df$max_T2_TROD_MCD,
  probs = 1 - alpha_TROD_MCD / 2
)

UCL_T_TROD_MCD_80 <- quantile(
  Resultados_bootstrap_df$max_T2_TROD_MCD_80,
  probs = 1 - alpha_TROD_MCD_80 / 2
)


UCL_T_MPCA <- quantile(
  Resultados_bootstrap_df$max_T2_MPCA,
  probs = 1 - alpha_MPCA / 2
)

UCL_T_MPCA_MCD <- quantile(
  Resultados_bootstrap_df$max_T2_MPCA_MCD,
  probs = 1 - alpha_MPCA_MCD / 2
)

UCL_T_MPCA_MCD_80 <- quantile(
  Resultados_bootstrap_df$max_T2_MPCA_MCD_80,
  probs = 1 - alpha_MPCA_MCD_80 / 2
)


UCL_T_MACRO <- quantile(
  Resultados_bootstrap_df$max_T2_MACRO,
  probs = 1 - alpha_MACRO / 2,
  na.rm = TRUE
)

UCL_T_MACRO_MCD <- quantile(
  Resultados_bootstrap_df$max_T2_MACRO_MCD,
  probs = 1 - alpha_MACRO_MCD / 2,
  na.rm = TRUE
)

UCL_T_MACRO_MCD_80 <- quantile(
  Resultados_bootstrap_df$max_T2_MACRO_MCD_80,
  probs = 1 - alpha_MACRO_MCD_80 / 2,
  na.rm = TRUE
)


# ----------------------------------------------------------
# Q charts
# ----------------------------------------------------------

UCL_Q_TROD <- quantile(
  Resultados_bootstrap_df$max_Q_TROD,
  probs = 1 - alpha_TROD / 2
)

UCL_Q_TROD_MCD <- quantile(
  Resultados_bootstrap_df$max_Q_TROD,
  probs = 1 - alpha_TROD_MCD / 2
)

UCL_Q_TROD_MCD_80 <- quantile(
  Resultados_bootstrap_df$max_Q_TROD,
  probs = 1 - alpha_TROD_MCD_80 / 2
)


UCL_Q_MPCA <- quantile(
  Resultados_bootstrap_df$max_Q_MPCA,
  probs = 1 - alpha_MPCA / 2
)

UCL_Q_MPCA_MCD <- quantile(
  Resultados_bootstrap_df$max_Q_MPCA,
  probs = 1 - alpha_MPCA_MCD / 2
)

UCL_Q_MPCA_MCD_80 <- quantile(
  Resultados_bootstrap_df$max_Q_MPCA,
  probs = 1 - alpha_MPCA_MCD_80 / 2
)


UCL_Q_MACRO <- quantile(
  Resultados_bootstrap_df$max_Q_MACRO,
  probs = 1 - alpha_MACRO / 2,
  na.rm = TRUE
)

UCL_Q_MACRO_MCD <- quantile(
  Resultados_bootstrap_df$max_Q_MACRO,
  probs = 1 - alpha_MACRO_MCD / 2,
  na.rm = TRUE
)

UCL_Q_MACRO_MCD_80 <- quantile(
  Resultados_bootstrap_df$max_Q_MACRO,
  probs = 1 - alpha_MACRO_MCD_80 / 2,
  na.rm = TRUE
)


# ==========================================================
# 9. In-control bootstrap signal rates
# ==========================================================

Resultados_bootstrap_df$SIGNAL_TROD <-
  as.integer(
    Resultados_bootstrap_df$max_T2_TROD > UCL_T_TROD |
      Resultados_bootstrap_df$max_Q_TROD > UCL_Q_TROD
  )


Resultados_bootstrap_df$SIGNAL_TROD_MCD <-
  as.integer(
    Resultados_bootstrap_df$max_T2_TROD_MCD > UCL_T_TROD_MCD |
      Resultados_bootstrap_df$max_Q_TROD > UCL_Q_TROD_MCD
  )


Resultados_bootstrap_df$SIGNAL_TROD_MCD_80 <-
  as.integer(
    Resultados_bootstrap_df$max_T2_TROD_MCD_80 >
      UCL_T_TROD_MCD_80 |
      Resultados_bootstrap_df$max_Q_TROD >
      UCL_Q_TROD_MCD_80
  )


Resultados_bootstrap_df$SIGNAL_MPCA <-
  as.integer(
    Resultados_bootstrap_df$max_T2_MPCA > UCL_T_MPCA |
      Resultados_bootstrap_df$max_Q_MPCA > UCL_Q_MPCA
  )


Resultados_bootstrap_df$SIGNAL_MPCA_MCD <-
  as.integer(
    Resultados_bootstrap_df$max_T2_MPCA_MCD >
      UCL_T_MPCA_MCD |
      Resultados_bootstrap_df$max_Q_MPCA >
      UCL_Q_MPCA_MCD
  )


Resultados_bootstrap_df$SIGNAL_MPCA_MCD_80 <-
  as.integer(
    Resultados_bootstrap_df$max_T2_MPCA_MCD_80 >
      UCL_T_MPCA_MCD_80 |
      Resultados_bootstrap_df$max_Q_MPCA >
      UCL_Q_MPCA_MCD_80
  )


Resultados_bootstrap_df$SIGNAL_MACRO <-
  as.integer(
    Resultados_bootstrap_df$max_T2_MACRO > UCL_T_MACRO |
      Resultados_bootstrap_df$max_Q_MACRO > UCL_Q_MACRO
  )


Resultados_bootstrap_df$SIGNAL_MACRO_MCD <-
  as.integer(
    Resultados_bootstrap_df$max_T2_MACRO_MCD >
      UCL_T_MACRO_MCD |
      Resultados_bootstrap_df$max_Q_MACRO >
      UCL_Q_MACRO_MCD
  )


Resultados_bootstrap_df$SIGNAL_MACRO_MCD_80 <-
  as.integer(
    Resultados_bootstrap_df$max_T2_MACRO_MCD_80 >
      UCL_T_MACRO_MCD_80 |
      Resultados_bootstrap_df$max_Q_MACRO >
      UCL_Q_MACRO_MCD_80
  )


mean_signal_rates <- c(
  TROD = mean(
    Resultados_bootstrap_df$SIGNAL_TROD
  ),
  TROD_MCD = mean(
    Resultados_bootstrap_df$SIGNAL_TROD_MCD
  ),
  TROD_MCD_80 = mean(
    Resultados_bootstrap_df$SIGNAL_TROD_MCD_80
  ),
  MPCA = mean(
    Resultados_bootstrap_df$SIGNAL_MPCA
  ),
  MPCA_MCD = mean(
    Resultados_bootstrap_df$SIGNAL_MPCA_MCD
  ),
  MPCA_MCD_80 = mean(
    Resultados_bootstrap_df$SIGNAL_MPCA_MCD_80
  ),
  MACRO = mean(
    Resultados_bootstrap_df$SIGNAL_MACRO,
    na.rm = TRUE
  ),
  MACRO_MCD = mean(
    Resultados_bootstrap_df$SIGNAL_MACRO_MCD,
    na.rm = TRUE
  ),
  MACRO_MCD_80 = mean(
    Resultados_bootstrap_df$SIGNAL_MACRO_MCD_80,
    na.rm = TRUE
  )
)

print(
  round(mean_signal_rates, 3)
)


# ==========================================================
# 10. Helper function for scenario evaluation
# ==========================================================

evaluate_scenario <- function(
    X_test,
    fault_indices,
    UCLs
) {

  results <- Real_test(X_test)

  n <- dim(X_test)[1]

  is_fault <- rep(FALSE, n)
  is_fault[fault_indices] <- TRUE


  signals <- list(

    TROD =
      results$T2_TROD > UCLs$TROD_T |
      results$Q_TROD > UCLs$TROD_Q,

    TROD_MCD =
      results$T2_TROD_MCD > UCLs$TROD_MCD_T |
      results$Q_TROD > UCLs$TROD_MCD_Q,

    TROD_MCD_80 =
      results$T2_TROD_MCD_80 >
        UCLs$TROD_MCD_80_T |
      results$Q_TROD >
        UCLs$TROD_MCD_80_Q,


    MPCA =
      results$T2_MPCA > UCLs$MPCA_T |
      results$Q_MPCA > UCLs$MPCA_Q,

    MPCA_MCD =
      results$T2_MPCA_MCD >
        UCLs$MPCA_MCD_T |
      results$Q_MPCA >
        UCLs$MPCA_MCD_Q,

    MPCA_MCD_80 =
      results$T2_MPCA_MCD_80 >
        UCLs$MPCA_MCD_80_T |
      results$Q_MPCA >
        UCLs$MPCA_MCD_80_Q,


    MACRO =
      results$T2_MACRO > UCLs$MACRO_T |
      results$Q_MACRO > UCLs$MACRO_Q,

    MACRO_MCD =
      results$T2_MACRO_MCD >
        UCLs$MACRO_MCD_T |
      results$Q_MACRO >
        UCLs$MACRO_MCD_Q,

    MACRO_MCD_80 =
      results$T2_MACRO_MCD_80 >
        UCLs$MACRO_MCD_80_T |
      results$Q_MACRO >
        UCLs$MACRO_MCD_80_Q
  )


  summary <- do.call(
    rbind,
    lapply(
      names(signals),
      function(method) {

        signal <- signals[[method]]

        data.frame(
          method = method,
          TP = sum(signal & is_fault),
          FP = sum(signal & !is_fault),
          FN = sum(!signal & is_fault),
          TN = sum(!signal & !is_fault)
        )
      }
    )
  )


  list(
    results = results,
    is_fault = is_fault,
    signals = signals,
    summary = summary
  )
}


# ==========================================================
# 11. UCL list
# ==========================================================

UCLs <- list(

  TROD_T = UCL_T_TROD,
  TROD_Q = UCL_Q_TROD,

  TROD_MCD_T = UCL_T_TROD_MCD,
  TROD_MCD_Q = UCL_Q_TROD_MCD,

  TROD_MCD_80_T = UCL_T_TROD_MCD_80,
  TROD_MCD_80_Q = UCL_Q_TROD_MCD_80,


  MPCA_T = UCL_T_MPCA,
  MPCA_Q = UCL_Q_MPCA,

  MPCA_MCD_T = UCL_T_MPCA_MCD,
  MPCA_MCD_Q = UCL_Q_MPCA_MCD,

  MPCA_MCD_80_T = UCL_T_MPCA_MCD_80,
  MPCA_MCD_80_Q = UCL_Q_MPCA_MCD_80,


  MACRO_T = UCL_T_MACRO,
  MACRO_Q = UCL_Q_MACRO,

  MACRO_MCD_T = UCL_T_MACRO_MCD,
  MACRO_MCD_Q = UCL_Q_MACRO_MCD,

  MACRO_MCD_80_T = UCL_T_MACRO_MCD_80,
  MACRO_MCD_80_Q = UCL_Q_MACRO_MCD_80
)


# ==========================================================
# 12. Case 1A: 21 shifted images
# ==========================================================

n_outliers <- 21

X_case1A <- replace_last_n_with_shifted(
  X_real,
  X_shifted,
  n_outliers
)

case1A <- evaluate_scenario(
  X_case1A,
  fault_indices = (n_images - n_outliers + 1):n_images,
  UCLs = UCLs
)

print(case1A$summary)


# ==========================================================
# 13. Case 1B: 52 shifted images
# ==========================================================

n_outliers <- 52

X_case1B <- replace_last_n_with_shifted(
  X_real,
  X_shifted,
  n_outliers
)

case1B <- evaluate_scenario(
  X_case1B,
  fault_indices = (n_images - n_outliers + 1):n_images,
  UCLs = UCLs
)

print(case1B$summary)


# ==========================================================
# 14. Case 1C: 84 shifted images at random locations
# ==========================================================

n_outliers <- 84

set.seed(508)

idx_case1C <- sample(
  seq_len(n_images),
  size = n_outliers,
  replace = FALSE
)

X_case1C <- replace_at_indices_with_broken(
  X_real,
  X_shifted,
  idx_case1C
)

case1C <- evaluate_scenario(
  X_case1C,
  fault_indices = idx_case1C,
  UCLs = UCLs
)

print(case1C$summary)


# ==========================================================
# 15. Case 2: contaminated images
# ==========================================================

set.seed(508)

idx_case2 <- seq(
  8,
  n_images,
  by = 14
)

X_case2 <- replace_at_indices_with_broken(
  X_real,
  X_contamination_strong,
  idx_case2
)

case2 <- evaluate_scenario(
  X_case2,
  fault_indices = idx_case2,
  UCLs = UCLs
)

print(case2$summary)


# ==========================================================
# 16. Save summary tables
# ==========================================================

write.xlsx(
  list(
    Case_1A = case1A$summary,
    Case_1B = case1B$summary,
    Case_1C = case1C$summary,
    Case_2 = case2$summary
  ),
  file = file.path(
    tables_dir,
    "real_data_detection_results.xlsx"
  ),
  overwrite = TRUE
)


# ==========================================================
# 17. Plotting functions
# ==========================================================

plot_control_chart <- function(
    results,
    statistic,
    UCL,
    fault_indices,
    output_file,
    y_label,
    legend_title = "Observation type"
) {

  plot_data <- results

  plot_data$observation_type <- "In-control"

  plot_data$observation_type[
    plot_data$Image %in% fault_indices
  ] <- "Outlier"


  pal <- c(
    "In-control" = "#1f77b4",
    "Outlier" = "#d62728"
  )


  ggplot(
    plot_data,
    aes(
      x = Image,
      y = .data[[statistic]]
    )
  ) +

    geom_line(
      color = "#89CFF0",
      linewidth = 0.8,
      alpha = 0.9
    ) +

    geom_point(
      aes(fill = observation_type),
      shape = 21,
      size = 3,
      color = "black",
      stroke = 0.35
    ) +

    geom_hline(
      yintercept = UCL,
      color = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +

    scale_fill_manual(
      values = pal,
      name = legend_title
    ) +

    labs(
      x = "Image",
      y = y_label
    ) +

    theme_bw(
      base_size = 13
    ) +

    theme(
      legend.position = "right",
      panel.grid.major =
        element_line(
          color = "grey85",
          linewidth = 0.35
        ),
      panel.grid.minor =
        element_line(
          color = "grey93",
          linewidth = 0.25
        ),
      panel.border =
        element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.7
        ),
      axis.title =
        element_text(size = 20),
      legend.text =
        element_text(size = 14),
      legend.title =
        element_text(size = 15)
    )
}


# ==========================================================
# 18. Figures for Case 1A
# ==========================================================

plot_case1A_T2 <- plot_control_chart(
  case1A$results,
  statistic = "T2_MACRO",
  UCL = UCL_T_MACRO,
  fault_indices =
    (n_images - 21 + 1):n_images,
  output_file =
    file.path(
      figures_dir,
      "MACRO_T2_C1LOW.png"
    ),
  y_label =
    expression(T[Macro]^2)
)

plot_case1A_Q <- plot_control_chart(
  case1A$results,
  statistic = "Q_MACRO",
  UCL = UCL_Q_MACRO,
  fault_indices =
    (n_images - 21 + 1):n_images,
  output_file =
    file.path(
      figures_dir,
      "MACRO_Q_C1LOW.png"
    ),
  y_label =
    expression(Q[Macro])
)


ggsave(
  file.path(
    figures_dir,
    "MACRO_T2_C1LOW.png"
  ),
  plot_case1A_T2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "MACRO_Q_C1LOW.png"
  ),
  plot_case1A_Q,
  width = 9,
  height = 5,
  dpi = 450
)


# ==========================================================
# 19. Figures for Case 1B
# ==========================================================

case1B_faults <-
  (n_images - 52 + 1):n_images


plot_case1B_macro_T2 <- plot_control_chart(
  case1B$results,
  "T2_MACRO",
  UCL_T_MACRO,
  case1B_faults,
  file.path(
    figures_dir,
    "MACRO_T2_C1MED.png"
  ),
  expression(T[Macro]^2)
)


plot_case1B_macro_Q <- plot_control_chart(
  case1B$results,
  "Q_MACRO",
  UCL_Q_MACRO,
  case1B_faults,
  file.path(
    figures_dir,
    "MACRO_Q_C1MED.png"
  ),
  expression(Q[Macro])
)


plot_case1B_trod_T2 <- plot_control_chart(
  case1B$results,
  "T2_TROD",
  UCL_T_TROD,
  case1B_faults,
  file.path(
    figures_dir,
    "TROD_T2_C1MED.png"
  ),
  expression(T[TROD]^2)
)


plot_case1B_trod_Q <- plot_control_chart(
  case1B$results,
  "Q_TROD",
  UCL_Q_TROD,
  case1B_faults,
  file.path(
    figures_dir,
    "TROD_Q_C1MED.png"
  ),
  expression(Q[TROD])
)


ggsave(
  file.path(
    figures_dir,
    "MACRO_T2_C1MED.png"
  ),
  plot_case1B_macro_T2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "MACRO_Q_C1MED.png"
  ),
  plot_case1B_macro_Q,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "TROD_T2_C1MED.png"
  ),
  plot_case1B_trod_T2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "TROD_Q_C1MED.png"
  ),
  plot_case1B_trod_Q,
  width = 9,
  height = 5,
  dpi = 450
)


# ==========================================================
# 20. Figures for Case 1C
# ==========================================================

plot_case1C_macro_T2 <- plot_control_chart(
  case1C$results,
  "T2_MACRO",
  UCL_T_MACRO,
  idx_case1C,
  file.path(
    figures_dir,
    "MACRO_T2_C1.png"
  ),
  expression(T[Macro]^2)
)


plot_case1C_macro_Q <- plot_control_chart(
  case1C$results,
  "Q_MACRO",
  UCL_Q_MACRO,
  idx_case1C,
  file.path(
    figures_dir,
    "MACRO_Q_C1.png"
  ),
  expression(Q[Macro])
)


plot_case1C_mpca_T2 <- plot_control_chart(
  case1C$results,
  "T2_MPCA",
  UCL_T_MPCA,
  idx_case1C,
  file.path(
    figures_dir,
    "MPCA_T2_C1.png"
  ),
  expression(T[MPCA]^2)
)


plot_case1C_mpca_Q <- plot_control_chart(
  case1C$results,
  "Q_MPCA",
  UCL_Q_MPCA,
  idx_case1C,
  file.path(
    figures_dir,
    "MPCA_Q_C1.png"
  ),
  expression(Q[MPCA])
)


ggsave(
  file.path(
    figures_dir,
    "MACRO_T2_C1.png"
  ),
  plot_case1C_macro_T2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "MACRO_Q_C1.png"
  ),
  plot_case1C_macro_Q,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "MPCA_T2_C1.png"
  ),
  plot_case1C_mpca_T2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "MPCA_Q_C1.png"
  ),
  plot_case1C_mpca_Q,
  width = 9,
  height = 5,
  dpi = 450
)


# ==========================================================
# 21. Figures for Case 2
# ==========================================================

plot_case2_trod_mcd_T2 <- plot_control_chart(
  case2$results,
  "T2_TROD_MCD_80",
  UCL_T_TROD_MCD_80,
  idx_case2,
  file.path(
    figures_dir,
    "TRODMCDF_T2_C2.png"
  ),
  expression(T[TROD-MCD[F]]^2)
)


plot_case2_trod_mcd_Q <- plot_control_chart(
  case2$results,
  "Q_TROD",
  UCL_Q_TROD_MCD_80,
  idx_case2,
  file.path(
    figures_dir,
    "TRODMCDF_Q_C2.png"
  ),
  expression(Q[TROD-MCD[F]])
)


plot_case2_trod_T2 <- plot_control_chart(
  case2$results,
  "T2_TROD",
  UCL_T_TROD,
  idx_case2,
  file.path(
    figures_dir,
    "TROD_T2_C2.png"
  ),
  expression(T[TROD]^2)
)


plot_case2_trod_Q <- plot_control_chart(
  case2$results,
  "Q_TROD",
  UCL_Q_TROD,
  idx_case2,
  file.path(
    figures_dir,
    "TROD_Q_C2.png"
  ),
  expression(Q[TROD])
)


ggsave(
  file.path(
    figures_dir,
    "TRODMCDF_T2_C2.png"
  ),
  plot_case2_trod_mcd_T2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "TRODMCDF_Q_C2.png"
  ),
  plot_case2_trod_mcd_Q,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "TROD_T2_C2.png"
  ),
  plot_case2_trod_T2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  file.path(
    figures_dir,
    "TROD_Q_C2.png"
  ),
  plot_case2_trod_Q,
  width = 9,
  height = 5,
  dpi = 450
)


# ==========================================================
# End of script
# ==========================================================





