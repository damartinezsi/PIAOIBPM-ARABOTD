# ==========================================================
# Lectura y redimensión de imágenes PNG (900x900 → 60x60)
# ==========================================================

# install.packages("magick")
library(magick)
library(grid)
library(abind)
library(rTensor)
library(pracma)     # para khatrirao
library(stats)      # para optim
library(cellWise)
library(robustbase)
library(rrcov)
library(RobStatTM)
library(ggplot2)
library(gridExtra)
library(openxlsx)

# ---------- Parámetros ----------
id_range <- 0:208


# Entrenamiento
dir_train_good <- "C:/Users/diego/Downloads/bottle-20251207T062318Z-3-001/bottle/train/good"
X_real <- load_tensor_images(dir_train_good)

# Prueba (cada categoría)
root_test <- "C:/Users/diego/Downloads/bottle-20251207T062318Z-3-001/bottle/test"
X_good          <- load_tensor_images(file.path(root_test, "good"))
X_broken_large  <- load_tensor_images(file.path(root_test, "broken_large"))
X_broken_small  <- load_tensor_images(file.path(root_test, "broken_small"))
X_contamination <- load_tensor_images(file.path(root_test, "contamination"))
X_contamination_strong <- load_tensor_images(file.path(root_test, "Contamination_strong"))


# Aplica a todo el tensor
X_shifted <- array(0, dim = dim(X_real))
for (i in 1:dim(X_real)[1]) {
  X_shifted[i,,,] <- shift_image(X_real[i,,,], up = 2, left = 2)
}

X_shifted2 <- array(0, dim = dim(X_real))
for (i in 1:dim(X_real)[1]) {
  X_shifted2[i,,,] <- shift_image(X_real[i,,,], up = 1, left = 1)
}

# Pequeño gráfico de las imágenes
plot(as.raster(X_real[1,,,]/255))
plot(as.raster(X_contamination_strong[2,,,]/255))
plot(as.raster(X_shifted[1,,,]/255))

# ==========================================================
# Estimación varianza retenida y gamma para MCD
# ==========================================================

library(parallel)

# Definir el número de clusters
num_cores <- detectCores() - 3 

# Crear el cluster
cl <- makeCluster(num_cores)

# Exportar los objetos necesarios al cluster
clusterExport(cl, varlist = c( "recon_trod_tensor",
                              "vecnorm", "makeDirections","unimcd",
                              "MacroPARAFAC2","Best_alpha_real","proj_qr",
                              "Real_images","X_real"), envir = environment())

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
#                Hallar valores de alpha y varianza explicada                  #
# ============================================================================ #

N_iter <-1000

# Establecer la misma semilla en cada nodo (para reproducibilidad)
clusterSetRNGStream(cl, 520)


# Correr simulaciones y ordenar resultados
results_ba <- parLapply(cl, 1:N_iter,
                        function(i) Best_alpha_real(i, X_IMAGES = X_real))
results_df_besta <- do.call(rbind, results_ba)


#Observar los mejores valores de alpha y la varianza retenida
apply(results_df_besta,2,mean)


# ============================================================================ #
#        Medir el tiempo de descomposición - Hallar los UCL                    #
# ============================================================================ #

clusterSetRNGStream(cl, 520)

results <- parLapply(cl, 1:N_iter,
                     function(i) Real_images(i, X_IMAGES = X_real))
Resultados_bootstrap_df <- do.call(rbind, results)

stopCluster(cl)

round(apply(Resultados_bootstrap_df[,17:25],2,mean),4)
# Nivel de significancia (Fase I, empírico 95%)
alpha_TROD      <- 0.051
alpha_TROD_MCD  <- 0.05
alpha_TROD_MCD_80  <- 0.05

alpha_MPCA      <- 0.05
alpha_MPCA_MCD  <- 0.05
alpha_MPCA_MCD_80  <- 0.051

alpha_MACRO     <- 0.05
alpha_MACRO_MCD <- 0.05
alpha_MACRO_MCD_80 <- 0.05


# ==========================================================
# === UCLs basados en percentiles empíricos (95%)
# ==========================================================

# Cartas T2 -------------------------------------------------
UCL_T_TROD      <- quantile(Resultados_bootstrap_df$max_T2_TROD,
                            probs = 1 - (alpha_TROD / 2))
UCL_T_TROD_MCD  <- quantile(Resultados_bootstrap_df$max_T2_TROD_MCD,
                            probs = 1 - (alpha_TROD_MCD / 2))
UCL_T_TROD_MCD_80  <- quantile(Resultados_bootstrap_df$max_T2_TROD_MCD_80,
                               probs = 1 - (alpha_TROD_MCD_80 / 2))

UCL_T_MPCA      <- quantile(Resultados_bootstrap_df$max_T2_MPCA,
                            probs = 1 - (alpha_MPCA / 2))
UCL_T_MPCA_MCD  <- quantile(Resultados_bootstrap_df$max_T2_MPCA_MCD,
                            probs = 1 - (alpha_MPCA_MCD / 2))
UCL_T_MPCA_MCD_80  <- quantile(Resultados_bootstrap_df$max_T2_MPCA_MCD_80,
                               probs = 1 - (alpha_MPCA_MCD_80 / 2))

UCL_T_MACRO     <- quantile(Resultados_bootstrap_df$max_T2_MACRO, 
                            probs = 1 - (alpha_MACRO / 2), na.rm = TRUE)
UCL_T_MACRO_MCD <- quantile(Resultados_bootstrap_df$max_T2_MACRO_MCD,
                            probs = 1 - (alpha_MACRO_MCD / 2), na.rm = TRUE)
UCL_T_MACRO_MCD_80 <- quantile(Resultados_bootstrap_df$max_T2_MACRO_MCD_80,
                               probs = 1 - (alpha_MACRO_MCD_80 / 2), na.rm = TRUE)

# Cartas Q --------------------------------------------------
UCL_Q_TROD      <- quantile(Resultados_bootstrap_df$max_Q_TROD,
                            probs = 1 - (alpha_TROD / 2))
UCL_Q_TROD_MCD  <- quantile(Resultados_bootstrap_df$max_Q_TROD,
                            probs = 1 - (alpha_TROD_MCD / 2))
UCL_Q_TROD_MCD_80  <- quantile(Resultados_bootstrap_df$max_Q_TROD,
                               probs = 1 - (alpha_TROD_MCD_80 / 2))

UCL_Q_MPCA      <- quantile(Resultados_bootstrap_df$max_Q_MPCA,
                            probs = 1 - (alpha_MPCA / 2))
UCL_Q_MPCA_MCD  <- quantile(Resultados_bootstrap_df$max_Q_MPCA,
                            probs = 1 - (alpha_MPCA_MCD / 2))
UCL_Q_MPCA_MCD_80  <- quantile(Resultados_bootstrap_df$max_Q_MPCA,
                               probs = 1 - (alpha_MPCA_MCD_80 / 2))

UCL_Q_MACRO     <- quantile(Resultados_bootstrap_df$max_Q_MACRO,
                            probs = 1 - (alpha_MACRO / 2), na.rm = TRUE)
UCL_Q_MACRO_MCD <- quantile(Resultados_bootstrap_df$max_Q_MACRO,
                            probs = 1 - (alpha_MACRO_MCD / 2), na.rm = TRUE)
UCL_Q_MACRO_MCD_80 <- quantile(Resultados_bootstrap_df$max_Q_MACRO,
                               probs = 1 - (alpha_MACRO_MCD_80 / 2), na.rm = TRUE)


# ==========================================================
# === Señales conjuntas T²–Q para cada método
# ==========================================================

# TROD ------------------------------------------------------
Resultados_bootstrap_df$SIGNAL_TROD      <- as.integer((Resultados_bootstrap_df$max_T2_TROD      > UCL_T_TROD)      | (Resultados_bootstrap_df$max_Q_TROD  > UCL_Q_TROD))
Resultados_bootstrap_df$SIGNAL_TROD_MCD  <- as.integer((Resultados_bootstrap_df$max_T2_TROD_MCD  > UCL_T_TROD_MCD)  | (Resultados_bootstrap_df$max_Q_TROD  > UCL_Q_TROD_MCD))
Resultados_bootstrap_df$SIGNAL_TROD_MCD_80  <- as.integer((Resultados_bootstrap_df$max_T2_TROD_MCD_80  > UCL_T_TROD_MCD_80)  | (Resultados_bootstrap_df$max_Q_TROD  > UCL_Q_TROD_MCD_80))

# MPCA ------------------------------------------------------
Resultados_bootstrap_df$SIGNAL_MPCA      <- as.integer((Resultados_bootstrap_df$max_T2_MPCA      > UCL_T_MPCA)      | (Resultados_bootstrap_df$max_Q_MPCA  > UCL_Q_MPCA))
Resultados_bootstrap_df$SIGNAL_MPCA_MCD  <- as.integer((Resultados_bootstrap_df$max_T2_MPCA_MCD  > UCL_T_MPCA_MCD)  | (Resultados_bootstrap_df$max_Q_MPCA  > UCL_Q_MPCA_MCD))
Resultados_bootstrap_df$SIGNAL_MPCA_MCD_80  <- as.integer((Resultados_bootstrap_df$max_T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80)  | (Resultados_bootstrap_df$max_Q_MPCA  > UCL_Q_MPCA_MCD_80))

# MACRO -----------------------------------------------------
Resultados_bootstrap_df$SIGNAL_MACRO     <- as.integer((Resultados_bootstrap_df$max_T2_MACRO     > UCL_T_MACRO)     | (Resultados_bootstrap_df$max_Q_MACRO > UCL_Q_MACRO))
Resultados_bootstrap_df$SIGNAL_MACRO_MCD <- as.integer((Resultados_bootstrap_df$max_T2_MACRO_MCD > UCL_T_MACRO_MCD) | (Resultados_bootstrap_df$max_Q_MACRO > UCL_Q_MACRO_MCD))
Resultados_bootstrap_df$SIGNAL_MACRO_MCD_80 <- as.integer((Resultados_bootstrap_df$max_T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) | (Resultados_bootstrap_df$max_Q_MACRO > UCL_Q_MACRO_MCD_80))


# ==========================================================
# === Tasas de señal promedio (esperadas cercanas a 0)
# ==========================================================

mean_vals <- c(
  TROD        = mean(Resultados_bootstrap_df$SIGNAL_TROD),
  TROD_MCD    = mean(Resultados_bootstrap_df$SIGNAL_TROD_MCD),
  TROD_MCD_80    = mean(Resultados_bootstrap_df$SIGNAL_TROD_MCD_80),
  MPCA        = mean(Resultados_bootstrap_df$SIGNAL_MPCA),
  MPCA_MCD    = mean(Resultados_bootstrap_df$SIGNAL_MPCA_MCD),
  MPCA_MCD_80    = mean(Resultados_bootstrap_df$SIGNAL_MPCA_MCD_80),
  MACRO       = mean(Resultados_bootstrap_df$SIGNAL_MACRO,     na.rm = TRUE),
  MACRO_MCD   = mean(Resultados_bootstrap_df$SIGNAL_MACRO_MCD, na.rm = TRUE),
  MACRO_MCD_80   = mean(Resultados_bootstrap_df$SIGNAL_MACRO_MCD_80, na.rm = TRUE)
)
print(round(mean_vals, 3))

############################################################################
#                   Caso 1A: 21 imágenes con cambio de centro 
############################################################################

n_outliers <- 21
set.seed(508)
X_defectuosas <- replace_last_n_with_shifted(X_real, X_shifted, n_outliers)
Resultados_Xprueba <- Real_test(X_defectuosas)
df <- Resultados_Xprueba
total_imgs <- dim(X_real)[1]

# Ubicación real de los outliers
is_fault <- rep(FALSE, total_imgs)
is_fault[(total_imgs - n_outliers + 1):total_imgs] <- TRUE


sig_TROD      <- (df$T2_TROD      > UCL_T_TROD)      | (df$Q_TROD  > UCL_Q_TROD)
sig_TROD_MCD  <- (df$T2_TROD_MCD  > UCL_T_TROD_MCD)  | (df$Q_TROD  > UCL_Q_TROD_MCD)
sig_TROD_MCD_80  <- (df$T2_TROD_MCD_80  > UCL_T_TROD_MCD_80)  | (df$Q_TROD  > UCL_Q_TROD_MCD_80)

sig_MPCA      <- (df$T2_MPCA      > UCL_T_MPCA)      | (df$Q_MPCA  > UCL_Q_MPCA)
sig_MPCA_MCD  <- (df$T2_MPCA_MCD  > UCL_T_MPCA_MCD)  | (df$Q_MPCA  > UCL_Q_MPCA_MCD)
sig_MPCA_MCD_80  <- (df$T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80)  | (df$Q_MPCA  > UCL_Q_MPCA_MCD_80)

sig_MACRO     <- (df$T2_MACRO     > UCL_T_MACRO)     | (df$Q_MACRO > UCL_Q_MACRO)
sig_MACRO_MCD <- (df$T2_MACRO_MCD > UCL_T_MACRO_MCD) | (df$Q_MACRO > UCL_Q_MACRO_MCD)
sig_MACRO_MCD_80 <- (df$T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) | (df$Q_MACRO > UCL_Q_MACRO_MCD_80)

# =========================
# TP / FP por carta (índices y conteos)
# =========================
TP_TROD      <- which(sig_TROD      &  is_fault); FP_TROD      <- which(sig_TROD      & !is_fault)
TP_TROD_MCD  <- which(sig_TROD_MCD  &  is_fault); FP_TROD_MCD  <- which(sig_TROD_MCD  & !is_fault)
TP_TROD_MCD_80  <- which(sig_TROD_MCD_80  &  is_fault); FP_TROD_MCD_80  <- which(sig_TROD_MCD_80  & !is_fault)

TP_MPCA      <- which(sig_MPCA      &  is_fault); FP_MPCA      <- which(sig_MPCA      & !is_fault)
TP_MPCA_MCD  <- which(sig_MPCA_MCD  &  is_fault); FP_MPCA_MCD  <- which(sig_MPCA_MCD  & !is_fault)
TP_MPCA_MCD_80  <- which(sig_MPCA_MCD_80  &  is_fault); FP_MPCA_MCD_80  <- which(sig_MPCA_MCD_80  & !is_fault)

TP_MACRO     <- which(sig_MACRO     &  is_fault); FP_MACRO     <- which(sig_MACRO     & !is_fault)
TP_MACRO_MCD <- which(sig_MACRO_MCD &  is_fault); FP_MACRO_MCD <- which(sig_MACRO_MCD & !is_fault)
TP_MACRO_MCD_80 <- which(sig_MACRO_MCD_80 &  is_fault); FP_MACRO_MCD_80 <- which(sig_MACRO_MCD_80 & !is_fault)

Resumen_TP_FP <- data.frame(
  carta = c("TROD","TROD_MCD","TROD_MCD_80","MPCA","MPCA_MCD","MPCA_MCD_80","MACRO","MACRO_MCD","MACRO_MCD_80"),
  TP = c(length(TP_TROD), length(TP_TROD_MCD), length(TP_TROD_MCD_80),
         length(TP_MPCA), length(TP_MPCA_MCD), length(TP_MPCA_MCD_80),
         length(TP_MACRO), length(TP_MACRO_MCD), length(TP_MACRO_MCD_80)),
  FP = c(length(FP_TROD), length(FP_TROD_MCD), length(FP_TROD_MCD_80),
         length(FP_MPCA), length(FP_MPCA_MCD), length(FP_MPCA_MCD_80),
         length(FP_MACRO), length(FP_MACRO_MCD), length(FP_MACRO_MCD_80))
)

print(Resumen_TP_FP)

df$tipo <- ifelse(is_fault==TRUE, "Shifted", "Original")

# ==========================================================
# CARTAS DE CONTROL - MACRO (T² y Q)
# ==========================================================

# Paleta
pal <- c("Original" = "#1f77b4", "Shifted" = "#d62728")

# === Carta T² (MACRO) ===
plot_T2_MACRO_C1A <- ggplot(df, aes(x = Image, y = T2_MACRO)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +          # línea azul clara fija
  geom_point(aes(fill = tipo),
             shape = 21, size = 3, color = "black", stroke = 0.35) +    # puntos con borde
  geom_hline(yintercept = UCL_T_MACRO, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 205,                     # al final del eje x
           y = UCL_T_MACRO * 0.95,                 # un poco por encima del UCL
           label = paste0("UCL = ", round(UCL_T_MACRO, 2)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_vline(xintercept = 188.5, color = "red", linetype = "dashed", linewidth = 1)+
  labs( x = "Image",
        y = expression(T[Macro]^2))+
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Carta Q (MACRO) ===
plot_Q_MACRO_C1A <- ggplot(df, aes(x = Image, y = Q_MACRO)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +
  geom_point(aes(fill = tipo),
             shape = 21, size = 3, color = "black", stroke = 0.35) +
  geom_hline(yintercept = UCL_Q_MACRO, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 205,                     
           y = UCL_Q_TROD * 0.96,                 
           label = paste0("UCL = ", round(UCL_Q_MACRO, 0)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_vline(xintercept = 188.5, color = "red", linetype = "dashed", linewidth = 1)+
  labs( x = "Image",
       y = expression(Q[Macro])) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Guardar los plots ===

ggsave(
  filename = "MACRO_T2_C1LOW.png",
  plot = plot_T2_MACRO_C1A,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  filename = "MACRO_Q_C1LOW.png",
  plot = plot_Q_MACRO_C1A,
  width = 9,
  height = 5,
  dpi = 450
)


############################################################################
#                   Caso 1B: 52 imágenes con cambio de centro 
############################################################################

n_outliers <- 52
set.seed(508)
X_defectuosas <- replace_last_n_with_shifted(X_real, X_shifted, n_outliers)
Resultados_Xprueba <- Real_test(X_defectuosas)
df <- Resultados_Xprueba
total_imgs <- dim(X_real)[1]

# Ubicación real de los outliers
is_fault <- rep(FALSE, total_imgs)
is_fault[(total_imgs - n_outliers + 1):total_imgs] <- TRUE


sig_TROD      <- (df$T2_TROD      > UCL_T_TROD)      | (df$Q_TROD  > UCL_Q_TROD)
sig_TROD_MCD  <- (df$T2_TROD_MCD  > UCL_T_TROD_MCD)  | (df$Q_TROD  > UCL_Q_TROD_MCD)
sig_TROD_MCD_80  <- (df$T2_TROD_MCD_80  > UCL_T_TROD_MCD_80)  | (df$Q_TROD  > UCL_Q_TROD_MCD_80)

sig_MPCA      <- (df$T2_MPCA      > UCL_T_MPCA)      | (df$Q_MPCA  > UCL_Q_MPCA)
sig_MPCA_MCD  <- (df$T2_MPCA_MCD  > UCL_T_MPCA_MCD)  | (df$Q_MPCA  > UCL_Q_MPCA_MCD)
sig_MPCA_MCD_80  <- (df$T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80)  | (df$Q_MPCA  > UCL_Q_MPCA_MCD_80)

sig_MACRO     <- (df$T2_MACRO     > UCL_T_MACRO)     | (df$Q_MACRO > UCL_Q_MACRO)
sig_MACRO_MCD <- (df$T2_MACRO_MCD > UCL_T_MACRO_MCD) | (df$Q_MACRO > UCL_Q_MACRO_MCD)
sig_MACRO_MCD_80 <- (df$T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) | (df$Q_MACRO > UCL_Q_MACRO_MCD_80)

# =========================
# TP / FP por carta (índices y conteos)
# =========================
TP_TROD      <- which(sig_TROD      &  is_fault); FP_TROD      <- which(sig_TROD      & !is_fault)
TP_TROD_MCD  <- which(sig_TROD_MCD  &  is_fault); FP_TROD_MCD  <- which(sig_TROD_MCD  & !is_fault)
TP_TROD_MCD_80  <- which(sig_TROD_MCD_80  &  is_fault); FP_TROD_MCD_80  <- which(sig_TROD_MCD_80  & !is_fault)

TP_MPCA      <- which(sig_MPCA      &  is_fault); FP_MPCA      <- which(sig_MPCA      & !is_fault)
TP_MPCA_MCD  <- which(sig_MPCA_MCD  &  is_fault); FP_MPCA_MCD  <- which(sig_MPCA_MCD  & !is_fault)
TP_MPCA_MCD_80  <- which(sig_MPCA_MCD_80  &  is_fault); FP_MPCA_MCD_80  <- which(sig_MPCA_MCD_80  & !is_fault)

TP_MACRO     <- which(sig_MACRO     &  is_fault); FP_MACRO     <- which(sig_MACRO     & !is_fault)
TP_MACRO_MCD <- which(sig_MACRO_MCD &  is_fault); FP_MACRO_MCD <- which(sig_MACRO_MCD & !is_fault)
TP_MACRO_MCD_80 <- which(sig_MACRO_MCD_80 &  is_fault); FP_MACRO_MCD_80 <- which(sig_MACRO_MCD_80 & !is_fault)

Resumen_TP_FP <- data.frame(
  carta = c("TROD","TROD_MCD","TROD_MCD_80","MPCA","MPCA_MCD","MPCA_MCD_80","MACRO","MACRO_MCD","MACRO_MCD_80"),
  TP = c(length(TP_TROD), length(TP_TROD_MCD), length(TP_TROD_MCD_80),
         length(TP_MPCA), length(TP_MPCA_MCD), length(TP_MPCA_MCD_80),
         length(TP_MACRO), length(TP_MACRO_MCD), length(TP_MACRO_MCD_80)),
  FP = c(length(FP_TROD), length(FP_TROD_MCD), length(FP_TROD_MCD_80),
         length(FP_MPCA), length(FP_MPCA_MCD), length(FP_MPCA_MCD_80),
         length(FP_MACRO), length(FP_MACRO_MCD), length(FP_MACRO_MCD_80))
)

print(Resumen_TP_FP)

df$tipo <- ifelse(is_fault==TRUE, "Shifted", "Original")

# ==========================================================
# CARTAS DE CONTROL - MACRO (T² y Q)
# ==========================================================

# Paleta
pal <- c("Original" = "#1f77b4", "Shifted" = "#d62728")

# === Carta T² (MACRO) ===
plot_T2_MACRO_C1 <- ggplot(df, aes(x = Image, y = T2_MACRO)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +          # línea azul clara fija
  geom_point(aes(fill = tipo),
             shape = 21, size = 3, color = "black", stroke = 0.35) +    # puntos con borde
  geom_hline(yintercept = UCL_T_MACRO, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 205,                     # al final del eje x
           y = UCL_T_MACRO * 0.95,                 # un poco por encima del UCL
           label = paste0("UCL = ", round(UCL_T_MACRO, 2)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_vline(xintercept = 157.5, color = "red", linetype = "dashed", linewidth = 1)+
  labs( x = "Image",
        y = expression(T[Macro]^2))+
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Carta Q (MACRO) ===
plot_Q_MACRO_C1 <- ggplot(df, aes(x = Image, y = Q_MACRO)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +
  geom_point(aes(fill = tipo),
             shape = 21, size = 3, color = "black", stroke = 0.35) +
  geom_hline(yintercept = UCL_Q_MACRO, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 205,                     
           y = UCL_Q_TROD * 0.96,                 
           label = paste0("UCL = ", round(UCL_Q_MACRO, 0)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_vline(xintercept = 157.5, color = "red", linetype = "dashed", linewidth = 1)+
  labs( x = "Image",
        y = expression(Q[Macro])) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Guardar los plots ===

ggsave(
  filename = "MACRO_T2_C1MED.png",
  plot = plot_T2_MACRO_C1,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  filename = "MACRO_Q_C1MED.png",
  plot = plot_Q_MACRO_C1,
  width = 9,
  height = 5,
  dpi = 450
)

# ==========================================================
# CARTAS DE CONTROL - TROD (T² y Q)
# ==========================================================

# === Carta T² (TROD) ===
plot_T2_TROD_C1 <- ggplot(df, aes(x = Image, y = T2_TROD)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +          # línea azul clara fija
  geom_point(aes(fill = tipo),
             shape = 22, size = 3, color = "black", stroke = 0.35) +    # puntos con borde
  geom_hline(yintercept = UCL_T_TROD, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 205,                     
           y = UCL_T_TROD * 0.95,                 
           label = paste0("UCL = ", round(UCL_T_TROD, 2)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_vline(xintercept = 157.5, color = "red", linetype = "dashed", linewidth = 1)+
  labs( x = "Image",
        y = expression(T[TROD]^2))+
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Carta Q (TROD) ===
plot_Q_TROD_C1 <- ggplot(df, aes(x = Image, y = Q_TROD)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +
  geom_point(aes(fill = tipo),
             shape = 22, size = 3, color = "black", stroke = 0.35) +
  geom_hline(yintercept = UCL_Q_TROD, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 215,                     
           y = UCL_Q_TROD * 0.97,                
           label = paste0("UCL = ", round(UCL_Q_TROD, 0)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_vline(xintercept = 157.5, color = "red", linetype = "dashed", linewidth = 1)+
  labs( x = "Image",
        y = expression(Q[TROD])) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )


# === Guardar los plots ===

ggsave(
  filename = "TROD_T2_C1MED.png",
  plot = plot_T2_TROD_C1,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  filename = "TROD_Q_C1MED.png",
  plot = plot_Q_TROD_C1,
  width = 9,
  height = 5,
  dpi = 450
)

############################################################################
#                   Caso 1C: 84 imágenes con cambio de centro 
############################################################################

n_outliers <- 84
set.seed(508)



idx_replace <-sample(1:209,n_outliers,FALSE)
X_defectuosas <- replace_at_indices_with_broken (X_real, X_shifted, 
                                                 idx_replace)
Resultados_Xprueba <- Real_test(X_defectuosas)
df <- Resultados_Xprueba
total_imgs <- dim(X_real)[1]

# Ubicación real de los outliers
is_fault <- rep(FALSE, total_imgs)
is_fault[idx_replace] <- TRUE


sig_TROD      <- (df$T2_TROD      > UCL_T_TROD)      | (df$Q_TROD  > UCL_Q_TROD)
sig_TROD_MCD  <- (df$T2_TROD_MCD  > UCL_T_TROD_MCD)  | (df$Q_TROD  > UCL_Q_TROD_MCD)
sig_TROD_MCD_80  <- (df$T2_TROD_MCD_80  > UCL_T_TROD_MCD_80)  | (df$Q_TROD  > UCL_Q_TROD_MCD_80)

sig_MPCA      <- (df$T2_MPCA      > UCL_T_MPCA)      | (df$Q_MPCA  > UCL_Q_MPCA)
sig_MPCA_MCD  <- (df$T2_MPCA_MCD  > UCL_T_MPCA_MCD)  | (df$Q_MPCA  > UCL_Q_MPCA_MCD)
sig_MPCA_MCD_80  <- (df$T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80)  | (df$Q_MPCA  > UCL_Q_MPCA_MCD_80)

sig_MACRO     <- (df$T2_MACRO     > UCL_T_MACRO)     | (df$Q_MACRO > UCL_Q_MACRO)
sig_MACRO_MCD <- (df$T2_MACRO_MCD > UCL_T_MACRO_MCD) | (df$Q_MACRO > UCL_Q_MACRO_MCD)
sig_MACRO_MCD_80 <- (df$T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) | (df$Q_MACRO > UCL_Q_MACRO_MCD_80)

# =========================
# TP / FP por carta (índices y conteos)
# =========================
TP_TROD      <- which(sig_TROD      &  is_fault); FP_TROD      <- which(sig_TROD      & !is_fault)
TP_TROD_MCD  <- which(sig_TROD_MCD  &  is_fault); FP_TROD_MCD  <- which(sig_TROD_MCD  & !is_fault)
TP_TROD_MCD_80  <- which(sig_TROD_MCD_80  &  is_fault); FP_TROD_MCD_80  <- which(sig_TROD_MCD_80  & !is_fault)

TP_MPCA      <- which(sig_MPCA      &  is_fault); FP_MPCA      <- which(sig_MPCA      & !is_fault)
TP_MPCA_MCD  <- which(sig_MPCA_MCD  &  is_fault); FP_MPCA_MCD  <- which(sig_MPCA_MCD  & !is_fault)
TP_MPCA_MCD_80  <- which(sig_MPCA_MCD_80  &  is_fault); FP_MPCA_MCD_80  <- which(sig_MPCA_MCD_80  & !is_fault)

TP_MACRO     <- which(sig_MACRO     &  is_fault); FP_MACRO     <- which(sig_MACRO     & !is_fault)
TP_MACRO_MCD <- which(sig_MACRO_MCD &  is_fault); FP_MACRO_MCD <- which(sig_MACRO_MCD & !is_fault)
TP_MACRO_MCD_80 <- which(sig_MACRO_MCD_80 &  is_fault); FP_MACRO_MCD_80 <- which(sig_MACRO_MCD_80 & !is_fault)

Resumen_TP_FP <- data.frame(
  carta = c("TROD","TROD_MCD","TROD_MCD_80","MPCA","MPCA_MCD","MPCA_MCD_80","MACRO","MACRO_MCD","MACRO_MCD_80"),
  TP = c(length(TP_TROD), length(TP_TROD_MCD), length(TP_TROD_MCD_80),
         length(TP_MPCA), length(TP_MPCA_MCD), length(TP_MPCA_MCD_80),
         length(TP_MACRO), length(TP_MACRO_MCD), length(TP_MACRO_MCD_80)),
  FP = c(length(FP_TROD), length(FP_TROD_MCD), length(FP_TROD_MCD_80),
         length(FP_MPCA), length(FP_MPCA_MCD), length(FP_MPCA_MCD_80),
         length(FP_MACRO), length(FP_MACRO_MCD), length(FP_MACRO_MCD_80))
)

print(Resumen_TP_FP)

df$tipo <- ifelse(is_fault==TRUE, "Shifted", "Original")

# ==========================================================
# CARTAS DE CONTROL - MACRO (T² y Q)
# ==========================================================

# Paleta
pal <- c("Original" = "#1f77b4", "Shifted" = "#d62728")

# === Carta T² (MACRO) ===
plot_T2_MACRO_C1 <- ggplot(df, aes(x = Image, y = T2_MACRO)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +          # línea azul clara fija
  geom_point(aes(fill = tipo),
             shape = 21, size = 3, color = "black", stroke = 0.35) +    # puntos con borde
  geom_hline(yintercept = UCL_T_MACRO, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 205,                     # al final del eje x
           y = UCL_T_MACRO * 0.95,                 # un poco por encima del UCL
           label = paste0("UCL = ", round(UCL_T_MACRO, 2)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  labs( x = "Image",
        y = expression(T[Macro]^2))+
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Carta Q (MACRO) ===
plot_Q_MACRO_C1 <- ggplot(df, aes(x = Image, y = Q_MACRO)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +
  geom_point(aes(fill = tipo),
             shape = 21, size = 3, color = "black", stroke = 0.35) +
  geom_hline(yintercept = UCL_Q_MACRO, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 205,                     
           y = UCL_Q_TROD * 0.96,                 
           label = paste0("UCL = ", round(UCL_Q_MACRO, 0)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  labs( x = "Image",
        y = expression(Q[Macro])) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Guardar los plots ===

ggsave(
  filename = "MACRO_T2_C1.png",
  plot = plot_T2_MACRO_C1,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  filename = "MACRO_Q_C1.png",
  plot = plot_Q_MACRO_C1,
  width = 9,
  height = 5,
  dpi = 450
)

# ==========================================================
# CARTAS DE CONTROL - MPCA (T² y Q)
# ==========================================================

# === Carta T² (MPCA) ===
plot_T2_MPCA_C1 <- ggplot(df, aes(x = Image, y = T2_MPCA)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +          # línea azul clara fija
  geom_point(aes(fill = tipo),
             shape = 22, size = 3, color = "black", stroke = 0.35) +    # puntos con borde
  geom_hline(yintercept = UCL_T_MPCA, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 205,                     
           y = UCL_T_MPCA * 0.95,                 
           label = paste0("UCL = ", round(UCL_T_MPCA, 2)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  labs( x = "Image",
        y = expression(T[MPCA]^2))+
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Carta Q (MPCA) ===
plot_Q_MPCA_C1 <- ggplot(df, aes(x = Image, y = Q_MPCA)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +
  geom_point(aes(fill = tipo),
             shape = 22, size = 3, color = "black", stroke = 0.35) +
  geom_hline(yintercept = UCL_Q_MPCA, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 215,                     
           y = UCL_Q_TROD * 0.97,                
           label = paste0("UCL = ", round(UCL_Q_MPCA, 0)),
           hjust = 1, vjust = 0, size = 5) +
  scale_fill_manual(values = pal, name = "Image version") +
  labs( x = "Image",
        y = expression(Q[MPCA])) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )


# === Guardar los plots ===

ggsave(
  filename = "MPCA_T2_C1.png",
  plot = plot_T2_MPCA_C1,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  filename = "MPCA_Q_C1.png",
  plot = plot_Q_MPCA_C1,
  width = 9,
  height = 5,
  dpi = 450
)
























########################################################################
#        Caso 2: imágenes con botellas contaminadas 
########################################################################


idx_replace <-seq(8,209,14)
set.seed(508)
X_defectuosas <- replace_at_indices_with_broken (X_real, X_contamination_strong, 
                                                 idx_replace)

Resultados_Xprueba <- Real_test(X_defectuosas)
df <- Resultados_Xprueba
total_imgs <- dim(X_real)[1]

# Ubicación real de los outliers
is_fault <- rep(FALSE, total_imgs)
is_fault[idx_replace] <- TRUE


sig_TROD      <- (df$T2_TROD      > UCL_T_TROD)      | (df$Q_TROD  > UCL_Q_TROD)
sig_TROD_MCD  <- (df$T2_TROD_MCD  > UCL_T_TROD_MCD)  | (df$Q_TROD  > UCL_Q_TROD_MCD)
sig_TROD_MCD_80  <- (df$T2_TROD_MCD_80  > UCL_T_TROD_MCD_80)  | (df$Q_TROD  > UCL_Q_TROD_MCD_80)

sig_MPCA      <- (df$T2_MPCA      > UCL_T_MPCA)      | (df$Q_MPCA  > UCL_Q_MPCA)
sig_MPCA_MCD  <- (df$T2_MPCA_MCD  > UCL_T_MPCA_MCD)  | (df$Q_MPCA  > UCL_Q_MPCA_MCD)
sig_MPCA_MCD_80  <- (df$T2_MPCA_MCD_80  > UCL_T_MPCA_MCD_80)  | (df$Q_MPCA  > UCL_Q_MPCA_MCD_80)

sig_MACRO     <- (df$T2_MACRO     > UCL_T_MACRO)     | (df$Q_MACRO > UCL_Q_MACRO)
sig_MACRO_MCD <- (df$T2_MACRO_MCD > UCL_T_MACRO_MCD) | (df$Q_MACRO > UCL_Q_MACRO_MCD)
sig_MACRO_MCD_80 <- (df$T2_MACRO_MCD_80 > UCL_T_MACRO_MCD_80) | (df$Q_MACRO > UCL_Q_MACRO_MCD_80)

# =========================
# TP / FP por carta (índices y conteos)
# =========================
TP_TROD      <- which(sig_TROD      &  is_fault); FP_TROD      <- which(sig_TROD      & !is_fault)
TP_TROD_MCD  <- which(sig_TROD_MCD  &  is_fault); FP_TROD_MCD  <- which(sig_TROD_MCD  & !is_fault)
TP_TROD_MCD_80  <- which(sig_TROD_MCD_80  &  is_fault); FP_TROD_MCD_80  <- which(sig_TROD_MCD_80  & !is_fault)

TP_MPCA      <- which(sig_MPCA      &  is_fault); FP_MPCA      <- which(sig_MPCA      & !is_fault)
TP_MPCA_MCD  <- which(sig_MPCA_MCD  &  is_fault); FP_MPCA_MCD  <- which(sig_MPCA_MCD  & !is_fault)
TP_MPCA_MCD_80  <- which(sig_MPCA_MCD_80  &  is_fault); FP_MPCA_MCD_80  <- which(sig_MPCA_MCD_80  & !is_fault)

TP_MACRO     <- which(sig_MACRO     &  is_fault); FP_MACRO     <- which(sig_MACRO     & !is_fault)
TP_MACRO_MCD <- which(sig_MACRO_MCD &  is_fault); FP_MACRO_MCD <- which(sig_MACRO_MCD & !is_fault)
TP_MACRO_MCD_80 <- which(sig_MACRO_MCD_80 &  is_fault); FP_MACRO_MCD_80 <- which(sig_MACRO_MCD_80 & !is_fault)

Resumen_TP_FP <- data.frame(
  carta = c("TROD","TROD_MCD","TROD_MCD_80","MPCA","MPCA_MCD","MPCA_MCD_80","MACRO","MACRO_MCD","MACRO_MCD_80"),
  TP = c(length(TP_TROD), length(TP_TROD_MCD), length(TP_TROD_MCD_80),
         length(TP_MPCA), length(TP_MPCA_MCD), length(TP_MPCA_MCD_80),
         length(TP_MACRO), length(TP_MACRO_MCD), length(TP_MACRO_MCD_80)),
  FP = c(length(FP_TROD), length(FP_TROD_MCD), length(FP_TROD_MCD_80),
         length(FP_MPCA), length(FP_MPCA_MCD), length(FP_MPCA_MCD_80),
         length(FP_MACRO), length(FP_MACRO_MCD), length(FP_MACRO_MCD_80))
)

print(Resumen_TP_FP)

df$tipo <- ifelse(is_fault==TRUE, "Contaminated", "Original")

# ==========================================================
# CARTAS DE CONTROL - TROD-MCD_F (T² y Q)
# ==========================================================

# Paleta
pal <- c("Original" = "#1f77b4", "Contaminated" = "#d62728")

# === Carta T² (TROD-MCD_F) ===
plot_T2_TRODMCDF_C2 <- ggplot(df, aes(x = Image, y = T2_TROD_MCD_80)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +          # línea azul clara fija
  geom_point(aes(fill = tipo),
             shape = 21, size = 3, color = "black", stroke = 0.35) +    # puntos con borde
  geom_hline(yintercept = UCL_T_TROD_MCD_80, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 30,                     # al final del eje x
           y = UCL_T_TROD_MCD_80 * 2,                 # un poco por encima del UCL
           label = paste0("UCL = ", round(UCL_T_TROD_MCD_80, 2)),
           hjust = 1, vjust = 0, size = 4) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_text(
    data = subset(df, Image %in% idx_replace),
    aes(label = Image),
    vjust = -0.9,            
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  )+
  labs( x = "Image",
        y = expression(T[TROD-MCD[F]]^2))+
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Carta Q (TROD-MCD_F) ===
plot_Q_TRODMCDF_C2 <- ggplot(df, aes(x = Image, y = Q_TROD)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +
  geom_point(aes(fill = tipo),
             shape = 21, size = 3, color = "black", stroke = 0.35) +
  geom_hline(yintercept = UCL_Q_TROD_MCD_80, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 100,                     
           y = UCL_Q_TROD_MCD_80 *1.01,                 
           label = paste0("UCL = ", round(UCL_Q_TROD_MCD_80, 0)),
           hjust = 1, vjust = 0, size = 4) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_text(
    data = subset(df, Image %in% idx_replace),
    aes(label = Image),
    vjust = -0.9,            
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  )+
  labs( x = "Image",
        y = expression(Q[TROD-MCD[F]])) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), 
    legend.title = element_text(size=15)
  )

# === Guardar los plots ===

ggsave(
  filename = "TRODMCDF_T2_C2.png",
  plot = plot_T2_TRODMCDF_C2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  filename = "TRODMCDF_Q_C2.png",
  plot = plot_Q_TRODMCDF_C2,
  width = 9,
  height = 5,
  dpi = 450
)

# ==========================================================
# CARTAS DE CONTROL - TROD (T² y Q)
# ==========================================================

# === Carta T² (TROD) ===
plot_T2_TROD_C2 <- ggplot(df, aes(x = Image, y = T2_TROD)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +          # línea azul clara fija
  geom_point(aes(fill = tipo),
             shape = 22, size = 3, color = "black", stroke = 0.35) +    # puntos con borde
  geom_hline(yintercept = UCL_T_TROD, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 210,                     
           y = UCL_T_TROD * 1.05,                 
           label = paste0("UCL = ", round(UCL_T_TROD, 2)),
           hjust = 1, vjust = 0, size = 4) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_text(
    data = subset(df, Image %in% idx_replace),
    aes(label = Image),
    vjust = -0.9,            
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  )+
  labs( x = "Image",
        y = expression(T[TROD]^2))+
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )

# === Carta Q (TROD) ===
plot_Q_TROD_C2 <- ggplot(df, aes(x = Image, y = Q_TROD)) +
  geom_line(color = "#89CFF0", linewidth = 0.8, alpha = 0.9) +
  geom_point(aes(fill = tipo),
             shape = 22, size = 3, color = "black", stroke = 0.35) +
  geom_hline(yintercept = UCL_Q_TROD, color = "black",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = 100,                     
           y = UCL_Q_TROD * 1.01,                
           label = paste0("UCL = ", round(UCL_Q_TROD, 0)),
           hjust = 1, vjust = 0, size = 4) +
  scale_fill_manual(values = pal, name = "Image version") +
  geom_text(
    data = subset(df, Image %in% idx_replace),
    aes(label = Image),
    vjust = -0.9,            
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  )+
  labs( x = "Image",
        y = expression(Q[TROD])) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title=element_text(size=20), legend.text = element_text(size=14), legend.title = element_text(size=15)
  )


# === Guardar los plots ===

ggsave(
  filename = "TROD_T2_C2.png",
  plot = plot_T2_TROD_C2,
  width = 9,
  height = 5,
  dpi = 450
)

ggsave(
  filename = "TROD_Q_C2.png",
  plot = plot_Q_TROD_C2,
  width = 9,
  height = 5,
  dpi = 450
)









