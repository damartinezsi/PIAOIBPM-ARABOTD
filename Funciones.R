# ============================================================================ #
#                   MacroPARAFAC para tensores de orden 4                      #
# ============================================================================ #

MacroPARAFAC2 <- function(X, ncomp = 6, MacroPARAFACpars = NULL, 
                          const = NULL){
  
  # The random seed is retained when leaving the function
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    seed.keep <- get(".Random.seed", envir = .GlobalEnv,
                     inherits = FALSE)
    on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
  }
  
  # Check input parameters:
  if (!is.array(X)) {
    stop("The input data must be an array.")
  }
  if (dim(X)[2] < 2 || dim(X)[3] < 2) {
    stop("The input data array must have at least 2 variables in each mode.")
  }
  if (is.null(MacroPARAFACpars)) {
    MacroPARAFACpars <- list()
  }
  if (!is.list(MacroPARAFACpars)) {
    stop("MacroPARAFACpars must be a list")
  }
  if (is.null(const)){
    const <- rep("uncons", 4)
  }
  # parameters for DDC
  if (!"DDCpars" %in% names(MacroPARAFACpars)) {
    MacroPARAFACpars$DDCpars <- list()
  }
  if (!"fracNA" %in% names(MacroPARAFACpars$DDCpars)) {
    MacroPARAFACpars$DDCpars$fracNA <-  0.5
  }
  if (!"numDiscrete" %in% names(MacroPARAFACpars$DDCpars)) {
    MacroPARAFACpars$DDCpars$numDiscrete <- 3
  }
  if (!"precScale" %in% names(MacroPARAFACpars$DDCpars)) {
    MacroPARAFACpars$DDCpars$precScale <- 1e-12
  }
  if (!"cleanNAfirst" %in% names(MacroPARAFACpars$DDCpars)) {
    MacroPARAFACpars$DDCpars$cleanNAfirst <- "rows"
  }
  if (!"tolProb" %in% names(MacroPARAFACpars$DDCpars)) {
    MacroPARAFACpars$DDCpars$tolProb <-  0.998
  }
  if (!"silent" %in% names(MacroPARAFACpars$DDCpars)) {
    MacroPARAFACpars$DDCpars$silent <- if (is.null(MacroPARAFACpars$silent)) {
      FALSE
    } else {
      MacroPARAFACpars$silent
    }
  }
  
  MacroPARAFACpars$DDCpars$fastDDC=TRUE
  
  # parameters for MacroPARAFAC
  if (!"alpha" %in% names(MacroPARAFACpars)) {
    MacroPARAFACpars$alpha <- 0.750
  }
  if (!"maxdir" %in% names(MacroPARAFACpars)) {
    MacroPARAFACpars$maxdir <- 250
  }
  if (!"distprob" %in% names(MacroPARAFACpars)) {
    MacroPARAFACpars$distprob <- 0.998
  }
  if (!"ctol" %in% names(MacroPARAFACpars)) {
    MacroPARAFACpars$ctol <- 1e-04
  }
  if (!"silent" %in% names(MacroPARAFACpars)) {
    MacroPARAFACpars$silent <- FALSE
  }
  
  # Put options in convenient order
  MacroPARAFACpars <- MacroPARAFACpars[c("DDCpars", "alpha", "maxdir",
                                         "distprob", "ctol", "silent")]
  # Input arguments
  DDCpars   <- MacroPARAFACpars$DDCpars
  alpha     <- MacroPARAFACpars$alpha
  maxdir    <- MacroPARAFACpars$maxdir
  distprob  <- MacroPARAFACpars$distprob
  ctol      <- MacroPARAFACpars$ctol
  silent    <- MacroPARAFACpars$silent
  
  if (alpha < 0.5 | alpha > 1)
    stop("Alpha is out of range: should be between 1/2 and 1.")
  
  # Check Dataset
  # needs to be done first, as some rows/columns may leave the analysis.
  X_ar <- X
  size <- dim(X_ar)
  J <- size[2]; K <- size[3];L<-size[4]
  X <- matrix(X_ar, ncol = J*K*L)
  
  checkOut <- cellWise::checkDataSet(X, DDCpars$fracNA, DDCpars$numDiscrete,
                                     DDCpars$precScale, DDCpars$silent, 
                                     DDCpars$cleanNAfirst)
  X <- checkOut$remX
  X <- as.matrix(X) 
  if (ncol(X) < J*K*L) stop("DDC has removed some columns of the unfolded array.")
  I <- nrow(X)
  X_ar <- array(X, c(I, J, K, L))
  
  h <- ceiling(alpha*(I+1)) 
  
  ##################################
  ##  Step 1: Detect deviating cells
  ##################################
  
  DDCpars$coreOnly <- TRUE # now only execute DDCcore, CheckDataSet already done
  resultDDC <- cellWise::DDC(X, DDCpars)
  DDCpars   <- resultDDC$DDCpars
  MacroPARAFACpars$DDCpars <- DDCpars 
  Ti <- resultDDC$Ti
  indrows1 <- which(abs(Ti) > sqrt(qchisq(DDCpars$tolProb, 1)))
  indrows2 <- intersect(order(abs(Ti), decreasing = T)[seq_len(I - h)],
                        indrows1)
  indcells <- resultDDC$indcells
  
  Xnai      <- X   # initialize NA-imputed data matrix
  Xci       <- X   # initialize cell-imputed data matrix
  
  indNA     <- which(is.na(Xci))          # positions of missing values
  indimp    <- unique(c(indcells, indNA)) # positions where values can be imputed
  
  Xci[indimp] <- resultDDC$Ximp[indimp]   # imputes all missings and cellwise outliers
  
  ############################
  # Step 2: Projection Pursuit
  ############################
  
  # Compute initial cell-imputed matrix
  
  XOimp <- Xci  # store imputed rowoutliers
  
  # Xind is TRUE for elements that are missing or cellwise (as matrix)
  Xind <- Xci; Xind[indimp] <- NA; Xind <- is.na(Xind)
  Xindna <- Xci
  
  # impute original NA values 
  Xnai[indNA] <- XOimp[indNA]
  
  # classical PCA to determine the rank of the data matrix
  XciSVD <- cellWise::truncPC(Xci)
  if (XciSVD$rank == 0) stop("All data points collapse!")
  
  rowcellind   <- tabulate(((indcells - 1) %% I) + 1, I)
  rowcellind[indrows2] <- J*K*L # pretends all cells flagged, so
  # the rows flagged as outlying by DDC will not be imputed.
  hrowsLowcell <- order(rowcellind)[seq_len(h)]
  indComb <- unique(hrowsLowcell,
                    which(rowcellind == rowcellind[hrowsLowcell[h]]))
  # these rows have the fewest cellwise outliers (includes ties)
  indComb <- setdiff(seq_len(I), indComb) # take complement
  
  # Use unimputed rowwise outliers for least outlying points
  # (original NA values are imputed)
  if (length(indComb) > 0) Xci[indComb, ] <- Xnai[indComb, ]
  
  # Find least outlying points by projections:
  alldir <- choose(I, 2)
  ndir   <- min(maxdir, alldir)
  all    <- (ndir == alldir)
  # calculate ndirect directions through two random data points
  D <- makeDirections(Xci, ndir, all = all)
  Dnorm <- vector(mode = "numeric", length = nrow(D))
  Dnorm <- apply(D, 1, vecnorm)
  Dnormr <- Dnorm[Dnorm > 1e-12]
  m <- length(Dnormr)
  D <- D[Dnorm > 1e-12, ]
  V <- diag(1 / Dnormr) %*% D
  # projected points in columns
  Y <- Xci %*% t(V)
  Z <- matrix(0, I, m)
  
  umcds <- cellWise::estLocScale(Y, nLocScale = DDCpars$nLocScale,
                                 type = "mcd", alpha = alpha, silent = TRUE)
  zeroscales <- which(umcds$scale < 1e-12)
  
  for (i in seq_len(length(zeroscales))) {
    umcdweights <- unimcd(Y[, zeroscales[i]], alpha = alpha)$weights
    if (robustbase::rankMM(Xci[umcdweights == 1, ]) == 1) {
      stop("At least ", sum(umcdweights),
           " observations are identical.")
    }
  }
  
  Z <- abs(scale(Y, umcds$loc, umcds$scale))
  
  H0 <- order(apply(Z, 1, max))
  H0 <- setdiff(H0, indrows2)
  H0 <- H0[seq_len(h)]
  indrows3 <- setdiff((seq_len(I)), H0) # rows not in H0
  
  #################################
  # Step 3: Initialisation loadings
  #################################
  
  Xind_ar <- array(Xind, c(I, J, K, L))
  Xnai_ar <- array(Xnai, c(I, J, K, L))
  
  # make new cell-imputed matrix
  # cell-imputed for the rows of H0, NA for the others
  Xci <- XOimp
  if (length(indrows3) > 0)  Xci[indrows3, ] <- Xnai[indrows3, ]
  Xci_ar <- array(Xci, c(I, J, K, L))
  
  # fully iterated PARAFAC on cell-imputed matrix (only H_0)
  # to find starting value for B and C loadings
  Xcih_ar <- Xci_ar[H0, , ,]
  Xcih.parafac <- multiway::parafac(Xcih_ar, nfac = ncomp, const = const, 
                                    ctol = ctol)
  B.H0init <- Xcih.parafac$B
  C.H0init <- Xcih.parafac$C
  D.H0init <- Xcih.parafac$D
  
  # update NAs and missing values in Xcih 
  krBCD.H0 <- rrcov3way::krp(rrcov3way::krp(D.H0init, C.H0init), B.H0init)
  Xhat.H0 <- tcrossprod(Xcih.parafac$A, krBCD.H0)
  Xhat.H0_ar <- array(Xhat.H0, c(h, J, K,L))
  Xindh_ar <- Xind_ar[H0, , , ] 
  Xcih_ar[Xindh_ar] <- Xhat.H0_ar[Xindh_ar]
  indimph  <- which(Xindh_ar) # indicates the cells from H0 to be imputed 
  
  ##############################
  # Step 4: Iterative estimation
  ##############################
  
  # PARAFAC with updating the NAs and missing values 
  Xcih.parafac.si <- multiway::parafac_4wayna(Xcih_ar, nfac = ncomp, 
                                              naid = indimph,
                                              const = const,
                                              Bstart = B.H0init, 
                                              Cstart = C.H0init,
                                              Dstart = D.H0init,
                                              ctol = ctol) 
  
  # scores and fitted values for all samples
  B.H0 <- Xcih.parafac.si$B
  C.H0 <- Xcih.parafac.si$C
  D.H0 <- Xcih.parafac.si$D
  krBCD.H0 <- rrcov3way::krp(rrcov3way::krp(D.H0init, C.H0init), B.H0init)
  invkrBCD.H0 <- pracma::pinv(t(krBCD.H0))
  # for H0 the scores are computed in parafac_4wayna
  # for the complement we use the DDC imputed matrix
  A <- XOimp %*% invkrBCD.H0
  Xhat <- tcrossprod(A, krBCD.H0)
  Xhat[H0, ] <- tcrossprod(Xcih.parafac.si$A, krBCD.H0)
  Xhat_ar <- array(Xhat, c(I, J, K, L))
  
  ## Impute cellwise and missing values in Xci
  Xci_ar[Xind_ar] <- Xhat_ar[Xind_ar]
  Xnai_ar[indNA] <- Xci_ar[indNA] ## update the missing values in Xnai
  
  Xfi_ar <- Xci_ar ## initialize the fully imputed Xfi
  if (length(indrows3) > 0) Xci_ar[indrows3, , , ] <- Xnai_ar[indrows3, , , ]
  
  Xci <- matrix(Xci_ar, ncol = J*K*L)
  Xfi <- matrix(Xfi_ar, ncol = J*K*L)
  
  ######################
  # Step 5: Re-weighting
  ######################
  
  ## new fitted values of Xci (since Xci has been updated)
  A <- Xci %*% invkrBCD.H0
  Xhat <- tcrossprod(A, krBCD.H0)
  
  RD <- apply(Xci - Xhat, 1, vecnorm)
  umcd <- unimcd(RD^(2/3), alpha = alpha)
  mRD <- umcd$loc; sRD <- umcd$scale
  cutoffRD <- sqrt(qnorm(distprob, mRD, sRD)^3)
  
  Hstar <- (RD <= cutoffRD)
  Hstar[indrows2] <- FALSE # removes outlying rows flagged by DDC
  Xindh_ar <- Xind_ar[Hstar, , , ] 
  indimph  <- which(Xindh_ar) # indicates the cells from Hstar to be imputed 
  
  Xrew.parafac.si <- multiway::parafac_4wayna(Xfi_ar[Hstar, , , ], 
                                              nfac = ncomp, 
                                              naid = indimph, 
                                              const = const,
                                              Bstart = B.H0, 
                                              Cstart = C.H0,
                                              Dstart = D.H0,
                                              ctol = ctol) 
  
  B.final <- Xrew.parafac.si$B
  C.final <- Xrew.parafac.si$C
  D.final <- Xrew.parafac.si$D
  krBCDfinal <- rrcov3way::krp(rrcov3way::krp(D.final, C.final), B.final)
  invkrBCDfinal <- pracma::pinv(t(krBCDfinal))
  
  # for Hstar the scores are computed in parafac_3wayna
  # for the complement we use Xci
  A <- Xci %*% invkrBCDfinal
  Xhat <- tcrossprod(A, krBCDfinal)
  Xhat[Hstar, ] <- tcrossprod(Xrew.parafac.si$A, krBCDfinal)
  Xhat_ar <- array(Xhat, c(I, J, K, L))
  
  ## Impute cellwise and missing values in Xfi
  Xfi_ar[Xind_ar] <- Xhat_ar[Xind_ar]
  Xnai_ar[indNA] <- Xfi_ar[indNA] ## update the missing values in Xnai
  Xci_ar <- Xfi_ar 
  
  # in Xci: do not update the rows not in Hstar 
  indrows4 <- which(!Hstar)
  if (length(indrows3) > 0) {
    Xci_ar[indrows4, , , ] <- Xnai_ar[indrows4, , , ]
  }
  
  Xfi <- matrix(Xfi_ar, ncol = J*K*L) # final fully imputed matrix
  
  ###########################################
  ## Step 6: Scores, distances and residuals
  ###########################################
  
  # recompute fitted values (cell-imputed no longer needed)
  A.final <- Xfi %*% invkrBCDfinal
  Xhat <- tcrossprod(A.final, krBCDfinal)
  Xhat_ar <- array(Xhat, c(I, J, K, L))
  
  # replace NAs such that their residual is zero (but do not change
  # the fitted values again)
  Xfi[indNA] <- Xhat[indNA]
  Xnai[indNA] <- Xfi[indNA]
  Xfi_ar <- array(Xfi, c(I, J, K, L))
  Xnai_ar <- array(Xnai, c(I, J, K, L)) 
  
  # scores distances
  SD <- sqrt(covMcd(A.final, nsamp = "deterministic", alpha = alpha)$mah)
  cutoffSD <- sqrt(qchisq(distprob, ncol(A.final)))
  highSD <- which(SD > cutoffSD)
  
  # standardised residuals and the imputed residuals 
  resid <- Xnai - Xhat
  RD <- apply(resid, 1, vecnorm)
  highRD <- which(RD > cutoffRD)
  
  resid <- resid + 0 * X
  residScale <- cellWise::estLocScale(resid,
                                      nLocScale = DDCpars$nLocScale,
                                      type = "1stepM",
                                      precScale = DDCpars$precScale,
                                      center = FALSE)$scale
  stdResid <- sweep(resid, 2, residScale, "/")
  indcells <- which(abs(stdResid) > sqrt(qchisq(DDCpars$tolProb, 1)))
  # stdResid_ar <- array(stdResid, c(I, J, K))
  
  resid.fi <- Xfi - Xhat
  RD.fi <- apply(resid.fi, 1, vecnorm)
  highRD.fi <- which(RD.fi > cutoffRD)
  
  Fullimp <- list(SD = SD, cutoffSD = cutoffSD,
                  RD = RD.fi, cutoffRD = cutoffRD,
                  highRD = highRD.fi)
  
  res <- list(A = A.final, B = B.final, C = C.final, D = D.final,
              Xhat = Xhat_ar,
              Xci = Xci_ar, Xnai = Xnai_ar, Xfi = Xfi_ar, 
              residScale = residScale,
              stdResid = stdResid, indcells = indcells, 
              RD = RD, cutoffRD = cutoffRD,
              SD = SD, cutoffSD = cutoffSD,
              highSD = highSD, highRD = highRD, 
              DDC = resultDDC, 
              Fullimp = Fullimp
  )
  
  return(res)
  
}

# ============================================================================ #
#                      Función para generar una elipse                         #
# ============================================================================ #

generate_ellipse_image <- function(center, color, radius_x, 
                                   radius_y, img_shape) {
  # Crear imagen negra
  image <- array(0, dim = img_shape)
  
  # Crear rejilla de coordenadas
  x <- matrix(rep(1:img_shape[2], each = img_shape[1]), nrow = img_shape[1])
  y <- matrix(rep(1:img_shape[1], img_shape[2]), nrow = img_shape[1])
  
  # Máscara de la elipse
  mask <- (((x - center[1]) / radius_x)^2 + ((y - center[2]) / radius_y)^2) <= 1
  
  # Aplicar color canal por canal
  for (c in 1:img_shape[3]) {
    channel <- image[,,c]
    channel[mask] <- color[c]
    image[,,c] <- channel
  }
  
  # Convertir a tipo double si se desea trabajar con valores reales
  image <- image
  
  return(image)
}

# ============================================================================ #
#                            Generar muestras                                  #
# ============================================================================ #

generate_samples<- function(n_samples, img_shape, d1, d2, d3, d4, 
                            missing_percentage) {
  images <- array(0, dim = c(n_samples, img_shape))
  
  for (i in 1:n_samples) {
    # Centro con desplazamiento
    center <- round(rnorm(2, mean = 30, sd = 1)) + d1
    
    # Color RGB modificado por d4
    color <- round(rnorm(3, 
                         mean = c(100 + d4 * 1, 20 + d4 * 0.04, 15 + d4 * 0.0225), 
                         sd = c(1, 0.04, 0.0225)))
    color <- pmin(pmax(color, 0), 255)
    
    # Radios con distorsiones d2 y d3
    radius_x <- 10 * (1 + d2) * (1 + d3)
    radius_y <- 10 * (1 + d2) / (1 + d3)
    
    # Imagen base
    image <- generate_ellipse_image(center, color, radius_x, radius_y, img_shape)
    
    # Total de píxeles
    total_pixels <- prod(img_shape)
    n_missing <- round(total_pixels * missing_percentage / 100)
    
    # Índices aleatorios a borrar (linealizados)
    missing_indices <- sample(total_pixels, n_missing)
    
    # Insertar NA y rearmar imagen
    image_vec <- as.vector(image)
    image_vec[missing_indices] <- NA
    image <- array(image_vec, dim = img_shape)
    
    # Guardar la imagen con valores NA
    images[i,,,] <- image
  }
  
  return(images)
}

# ============================================================================ #
#                         Reconstruir tensor con TROD                          #
# ============================================================================ #

recon_trod_tensor <- function(lambda, factors) {

  U2 <- factors[[2]]
  U3 <- factors[[3]]
  U4 <- factors[[4]]
  
  # Khatri-Rao product de C y B, resultado (J*K x R)
  U43_kr <- khatri_rao(U4, U3)
  
  # Multiplicación para obtener cada frontal: A %*% diag(lambda) %*% t(CB_kr)
  unfolded_tensor <- U2 %*% diag(as.vector(lambda)) %*% t(U43_kr)  # Resultado: (I x JK)
  
  # Finalmente, volver al arreglo 3D
  dims <- c(nrow(U2), nrow(U3), nrow(U4))
  array(unfolded_tensor, dim = dims)
}

# ============================================================================ #
#                      Funciones auxiliares de MacroPARAFAC                    #
# ============================================================================ #

vecnorm <- function(x, p = 2) {
  sum(x ^ p) ^ (1 / p) 
}

makeDirections <- function(data, ndirect, all = TRUE){ # extracts directions between pairs of points
  #
  # We first make two functions:
  
  uniran <- function(seed = 0) {
    seed <- floor(seed * 5761) + 999
    quot <- floor(seed / 65536)
    seed <- floor(seed) - floor(quot * 65536)
    random <- seed / 65536
    list(seed = seed, random = random)
  }
  
  randomset <- function(n, k, seed) {
    ranset <- vector(mode = "numeric", length = k)
    for (j in seq_len(k)) {
      r <- uniran(seed)
      seed <- r$seed
      num <- floor(r$random * n) + 1
      if (j > 1) {
        while (any(ranset == num)) {
          r <- uniran(seed)
          seed <- r$seed
          num <- floor(r$random * n) + 1
        }
      }
      ranset[j] <- num
    }
    ans <- list()
    ans$seed <- seed
    ans$ranset <- ranset
    ans
  }
  
  # Now the actual computation:
  if (all) {
    cc <- utils::combn(nrow(data), 2)
    B2 <- data[cc[1, ], ] - data[cc[2, ], ]
  } else {
    n <- nrow(data)
    p <- ncol(data)
    r <- 1
    B2 <- matrix(0, ndirect, p)
    seed <- 0
    while (r <= ndirect) {
      sseed <- randomset(n, 2, seed)
      seed  <- sseed$seed
      B2[r, ] <- data[sseed$ranset[1], ] - data[sseed$ranset[2], ]
      r <- r + 1
    }
  }
  return(B2)
}

unimcd <- function (y, alpha = NULL, center = TRUE) {
  if (is.null(alpha)) {
    alpha <- 0.5
  }
  center <- as.numeric(center)
  res <- tryCatch(.Call("_cellWise_unimcd_cpp", y, alpha, center, 
                        PACKAGE = "cellWise"), `std::range_error` = function(e) {
                          conditionMessage(e)
                        })
  return(list(loc = res$loc, scale = res$scale, weights = drop(res$weights)))
}


# ============================================================================ #
#                 Buscar mejor alpha y Criterio de selección                   #
# ============================================================================ #

Best_alpha <- function(i, delta = 0, n_samples = 80, prop_in_control = 0.9,
                       d1delta=0, d2delta=0, d3delta=0, d4delta=0) {
  n_in_control <- round(n_samples * prop_in_control)
  n_out_control <- n_samples - n_in_control
  
  img_shape = c(60, 60, 3)
  rank_cp = 3
  rank_macro = 7

  tensor_in_control <- generate_samples(
    n_samples = n_in_control,
    img_shape = img_shape,
    d1 = c(0, 0),
    d2 = 0,
    d3 = 0,
    d4 = 0,
    missing_percentage = 0
  )
  
  # Generar muestras fuera de control (con delta)
  tensor_out_control <- generate_samples(
    n_samples = n_out_control,
    img_shape = img_shape,
    d1 = c(d1delta, d1delta),
    d2 = d2delta,
    d3 = d3delta,
    d4 = d4delta,
    missing_percentage = 0
  )
  # Unir los datos
  tensor_data <- abind::abind(tensor_in_control, tensor_out_control, along = 1)
  
  # Reordenar las imágenes
  perm_idx <- sample(1:dim(tensor_data)[1])
  tensor_data <- tensor_data[perm_idx,,,]
  
  X <- tensor_data + array(rnorm(length(tensor_data), mean = 0, sd = 0.01), dim = dim(tensor_data))
  
  ##### TROD - CP decomposition
  InicioCP <- Sys.time()
  cp_result <- cp(as.tensor(X), num_components = rank_cp)
  FinalCP <- Sys.time()
  TimeCP <- FinalCP - InicioCP
  
  normU1 <- apply(cp_result$U[[1]], 2, vecnorm)
  normU2 <- apply(cp_result$U[[2]], 2, vecnorm)
  normU3 <- apply(cp_result$U[[3]], 2, vecnorm)
  normU4 <- apply(cp_result$U[[4]], 2, vecnorm)
  
  U1CP <- sweep(cp_result$U[[1]], 2, normU1, "/")
  U2CP <- sweep(cp_result$U[[2]], 2, normU2, "/")
  U3CP <- sweep(cp_result$U[[3]], 2, normU3, "/")
  U4CP <- sweep(cp_result$U[[4]], 2, normU4, "/")
  
  UpdateLambdasCP <- cp_result$lambdas * normU1 * normU2 * normU3 * normU4
  
  lambdas_trod <- matrix(0, nrow = n_samples, ncol = rank_cp)
  for (ii in 1:n_samples) {
    for (j in 1:rank_cp) {
      lambdas_trod[ii, j] <- UpdateLambdasCP[j] * U1CP[ii, j]
    }
  }
  
  ##### MPCA
  Xmpca <- aperm(X, c(4, 2, 3, 1))
  InicioMPCA <- Sys.time()
  MPCA_result <- mpca(as.tensor(Xmpca), ranks = c(1, 2, 3))
  FinalMPCA <- Sys.time()
  TimeMPCA <- FinalMPCA - InicioMPCA
  
  ##### MACRO-PARAFAC
  InicioMACRO <- Sys.time()
  
  MACRO_result <-try(MacroPARAFAC2(X, ncomp = rank_macro))
  
  FinalMACRO <- Sys.time()
  TimeMACRO <- FinalMACRO - InicioMACRO
  
  if (inherits(MACRO_result, "try-error")) {
    MACRO_result <- NULL 
  }
  
  if (!is.null(MACRO_result)) {
    normaA <- apply(MACRO_result$A, 2, vecnorm)
    normaB <- apply(MACRO_result$B, 2, vecnorm)
    normaC <- apply(MACRO_result$C, 2, vecnorm)
    normaD <- apply(MACRO_result$D, 2, vecnorm)
    
    A_norm <- sweep(MACRO_result$A, 2, normaA, "/")
    B_norm <- sweep(MACRO_result$B, 2, normaB, "/")
    C_norm <- sweep(MACRO_result$C, 2, normaC, "/")
    D_norm <- sweep(MACRO_result$D, 2, normaD, "/")
    
    UpdateLambdasMACRO <- normaA * normaB * normaC * normaD
    
    lambdas_MACRO <- matrix(0, nrow = n_samples, ncol = rank_macro)
    for (ii in 1:n_samples) {
      for (j in 1:rank_macro) {
        lambdas_MACRO[ii, j] <- UpdateLambdasMACRO[j] * A_norm[ii, j]
      }
    }
  }
  
  ##### Inicialización para reconstrucciones
  num_elements <- prod(img_shape)
  theta_TROD <- matrix(0, rank_cp, n_samples)
  theta_MPCA <- matrix(0, prod(dim(MPCA_result$Z_ext)[-4]), n_samples)
  theta_MACRO <- matrix(0, rank_macro, n_samples)
  
  for (ii in 1:n_samples) {
    y <- X[ii,,,]
    y_tensor <- as.tensor(y)
    
    # TROD
    z_trod <- recon_trod_tensor(lambdas_trod[ii, ], list(U1CP, U2CP, U3CP, U4CP))
    theta_TROD[, ii] <- lambdas_trod[ii, ]
    
    # MPCA
    z_mpca <- ttm(ttm(ttm(y_tensor, t(MPCA_result$U[[2]]), m = 1),
                      t(MPCA_result$U[[3]]), m = 2),
                  t(MPCA_result$U[[1]]), m = 3)
    
    y_recon_mpca <- ttm(ttm(ttm(z_mpca, MPCA_result$U[[2]], m = 1),
                            MPCA_result$U[[3]], m = 2),
                        MPCA_result$U[[1]], m = 3)
    
    theta_MPCA[, ii] <- as.vector(z_mpca@data)
    
    # MACRO
    if (!is.null(MACRO_result)) {z_macro <- recon_trod_tensor(lambdas_MACRO[ii, ], list(A_norm, B_norm, C_norm, D_norm))
    theta_MACRO[, ii] <- lambdas_MACRO[ii, ]
    }
  }
  
  theta_TROD <- t(theta_TROD)
  
  theta_MPCA <- t(theta_MPCA)
  
  if (!is.null(MACRO_result)) {
    theta_MACRO <- t(theta_MACRO)
  }
  
  # Buscar los mejores alpha para MCD -----------------------------------------
  
  result_TROD  <- tryCatch(bootstrap_mcd(theta_TROD,  seq(0.5,0.975,by=0.025), B=50, classifier="MD"),
                           error = function(e) list(best_alpha = NA_real_))
  
  result_MPCA  <- tryCatch(bootstrap_mcd(theta_MPCA,  seq(0.5,0.975,by=0.025), B=50, classifier="MD"),
                           error = function(e) list(best_alpha = NA_real_))
  
  result_MACRO <- tryCatch(bootstrap_mcd(theta_MACRO, seq(0.5,0.975,by=0.025), B=50, classifier="MD"),
                           error = function(e) list(best_alpha = NA_real_))
  
  
  # ---------  Criterio de Yan & Paynabar ---------------------------------
  
  X_tensor <- as.tensor(X)
  ## -------- TROD --------
  P2_TROD <- proj_qr(U2CP)
  P3_TROD <- proj_qr(U3CP)
  P4_TROD <- proj_qr(U4CP)
  
  X_TEST_TROD <- ttm(ttm(ttm(X_tensor, P2_TROD, 2),
                         P3_TROD, 3),
                     P4_TROD, 4)
  
  CRIT_TROD <- (rTensor::fnorm(X_TEST_TROD) / rTensor::fnorm(X_tensor))^2
  
  
  ## -------- MPCA --------
  P2_MPCA <- proj_qr(MPCA_result$U[[2]])
  P3_MPCA <- proj_qr(MPCA_result$U[[3]])
  P4_MPCA <- proj_qr(MPCA_result$U[[1]])
  
  X_TEST_MPCA <- ttm(ttm(ttm(X_tensor, P2_MPCA, 2),
                         P3_MPCA, 3),
                     P4_MPCA, 4)
  CRIT_MPCA <- (rTensor::fnorm(X_TEST_MPCA) / rTensor::fnorm(X_tensor))^2
  
  
  ## -------- MACRO --------
  P2_MACRO <- proj_qr(B_norm)
  P3_MACRO <- proj_qr(C_norm)
  P4_MACRO <- diag(3)  
  
  X_TEST_MACRO <- ttm(ttm(ttm(X_tensor, P2_MACRO, 2),
                          P3_MACRO, 3),
                      P4_MACRO, 4)
  CRIT_MACRO <- (rTensor::fnorm(X_TEST_MACRO) / rTensor::fnorm(X_tensor))^2
  

  return(data.frame(
    # Valores de las carta T2
    alpha_T2_TROD_mcd = result_TROD$best_alpha,
    alpha_T2_MPCA_mcd = result_MPCA$best_alpha,
    alpha_T2_MACRO_mcd = result_MACRO$best_alpha,
    Crit_TROD=CRIT_TROD,
    Crit_MPCA=CRIT_MPCA,
    Crit_MACRO=CRIT_MACRO
  ))
}

proj_qr <- function(U) {
  qU <- qr(U)
  Q  <- qr.Q(qU)      # columnas ortonormales que generan el mismo subespacio
  Q %*% t(Q)          # proyector QQ'
}

# ============================================================================ #
# ---------------------------------------------------------------------------- #
#                        Aplicación con datos reales                           #
# ---------------------------------------------------------------------------- #
# ============================================================================ #

load_tensor_images <- function(dir_path, id_range = 0:208, h = 100, w = 100) {
  # Usa nombre corto si Windows
  dir_path <- tryCatch(utils::shortPathName(dir_path), error = function(e) dir_path)
  
  # Intenta con 000.png...208.png, si no existen, usa list.files
  candidate_files <- file.path(dir_path, sprintf("%03d.png", id_range))
  existing_files <- candidate_files[file.exists(candidate_files)]
  
  if (length(existing_files) == 0L) {
    existing_files <- list.files(dir_path, pattern = "\\.png$", full.names = TRUE, ignore.case = TRUE)
    existing_files <- existing_files[order(tolower(basename(existing_files)))]
  }
  if (length(existing_files) == 0L) stop("No se encontraron imágenes PNG en la carpeta indicada.")
  
  # Función auxiliar
  load_and_resize <- function(fp) {
    img <- image_read(fp)
    img <- image_convert(img, colorspace = "RGB")
    img <- image_resize(img, paste0(w, "x", h, "!"))
    arr <- image_data(img, channels = "rgb")
    arr <- aperm(arr, c(3, 2, 1))
    storage.mode(arr) <- "integer"
    array(as.integer(arr), dim = dim(arr))
  }
  
  # Prealocar y cargar
  N <- length(existing_files)
  probe <- load_and_resize(existing_files[1])
  tensor <- array(NA_integer_, dim = c(N, dim(probe)))
  tensor[1,,,] <- probe
  
  if (N > 1) {
    for (i in 2:N) {
      tensor[i,,,] <- load_and_resize(existing_files[i])
      if (i %% 20 == 0 || i == N) cat(sprintf("[%s] %d/%d\n", basename(dir_path), i, N))
    }
  }
  tensor
}

# ============================================================================ #
#                    Funcion para hallar var con bootstrap                     #
# ============================================================================ #

Best_alpha_real <- function(i, X_IMAGES, n_samples = 209) {
  
  img_shape = c(100, 100, 3)
  rank_cp =4 
  rank_macro = 4
  core_mpca=c(1,4,4)
  
  idx_sample <- sample.int(209, n_samples, replace = TRUE)
  X_sample   <- X_IMAGES[idx_sample,,, , drop = FALSE]
  X <- X_sample + array(rnorm(length(X_IMAGES), mean = 0, sd = 0.01), dim = dim(X_IMAGES))
  
  ##### TROD - CP decomposition
  InicioCP <- Sys.time()
  cp_result <- cp(as.tensor(X), num_components = rank_cp)
  FinalCP <- Sys.time()
  TimeCP <- FinalCP - InicioCP
  
  normU1 <- apply(cp_result$U[[1]], 2, vecnorm)
  normU2 <- apply(cp_result$U[[2]], 2, vecnorm)
  normU3 <- apply(cp_result$U[[3]], 2, vecnorm)
  normU4 <- apply(cp_result$U[[4]], 2, vecnorm)
  
  U1CP <- sweep(cp_result$U[[1]], 2, normU1, "/")
  U2CP <- sweep(cp_result$U[[2]], 2, normU2, "/")
  U3CP <- sweep(cp_result$U[[3]], 2, normU3, "/")
  U4CP <- sweep(cp_result$U[[4]], 2, normU4, "/")
  
  UpdateLambdasCP <- cp_result$lambdas * normU1 * normU2 * normU3 * normU4
  
  lambdas_trod <- matrix(0, nrow = n_samples, ncol = rank_cp)
  for (ii in 1:n_samples) {
    for (j in 1:rank_cp) {
      lambdas_trod[ii, j] <- UpdateLambdasCP[j] * U1CP[ii, j]
    }
  }
  
  ##### MPCA
  Xmpca <- aperm(X, c(4, 2, 3, 1))
  InicioMPCA <- Sys.time()
  MPCA_result <- mpca(as.tensor(Xmpca), ranks = core_mpca)
  FinalMPCA <- Sys.time()
  TimeMPCA <- FinalMPCA - InicioMPCA
  
  ##### MACRO-PARAFAC
  InicioMACRO <- Sys.time()
  
  MACRO_result <-try(MacroPARAFAC2(X, ncomp = rank_macro))
  
  FinalMACRO <- Sys.time()
  TimeMACRO <- FinalMACRO - InicioMACRO
  
  if (inherits(MACRO_result, "try-error")) {
    MACRO_result <- NULL 
  }
  
  if (!is.null(MACRO_result)) {
    normaA <- apply(MACRO_result$A, 2, vecnorm)
    normaB <- apply(MACRO_result$B, 2, vecnorm)
    normaC <- apply(MACRO_result$C, 2, vecnorm)
    normaD <- apply(MACRO_result$D, 2, vecnorm)
    
    A_norm <- sweep(MACRO_result$A, 2, normaA, "/")
    B_norm <- sweep(MACRO_result$B, 2, normaB, "/")
    C_norm <- sweep(MACRO_result$C, 2, normaC, "/")
    D_norm <- sweep(MACRO_result$D, 2, normaD, "/")
    
    UpdateLambdasMACRO <- normaA * normaB * normaC * normaD
    
    lambdas_MACRO <- matrix(0, nrow = n_samples, ncol = rank_macro)
    for (ii in 1:n_samples) {
      for (j in 1:rank_macro) {
        lambdas_MACRO[ii, j] <- UpdateLambdasMACRO[j] * A_norm[ii, j]
      }
    }
  }
  
  ##### Inicialización para reconstrucciones
  num_elements <- prod(img_shape)
  theta_TROD <- matrix(0, rank_cp, n_samples)
  theta_MPCA <- matrix(0, prod(dim(MPCA_result$Z_ext)[-4]), n_samples)
  theta_MACRO <- matrix(0, rank_macro, n_samples)
  
  for (ii in 1:n_samples) {
    y <- X[ii,,,]
    y_tensor <- as.tensor(y)
    
    # TROD
    z_trod <- recon_trod_tensor(lambdas_trod[ii, ], list(U1CP, U2CP, U3CP, U4CP))
    theta_TROD[, ii] <- lambdas_trod[ii, ]
    
    # MPCA
    z_mpca <- ttm(ttm(ttm(y_tensor, t(MPCA_result$U[[2]]), m = 1),
                      t(MPCA_result$U[[3]]), m = 2),
                  t(MPCA_result$U[[1]]), m = 3)
    
    y_recon_mpca <- ttm(ttm(ttm(z_mpca, MPCA_result$U[[2]], m = 1),
                            MPCA_result$U[[3]], m = 2),
                        MPCA_result$U[[1]], m = 3)
    
    theta_MPCA[, ii] <- as.vector(z_mpca@data)
    
    # MACRO
    if (!is.null(MACRO_result)) {z_macro <- recon_trod_tensor(lambdas_MACRO[ii, ], list(A_norm, B_norm, C_norm, D_norm))
    theta_MACRO[, ii] <- lambdas_MACRO[ii, ]
    }
  }
  
  theta_TROD <- t(theta_TROD)
  
  theta_MPCA <- t(theta_MPCA)
  
  if (!is.null(MACRO_result)) {
    theta_MACRO <- t(theta_MACRO)
  }
  
  # Buscar los mejores alpha para MCD -----------------------------------------
  
  result_TROD  <- tryCatch(bootstrap_mcd(theta_TROD,  seq(0.5,0.975,by=0.025), B=50, classifier="MD"),
                           error = function(e) list(best_alpha = NA_real_))
  
  result_MPCA  <- tryCatch(bootstrap_mcd(theta_MPCA,  seq(0.5,0.975,by=0.025), B=50, classifier="MD"),
                           error = function(e) list(best_alpha = NA_real_))
  
  result_MACRO <- tryCatch(bootstrap_mcd(theta_MACRO, seq(0.5,0.975,by=0.025), B=50, classifier="MD"),
                           error = function(e) list(best_alpha = NA_real_))
  
  
  # ---------  Criterio de Yan & Paynabar ---------------------------------
  
  X_tensor <- as.tensor(X)
  ## -------- TROD --------
  P2_TROD <- proj_qr(U2CP)
  P3_TROD <- proj_qr(U3CP)
  P4_TROD <- proj_qr(U4CP)
  
  X_TEST_TROD <- ttm(ttm(ttm(X_tensor, P2_TROD, 2),
                         P3_TROD, 3),
                     P4_TROD, 4)
  
  CRIT_TROD <- (rTensor::fnorm(X_TEST_TROD) / rTensor::fnorm(X_tensor))^2
  
  
  ## -------- MPCA --------
  P2_MPCA <- proj_qr(MPCA_result$U[[2]])
  P3_MPCA <- proj_qr(MPCA_result$U[[3]])
  P4_MPCA <- proj_qr(MPCA_result$U[[1]])
  
  X_TEST_MPCA <- ttm(ttm(ttm(X_tensor, P2_MPCA, 2),
                         P3_MPCA, 3),
                     P4_MPCA, 4)
  CRIT_MPCA <- (rTensor::fnorm(X_TEST_MPCA) / rTensor::fnorm(X_tensor))^2
  
  
  ## -------- MACRO --------
  P2_MACRO <- proj_qr(B_norm)
  P3_MACRO <- proj_qr(C_norm)
  P4_MACRO <- diag(3)  
  
  X_TEST_MACRO <- ttm(ttm(ttm(X_tensor, P2_MACRO, 2),
                          P3_MACRO, 3),
                      P4_MACRO, 4)
  CRIT_MACRO <- (rTensor::fnorm(X_TEST_MACRO) / rTensor::fnorm(X_tensor))^2
  
  
  return(data.frame(
    # Valores de las carta T2
    alpha_T2_TROD_mcd = result_TROD$best_alpha,
    alpha_T2_MPCA_mcd = result_MPCA$best_alpha,
    alpha_T2_MACRO_mcd = result_MACRO$best_alpha,
    Crit_TROD=CRIT_TROD,
    Crit_MPCA=CRIT_MPCA,
    Crit_MACRO=CRIT_MACRO
  ))
}

# ============================================================================ #
#                    Funcion para hallar UCL con bootstrap                     #
# ============================================================================ #

Real_images    <- function(i, X_IMAGES, n_samples = 209) {
  
  img_shape = c(100, 100, 3)
  
  idx_sample <- sample.int(209, n_samples, replace = TRUE)
  X_sample   <- X_IMAGES[idx_sample,,, , drop = FALSE] 
  
  X <- X_sample + array(rnorm(length(X_sample), mean = 0, sd = 0.01), dim = dim(X_sample))
  
  ##### TROD - CP decomposition
  InicioCP <- Sys.time()
  cp_result <- cp(as.tensor(X), num_components = 4)
  FinalCP <- Sys.time()
  TimeCP <- FinalCP - InicioCP
  
  normU1 <- apply(cp_result$U[[1]], 2, vecnorm)
  normU2 <- apply(cp_result$U[[2]], 2, vecnorm)
  normU3 <- apply(cp_result$U[[3]], 2, vecnorm)
  normU4 <- apply(cp_result$U[[4]], 2, vecnorm)
  
  U1CP <- sweep(cp_result$U[[1]], 2, normU1, "/")
  U2CP <- sweep(cp_result$U[[2]], 2, normU2, "/")
  U3CP <- sweep(cp_result$U[[3]], 2, normU3, "/")
  U4CP <- sweep(cp_result$U[[4]], 2, normU4, "/")
  
  UpdateLambdasCP <- cp_result$lambdas * normU1 * normU2 * normU3 * normU4
  
  lambdas_trod <- matrix(0, nrow = n_samples, ncol = 4)
  for (ii in 1:n_samples) {
    for (j in 1:4) {
      lambdas_trod[ii, j] <- UpdateLambdasCP[j] * U1CP[ii, j]
    }
  }
  
  ##### MPCA
  Xmpca <- aperm(X, c(4, 2, 3, 1))
  InicioMPCA <- Sys.time()
  MPCA_result <- mpca(as.tensor(Xmpca), ranks = c(1,4,4))
  FinalMPCA <- Sys.time()
  TimeMPCA <- FinalMPCA - InicioMPCA
  
  ##### MACRO-PARAFAC
  InicioMACRO <- Sys.time()
  MACRO_result <-try(MacroPARAFAC2(X, ncomp = 4))
  FinalMACRO <- Sys.time()
  TimeMACRO <- FinalMACRO - InicioMACRO
  
  if (inherits(MACRO_result, "try-error")) {
    MACRO_result <- NULL 
  }
  
  if (!is.null(MACRO_result)) {
    normaA <- apply(MACRO_result$A, 2, vecnorm)
    normaB <- apply(MACRO_result$B, 2, vecnorm)
    normaC <- apply(MACRO_result$C, 2, vecnorm)
    normaD <- apply(MACRO_result$D, 2, vecnorm)
    
    A_norm <- sweep(MACRO_result$A, 2, normaA, "/")
    B_norm <- sweep(MACRO_result$B, 2, normaB, "/")
    C_norm <- sweep(MACRO_result$C, 2, normaC, "/")
    D_norm <- sweep(MACRO_result$D, 2, normaD, "/")
    
    UpdateLambdasMACRO <- normaA * normaB * normaC * normaD
    
    lambdas_MACRO <- matrix(0, nrow = n_samples, ncol = 4)
    for (ii in 1:n_samples) {
      for (j in 1:4) {
        lambdas_MACRO[ii, j] <- UpdateLambdasMACRO[j] * A_norm[ii, j]
      }
    }
  }
  
  ##### Inicialización para reconstrucciones
  num_elements <- prod(img_shape)
  theta_TROD <- matrix(0, 4, n_samples)
  theta_MPCA <- matrix(0, prod(dim(MPCA_result$Z_ext)[-4]), n_samples)
  theta_MACRO <- matrix(0, 4, n_samples)
  e_TROD <- matrix(0, num_elements, n_samples)
  e_MPCA <- matrix(0, num_elements, n_samples)
  e_MACRO <- matrix(0, num_elements, n_samples)
  
  for (ii in 1:n_samples) {
    y <- X[ii,,,]
    y_tensor <- as.tensor(y)
    
    # TROD
    z_trod <- recon_trod_tensor(lambdas_trod[ii, ], list(U1CP, U2CP, U3CP, U4CP))
    theta_TROD[, ii] <- lambdas_trod[ii, ]
    e_TROD[, ii] <- as.vector(y) - as.vector(z_trod)
    
    # MPCA
    z_mpca <- ttm(ttm(ttm(y_tensor, t(MPCA_result$U[[2]]), m = 1),
                      t(MPCA_result$U[[3]]), m = 2),
                  t(MPCA_result$U[[1]]), m = 3)
    
    y_recon_mpca <- ttm(ttm(ttm(z_mpca, MPCA_result$U[[2]], m = 1),
                            MPCA_result$U[[3]], m = 2),
                        MPCA_result$U[[1]], m = 3)
    
    theta_MPCA[, ii] <- as.vector(z_mpca@data)
    e_MPCA[, ii] <- as.vector(y) - as.vector(y_recon_mpca@data)
    
    # MACRO
    if (!is.null(MACRO_result)) {z_macro <- recon_trod_tensor(lambdas_MACRO[ii, ], list(A_norm, B_norm, C_norm, D_norm))
    theta_MACRO[, ii] <- lambdas_MACRO[ii, ]
    e_MACRO[, ii] <- as.vector(y) - as.vector(z_macro)
    }
  }
  
  theta_TROD <- t(theta_TROD)
  
  theta_MPCA <- t(theta_MPCA)
  
  if (!is.null(MACRO_result)) {
    theta_MACRO <- t(theta_MACRO)
  }
  
  # T2 clásico ################################################################
  
  # Para TROD ----------------------------------------------------------------
  Inicio_TROD_T2 <- Sys.time()
  T2_TROD <- tryCatch(
    apply(theta_TROD, 1, function(x) mahalanobis(x, colMeans(theta_TROD), cov(theta_TROD))),
    error = function(e) rep(NA_real_, nrow(theta_TROD))
  )
  Final_TROD_T2 <- Sys.time()
  Time_TROD_T2 <- Final_TROD_T2 - Inicio_TROD_T2
  
  # Para MPCA -----------------------------------------------------------------
  Inicio_MPCA_T2 <- Sys.time()
  T2_MPCA <- tryCatch(
    apply(theta_MPCA, 1, function(x) mahalanobis(x, colMeans(theta_MPCA), cov(theta_MPCA))),
    error = function(e) rep(NA_real_, nrow(theta_MPCA))
  )
  Final_MPCA_T2 <- Sys.time()
  Time_MPCA_T2 <- Final_MPCA_T2 - Inicio_MPCA_T2
  
  # Para MACROparafac --------------------------------------------------------  
  Inicio_MACRO_T2 <- Sys.time()
  if (!is.null(MACRO_result)) {
    T2_MACRO <- tryCatch(
      apply(theta_MACRO, 1, function(x) mahalanobis(x, colMeans(theta_MACRO), cov(theta_MACRO))),
      error = function(e) rep(NA_real_, nrow(theta_MACRO))
    )
  } else {
    T2_MACRO <- NA
  }
  Final_MACRO_T2 <- Sys.time()
  Time_MACRO_T2 <- Final_MACRO_T2 - Inicio_MACRO_T2
  
  # T2 robustos (MCD) ###################################################
  
  # Para TROD-----------------------------------------------------------------
  Inicio_TROD_MCD_T2 <- Sys.time()
  trod_mcd <- tryCatch(CovMcd(theta_TROD, alpha = 0.7679500),
                       error = function(e) NULL)
  T2_TROD_MCD <- tryCatch(
    if (is.null(trod_mcd)) rep(NA_real_, nrow(theta_TROD)) else
      mahalanobis(theta_TROD, trod_mcd$center, trod_mcd$cov),
    error = function(e) rep(NA_real_, nrow(theta_TROD))
  )
  Final_TROD_MCD_T2 <- Sys.time()
  Time_TROD_MCD_T2 <- Final_TROD_MCD_T2 - Inicio_TROD_MCD_T2
  
  # 80 %-----------------------------------------------------------------
  Inicio_TROD_MCD_T2_80 <- Sys.time()
  trod_mcd_80 <- tryCatch(CovMcd(theta_TROD, alpha = 0.8),
                       error = function(e) NULL)
  T2_TROD_MCD_80 <- tryCatch(
    if (is.null(trod_mcd_80)) rep(NA_real_, nrow(theta_TROD)) else
      mahalanobis(theta_TROD, trod_mcd_80$center, trod_mcd_80$cov),
    error = function(e) rep(NA_real_, nrow(theta_TROD))
  )
  Final_TROD_MCD_T2_80 <- Sys.time()
  Time_TROD_MCD_T2_80 <- Final_TROD_MCD_T2_80 - Inicio_TROD_MCD_T2_80
  
  # Para MPCA-----------------------------------------------------------------
  Inicio_MPCA_MCD_T2 <- Sys.time()
  MPCA_mcd <- tryCatch(CovMcd(theta_MPCA, alpha = 0.7357250),
                       error = function(e) NULL)
  T2_MPCA_MCD <- tryCatch(
    if (is.null(MPCA_mcd)) rep(NA_real_, nrow(theta_MPCA)) else
      mahalanobis(theta_MPCA, MPCA_mcd$center, MPCA_mcd$cov),
    error = function(e) rep(NA_real_, nrow(theta_MPCA))
  )
  Final_MPCA_MCD_T2 <- Sys.time()
  Time_MPCA_MCD_T2 <- Final_MPCA_MCD_T2 - Inicio_MPCA_MCD_T2
  
  # MCD 80 ----------------------------------------------------------------------
  Inicio_MPCA_MCD_T2_80 <- Sys.time()
  MPCA_mcd_80 <- tryCatch(CovMcd(theta_MPCA, alpha = 0.8),
                          error = function(e) NULL)
  T2_MPCA_MCD_80 <- tryCatch(
    if (is.null(MPCA_mcd_80)) rep(NA_real_, nrow(theta_MPCA)) else
      mahalanobis(theta_MPCA, MPCA_mcd_80$center, MPCA_mcd_80$cov),
    error = function(e) rep(NA_real_, nrow(theta_MPCA))
  )
  Final_MPCA_MCD_T2_80 <- Sys.time()
  Time_MPCA_MCD_T2_80 <- Final_MPCA_MCD_T2_80 - Inicio_MPCA_MCD_T2_80
  
  # Para MACROparafac----------------------------------------------------------
  Inicio_MACRO_MCD_T2 <- Sys.time()
  MACRO_mcd <- tryCatch(CovMcd(theta_MACRO, alpha = 0.8378250),
                        error = function(e) NULL)
  T2_MACRO_MCD <- tryCatch(
    if (is.null(MACRO_mcd)) rep(NA_real_, nrow(theta_MACRO)) else
      mahalanobis(theta_MACRO, MACRO_mcd$center, MACRO_mcd$cov),
    error = function(e) rep(NA_real_, nrow(theta_MACRO))
  )
  Final_MACRO_MCD_T2 <- Sys.time()
  Time_MACRO_MCD_T2 <- Final_MACRO_MCD_T2 - Inicio_MACRO_MCD_T2
  
  # MCD 80 ----------------------------------------------------------------------
  Inicio_MACRO_MCD_T2_80 <- Sys.time()
  MACRO_mcd_80 <- tryCatch(CovMcd(theta_MACRO, alpha = 0.8),
                           error = function(e) NULL)
  
  T2_MACRO_MCD_80 <- tryCatch(
    if (is.null(MACRO_mcd_80)) rep(NA_real_, nrow(theta_MACRO)) else
      mahalanobis(theta_MACRO, MACRO_mcd_80$center, MACRO_mcd_80$cov),
    error = function(e) rep(NA_real_, nrow(theta_MACRO))
  )
  Final_MACRO_MCD_T2_80 <- Sys.time()
  Time_MACRO_MCD_T2_80 <- Final_MACRO_MCD_T2_80 - Inicio_MACRO_MCD_T2_80
  
  # Q ##########################################################################
  
  Inicio_TROD_Q <- Sys.time()
  Q_TROD <- colSums(e_TROD^2)
  Final_TROD_Q <- Sys.time()
  Time_TROD_Q <- Final_TROD_Q - Inicio_TROD_Q
  
  Inicio_MPCA_Q <- Sys.time()
  Q_MPCA <- colSums(e_MPCA^2)
  Final_MPCA_Q <- Sys.time()
  Time_MPCA_Q <- Final_MPCA_Q - Inicio_MPCA_Q
  
  Inicio_MACRO_Q <- Sys.time()
  Q_MACRO <- if (!is.null(MACRO_result)) colSums(e_MACRO^2) else NA
  Final_MACRO_Q <- Sys.time()
  Time_MACRO_Q <- Final_MACRO_Q - Inicio_MACRO_Q
  
  #cat("Fin de la iteración:", i, "- delta:", delta_sim,"\n")
  cat("Fin de la iteración:", i,"\n")
  return(data.frame(
    iter = i,
    # Tiempos de las descomposiciones
    time_TROD_des = as.numeric(TimeCP, units = "secs"),
    time_MPCA_des = as.numeric(TimeMPCA, units = "secs"),
    time_MACRO_des = as.numeric(TimeMACRO, units = "secs"),
    
    # Valores de las carta T2
    max_T2_TROD = max(T2_TROD),
    max_T2_MPCA = max(T2_MPCA),
    max_T2_MACRO = if (!is.null(T2_MACRO)) max(T2_MACRO) else NA,
    
    # Valores de las carta T2 Robustas
    max_T2_TROD_MCD = max(T2_TROD_MCD),
    max_T2_TROD_MCD_80 = max(T2_TROD_MCD_80),
    max_T2_MPCA_MCD = max(T2_MPCA_MCD),
    max_T2_MPCA_MCD_80 = max(T2_MPCA_MCD_80),
    max_T2_MACRO_MCD = max(T2_MACRO_MCD),
    max_T2_MACRO_MCD_80 = max(T2_MACRO_MCD_80),
    
    # Valores de las cartas Q
    max_Q_TROD = max(Q_TROD),
    max_Q_MPCA = max(Q_MPCA),
    max_Q_MACRO = if (!is.null(Q_MACRO)) max(Q_MACRO) else NA,
    
    # Tiempos de las cartas 
    time_TROD =as.numeric(TimeCP+ Time_TROD_T2+Time_TROD_Q, units = "secs"),
    time_TROD_MCD = as.numeric(TimeCP+Time_TROD_MCD_T2+Time_TROD_Q, units = "secs"),
    time_TROD_MCD_80 = as.numeric(TimeCP+Time_TROD_MCD_T2_80+Time_TROD_Q, units = "secs") ,
    
    time_MPCA = as.numeric(TimeMPCA+Time_MPCA_T2+Time_MPCA_Q, units = "secs") ,
    time_MPCA_MCD = as.numeric(TimeMPCA+Time_MPCA_MCD_T2+Time_MPCA_Q, units = "secs") ,
    time_MPCA_MCD_80 = as.numeric(TimeMPCA+Time_MPCA_MCD_T2_80+Time_MPCA_Q, units = "secs") ,
    
    time_MACRO =as.numeric(TimeMACRO+Time_MACRO_T2+Time_MACRO_Q, units = "secs") ,
    time_MACRO_MCD =as.numeric(TimeMACRO+Time_MACRO_MCD_T2+Time_MACRO_Q, units = "secs") ,
    time_MACRO_MCD_80 = as.numeric(TimeMACRO+Time_MACRO_MCD_T2_80+Time_MACRO_Q, units = "secs")
  ))
}


# ============================================================================ #
#                           Desplazar las imágenes                             #
# ============================================================================ #

shift_image <- function(img, up = 3, left = 3) {
  h <- dim(img)[1]
  w <- dim(img)[2]
  c <- dim(img)[3]
  
  # Creamos la imagen desplazada con fondo blanco (255)
  out <- array(255, dim = c(h, w, c))
  
  # Copiamos los píxeles desplazados
  out[1:(h - up), 1:(w - left), ] <- img[(1 + up):h, (1 + left):w, ]
  return(out)
}

# ============================================================================ #
#                         Reemplazar con desplazadas                           #
# ============================================================================ #

replace_last_n_with_shifted <- function(X_base, X_shifted, n_replace) {
  X_modified <- X_base
  total <- dim(X_base)[1]
  idx_replace <- (total - n_replace + 1):total
  X_modified[idx_replace,,,] <- X_shifted[idx_replace,,,]
  return(X_modified)
}


# ============================================================================ #
#                    Poner a prueba con imágenes reales                        #
# ============================================================================ #

Real_test <- function(X_IMAGES=X_defectuosas) {
  
  img_shape = c(100, 100, 3)
  
  X <- X_IMAGES + array(rnorm(length(X_IMAGES), mean = 0, sd = 0.01), dim = dim(X_IMAGES))
  
  ##### TROD - CP decomposition
  InicioCP <- Sys.time()
  cp_result <- cp(as.tensor(X), num_components = 4)
  FinalCP <- Sys.time()
  TimeCP <- FinalCP - InicioCP
  
  normU1 <- apply(cp_result$U[[1]], 2, vecnorm)
  normU2 <- apply(cp_result$U[[2]], 2, vecnorm)
  normU3 <- apply(cp_result$U[[3]], 2, vecnorm)
  normU4 <- apply(cp_result$U[[4]], 2, vecnorm)
  
  U1CP <- sweep(cp_result$U[[1]], 2, normU1, "/")
  U2CP <- sweep(cp_result$U[[2]], 2, normU2, "/")
  U3CP <- sweep(cp_result$U[[3]], 2, normU3, "/")
  U4CP <- sweep(cp_result$U[[4]], 2, normU4, "/")
  
  UpdateLambdasCP <- cp_result$lambdas * normU1 * normU2 * normU3 * normU4
  
  lambdas_trod <- matrix(0, nrow = n_samples, ncol = 4)
  for (ii in 1:n_samples) {
    for (j in 1:4) {
      lambdas_trod[ii, j] <- UpdateLambdasCP[j] * U1CP[ii, j]
    }
  }
  
  ##### MPCA
  Xmpca <- aperm(X, c(4, 2, 3, 1))
  InicioMPCA <- Sys.time()
  MPCA_result <- mpca(as.tensor(Xmpca), ranks = c(1,4,4))
  FinalMPCA <- Sys.time()
  TimeMPCA <- FinalMPCA - InicioMPCA
  
  ##### MACRO-PARAFAC
  InicioMACRO <- Sys.time()
  MACRO_result <-try(MacroPARAFAC2(X, ncomp = 4))
  FinalMACRO <- Sys.time()
  TimeMACRO <- FinalMACRO - InicioMACRO
  
  if (inherits(MACRO_result, "try-error")) {
    MACRO_result <- NULL 
  }
  
  if (!is.null(MACRO_result)) {
    normaA <- apply(MACRO_result$A, 2, vecnorm)
    normaB <- apply(MACRO_result$B, 2, vecnorm)
    normaC <- apply(MACRO_result$C, 2, vecnorm)
    normaD <- apply(MACRO_result$D, 2, vecnorm)
    
    A_norm <- sweep(MACRO_result$A, 2, normaA, "/")
    B_norm <- sweep(MACRO_result$B, 2, normaB, "/")
    C_norm <- sweep(MACRO_result$C, 2, normaC, "/")
    D_norm <- sweep(MACRO_result$D, 2, normaD, "/")
    
    UpdateLambdasMACRO <- normaA * normaB * normaC * normaD
    
    lambdas_MACRO <- matrix(0, nrow = n_samples, ncol = 4)
    for (ii in 1:n_samples) {
      for (j in 1:4) {
        lambdas_MACRO[ii, j] <- UpdateLambdasMACRO[j] * A_norm[ii, j]
      }
    }
  }
  
  ##### Inicialización para reconstrucciones
  num_elements <- prod(img_shape)
  theta_TROD <- matrix(0, 4, n_samples)
  theta_MPCA <- matrix(0, prod(dim(MPCA_result$Z_ext)[-4]), n_samples)
  theta_MACRO <- matrix(0, 4, n_samples)
  e_TROD <- matrix(0, num_elements, n_samples)
  e_MPCA <- matrix(0, num_elements, n_samples)
  e_MACRO <- matrix(0, num_elements, n_samples)
  
  for (ii in 1:n_samples) {
    y <- X[ii,,,]
    y_tensor <- as.tensor(y)
    
    # TROD
    z_trod <- recon_trod_tensor(lambdas_trod[ii, ], list(U1CP, U2CP, U3CP, U4CP))
    theta_TROD[, ii] <- lambdas_trod[ii, ]
    e_TROD[, ii] <- as.vector(y) - as.vector(z_trod)
    
    # MPCA
    z_mpca <- ttm(ttm(ttm(y_tensor, t(MPCA_result$U[[2]]), m = 1),
                      t(MPCA_result$U[[3]]), m = 2),
                  t(MPCA_result$U[[1]]), m = 3)
    
    y_recon_mpca <- ttm(ttm(ttm(z_mpca, MPCA_result$U[[2]], m = 1),
                            MPCA_result$U[[3]], m = 2),
                        MPCA_result$U[[1]], m = 3)
    
    theta_MPCA[, ii] <- as.vector(z_mpca@data)
    e_MPCA[, ii] <- as.vector(y) - as.vector(y_recon_mpca@data)
    
    # MACRO
    if (!is.null(MACRO_result)) {z_macro <- recon_trod_tensor(lambdas_MACRO[ii, ], list(A_norm, B_norm, C_norm, D_norm))
    theta_MACRO[, ii] <- lambdas_MACRO[ii, ]
    e_MACRO[, ii] <- as.vector(y) - as.vector(z_macro)
    }
  }
  
  theta_TROD <- t(theta_TROD)
  
  theta_MPCA <- t(theta_MPCA)
  
  if (!is.null(MACRO_result)) {
    theta_MACRO <- t(theta_MACRO)
  }
  
  # T2 clásico ################################################################
  
  # Para TROD ----------------------------------------------------------------
  Inicio_TROD_T2 <- Sys.time()
  T2_TROD <- tryCatch(
    apply(theta_TROD, 1, function(x) mahalanobis(x, colMeans(theta_TROD), cov(theta_TROD))),
    error = function(e) rep(NA_real_, nrow(theta_TROD))
  )
  Final_TROD_T2 <- Sys.time()
  Time_TROD_T2 <- Final_TROD_T2 - Inicio_TROD_T2
  
  # Para MPCA -----------------------------------------------------------------
  Inicio_MPCA_T2 <- Sys.time()
  T2_MPCA <- tryCatch(
    apply(theta_MPCA, 1, function(x) mahalanobis(x, colMeans(theta_MPCA), cov(theta_MPCA))),
    error = function(e) rep(NA_real_, nrow(theta_MPCA))
  )
  Final_MPCA_T2 <- Sys.time()
  Time_MPCA_T2 <- Final_MPCA_T2 - Inicio_MPCA_T2
  
  # Para MACROparafac --------------------------------------------------------  
  Inicio_MACRO_T2 <- Sys.time()
  if (!is.null(MACRO_result)) {
    T2_MACRO <- tryCatch(
      apply(theta_MACRO, 1, function(x) mahalanobis(x, colMeans(theta_MACRO), cov(theta_MACRO))),
      error = function(e) rep(NA_real_, nrow(theta_MACRO))
    )
  } else {
    T2_MACRO <- NA
  }
  Final_MACRO_T2 <- Sys.time()
  Time_MACRO_T2 <- Final_MACRO_T2 - Inicio_MACRO_T2
  
  # T2 robustos (MCD) ###################################################
  
  # Para TROD-----------------------------------------------------------------
  Inicio_TROD_MCD_T2 <- Sys.time()
  trod_mcd <- tryCatch(CovMcd(theta_TROD, alpha = 0.7679500),
                       error = function(e) NULL)
  T2_TROD_MCD <- tryCatch(
    if (is.null(trod_mcd)) rep(NA_real_, nrow(theta_TROD)) else
      mahalanobis(theta_TROD, trod_mcd$center, trod_mcd$cov),
    error = function(e) rep(NA_real_, nrow(theta_TROD))
  )
  Final_TROD_MCD_T2 <- Sys.time()
  Time_TROD_MCD_T2 <- Final_TROD_MCD_T2 - Inicio_TROD_MCD_T2
  
  # 80 %-----------------------------------------------------------------
  Inicio_TROD_MCD_T2_80 <- Sys.time()
  trod_mcd_80 <- tryCatch(CovMcd(theta_TROD, alpha = 0.8),
                          error = function(e) NULL)
  T2_TROD_MCD_80 <- tryCatch(
    if (is.null(trod_mcd_80)) rep(NA_real_, nrow(theta_TROD)) else
      mahalanobis(theta_TROD, trod_mcd_80$center, trod_mcd_80$cov),
    error = function(e) rep(NA_real_, nrow(theta_TROD))
  )
  Final_TROD_MCD_T2_80 <- Sys.time()
  Time_TROD_MCD_T2_80 <- Final_TROD_MCD_T2_80 - Inicio_TROD_MCD_T2_80
  
  # Para MPCA-----------------------------------------------------------------
  Inicio_MPCA_MCD_T2 <- Sys.time()
  MPCA_mcd <- tryCatch(CovMcd(theta_MPCA, alpha = 0.7357250),
                       error = function(e) NULL)
  T2_MPCA_MCD <- tryCatch(
    if (is.null(MPCA_mcd)) rep(NA_real_, nrow(theta_MPCA)) else
      mahalanobis(theta_MPCA, MPCA_mcd$center, MPCA_mcd$cov),
    error = function(e) rep(NA_real_, nrow(theta_MPCA))
  )
  Final_MPCA_MCD_T2 <- Sys.time()
  Time_MPCA_MCD_T2 <- Final_MPCA_MCD_T2 - Inicio_MPCA_MCD_T2
  
  # MCD 80 ----------------------------------------------------------------------
  Inicio_MPCA_MCD_T2_80 <- Sys.time()
  MPCA_mcd_80 <- tryCatch(CovMcd(theta_MPCA, alpha = 0.8),
                          error = function(e) NULL)
  T2_MPCA_MCD_80 <- tryCatch(
    if (is.null(MPCA_mcd_80)) rep(NA_real_, nrow(theta_MPCA)) else
      mahalanobis(theta_MPCA, MPCA_mcd_80$center, MPCA_mcd_80$cov),
    error = function(e) rep(NA_real_, nrow(theta_MPCA))
  )
  Final_MPCA_MCD_T2_80 <- Sys.time()
  Time_MPCA_MCD_T2_80 <- Final_MPCA_MCD_T2_80 - Inicio_MPCA_MCD_T2_80
  
  # Para MACROparafac----------------------------------------------------------
  Inicio_MACRO_MCD_T2 <- Sys.time()
  MACRO_mcd <- tryCatch(CovMcd(theta_MACRO, alpha = 0.8378250),
                        error = function(e) NULL)
  T2_MACRO_MCD <- tryCatch(
    if (is.null(MACRO_mcd)) rep(NA_real_, nrow(theta_MACRO)) else
      mahalanobis(theta_MACRO, MACRO_mcd$center, MACRO_mcd$cov),
    error = function(e) rep(NA_real_, nrow(theta_MACRO))
  )
  Final_MACRO_MCD_T2 <- Sys.time()
  Time_MACRO_MCD_T2 <- Final_MACRO_MCD_T2 - Inicio_MACRO_MCD_T2
  
  # MCD 80 ----------------------------------------------------------------------
  Inicio_MACRO_MCD_T2_80 <- Sys.time()
  MACRO_mcd_80 <- tryCatch(CovMcd(theta_MACRO, alpha = 0.8),
                           error = function(e) NULL)
  
  T2_MACRO_MCD_80 <- tryCatch(
    if (is.null(MACRO_mcd_80)) rep(NA_real_, nrow(theta_MACRO)) else
      mahalanobis(theta_MACRO, MACRO_mcd_80$center, MACRO_mcd_80$cov),
    error = function(e) rep(NA_real_, nrow(theta_MACRO))
  )
  Final_MACRO_MCD_T2_80 <- Sys.time()
  Time_MACRO_MCD_T2_80 <- Final_MACRO_MCD_T2_80 - Inicio_MACRO_MCD_T2_80
  
  # Q ##########################################################################
  
  Inicio_TROD_Q <- Sys.time()
  Q_TROD <- colSums(e_TROD^2)
  Final_TROD_Q <- Sys.time()
  Time_TROD_Q <- Final_TROD_Q - Inicio_TROD_Q
  
  Inicio_MPCA_Q <- Sys.time()
  Q_MPCA <- colSums(e_MPCA^2)
  Final_MPCA_Q <- Sys.time()
  Time_MPCA_Q <- Final_MPCA_Q - Inicio_MPCA_Q
  
  Inicio_MACRO_Q <- Sys.time()
  Q_MACRO <- if (!is.null(MACRO_result)) colSums(e_MACRO^2) else NA
  Final_MACRO_Q <- Sys.time()
  Time_MACRO_Q <- Final_MACRO_Q - Inicio_MACRO_Q
  
  return(data.frame(
    Image = 1:209,
    # Valores de las carta T2
    T2_TROD = (T2_TROD),
    T2_MPCA = (T2_MPCA),
    T2_MACRO = (T2_MACRO),
    
    # Valores de las carta T2 Robustas
    T2_TROD_MCD = (T2_TROD_MCD),
    T2_TROD_MCD_80 = (T2_TROD_MCD_80),
    T2_MPCA_MCD = (T2_MPCA_MCD),
    T2_MPCA_MCD_80 = (T2_MPCA_MCD_80),
    T2_MACRO_MCD = (T2_MACRO_MCD),
    T2_MACRO_MCD_80 = (T2_MACRO_MCD_80),
    
    # Valores de las cartas Q
    Q_TROD = (Q_TROD),
    Q_MPCA = (Q_MPCA),
    Q_MACRO =(Q_MACRO)
  ))
}

# ============================================================================ #
#                           Reemplazar cada k img                              #
# ============================================================================ #


replace_at_indices_with_broken <- function(X_base, X_broken, idx_replace) {
  X_mod <- X_base
  total <- dim(X_base)[1]
  n_broken <- dim(X_broken)[1]

  k <- length(idx_replace)
  idx_broken <- 1:k
  
  X_mod[idx_replace,,,] <- X_broken[idx_broken,,,]
  return(X_mod)
}

