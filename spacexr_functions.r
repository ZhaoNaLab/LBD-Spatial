run.RCTD.replicates.local <- function (RCTD.replicates, doublet_mode = "doublet") 
{
  if (!(doublet_mode %in% c("doublet", "multi", "full"))) 
    stop(paste0("run.RCTD.replicates: doublet_mode=", doublet_mode, 
                " is not a valid choice. Please set doublet_mode=doublet, multi, or full."))
  for (i in 1:length(RCTD.replicates@RCTD.reps)) {
    message(paste("run.RCTD.replicates: running RCTD for replicate", 
                  i))
    RCTD.replicates@RCTD.reps[[i]] <- run.RCTD.local(RCTD.replicates@RCTD.reps[[i]], 
                                               doublet_mode = doublet_mode)
  }
  return(RCTD.replicates)
}

run.RCTD.local <- function (RCTD, doublet_mode = "doublet") 
{
  if (!(doublet_mode %in% c("doublet", "multi", "full"))) 
    stop(paste0("run.RCTD: doublet_mode=", doublet_mode, 
                " is not a valid choice. Please set doublet_mode=doublet, multi, or full."))
  RCTD@config$RCTDmode <- doublet_mode
  RCTD <- fitBulk(RCTD)
  RCTD <- choose_sigma_c.local(RCTD)
  RCTD <- fitPixels.local(RCTD, doublet_mode = doublet_mode)
}

choose_sigma_c.local <- function (RCTD) 
{
  puck = RCTD@spatialRNA
  MIN_UMI = RCTD@config$UMI_min_sigma
  sigma = 100
  Q1 <- readRDS(system.file("extdata", "Qmat/Q_mat_1.rds", 
                            package = "spacexr"))
  Q2 <- readRDS(system.file("extdata", "Qmat/Q_mat_2.rds", 
                            package = "spacexr"))
  Q3 <- readRDS(system.file("extdata", "Qmat/Q_mat_3.rds", 
                            package = "spacexr"))
  Q4 <- readRDS(system.file("extdata", "Qmat/Q_mat_4.rds", 
                            package = "spacexr"))
  Q5 <- readRDS(system.file("extdata", "Qmat/Q_mat_5.rds", 
                            package = "spacexr"))
  Q_mat_all <- c(Q1, Q2, Q3, Q4, Q5)
  sigma_vals <- names(Q_mat_all)
  X_vals <- readRDS(system.file("extdata", "Qmat/X_vals.rds", 
                                package = "spacexr"))
  N_fit = min(RCTD@config$N_fit, sum(puck@nUMI > MIN_UMI))
  if (N_fit == 0) {
    stop(paste("choose_sigma_c determined a N_fit of 0! This is probably due to unusually low UMI counts per bead in your dataset. Try decreasing the parameter UMI_min_sigma. It currently is", 
               MIN_UMI, "but none of the beads had counts larger than that."))
  }
  fit_ind = sample(names(puck@nUMI[puck@nUMI > MIN_UMI]), N_fit)
  beads = t(as.matrix(puck@counts[RCTD@internal_vars$gene_list_reg, 
                                  fit_ind]))
  message(paste("chooseSigma: using initial Q_mat with sigma = ", 
                sigma/100))
  for (iter in 1:RCTD@config$N_epoch) {
    set_likelihood_vars(Q_mat_all[[as.character(sigma)]], 
                        X_vals)
    results = decompose_batch.local(puck@nUMI[fit_ind], RCTD@cell_type_info$renorm[[1]], 
                              beads, RCTD@internal_vars$gene_list_reg, constrain = F, 
                              max_cores = RCTD@config$max_cores)
    weights = Matrix(0, nrow = N_fit, ncol = RCTD@cell_type_info$renorm[[3]])
    rownames(weights) = fit_ind
    colnames(weights) = RCTD@cell_type_info$renorm[[2]]
    for (i in 1:N_fit) weights[i, ] = results[[i]]$weights
    prediction <- sweep(as.matrix(RCTD@cell_type_info$renorm[[1]][RCTD@internal_vars$gene_list_reg, 
    ]) %*% t(as.matrix(weights)), 2, puck@nUMI[fit_ind], 
    "*")
    message(paste("Likelihood value:", calc_log_l_vec(as.vector(prediction), 
                                                      as.vector(t(beads)))))
    sigma_prev <- sigma
    sigma <- chooseSigma(prediction, t(beads), Q_mat_all, 
                         X_vals, sigma)
    message(paste("Sigma value: ", sigma/100))
    if (sigma == sigma_prev) 
      break
  }
  RCTD@internal_vars$sigma <- sigma/100
  RCTD@internal_vars$Q_mat <- Q_mat_all[[as.character(sigma)]]
  RCTD@internal_vars$X_vals <- X_vals
  return(RCTD)
}

fitPixels.local <- function(RCTD, doublet_mode = "doublet") {
  RCTD@internal_vars$cell_types_assigned <- TRUE
  RCTD@config$RCTDmode <- doublet_mode
  set_likelihood_vars(RCTD@internal_vars$Q_mat, RCTD@internal_vars$X_vals)
  cell_type_info <- RCTD@cell_type_info$renorm
  if(doublet_mode == "doublet") {
    results = process_beads_batch(cell_type_info, RCTD@internal_vars$gene_list_reg, RCTD@spatialRNA, class_df = RCTD@internal_vars$class_df,
                                  constrain = F, MAX_CORES = RCTD@config$max_cores, MIN.CHANGE = RCTD@config$MIN_CHANGE_REG,
                                  CONFIDENCE_THRESHOLD = RCTD@config$CONFIDENCE_THRESHOLD, DOUBLET_THRESHOLD = RCTD@config$DOUBLET_THRESHOLD)
    return(gather_results(RCTD, results))
  } else if(doublet_mode == "full") {
    beads = t(as.matrix(RCTD@spatialRNA@counts[RCTD@internal_vars$gene_list_reg,]))
    results = decompose_batch.local(RCTD@spatialRNA@nUMI, cell_type_info[[1]], beads, RCTD@internal_vars$gene_list_reg, constrain = F,
                              max_cores = RCTD@config$max_cores, MIN.CHANGE = RCTD@config$MIN_CHANGE_REG)
    weights = Matrix(0, nrow = length(results), ncol = RCTD@cell_type_info$renorm[[3]])
    rownames(weights) = colnames(RCTD@spatialRNA@counts); colnames(weights) = RCTD@cell_type_info$renorm[[2]];
    for(i in 1:dim(weights)[1])
      weights[i,] = results[[i]]$weights
    RCTD@results <- list(weights = weights)
    return(RCTD)
  } else if(doublet_mode == "multi") {
    RCTD@results = process_beads_multi(cell_type_info, RCTD@internal_vars$gene_list_reg, RCTD@spatialRNA, class_df = RCTD@internal_vars$class_df,
                                       constrain = F, MAX_CORES = RCTD@config$max_cores,
                                       MIN.CHANGE = RCTD@config$MIN_CHANGE_REG, MAX.TYPES = RCTD@config$MAX_MULTI_TYPES,
                                       CONFIDENCE_THRESHOLD = RCTD@config$CONFIDENCE_THRESHOLD, DOUBLET_THRESHOLD = RCTD@config$DOUBLET_THRESHOLD)
    return(RCTD)
  } else {
    stop(paste0("fitPixels: doublet_mode=",doublet_mode, " is not a valid choice. Please set doublet_mode=doublet, multi, or full."))
  }
}

decompose_batch.local <- function(nUMI, cell_type_means, beads, gene_list, constrain = T, OLS = F, max_cores = 8, MIN.CHANGE = 0.001) {
  #out_file = "logs/decompose_batch_log.txt"
  #if (file.exists(out_file))
  #  file.remove(out_file)
  if(max_cores > 1) {
    numCores = parallel::detectCores()
    if(parallel::detectCores() > max_cores)
      numCores <- max_cores
    cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="")
    doParallel::registerDoParallel(cl)
    environ = c('decompose_full','solveIRWLS.weights',
                'solveOLS','solveWLS', 'Q_mat', 'K_val','X_vals', 'SQ_mat')
    #for(i in 1:100) {
    weights <- foreach::foreach(i = 1:(dim(beads)[1]), .packages = c("quadprog"), .export = environ) %dopar% {
      #if(i %% 100 == 0)
      #  cat(paste0("Finished sample: ",i,"\n"), file=out_file, append=TRUE)
      source("/research/bsi/projects/PI/tertiary/Kostallari_Enis_m139122/s307685.spatial_mouse_liver/processing/scripts/spacexr_functions.r")
      assign("Q_mat",Q_mat, envir = globalenv()); assign("X_vals",X_vals, envir = globalenv())
      assign("K_val",K_val, envir = globalenv()); assign("SQ_mat",SQ_mat, envir = globalenv());
      decompose_full(data.matrix(cell_type_means[gene_list,]*nUMI[i]), nUMI[i], beads[i,], constrain = constrain, OLS = OLS, MIN_CHANGE = MIN.CHANGE)
    }
    parallel::stopCluster(cl)
  } else {
    weights <- list()
    for(i in 1:(dim(beads)[1])) {
      weights[[i]] <- decompose_full(data.matrix(cell_type_means[gene_list,]*nUMI[i]), nUMI[i], beads[i,], constrain = constrain, OLS = OLS, MIN_CHANGE = MIN.CHANGE)
    }
  }
  return(weights)
}

decompose_full <- function(cell_type_profiles, nUMI, bead, constrain = TRUE, OLS = FALSE, verbose = F, n.iter = 50, MIN_CHANGE = 0.001, bulk_mode = F) {
  results = solveIRWLS.weights(cell_type_profiles,bead,nUMI,OLS = OLS, constrain = constrain,
                               verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, bulk_mode = bulk_mode)
  return(results)
}

solveOLS<-function(S,B, solution, constrain = T){
  D<-t(S)%*%S
  d<-t(S)%*%B
  norm_factor <- norm(D,"2")
  D <- D / norm_factor
  d <- d / norm_factor
  epsilon <- 1e-7; D <- D + epsilon * diag(length(d))
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  if(constrain) {
    A_const = t(rbind(1,A))
    b_const <-c(1 - sum(solution),bzero)
    solution <- quadprog::solve.QP(D,d,A_const,b_const,meq=1)$solution
  } else {
    solution <- quadprog::solve.QP(D,d,A,bzero,meq=0)$solution
  }
  names(solution)<-colnames(S)
  return(solution)
}

calc_log_l_vec <- function(lambda, Y, return_vec = FALSE) {
  log_l_vec <- -calc_Q_all(Y,lambda)$d0_vec
  if(return_vec)
    return(log_l_vec)
  return(sum(log_l_vec))
}


chooseSigma <- function(prediction, counts, Q_mat_all, X_vals, sigma) {
  X = as.vector(prediction); X = pmax(X, 1e-4)
  Y = as.vector(counts); num_sample = min(1000000, length(X)) #300000
  use_ind = sample(1:length(X), num_sample)
  X = X[use_ind]; Y = Y[use_ind]
  mult_fac_vec <- (8:12)/10; sigma_ind <- c(10:70, (36:100)*2)
  si <- which(sigma_ind == round(sigma))
  sigma_ind <- sigma_ind[max(1,si - 8):min(si+8,length(sigma_ind))]
  score_vec <- numeric(length(sigma_ind))
  for(i in 1:length(sigma_ind)) {
    sigma <- sigma_ind[i]
    set_likelihood_vars(Q_mat_all[[as.character(sigma)]],X_vals)
    best_val <- calc_log_l_vec(X*mult_fac_vec[1], Y)
    for(mult_fac in mult_fac_vec[2:length(mult_fac_vec)])
      best_val <- min(best_val, calc_log_l_vec(X*mult_fac, Y))
    score_vec[i] <- best_val
    #score_vec[i] <- calc_log_l_vec(X, Y)
  }
  sigma = sigma_ind[which.min(score_vec)]
  return(sigma)
}


process_bead_doublet <- function(cell_type_info, gene_list, UMI_tot, bead, class_df = NULL, constrain = T, verbose = F,
                                 MIN.CHANGE = 0.001, CONFIDENCE_THRESHOLD = 10, DOUBLET_THRESHOLD = 25) {
  cell_type_profiles <- cell_type_info[[1]][gene_list,]
  cell_type_profiles = cell_type_profiles * UMI_tot
  cell_type_profiles = data.matrix(cell_type_profiles)
  QL_score_cutoff = CONFIDENCE_THRESHOLD; doublet_like_cutoff = DOUBLET_THRESHOLD
  results_all = decompose_full(cell_type_profiles, UMI_tot, bead, constrain = constrain, verbose = verbose, MIN_CHANGE = MIN.CHANGE)
  all_weights <- results_all$weights
  conv_all <- results_all$converged
  initial_weight_thresh = 0.01; cell_type_names = cell_type_info[[2]]
  candidates <- names(which(all_weights > initial_weight_thresh))
  if(length(candidates) == 0)
    candidates = cell_type_info[[2]][1:min(3,cell_type_info[[3]])]
  if(length(candidates) == 1)
    if(candidates[1] == cell_type_info[[2]][1])
      candidates = c(candidates, cell_type_info[[2]][2])
  else
    candidates = c(candidates, cell_type_info[[2]][1])
  score_mat = Matrix(0, nrow = length(candidates), ncol = length(candidates))
  rownames(score_mat) = candidates; colnames(score_mat) = candidates
  min_score = 0
  first_type = NULL; second_type = NULL
  first_class = F; second_class = F #indicates whether the first (resp second) refers to a class rather than a type
  for(i in 1:(length(candidates)-1)) {
    type1 = candidates[i]
    for(j in (i+1):length(candidates)) {
      type2 = candidates[j]
      score = decompose_sparse(cell_type_profiles, UMI_tot, bead, type1, type2, score_mode = T, constrain = constrain, verbose = verbose, MIN.CHANGE = MIN.CHANGE)
      score_mat[i,j] = score; score_mat[j,i] = score
      if(is.null(second_type) || score < min_score) {
        first_type <- type1; second_type <- type2
        min_score = score
      }
    }
  }
  type1_pres = check_pairs_type(cell_type_profiles, bead, UMI_tot, score_mat, min_score, first_type, class_df, QL_score_cutoff, constrain, MIN.CHANGE = MIN.CHANGE)
  type2_pres = check_pairs_type(cell_type_profiles, bead, UMI_tot, score_mat, min_score, second_type, class_df, QL_score_cutoff, constrain, MIN.CHANGE = MIN.CHANGE)
  if(!type1_pres$all_pairs_class && !type2_pres$all_pairs_class) {
    spot_class <- "reject"
    singlet_score = min_score + 2 * doublet_like_cutoff #arbitrary
  }
  else if(type1_pres$all_pairs_class && !type2_pres$all_pairs_class) {
    first_class <- !type1_pres$all_pairs
    singlet_score = type1_pres$singlet_score
    spot_class = "doublet_uncertain"
  } else if(!type1_pres$all_pairs_class && type2_pres$all_pairs_class) {
    first_class <- !type2_pres$all_pairs
    singlet_score = type2_pres$singlet_score
    temp = first_type; first_type = second_type; second_type = temp
    spot_class = "doublet_uncertain"
  } else {
    spot_class = "doublet_certain"
    singlet_score = min(type1_pres$singlet_score, type2_pres$singlet_score)
    first_class <- !type1_pres$all_pairs; second_class <- !type2_pres$all_pairs
    if(type2_pres$singlet_score < type1_pres$singlet_score) {
      temp = first_type; first_type = second_type; second_type = temp
      first_class <- !type2_pres$all_pairs; second_class <- !type1_pres$all_pairs
    }
  }
  if(singlet_score - min_score < doublet_like_cutoff)
    spot_class = "singlet"
  doublet_results = decompose_sparse(cell_type_profiles, UMI_tot, bead, first_type, second_type, constrain = constrain, MIN.CHANGE = MIN.CHANGE)
  doublet_weights = doublet_results$weights; conv_doublet = doublet_results$converged
  spot_class <- factor(spot_class, c("reject", "singlet", "doublet_certain", "doublet_uncertain"))
  return(list(all_weights = all_weights, spot_class = spot_class, first_type = first_type, second_type = second_type,
              doublet_weights = doublet_weights, min_score = min_score, singlet_score = singlet_score,
              conv_all = conv_all, conv_doublet = conv_doublet, score_mat = score_mat,
              first_class = first_class, second_class = second_class))
}

solveIRWLS.weights <-function(S,B,nUMI, OLS=FALSE, constrain = TRUE, verbose = FALSE,
                              n.iter = 50, MIN_CHANGE = .001, bulk_mode = F, solution = NULL){
  if(!bulk_mode)
    B[B > K_val] <- K_val
  solution <- numeric(dim(S)[2])
  solution[] <- 1/length(solution)
  if(OLS) {
    solution<-solveOLS(S,B, solution, constrain = constrain) #first solve OLS, use this solution to find a starting point for the weights
    return(list(weights = solution, converged = T))
  }
  #solution <- runif(length(solution))*2 / length(solution) # random initialization
  names(solution) <- colnames(S)
  
  S_mat <<- matrix(0,nrow = dim(S)[1],ncol = dim(S)[2]*(dim(S)[2] + 1)/2)
  counter = 1
  for(i in 1:dim(S)[2])
    for(j in i:dim(S)[2]) {
      S_mat[,counter] <<- S[,i] * S[,j] # depends on n^2
      counter <- counter + 1
    }
  
  iterations<-0 #now use dampened WLS, iterate weights until convergence
  changes<-c()
  change<-1;
  while(change > MIN_CHANGE && iterations<n.iter){
    new_solution<-solveWLS(S,B,solution, nUMI,constrain=constrain, bulk_mode = bulk_mode)
    change<-norm(as.matrix(new_solution-solution))
    if(verbose) {
      print(paste("Change:",change))
      print(solution)
    }
    solution <- new_solution
    iterations<-iterations+1
  }
  return(list(weights = solution, converged = (change <= MIN_CHANGE)))
}


solveWLS<-function(S,B,initialSol, nUMI, bulk_mode = F, constrain = F){
  solution<-pmax(initialSol,0)
  prediction = abs(S%*%solution)
  threshold = max(1e-4, nUMI * 1e-7)
  prediction[prediction < threshold] <- threshold
  gene_list = rownames(S)
  derivatives <- get_der_fast(S, B, gene_list, prediction, bulk_mode = bulk_mode)
  d_vec <- -derivatives$grad
  D_mat <- psd(derivatives$hess)
  norm_factor <- norm(D_mat,"2")
  D_mat <- D_mat / norm_factor
  d_vec <- d_vec / norm_factor
  epsilon <- 1e-7; D_mat <- D_mat + epsilon * diag(length(d_vec))
  A<-cbind(diag(dim(S)[2]))
  bzero<- (-solution)
  alpha = 0.3
  if(constrain) {
    A_const = t(rbind(1,A))
    b_const <-c(1 - sum(solution),bzero)
    solution <- solution + alpha*quadprog::solve.QP(D_mat,d_vec,A_const,b_const,meq=1)$solution
  } else {
    solution <- solution + alpha*quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution
  }
  names(solution)<-colnames(S)
  return(solution)
}


get_der_fast <- function(S, B, gene_list, prediction, bulk_mode = F) {
  if(bulk_mode) {
    #d1_vec <- -t(log(prediction) - log(B))
    #d2_vec <- -t(1/prediction)
    d1_vec <- -2*t((log(prediction) - log(B))/prediction)
    d2_vec <- -2*t((1 - log(prediction) + log(B))/prediction^2)
  } else {
    d1_d2 <- get_d1_d2(B, prediction)
    d1_vec <- d1_d2$d1_vec
    d2_vec <- d1_d2$d2_vec
  }
  grad = -d1_vec %*% S;
  hess_c <- -d2_vec %*% S_mat
  hess <- matrix(0,nrow = dim(S)[2], ncol = dim(S)[2])
  counter = 1
  for(i in 1:dim(S)[2]) {
    l <- dim(S)[2] - i
    hess[i,i:dim(S)[2]] <- hess_c[counter:(counter+l)]
    hess[i,i] <- hess[i,i] / 2
    counter <- counter + l + 1
  }
  hess <- hess + t(hess)
  return(list(grad=grad, hess=hess))
}

#return positive semidefinite part
psd <- function(H, epsilon = 1e-3) {
  eig <- eigen(H);
  if(length(H) == 1)
    P <- eig$vectors %*% pmax(eig$values,epsilon) %*% t(eig$vectors)
  else
    P <- eig$vectors %*% diag(pmax(eig$values,epsilon)) %*% t(eig$vectors)
  return(P)
}


get_d1_d2 <- function(Y, lambda) {
  d_all <- calc_Q_all(Y, lambda)
  return(list(d1_vec = d_all$d1_vec, d2_vec = d_all$d2_vec))
}


calc_Q_all <- function(Y, lambda) {
  Y[Y > K_val] <- K_val
  epsilon <- 1e-4; X_max <- max(X_vals); delta <- 1e-6
  lambda <- pmin(pmax(epsilon, lambda),X_max - epsilon)
  
  l <- floor((lambda/delta)^(1/2))
  m <- pmin(l - 9,40) + pmax(ceiling(sqrt(pmax(l-48.7499,0)*4))-2,0)
  ti1 <- X_vals[m]; ti <- X_vals[m+1]; hi <- ti - ti1
  #Q0 <- cbind(Y+1, m); Q1 <- cbind(Y+1, m+1)
  Q0 <- cbind(Y+1, m); Q1 <- Q0; Q1[,2] <- Q1[,2] + 1
  fti1 <- Q_mat[Q0]; fti <- Q_mat[Q1]
  zi1 <- SQ_mat[Q0]; zi <- SQ_mat[Q1]
  diff1 <- lambda - ti1; diff2 <- ti - lambda
  diff3 <- fti/hi-zi*hi/6; diff4 <- fti1/hi-zi1*hi/6
  zdi <- zi / hi; zdi1 <- zi1 / hi
  # cubic spline interpolation
  d0_vec <- zdi*(diff1)^3/6 + zdi1*(diff2)^3/6 + diff3*diff1 + diff4*diff2
  d1_vec <- zdi*(diff1)^2/2 - zdi1*(diff2)^2/2 + diff3 - diff4
  d2_vec <- zdi*(diff1) + zdi1*(diff2)
  return(list(d0_vec = d0_vec, d1_vec = d1_vec, d2_vec = d2_vec))
}


run.CSIDE.replicates.local <- function (RCTD.replicates, cell_types, explanatory.variable.replicates = NULL, 
                                        X.replicates = NULL, cell_type_threshold = 125, gene_threshold = 5e-05, 
                                        doublet_mode = T, weight_threshold = NULL, sigma_gene = T, 
                                        PRECISION.THRESHOLD = 0.05, cell_types_present = NULL, fdr = 0.01, 
                                        population_de = F, replicate_index = NULL, normalize_expr = F, 
                                        test_genes_sig_individual = F, de_mode = "single", df = 15, 
                                        barcodes = NULL, log_fc_thresh = 0.4, test_error = F, params_to_test = NULL, 
                                        test_mode = "individual") 
{
  if (!(de_mode %in% c("single", "nonparam", "general"))) 
    stop("run.CISDE.replicates: de_mode must be set to \"single\", \"general\", or \"nonparam\".")
  if (is.null(cell_types)) 
    stop("run.CSIDE.replicates: cell_types must not be null.")
  if (is.null(replicate_index)) 
    replicate_index <- 1:length(RCTD.replicates@RCTD.reps)
  if (any(!(replicate_index %in% 1:length(RCTD.replicates@RCTD.reps)))) 
    stop("run.CSIDE.replicates: replicate_index must be a subest of 1:N_replicates")
  if (test_error) 
    warning("run.CSIDE.replicates: test_error is TRUE, so this run will just test C-SIDE for errors without running C-SIDE.")
  for (i in replicate_index) {
    if (test_error) 
      message(paste("run.CSIDE.replicates: testing CSIDE for errors for replicate", 
                    i))
    else message(paste("run.CSIDE.replicates: running CSIDE for replicate", 
                       i))
    if (de_mode == "single") {
      if (is.null(explanatory.variable.replicates)) 
        stop("run.CSIDE.replicates: if de_mode = single, explanatory.variable.replicates cannot be null.")
      if (class(explanatory.variable.replicates) != "list") 
        stop("run.CSIDE.replicates: explanatory.variable.replicates must be a list of explanatory variable vectors for each replicate.")
      if (length(RCTD.replicates@RCTD.reps) != length(explanatory.variable.replicates)) 
        stop("create.RCTD.replicates: length(explanatory.variable.replicates) is not equal to the number of RCTD replicates, as required.")
      RCTD.replicates@RCTD.reps[[i]] <- run.CSIDE.single(RCTD.replicates@RCTD.reps[[i]], 
                                                         explanatory.variable.replicates[[i]], cell_types = cell_types, 
                                                         cell_type_threshold = cell_type_threshold, gene_threshold = gene_threshold, 
                                                         doublet_mode = doublet_mode, weight_threshold = weight_threshold, 
                                                         sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD, 
                                                         test_genes_sig = test_genes_sig_individual, cell_types_present = cell_types_present, 
                                                         fdr = fdr, log_fc_thresh = log_fc_thresh, test_error = test_error)
    }
    else if (de_mode == "nonparam") {
      RCTD.replicates@RCTD.reps[[i]] <- run.CSIDE.nonparam(RCTD.replicates@RCTD.reps[[i]], 
                                                           df = df, barcodes = barcodes, cell_types = cell_types, 
                                                           cell_type_threshold = cell_type_threshold, gene_threshold = gene_threshold, 
                                                           doublet_mode = doublet_mode, weight_threshold = weight_threshold, 
                                                           test_genes_sig = test_genes_sig_individual, sigma_gene = sigma_gene, 
                                                           PRECISION.THRESHOLD = PRECISION.THRESHOLD, cell_types_present = cell_types_present, 
                                                           fdr = fdr, test_error = test_error)
    }
    else {
      if (is.null(X.replicates)) 
        stop("run.CSIDE.replicates: if de_mode = single, X.replicates cannot be null.")
      if (class(X.replicates) != "list") 
        stop("run.CSIDE.replicates: X.replicates must be a list of design matrices for each replicate.")
      if (length(RCTD.replicates@RCTD.reps) != length(X.replicates)) 
        stop("run.CSIDE.replicates: length(X.replicates) is not equal to the number of RCTD replicates, as required.")
      X <- X.replicates[[i]]
      if (length(setdiff(rownames(X), rownames(RCTD.replicates@RCTD.reps[[i]]@results$weights))) > 
          0) 
        warning("run.CSIDE.replicates: some elements of rownames(X.replicates) do not appear in myRCTD object (myRCTD@results$weights) for this replicate, but they are required to be a subset.")
      RCTD.replicates@RCTD.reps[[i]] <- run.CSIDE(RCTD.replicates@RCTD.reps[[i]], 
                                                  X, rownames(X), cell_types = cell_types, cell_type_threshold = cell_type_threshold, 
                                                  gene_threshold = gene_threshold, doublet_mode = doublet_mode, 
                                                  test_genes_sig = test_genes_sig_individual, weight_threshold = weight_threshold, 
                                                  sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD, 
                                                  cell_types_present = cell_types_present, fdr = fdr, 
                                                  log_fc_thresh = log_fc_thresh, test_error = test_error, 
                                                  params_to_test = params_to_test)
    }
  }
  if (population_de) 
    RCTD.replicates <- CSIDE.population.inference(RCTD.replicates, 
                                                  log_fc_thresh = log_fc_thresh, fdr = fdr)
  return(RCTD.replicates)
}

run.CSIDE <- function(myRCTD, X, barcodes, cell_types = NULL, gene_threshold = 5e-5, cell_type_threshold = 125,
                      doublet_mode = T, test_mode = 'individual', weight_threshold = NULL,
                      sigma_gene = T, PRECISION.THRESHOLD = 0.05, cell_types_present = NULL,
                      test_genes_sig = T, fdr = .01, cell_type_specific = NULL,
                      params_to_test = NULL, normalize_expr = F, logs=F, log_fc_thresh = 0.4,
                      cell_type_filter = NULL, test_error = F, fdr_method = 'BH') {
  X <- check_designmatrix(X, 'run.CSIDE', require_2d = TRUE)
  if(is.null(cell_type_specific))
    cell_type_specific <- !logical(dim(X)[2])
  check_cell_type_specific(cell_type_specific, dim(X)[2])
  X1 <- X[,!cell_type_specific];
  if(any(!cell_type_specific))
    X2 <- X[,cell_type_specific]
  else
    X2 <- X
  return(run.CSIDE.general(myRCTD, X1, X2, barcodes, cell_types, cell_type_threshold = cell_type_threshold,
                           gene_threshold = gene_threshold,
                           doublet_mode = doublet_mode, test_mode = test_mode, weight_threshold = weight_threshold,
                           sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD, params_to_test = params_to_test,
                           cell_types_present = cell_types_present, test_genes_sig = test_genes_sig,
                           fdr = fdr, normalize_expr = normalize_expr, logs=logs,
                           cell_type_filter = cell_type_filter, log_fc_thresh = log_fc_thresh, test_error = test_error, fdr_method = fdr_method))
}

run.CSIDE.general <- function(myRCTD, X1, X2, barcodes, cell_types = NULL, gene_threshold = 5e-5, cell_type_threshold = 125,
                              doublet_mode = T, test_mode = 'individual', weight_threshold = NULL,
                              sigma_gene = T, PRECISION.THRESHOLD = 0.05, cell_types_present = NULL,
                              test_genes_sig = T, fdr = .01, params_to_test = NULL, normalize_expr = F,
                              logs=F, cell_type_filter = NULL, log_fc_thresh = 0.4, test_error = FALSE, fdr_method = 'BH') {
  if(gene_threshold == .01 || fdr == 0.25 || cell_type_threshold <= 10 ||
     (!is.null(weight_threshold) && weight_threshold == 0.1))
    warning('run.CSIDE.general: some parameters are set to the CSIDE vignette values, which are intended for testing but not proper execution. For more accurate results, consider using the default parameters to this function.')
  if(doublet_mode && myRCTD@config$RCTDmode != 'doublet')
    stop('run.CSIDE.general: attempted to run CSIDE in doublet mode, but RCTD was not run in doublet mode. Please run CSIDE in full mode (doublet_mode = F) or run RCTD in doublet mode.')
  if(!any("cell_types_assigned" %in% names(myRCTD@internal_vars)) || !myRCTD@internal_vars$cell_types_assigned)
    stop('run.CSIDE.general: cannot run CSIDE unless cell types have been assigned. If cell types have been assigned, you may run "myRCTD <- set_cell_types_assigned(myRCTD)".')
  if((myRCTD@config$doublet_mode != 'multi') && (length(setdiff(barcodes,rownames(myRCTD@results$weights))) > 0)) {
    warning('run.CSIDE.general: some elements of barcodes do not appear in myRCTD object (myRCTD@results$weights), but they are required to be a subset. Downsampling barcodes to the intersection of the two sets.')
    barcodes <- intersect(barcodes,rownames(myRCTD@results$weights))
  }
  cell_type_info <- myRCTD@cell_type_info$info
  if(doublet_mode) {
    my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
    thresh <- 0.999
  } else if(myRCTD@config$doublet_mode == "multi") {
    my_beta <- get_beta_multi(barcodes, cell_type_info[[2]], myRCTD@results, myRCTD@spatialRNA@coords)
    thresh <- 0.999
  } else {
    my_beta <- as.matrix(sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), '/'))
    thresh <- 0.8
  }
  if(!is.null(weight_threshold))
    thresh <- weight_threshold
  cell_types <- choose_cell_types(myRCTD, barcodes, doublet_mode, cell_type_threshold, cell_types,
                                  my_beta, thresh, cell_type_filter)
  if(length(cell_types) < 1)
    stop('run.CSIDE.general: zero cell types remain. Cannot run CSIDE with zero cell types.')
  message(paste0("run.CSIDE.general: running CSIDE with cell types ",paste(cell_types, collapse = ', ')))
  X1 <- check_designmatrix(X1, 'run.CSIDE.general')
  X2 <- check_designmatrix(X2, 'run.CSIDE.general', require_2d = TRUE)
  if(!(test_mode %in% c('individual', 'categorical')))
    stop(c('run.CSIDE.general: not valid test_mode = ',test_mode,'. Please set test_mode = "categorical" or "individual".'))
  if(is.null(params_to_test))
    if(test_mode == 'individual')
      params_to_test <- min(2, dim(X2)[2])
  else
    params_to_test <- 1:dim(X2)[2]
  if(normalize_expr && (test_mode != 'individual' || length(params_to_test) > 1))
    stop('run.CSIDE.general: Setting normalize_expr = TRUE is only valid for testing single parameters with test_mode = individual')
  message(paste0("run.CSIDE.general: configure params_to_test = ",
                 paste(paste0(params_to_test, ', ', collapse = ""))))
  if(any(!(params_to_test %in% 1:dim(X2)[2])))
    stop(c('run.CSIDE.general: params_to_test must be a vector of integers from 1 to dim(X2)[2] = ', dim(X2)[2],
           'please make sure that tested parameters are in the required range.'))
  if(test_mode == 'categorical' && any(!(X2[,params_to_test] %in% c(0,1))))
    stop(c('run.CSIDE.general: for test_mode = categorical, colums params_to_test, ',params_to_test,', must have values 0 or 1.'))
  if(is.null(cell_types_present))
    cell_types_present <- cell_types
  if(any(!(barcodes %in% rownames(X1))) || any(!(barcodes %in% rownames(X2))))
    stop('run.CSIDE.general: some barcodes do not appear in the rownames of X1 or X2.')
  puck = myRCTD@originalSpatialRNA
  gene_list_tot <- filter_genes(puck, threshold = gene_threshold)
  if(length(gene_list_tot) == 0)
    stop('run.CSIDE.general: no genes past threshold. Please consider lowering gene_threshold.')
  if(length(intersect(gene_list_tot,rownames(cell_type_info[[1]]))) == 0)
    stop('run.CSIDE.general: no genes that past threshold were contained in the single cell reference. Please lower gene threshold or ensure that there is agreement between the single cell reference genes and the SpatialRNA genes.')
  nUMI <- puck@nUMI[barcodes]
  res <- filter_barcodes_cell_types(barcodes, cell_types, my_beta, thresh = thresh)
  if(test_error)
    return(myRCTD)
  barcodes <- res$barcodes; my_beta <- res$my_beta
  sigma_init <- as.character(100*myRCTD@internal_vars$sigma)
  if(sigma_gene) {
    set_global_Q_all()
    sigma_set <- sigma_init
    set_likelihood_vars(Q_mat_all[[sigma_init]], X_vals, sigma = sigma_set)
  } else {
    set_likelihood_vars_sigma(sigma_init)
  }
  gene_fits <- get_de_gene_fits(X1[barcodes, , drop = FALSE],X2[barcodes, , drop = FALSE],my_beta, nUMI[barcodes], gene_list_tot,
                                cell_types, restrict_puck(puck, barcodes), barcodes, sigma_init,
                                test_mode, numCores = myRCTD@config$max_cores, sigma_gene = sigma_gene,
                                PRECISION.THRESHOLD = PRECISION.THRESHOLD, params_to_test = params_to_test,
                                logs=logs)
  if(normalize_expr)
    myRCTD <- normalize_de_estimates(myRCTD, normalize_expr = normalize_expr,
                                     param_position = params_to_test)
  if(test_genes_sig) {
    both_gene_list <- get_sig_genes(puck, myRCTD, gene_list_tot, cell_types, my_beta, barcodes, nUMI,
                                    gene_fits, cell_types_present, X2, test_mode, fdr = fdr,
                                    params_to_test = params_to_test, normalize_expr = normalize_expr,
                                    log_fc_thresh = log_fc_thresh, fdr_method = fdr_method)
    sig_gene_list <- both_gene_list$sig_gene_list; all_gene_list <- both_gene_list$all_gene_list
  } else {
    sig_gene_list <- NULL
    all_gene_list <- NULL
  }
  myRCTD@internal_vars_de <- list(barcodes = barcodes, cell_types = cell_types, doublet_mode = doublet_mode,
                                  cell_types_present = cell_types_present,
                                  my_beta = my_beta, X1 = X1, X2 = X2,
                                  test_mode = test_mode, params_to_test = params_to_test)
  myRCTD@de_results <- list(gene_fits = gene_fits, sig_gene_list = sig_gene_list, all_gene_list = all_gene_list)
  return(myRCTD)
}


get_de_gene_fits <- function(X1,X2,my_beta, nUMI, gene_list, cell_types, puck, barcodes, sigma_init, test_mode,
                             numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.05, params_to_test = 2, logs=F) {
  results_list <- fit_de_genes(X1,X2,my_beta, nUMI, gene_list, puck, barcodes,
                               sigma_init, test_mode, numCores = numCores,
                               sigma_gene = sigma_gene,
                               PRECISION.THRESHOLD = PRECISION.THRESHOLD, logs = logs)
  N_genes <- length(results_list)
  intercept_val <- matrix(0,nrow = N_genes, ncol = length(cell_types))
  mean_val <- matrix(0,nrow = N_genes, ncol = length(cell_types))
  all_vals <- array(0, dim = c(N_genes, dim(X2)[2],length(cell_types)))
  dimnames(all_vals)[[1]] <- gene_list
  dimnames(all_vals)[[3]] <- cell_types
  con_val <- logical(N_genes)
  ll_val <- numeric(N_genes)
  n_val <- numeric(N_genes)
  sigma_g <- numeric(N_genes)
  names(sigma_g) <- gene_list
  I_val <- list()
  names(n_val) <- gene_list
  names(con_val) <- gene_list
  names(ll_val) <- gene_list
  rownames(mean_val) <- gene_list; colnames(mean_val) <- cell_types
  rownames(intercept_val) <- gene_list; colnames(intercept_val) <- cell_types
  d_vals <- matrix(0,nrow=N_genes,ncol=dim(X2)[2]*length(cell_types))
  s_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  precision_mat <- matrix(0, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  con_all <- matrix(FALSE, nrow = N_genes, ncol = dim(X2)[2]*length(cell_types))
  con_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cell_types))
  error_mat <- matrix(FALSE, nrow = N_genes, ncol = length(cell_types))
  rownames(precision_mat) <- gene_list; rownames(con_all) <- gene_list
  rownames(s_mat) <- gene_list; rownames(con_mat) <- gene_list; rownames(error_mat) <- gene_list
  colnames(s_mat) <- get_param_names(X1,X2, cell_types)
  colnames(precision_mat) <- get_param_names(X1,X2, cell_types)
  colnames(con_all) <- get_param_names(X1,X2, cell_types)
  colnames(con_mat) <- cell_types
  colnames(error_mat) <- cell_types
  rownames(d_vals) <- gene_list
  for(i in 1:N_genes) {
    sigma_g[i] <- results_list[[i]]$sigma_s_best
    res <- results_list[[i]]$res
    d_vals[i,] <- res$d
    mean_val[i,] <- res$alpha2[params_to_test[1],]
    intercept_val[i,] <- res$alpha2[1,]
    all_vals[i, ,] <- res$alpha2
    con_val[i] <- res$converged
    precision_mat[i,] <- res$precision
    ll_val[i] <- res$log_l
    n_val[i] <- res$n.iter
    I_val[[i]] <- res$I
    s_mat[i,] <- sqrt(diag(I_val[[i]]))
    con_mat[i,] <- res$converged_vec
    con_all[i,] <- res$precision < PRECISION.THRESHOLD
    error_mat[i,] <- res$error_vec
  }
  return(list(mean_val = mean_val, con_val = con_val, ll_val = ll_val, I_val = I_val, s_mat = s_mat,
              n.iter = n_val,d_vals = d_vals, intercept_val = intercept_val, all_vals = all_vals,
              precision_mat = precision_mat, sigma_g = sigma_g, con_mat = con_mat, con_all = con_all, error_mat = error_mat))
}

fit_de_genes <- function(X1,X2,my_beta, nUMI, gene_list, puck, barcodes, sigma_init, test_mode, numCores = 4, sigma_gene = T, PRECISION.THRESHOLD = 0.05,
                         logs=F) {
  results_list <- list()
  if(numCores == 1) {
    for(i in 1:length(gene_list)) {
      message(i)
      gene <- gene_list[i]
      Y <- puck@counts[gene, barcodes]
      results_list[[i]] <- estimate_gene_wrapper(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = sigma_gene, PRECISION.THRESHOLD = PRECISION.THRESHOLD)
    }
  } else {
    cl <- parallel::makeCluster(numCores,setup_strategy = "sequential",outfile="") #makeForkCluster
    doParallel::registerDoParallel(cl)
    environ = c('estimate_effects_trust', 'solveIRWLS.effects_trust', 'K_val','X_vals',
                'calc_log_l_vec', 'get_d1_d2', 'calc_Q_all','psd','construct_hess_fast',
                'choose_sigma_gene', 'estimate_gene_wrapper', 'check_converged_vec', 'calc_log_l_vec_fast')
    if(sigma_gene)
      environ <- c(environ, 'Q_mat_all', 'SQ_mat_all')
    else
      environ <- c(environ, 'Q_mat', 'SQ_mat')
    if (logs) {
      out_file = "logs/de_log.txt"
      if(!dir.exists('logs'))
        dir.create('logs')
      if(file.exists(out_file))
        file.remove(out_file)
    }
    results_list <- foreach::foreach(i = 1:length(gene_list), .packages = c("quadprog", "Rfast"), .export = environ) %dopar% {
      if (logs) {
        if(i %% 1 == 0) { ##10
          cat(paste0("Testing sample: ",i," gene ", gene_list[i],"\n"), file=out_file, append=TRUE)
        }
      }
      assign("X_vals",X_vals, envir = globalenv()); assign("K_val",K_val, envir = globalenv());
      if(sigma_gene) {
        assign("Q_mat_all",Q_mat_all, envir = globalenv());
        assign("SQ_mat_all",SQ_mat_all, envir = globalenv());
      } else {
        assign("Q_mat",Q_mat, envir = globalenv()); assign("SQ_mat",SQ_mat, envir = globalenv())
      }
      gene <- gene_list[i]
      Y <- puck@counts[gene, barcodes]
      res <- estimate_gene_wrapper(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3, sigma_gene = sigma_gene)
    }
    parallel::stopCluster(cl)
  }
  return(results_list)
}




