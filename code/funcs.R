## enrichment functions

create_design <- function(y, obs_id){
  if(is.null(obs_id)) {
    #for unpaired datasets
    design <- model.matrix(~y) 
  } else{
    design <- model.matrix(~y+obs_id) 
  }
  return(design)
}

## FGSEA
cfgsea <- function(eset, y, obs_id, db){
  # differential expression analysis
  design <- create_design(y, obs_id)
  
  fit <- limma::eBayes(limma::lmFit(eset, design))
  top <- limma::topTable(fit, coef = 2, n = nrow(fit))
  
  # fgsea
  stats = top$t
  names(stats) <- rownames(top)
  run_fgsea <- fgsea::fgsea(pathways = db, 
                            stats    = stats,
                            minSize  = min(sapply(db, length)),
                            maxSize  = max(sapply(db, length)))
  result <- run_fgsea$pval
  names(result) <- run_fgsea$pathway
  result[order(result)]
}

cfgsea_p = function(eset, y, obs_id, dbs, cores){
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", cores), type = "SOCK")
  parallel::clusterExport(cl, varlist = c("eset", "y", "obs_id", "cfgsea", "create_design"), 
                          envir = environment())
  result <- parallel::parLapply(cl, dbs, function(db, eset, y, obs_id) {
    cfgsea(eset, y, obs_id, db)
  }, eset, y, obs_id)
  parallel::stopCluster(cl)
  return(result)
}

# GSVA
cgsva = function(eset, y, obs_id, db){
  gsva_eset <- GSVA::gsva(eset, db, verbose=FALSE)
  design <- create_design(y, obs_id)
  fit <- limma::eBayes(limma::lmFit(gsva_eset, design))
  top <- limma::topTable(fit, coef = 2, n = nrow(fit))
  pval <- top$P.Value
  names(pval) <- rownames(top)
  pval[order(pval)]
}

cgsva_p = function(eset, y, obs_id, dbs, cores){
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", cores), type = "SOCK")
  parallel::clusterExport(cl, varlist = c("eset", "y", "obs_id", "cgsva", "create_design"), 
                          envir = environment())
  result = parallel::parLapply(cl, dbs, function(db, eset, y, obs_id) {
    cgsva(eset, y, obs_id, db)
  }, eset, y, obs_id)
  parallel::stopCluster(cl)
  return(result)
}

# PCA
cpca = function(eset, y, obs_id, db){
  pc <- t(sapply(db, function(ligand){
    genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
    prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
  }))
  design <- create_design(y, obs_id)
  fit <- limma::eBayes(limma::lmFit(pc, design))
  top <- limma::topTable(fit, coef = 2, n = nrow(fit))
  pval <- top$P.Value
  names(pval) <- rownames(top)
  pval[order(pval)]
}

cpca_p = function(eset, y, obs_id, dbs, cores){
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", cores), type = "SOCK")
  parallel::clusterExport(cl, varlist = c("eset", "y", "obs_id", "cpca", "create_design"), 
                          envir = environment())
  result <- parallel::parLapply(cl, dbs, function(db, eset, y, obs_id) {
    cpca(eset, y, obs_id, db)
  }, eset, y, obs_id)
  parallel::stopCluster(cl)
  return(result)
}


# PCR
cpcr = function(eset, y, obs_id, db){
  pcs <- sapply(db, function(ligand){
    genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
    prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
  })
  fit <- mixOmics::plsda(pcs, y)
  coef <- abs(mixOmics::selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(mixOmics::selectVar(fit, comp=1)$value)
  coef[order(coef, decreasing = TRUE)]
}

cpcr_p = function(eset, y, obs_id, dbs, cores){
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", cores), type = "SOCK")
  parallel::clusterExport(cl, varlist = c("eset", "y", "obs_id", "cpcr"), 
                          envir = environment())
  result <- parallel::parLapply(cl, dbs, function(db, eset, y, obs_id) {
    cpcr(eset, y, obs_id, db)
  }, eset, y, obs_id)
  parallel::stopCluster(cl)
  return(result)
}


# GSVAR
cgsvar = function(eset, y, obs_id, db){
  gsva_eset <- t(GSVA::gsva(eset, db, verbose=FALSE))
  fit <- mixOmics::plsda(gsva_eset, y)
  coef <- abs(mixOmics::selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(mixOmics::selectVar(fit, comp=1)$value)
  coef[order(coef, decreasing = TRUE)]
}

cgsvar_p = function(eset, y, obs_id, dbs, cores){
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", cores), type = "SOCK")
  parallel::clusterExport(cl, varlist = c("eset", "y", "obs_id", "cgsvar"), 
                          envir = environment())
  result <- parallel::parLapply(cl, dbs, function(db, eset, y, obs_id) {
    cgsvar(eset, y, obs_id, db)
  }, eset, y, obs_id)
  parallel::stopCluster(cl)
  return(result)
}
