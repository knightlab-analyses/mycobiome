# NOTE: This function was adapted nearly verbatim from the Bioconductor R pvca package:
# https://www.bioconductor.org/packages/release/bioc/html/pvca.html
#' A function for principal variance component analysis
#'
#' The function is written based on the 'pvcaBatchAssess' function of the PVCA R package
#' and slightly changed to make it more efficient and flexible for sequencing read counts data.
#' (http://watson.nci.nih.gov/bioc_mirror/packages/release/bioc/manuals/pvca/man/pvca.pdf)
#'
#' @param counts The Normalized(e.g. TMM)/ log-transformed reads count matrix from sequencing data (row:gene/feature, col:sample)
#' @param meta  The Meta data matrix containing predictor variables (row:sample, col:predictor)
#' @param threshold The proportion of the variation in read counts explained by top k PCs. This value determines the number of top PCs to be used in pvca.
#' @param inter TRUE/FALSE - include/do not include pairwise interactions of predictors
#'
#' @return std.prop.val The vector of proportions of variation explained by each predictor.
#'
#' @export
#'

PVCA <- function(counts, meta, threshold, inter){
  
  counts.center <- t(apply(counts, 1, scale, center=TRUE, scale=FALSE))
  cor.counts <- cor(counts.center)
  dim(cor.counts)
  eigen.counts <- eigen(cor.counts)
  eigen.mat <- eigen.counts$vectors
  eigen.val <- eigen.counts$values
  n.eigen <- length(eigen.val)
  eigen.val.sum <- sum(eigen.val)
  percents.pcs <- eigen.val/eigen.val.sum
  meta <- as.data.frame(meta)
  
  all <- 0
  npc.in <- 0
  for(i in 1:n.eigen){
    all <- all + percents.pcs[i]
    npc.in <- npc.in + 1
    if(all > threshold){break}
  }
  if (npc.in < 3) {npc <- 3}
  
  pred.list <- colnames(meta)
  meta <- droplevels(meta)
  
  n.preds <- ncol(meta) + 1
  if(inter) {n.preds <- n.preds + choose(ncol(meta),2)}
  
  ran.pred.list <- c()
  for(i in 1:ncol(meta)){
    ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i],")"))
  }
  ##interactions
  if(inter){
    for(i in 1:(ncol(meta)-1)){
      for(j in (i+1):ncol(meta)){
        ran.pred.list <- c(ran.pred.list, paste0("(1|", pred.list[i], ":", pred.list[j], ")"))
        pred.list <- c(pred.list, paste0(pred.list[i], ":", pred.list[j]))
      }
    }
  }
  formula <- paste(ran.pred.list, collapse = " + ")
  formula <- paste("pc", formula, sep=" ~ ")
  ran.var.mat <- NULL
  for(i in 1:npc.in){
    dat <- cbind(eigen.mat[,i],meta)
    colnames(dat) <- c("pc",colnames(meta))
    Rm1ML <- lme4::lmer(formula, dat, REML = TRUE, verbose = FALSE, na.action = na.omit)
    var.vec <- unlist(VarCorr(Rm1ML))
    ran.var.mat <- rbind(ran.var.mat, c(var.vec[pred.list], resid = sigma(Rm1ML)^2))
  }
  ran.var.mat.std <- ran.var.mat/rowSums(ran.var.mat)
  wgt.vec <- eigen.val/eigen.val.sum
  prop.var <- colSums(ran.var.mat.std*wgt.vec[1:npc.in])
  std.prop.var <- prop.var/sum(prop.var)
  std.prop.var
}