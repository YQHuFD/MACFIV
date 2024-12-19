library(tidyverse) 
library(rootSolve)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(foreach)
library(doParallel)
library(glmnet)

setwd("C:/Users/17317/Desktop/NCI_MASIRAL_code/real_analysis")
ARIC = read.csv("ARIC_pheno.csv")
ARIC = (ARIC %>% drop_na(bmi_baseline,dbp_baseline))
ARIC.male = ARIC
#signal = read.table("sbp.1e-3.txt") 

sample.id.old = as.vector(unlist(read.table("ID2.txt"))) #All the individuals' ID
#sample.id.new = sample.id.old[sample.id.old %in% ARIC[,2]]
sample.id.new.male = sample.id.old[sample.id.old %in% ARIC.male[,2]]

valid.sample.number = length(sample.id.new.male)

computation.id = sample.id.new.male

#genotype.orig = NULL
#for(i in 1:22)
#{
#  path = paste0("chr", i, ".DS.txt")
#  genotype.orig = cbind(genotype.orig, t(read.table(file = path)))
#} 
genotype.orig = t(read.table(file = 'C:/Users/17317/Desktop/NCI_MASIRAL_code/real_analysis/bmi_snp.txt'))


geno.name = genotype.orig[1,] #name of snps
row = dim(genotype.orig)[1] #number of snps
genotype = data.frame(id = c(sample.id.old), geno.name = genotype.orig[2:row,]) # integrate samples and snps 
names(genotype) = c("id", geno.name)# modify its name
genotype = genotype[genotype[,1] %in% sample.id.new.male,]


genotype.computation = genotype[genotype[,1] %in% computation.id,]

ARIC.computation = ARIC[ARIC[,2]%in% computation.id, ]

pheno_x.computation = ARIC.computation[, colnames(ARIC.computation) %in% c("bmi_baseline")]
pheno_x.inveranktrans.computation = qnorm((rank(pheno_x.computation, na.last="keep")-0.5)/sum(!is.na(pheno_x.computation)))
pheno_y.computation = ARIC.computation[, colnames(ARIC.computation) %in% c("dbp_baseline")] #phenotype y

tmp1 = as.data.frame(genotype.computation[,2:dim(genotype.computation)[2]])
tmp1 = apply(tmp1,2,function(x) as.numeric(as.character(x)))
tmp2 = t(apply(tmp1,1,function(x)x-apply(tmp1,2,mean)))

Z <- unname(tmp2)
X <- as.matrix((pheno_x.computation - mean(pheno_x.computation))/sd(pheno_x.computation))
Y <- as.matrix((pheno_y.computation - mean(pheno_y.computation))/sd(pheno_y.computation))


pearson <- rep(0,ncol(Z))

for (j in 1:ncol(Z)){
  pearson[j] <- cor(X,Z[ ,j])
}


index <- 1:length(pearson)      
data <- data.frame(Index = index, Coefficient = pearson)

p <- ggplot(data, aes(x = Index, y = Coefficient)) +
  geom_point(color = "blue", alpha = 0.8, size = 1.8) + 
  labs(
    title = "Correlation between BMI and Potential Instruments in Dataset",
    x = "Index of Potential Instruments",
    y = "Correlation with BMI"
  ) +
  scale_x_continuous(
    breaks = seq(0, max(index), by = 20),  
    limits = c(0, max(index))              
  ) +
  scale_y_continuous(
    limits = c(-0.1, 0.1),                 
    breaks = seq(-0.1, 0.1, by = 0.05)     
  ) +
  #theme_classic(base_size = 12) +          
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
    axis.text.x = element_text(size = 10),                             
    axis.text.y = element_text(size = 10),                             
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),   
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1)
  )




gradient <- function(Y, B, Z, CF, theta) {
  beta <- theta[1:ncol(B)]
  alpha <- theta[(ncol(B) + 1):(ncol(B) + ncol(Z))]
  rho <- theta[length(theta)]
  
  residual <- Y - B %*% beta - Z %*% alpha - rho * CF
  
  grad_beta <- -t(B) %*% residual
  grad_alpha <- -t(Z) %*% residual + SCAD_derivative(alpha)
  grad_rho <- -sum(CF * residual)
  
  return(c(grad_beta, grad_alpha, grad_rho))
}


hessian_sample <- function(B, Z, CF, alpha) {
  H_beta_beta <- t(B) %*% B
  H_alpha_alpha <- t(Z) %*% Z + SCAD_second_derivative(alpha)
  H_rho_rho <- sum(CF^2)
  
  H_beta_alpha <- t(B) %*% Z
  H_beta_rho <- t(B) %*% CF
  H_alpha_rho <- t(Z) %*% CF
  
  H <- rbind(
    cbind(H_beta_beta, H_beta_alpha, H_beta_rho),
    cbind(t(H_beta_alpha), H_alpha_alpha, H_alpha_rho),
    cbind(t(H_beta_rho), t(H_alpha_rho), H_rho_rho)
  )
  
  return(H)
}


SCAD_derivative <- function(alpha, lambda = best_lambda, a = 3.7) {

  grad <- ifelse(abs(alpha) <= lambda, lambda, 
                 ifelse(abs(alpha) <= a * lambda, (a * lambda - abs(alpha)) / (a - 1), 0))
  return(sign(alpha) * grad)
}


SCAD_second_derivative <- function(alpha, lambda = best_lambda, a = 3.7) {

  hess <- ifelse(abs(alpha) <= lambda, 0, 
                 ifelse(abs(alpha) <= a * lambda, -1 / (a - 1), 0))
  return(diag(hess))
}


compute_H_V <- function(Y, B, Z, CF, theta_hat) {
  n <- length(Y)
  
  # Hessian H 
  H <- matrix(0, nrow = length(theta_hat), ncol = length(theta_hat))  
  for (i in 1:n) {
    H_i <- hessian_sample(B[i, , drop = FALSE], Z[i, , drop = FALSE], CF[i], theta_hat[(ncol(B) + 1):(ncol(B) + ncol(Z))])
    H <- H + H_i 
  }
  
  # Covariance V 
  gradients <- sapply(1:n, function(i) {
    gradient(Y[i, drop = FALSE], B[i, , drop = FALSE], Z[i, , drop = FALSE ], CF[i], theta_hat)
  })
  
  H <- H / n
  V <- cov(t(gradients))
  
  return(list(H = H, V = V))
}


compute_U <- function(Y, B, Z, CF, theta_hat, beta_length) {
  HV <- compute_H_V(Y, B, Z, CF, theta_hat)
  H <- HV$H
  V <- HV$V
  
  H_inv <- ginv(H)
  U <- H_inv %*% V %*% H_inv
  
  U_beta_beta <- U[1:beta_length, 1:beta_length]
  return(U_beta_beta)
}


n <- nrow(Z)
p <- ncol(Z)

B <- B_spline_basis
theta_hat <- model$beta
U <- compute_U(Y, B, Z, CF, theta_hat, ncol(B))

sd_g <- sapply(1:length(X), function(i) {
  B <- B_spline_basis[i, , drop = FALSE] 
  sqrt(B %*% U %*% t(B) / n)  
})


z_critical <- qnorm(0.975)
ci_lower <- f_hat - z_critical * sd_g
ci_upper <- f_hat + z_critical * sd_g
ci_f <- data.frame(
  X = X,
  f_hat = f_hat,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)


ci_f_filtered <- ci_f %>% filter(X >= -1.5 & X <= 2.5)

max_point_filtered <- ci_f_filtered[which.max(ci_f_filtered$f_hat), ]

label_text_with_math <- expression(BMI ~ "\u2248" ~ 33.41 ~ "kg/m"^2)



p <- ggplot(ci_f_filtered, aes(x = X)) +

  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "gray", alpha = 0.5)+

  geom_line(aes(y = f_hat), color = "black", size = 0.5) +

  geom_point(data = max_point_filtered, aes(x = X, y = f_hat), color = "red", size = 3) +

  geom_text(
    data = max_point_filtered,
    aes(x = X, y = f_hat),
    label = label_text_with_math,
    vjust = -1.5,  
    hjust = 0.5,   
    color = "black", size = 4  
  ) +
  labs(x = expression(BMI ~ (standardized)), y = expression(DBP ~ (standardized)) ) +

  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),                 
    panel.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
    axis.line = element_blank(),                 
    axis.ticks = element_line(color = "black", size = 0.5),  
    axis.title = element_text(face = "bold", size = 14),     
    axis.text = element_text(color = "black", size = 12)     
  )

