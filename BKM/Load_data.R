loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# Create the matrix ####
BKM_HMM_Resultsparams <- loadRData("Results/BKM_HMM_Resultsparams.RData")
BKM_DA_Resultsparams <- loadRData("Results/BKM_DA_Resultsparams.RData")
BKM_DA_scaled_Resultsparams <- loadRData("Results/BKM_DA_scaled_Resultsparams.RData")

BKM_HMM_Resultsparams <- BKM_HMM_Resultsparams[,-c(1:4)]

BKM_Resultsparams = cbind(BKM_DA_Resultsparams, BKM_DA_scaled_Resultsparams)

BKM_Resultsparams = cbind(BKM_Resultsparams, BKM_HMM_Resultsparams)
 
latex_rnames <- c("$\\sigma_{y}^{2}$","$\\alpha_{\\rho}$","$\\beta_{\\rho}$","$\\alpha_{\\lambda}$",
                  "$\\beta_{lambda}$","$\\alpha_{1}$","$\\beta_{1}$","$\\alpha_{a}$","$\\beta_{a}$", 
                  "$N_{3}$","$N_{13}$","$N_{23}$","$N_{33}$")


dimnames(BKM_Resultsparams)[[1]] <- latex_rnames

# Print to LATEX ####
sink("Results/BKM_DA_HMM_params.tex")
cat("\\begin{landscape} \n")
cat("\\begin{table}\n")
cat("\\begin{tabular}{l| rr rr |rr rr rr }\n")

cat("&") 
for (ii in 1:ncol(BKM_Resultsparams)){
  if (ii == ncol(BKM_Resultsparams)){
    cat(dimnames(BKM_Resultsparams)[[2]][ii]," \\\\ \\hline \n ")
  } else {
    cat(dimnames(BKM_Resultsparams)[[2]][ii]," & ")
  }
}
for (jj in 1:nrow(BKM_Resultsparams)){
  cat(dimnames(BKM_Resultsparams)[[1]][jj]," & ")
  for (ii in 1:ncol(BKM_Resultsparams)){
    if (ii == ncol(BKM_Resultsparams)){
      cat(sprintf("%6.4f",BKM_Resultsparams[jj,ii]), " \\\\ \n ")
    } else {
      cat(sprintf("%6.4f",BKM_Resultsparams[jj,ii]), " & ")
    }
  }
}

cat("\\hline \\hline\n")
cat("\\end{tabular}\n")
cat("\\caption{Results for HMM in DA (3 chains, 2000 \\texttt{iter},  100 \\texttt{adapt}, 500 burn-in) 
for the scaled (by 10) lapwings data, against the results for the full DA 
with the original and the scaled (by 10) data (both : 1 chain, 2000 \\texttt{iter}, 100 \\texttt{adapt}, 500 burn-in). 
($Na$ starts at $t=3$ because of different timing of the components of the integrated model.)}\n")
cat("\\end{table}\n")
cat("\\end{landscape} \n")
sink()


# Time ####
Time <- matrix(NaN, nrow = 3, ncol = 2,dimnames=list(c("DA","DA scaled","HMM"),c("Init","Sample")))
load("Results/BKM_model_linux.RData")
Time[1,1] <- time_init[3]
load("Results/BKM_model_scaled_linux.RData")
Time[2,1] <- time_init[3]
load("Results/BKM_HMM_model_ada100_linux.RData")
Time[3,1] <- time_HMM_init[3]
load("Results/BKM_iter2000_ada100_linux.RData")
Time[1,2] <- time_sample[3]
load("Results/BKM_iter2000_ada100_scaled_linux.RData")
Time[2,2] <- time_sample[3]
load("Results/BKM_HMM_try_iter2000_ada100_linux.RData")
Time[3,2] <- time_HMM_sample[3]



sink("Results/BKM_DA_HMM_time.tex")
cat("\\begin{table}\n")
cat("\\begin{tabular}{l| rr}\n")

cat("&") 
for (ii in 1:ncol(Time)){
  if (ii == ncol(Time)){
    cat(dimnames(Time)[[2]][ii]," \\\\ \\hline \n ")
  } else {
    cat(dimnames(Time)[[2]][ii]," & ")
  }
}
for (jj in 1:nrow(Time)){
  cat(dimnames(Time)[[1]][jj]," & ")
  for (ii in 1:ncol(Time)){
    if (ii == ncol(Time)){
      cat(sprintf("%6.4f",Time[jj,ii]), " \\\\ \n ")
    } else {
      cat(sprintf("%6.4f",Time[jj,ii]), " & ")
    }
  }
}

cat("\\hline \\hline\n")
cat("\\end{tabular}\n")
cat("\\caption{Time comparison (\\texttt{JAGS} model initialisation and \\texttt{coda} sampling) 
    for HMM in DA and the full DA with the original and the scaled data.}\n")
cat("\\end{table}\n")
sink()
