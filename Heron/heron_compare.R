load("Results/Heron_HMMnorm_ESS.RData")
load("Results/Heron_DA_ESS.RData")


# names_all <- c("DA","Q1=10, Q3=5","Q1=50, Q3=40","Q1=100, Q3=70") 
names_all <- c("DA","DA",
  "Q1=10, Q3=5","Q1=10, Q3=5",
  "Q1=50, Q3=40","Q1=50, Q3=40",
  "Q1=100, Q3=70","Q1=100, Q3=70")

myColors <- c("black","gray",
              "blue","skyblue3",
              "green","limegreen",
              "red","magenta")

ind = c(1,3,5,7)

# ESS ####
# BARPLOT ESS X2
plotname = "Figures/heron_COMP_ESS_X2.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

barplot(t(cbind(ESS_DA_X2[2:72,], ESS_HMM_X2[2:72,])[,ind]), beside= TRUE, 
        legend =  names_all[ind], col = myColors[ind])
dev.off()

# BARPLOT ESS X4
plotname = "Figures/heron_COMP_ESS_X4.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

barplot(t(cbind(ESS_DA_X4[2:72,], ESS_HMM_X4[2:72,])[,ind]), beside= TRUE, 
        legend =  names_all[ind], col = myColors[ind])
dev.off()

# BARPLOT ESS param
plotname = "Figures/heron_COMP_ESS_param.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

barplot(t(cbind(ESS_DA_param, ESS_HMM_param)[,ind]), beside= TRUE, 
        legend =  names_all[ind], col = myColors[ind])
dev.off()

# ESS PER SEC ####
# BARPLOT ESS PER SEC X2
plotname = "Figures/heron_COMP_ESS_persec_X2.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

barplot(t(cbind(ESS_persec_DA_X2[2:72,], ESS_persec_HMM_X2[2:72,])[,ind]), beside= TRUE, 
        legend =  names_all[ind], col = myColors[ind])
dev.off()

# BARPLOT ESS PER SEC X4
plotname = "Figures/heron_COMP_ESS_persec_X4.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

barplot(t(cbind(ESS_persec_DA_X4[2:72,], ESS_persec_HMM_X4[2:72,])[,ind]), beside= TRUE, 
        legend =  names_all[ind], col = myColors[ind])
dev.off()

# BARPLOT ESS PER SEC param
plotname = "Figures/heron_COMP_ESS_persec_param.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

barplot(t(cbind(ESS_persec_DA_param, ESS_persec_HMM_param)[,ind]), beside= TRUE, 
        legend =  names_all[ind], col = myColors[ind])
dev.off()
