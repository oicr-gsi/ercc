# =============================================
#
# Plots for 
#
# 1. dose response curve
#
# 2. Fold changes analysis
#
# To be used in ERCC workflow
#
# Usage:
#
# Rscript ercc_plots.R data[RPKM table] contr[File with ERCC mix info] \
#         prefix[prefix for image] type[DoseResponse or FoldDifference] \
#         MIX1[id of sample with Mix1] MIX2[id of sample with Mix2] \
#         WIDTH[in pixels, optional. default=650] \
#         HEIGHT[in pixels, optional. default=650] \
#         
# =============================================
# /data/Projects/RKT1/ERCC_Controls_Analysis_v2.txt
library(RColorBrewer)

USAGE = list("Rscript ercc.R, order of arguments:\n","1. [Table with RPKMs]\n","2. [File with ERCC mix info]\n",
"3. [prefix for image]\n","4. [DoseResponse or FoldDifference]\n","5. [id of sample with Mix1]\n",
"6. [id of sample with Mix2]\n","7. [width in pixels, optional. default=650]\n",
"8. [height in pixels, optional. default=650]\n")

cmd_args = commandArgs(trailingOnly = TRUE);

if (length(cmd_args) < 6) {
  for (l in USAGE) {
    cat(l)
  }
  quit()
}

DATA<-read.table(cmd_args[1],header=T,check.names=FALSE)
CONTR<-read.table(cmd_args[2], header=TRUE, sep="\t", stringsAsFactors = FALSE)
par(mfrow=c(1,1))
PREFIX = cmd_args[3]
TYPE   = cmd_args[4] # DoseResponse or FoldChange
MIX1   = cmd_args[5]
MIX2   = cmd_args[6]
WIDTH  = 650
HEIGHT = 650

print(paste(WIDTH,HEIGHT,sep="\t"))
# Image dimensions
if (length(cmd_args) == 8) {
 WIDTH  = as.integer(cmd_args[7])
 HEIGHT = as.integer(cmd_args[8])
}
print(paste(WIDTH,HEIGHT,sep="\t"))
TITLE = paste(PREFIX,TYPE,sep="_") 
  
CONTR<-read.table("/data/Projects/RKT1/ERCC_Controls_Analysis_v2.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE)
par(mfrow=c(1,1))

# =======================
#  Dose Response Plotting
# =======================

plotDoseResponse<-function(DATA,mix1,mix2,MAIN) {
  LEGEND = c("Mix 1","Mix 2")
  COL=c("steelblue","salmon")
  if (!any(colnames(DATA) == mix1)) {
    LEGEND = c("Mix 2")
    COL = c("salmon")
  } else if (!any(colnames(DATA) == mix2)) {
    LEGEND = c("Mix 1")
    COL = c("steelblue")
  }

  plot(NULL,
       ylim = c(-3,12),
       xlim = c(-5,25),
       main = MAIN,
       cex.axis = 1.2,
       cex.lab  = 1.2,
       xlab="Log2 ERCC Spike-in amount, amol/ul", 
       ylab="Log2 Normalized ERCC Counts")
  legend(20,2, legend=LEGEND, col=COL, pch = 16)
  P<-NULL
  for (f in c(mix1)) {
    if (any(colnames(DATA) == f)) {
      POINTS<-DATA[,c(1,2,which(colnames(DATA) == f))]
      POINTS$Am = 0
      for (e in DATA$Feature) {
        POINTS$Am[which(POINTS$Feature == e)] = CONTR$concentration.in.Mix.1..attomoles.ul.[which(CONTR$ERCC.ID == e)]
      }
      POINTS$AMLOG<-log2(POINTS$Am)
      POINTS$VCLOG<-log2(POINTS[,which(colnames(POINTS) == f)])
      points(POINTS[,5:6], col = "steelblue")
      P = POINTS[,c(1,5:6)]
    }
  }
  
  for (f in c(mix2)) {
    if (any(colnames(DATA) == f)) {
      POINTS<-DATA[,c(1,2,which(colnames(DATA) == f))]
      POINTS$Am = 0
      for (e in DATA$Feature) {
        POINTS$Am[which(POINTS$Feature == e)] = CONTR$concentration.in.Mix.2..attomoles.ul.[which(CONTR$ERCC.ID == e)]
      }
      POINTS$AMLOG<-log2(POINTS$Am)
      POINTS$VCLOG<-log2(POINTS[,which(colnames(POINTS) == f)])
      points(POINTS[,5:6], col = "salmon")
      if (!is.null(P)) {
        P<-rbind(P,POINTS[,c(1,5:6)])
      } else {
        P<-POINTS[,c(1,5:6)]
      }
    }
  }
  
  
  P[which(P$AMLOG==-Inf),2] = NA
  P[which(P$VCLOG==-Inf),3] = NA
  L<-with(P, lm(VCLOG~AMLOG))
  
  abline(L$coefficients[1], L$coefficients[2] , col="red", lty = 2) # regression line (y~x)
  # 1 RPKM line:
  abline(h = 0, col="gray60", lty = 2)
  R2<-cor(P$AMLOG,y=P$VCLOG, method = "pearson", use = "complete.obs")^2
  text(-5, 1, "1RPKM", cex = 1.1, pos = 4)
  text(-5, 9, paste("n = ", length(unique(P$Feature))), cex = 1.1, pos = 4)
  text(-5, 8, paste("R",toupper("2")," = ",round(R2,4),sep=""), cex = 1.1, pos = 4)
  text(-5, 7, paste(labels="Slope = ",round(L$coefficients[2],4)), cex = 1.1, pos = 4)
  
}


# ==============================================================================
#
#    Fold Change Response Analysis, need two samples here
#
# ==============================================================================
plotFoldChange<-function(DATA,mix1,mix2,MAIN="Fold Change Response") {
  
  plot(NULL, 
       ylim = c(-4,4),
       xlim = c(-4,4),
       main = MAIN,
       cex.axis = 1.2,
       cex.lab  = 1.2,
       xlab="Expected Log2 Ratio", 
       ylab="Observed Log2 Ratio")
  
  GROUPCOLS<-c(A="salmon",B="green",C="steelblue",D="darkmagenta")
  for (f in 1:length(mix1)) {
    
    POINTS<-DATA[,c(1,which(colnames(DATA) == mix1[f]),which(colnames(DATA) == mix2[f]))]
    POINTS$ExpFoldLog2 = 0
    POINTS$FoldLog2    = log2(POINTS[,2]/POINTS[,3])
    POINTS$Col = "black"
    
    for (e in DATA$Feature) {
      POINTS$ExpFoldLog2[which(POINTS$Feature == e)] = CONTR$log2.fold.change.[which(CONTR$ERCC.ID == e)]
      POINTS$Col[which(POINTS$Feature == e)] = GROUPCOLS[CONTR$subgroup[which(CONTR$ERCC.ID == e)]]
    }
    points(POINTS[,4:5], col = POINTS$Col, pch = 19)
  }
  
  # Indexes are very important! 
  POINTS[which(POINTS$FoldLog2==-Inf),5] = NA
  POINTS[which(POINTS$ExpFoldLog2==-Inf),4] = NA
  POINTS[which(POINTS$FoldLog2==Inf),5] = NA
  POINTS[which(POINTS$ExpFoldLog2==Inf),4] = NA
  POINTS[which(POINTS$FoldLog2=="NaN"),5] = NA
  POINTS[which(POINTS$ExpFoldLog2=="NaN"),4] = NA
  
  L<-with(POINTS, lm(FoldLog2~ExpFoldLog2, na.action = na.exclude))
  Icept = L$coefficients[1]
  # regression line (y~x)
  abline(L$coefficients[1], L$coefficients[2] , col="red", lty = 2) 
  # 1 RPKM line:
  abline(h = 0, col="gray60", lty = 2)
  abline(v = 0, col="gray60", lty = 2)
  R2<-cor(POINTS$FoldLog2,y=POINTS$ExpFoldLog2, method = "pearson", use = "complete.obs")^2
  
  
  text(1.0, -1.5, paste("n = ",length(unique(POINTS$Feature))), cex = 1.1, pos = 4)
  text(1.0, -2, bquote(R^2 == .(round(R2, 4))), cex = 1.1, pos = 4)
  text(1.0, -2.5, paste(labels="Slope = ",round(L$coefficients[2],4)), cex = 1.1, pos = 4)
  text(1.0, -3, paste(labels="Icept = ",round(L$coefficients[1],4)), cex = 1.1, pos = 4)
  
  legend(-4, 4, legend = names(GROUPCOLS), col = GROUPCOLS, pch = 19)
  
}

if (TYPE == "DoseResponse") {
  png(filename=paste0(TITLE,".png"),width=WIDTH,height=HEIGHT,units="px",pointsize=15,bg="white",type='cairo')
  plotDoseResponse(DATA,MIX1,MIX2,TITLE)
  dev.off()
  
} else if (TYPE == "FoldDifference") {
  png(filename=paste0(TITLE,".png"),width=WIDTH,height=HEIGHT,units="px",pointsize=15,bg="white",type='cairo')
  plotFoldChange(DATA,MIX1,MIX2,TITLE)
  dev.off()
}
