###########################################################################
# 001.PLOTS.R
#
# Author: Henrik Bengtsson
###########################################################################
plotMultiArrayCACB <- function(cacb, pch=19, cex=2, ..., Clim=c(-0.2,3.4), xlab=expression(C[A]), ylab=expression(C[B]), dataSet=NULL, chipType=NULL, unitName=NULL, tagsT=NULL) {
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  for (ss in 1:2) {
    axis(side=ss, at=0:round(max(Clim)+1L));
  }
  box();
  lines(x=c(0,2), y=c(2,0), lty=2, lwd=2, col="#999999");
  lines(x=c(0,Clim[2]), y=c(0,Clim[2]), lty=2, lwd=2, col="#999999");
  lines(x=c(0,0), y=c(0,Clim[2]), lty=2, lwd=2, col="#999999");
  lines(x=c(0,Clim[2]), y=c(0,0), lty=2, lwd=2, col="#999999");

  if (!is.null(unitName)) {
    stext(side=3, pos=0, sprintf("Unit: %s", unitName), cex=1.5);
  }

  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }

  if (!is.null(dataSet) & !is.null(chipType)) {
    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }

  stext(side=4, pos=1, sprintf("n=%d", nrow(cacb)), cex=1.5);

  points(cacb, pch=pch, cex=cex, ...);
} # plotMultiArrayCACB()


##########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Added option "auto" for argument 'reference' of extractSignals().
# o Created.
###########################################################################
