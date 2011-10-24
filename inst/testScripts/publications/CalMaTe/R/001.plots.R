###########################################################################
# 001.PLOTS.R
#
# Author: Henrik Bengtsson
###########################################################################

plotBvsB <- function(beta, col=NULL, pch=19, Blim=c(-0.2,+1.2), xlab="Normal BAF", ylab="Tumor BAF", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Blim, ylim=Blim, xlab=xlab, ylab=ylab, axes=FALSE);
  for (ss in 1:2) {
    axis(side=ss, at=c(0,1/2,1), label=c("0","1/2","1"));
  }
  box();
  abline(a=0, b=1, lty=3, lwd=2);

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  stext(side=4, pos=0, line=0, sprintf("n=%d", nrow(beta)), cex=1.5);

  points(beta, pch=pch, col=col, cex=1);

  if (is.element("dens", anTags)) {
    for (dd in 1:2) {
      d <- density(beta[,dd], adjust=0.4, from=Blim[1], to=Blim[2], na.rm=TRUE);
      draw(d, lwd=3, col="#666666", side=dd, 
                       height=0.05, scale="relative", xpd=FALSE);
    }
  }
} # plotBvsB()



plotCvsB <- function(ascn, col=NULL, pch=19, Blim=c(-0.2,+1.2), Clim=c(-0.4,4), xlab="BAF", ylab="Total CN", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
#  par(mar=c(4,5,3,3)+0.1, cex=0.8, cex.lab=2.4, cex.axis=2.2);
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Blim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1, at=c(0,1/2,1), label=c("0","1/2","1"));
  axis(side=2, at=0:round(max(Clim)+1L));
  box();

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  stext(side=4, pos=0, line=0, sprintf("n=%d", nrow(ascn)), cex=1.5);

  bc <- ascn;
  bc[,2] <- ascn[,1] + ascn[,2];
  bc[,1] <- ascn[,2] / bc[,2];
  points(bc, pch=pch, col=col, cex=1);

  if (is.element("dens", anTags)) {
    for (dd in 1:2) {
      d <- density(bc[,dd], adjust=0.4, from=Clim[1], to=Clim[2], na.rm=TRUE);
      draw(d, lwd=3, col="#666666", side=dd, 
                       height=0.05, scale="relative", xpd=FALSE);
    }
  }
} # plotCvsB()



plotASCN <- function(ascn, col=NULL, pch=19, Clim=c(-0.4,4), xlab=expression(C[A]), ylab=expression(C[B]), dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
#  par(mar=c(4,5,3,3)+0.1, cex=0.8, cex.lab=2.4, cex.axis=2.2);
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  for (ss in 1:2) {
    axis(side=ss, at=0:round(max(Clim)+1L));
  }
  box();
  abline(a=0, b=1, lty=3, lwd=2);

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }

  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  stext(side=4, pos=0, line=0, sprintf("n=%d", nrow(ascn)), cex=1.5);


  points(ascn, pch=pch, col=col, cex=1);

  if (is.element("dens", anTags)) {
    for (dd in 1:2) {
      d <- density(ascn[,dd], adjust=0.4, from=Clim[1], to=Clim[2], na.rm=TRUE);
      draw(d, lwd=3, col="#666666", side=dd, 
                       height=0.05, scale="relative", xpd=FALSE);
    }
  }
} # plotASCN()



##########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Added option "auto" for argument 'reference' of extractSignals().
# o Created.
###########################################################################
