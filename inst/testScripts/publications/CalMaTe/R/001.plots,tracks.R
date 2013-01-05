###########################################################################
# 001.PLOTS.R
#
# Author: Henrik Bengtsson
###########################################################################

plotTrackC1C2 <- function(ascn, x, muN, smooth=NULL, pch=".", cex=1, lwd=3, xlim=NULL, Clim=c(-0.4,3.4), xlab="Position (Mb)", ylab=expression(list(C[1],C[2])), dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
##   if (is.null(xlim)) {
##     xlim <- range(x, na.rm=TRUE);
##     dx <- diff(xlim);
##     xlim[1] <- xlim[1] - 0.05*dx;
##     xlim[2] <- xlim[2] + 0.05*dx;
##   }

#  par(mar=c(4,5,3,3)+0.1, cex=0.8, cex.lab=2.4, cex.axis=2.2);
#  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  par(mar=c(2.5,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
#  plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  plot(x, ascn[,1], type="n", xlim=xlim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=seq(from=0, to=as.integer(Clim[2]+1L)));
  box();

  if (!is.null(segName)) {
    if (segName == "mix") {
      segName <- NULL;
    }
  }

  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segTag)) {
    stext(side=3, pos=0, paste(c(segName, segTag), collapse=" @ "), cex=2.0);
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

  hets <- which(muN == 1/2);
  xH <- x[hets];
  c1c2 <- ascn[hets,];
  rr <- which(c1c2[,2] < c1c2[,1]);
  c1c2[rr,] <- c1c2[rr,2:1];

  c1 <- c1c2[,1];
  c2 <- c1c2[,2];

  colC1C2 <- c(blue="#0000ff", red="#ff0000");
  colC1C2b <- sprintf("%s66", colC1C2);
  points(xH, c1, pch=pch, cex=cex, col=colC1C2b[1]);
  points(xH, c2, pch=pch, cex=cex, col=colC1C2b[2]);

  if (!is.null(segName)) {
    abline(h=median(c1,na.rm=TRUE), lwd=lwd, col=colC1C2[1]);
    abline(h=median(c2,na.rm=TRUE), lwd=lwd, col=colC1C2[2]);

    if (is.element("dens", anTags)) {
      d <- density(c1, adjust=0.8, na.rm=TRUE);
      draw(d, lwd=3, col="blue", side=2, height=0.05, scale="relative", xpd=FALSE);
  
      d <- density(c2, adjust=0.8, na.rm=TRUE);
      draw(d, lwd=3, col="red", side=2, height=0.05, scale="relative", xpd=FALSE);
    }
  }

  by <- NULL;
  if (is.element("1Mb", c(smooth, anTags))) {
    by <- 1;
  } else if (is.element("2Mb", c(smooth, anTags))) {
    by <- 2;
  } else if (is.element("5Mb", c(smooth, anTags))) {
    by <- 5;
  }

  if (!is.null(by)) {
    c1s <- binnedSmoothing(c1, xH, by=by);
    c2s <- binnedSmoothing(c2, xH, by=by);
    xOut <- attr(c1s, "xOut");
    lwd <- 4;
    lines(xOut, c1s, lwd=lwd+4, col="white");
    lines(xOut, c2s, lwd=lwd+4, col="white");
    lines(xOut, c1s, lwd=lwd, col="blue");
    lines(xOut, c2s, lwd=lwd, col="red");
  }
} # plotTrackC1C2()


plotTrackTCN <- function(C, x, col=NULL, pch=".", cex=1, xlim=NULL, Clim=c(-0.4,5.4), ylab="Total CN", xlab="Position (Mb)", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(2.5,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
  plot(x, C, pch=pch, cex=cex, col=col, xlim=xlim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=seq(from=0, to=as.integer(Clim[2]+1L)));
  box();

  if (!is.null(segName)) {
    if (segName == "mix") {
      segName <- NULL;
    }
  }
  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segTag)) {
    stext(side=3, pos=0, paste(c(segName, segTag), collapse=" @ "), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName, cex=1.5, line=-0.2);
  }
} # plotTrackTCN()



plotTrackBAF <- function(B, x, col=NULL, pch=".", cex=1, xlim=NULL, Blim=c(-0.2,+1.2), ylab="BAF", xlab="Position (Mb)", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(2.5,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
  plot(x, B, pch=pch, cex=cex, col=col, xlim=xlim, ylim=Blim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=c(0, 1/2, 1), label=c(0, "1/2", 1));
  box();

  if (!is.null(segName)) {
    if (segName == "mix") {
      segName <- NULL;
    }
  }
  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segTag)) {
    stext(side=3, pos=0, paste(c(segName, segTag), collapse=" @ "), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName, cex=1.5, line=-0.2);
  }
} # plotTrackBAF()



plotTrackDH <- function(B, x, muN, col=NULL, pch=".", cex=1, xlim=NULL, Blim=c(0,+1.2), ylab="DH", xlab="Position (Mb)", dataSet=NULL, chipType=NULL, sampleName=NULL, chrTag=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(2.5,2.5,1.1,1)+0.1, tcl=-0.3, mgp=c(1.4,0.4,0), cex=2);
  DH <- 2*abs(B - 1/2);
  DH[muN != 1/2] <- NA;
  plot(x, DH, pch=pch, cex=cex, col=col, xlim=xlim, ylim=Blim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=c(0, 1), label=c(0, 1));
  box();

  if (!is.null(segName)) {
    if (segName == "mix") {
      segName <- NULL;
    }
  }
  # Text annotations
  if (!is.null(chrTag)) {
    stext(side=3, pos=0, chrTag,  cex=1.5);
  } else if (!is.null(segTag)) {
    stext(side=3, pos=0, paste(c(segName, segTag), collapse=" @ "), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName, cex=1.5, line=-0.2);
  }
} # plotTrackDH()



##########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Added option "auto" for argument 'reference' of extractSignals().
# o Created.
###########################################################################
