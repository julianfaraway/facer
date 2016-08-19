#' Read in a face file
#'
#' @param fname The name of the function
#' @param mfconfig motion file configuration
#' @param numskip The number of lines to skip at the beginning of the motion file
#' @return a data matrix of the motion file
#'
#' @export
read_face_file <- function(fname,mfconfig,numskip=5){
  sm <- read.table(fname,skip=numskip,sep="\t",colClasses=rep("numeric",3*mfconfig$nmarkers+2),
                   flush=TRUE,fill=TRUE)[,3:(3*mfconfig$nmarkers+2)]
  sm = data.matrix(sm)
  shapex <- array(NA, dim = c(mfconfig$nmarkers, 3, mfconfig$nframes))

  for(i in 1:mfconfig$nframes){
    shapex[,,i] <- matrix(sm[i,],nrow=mfconfig$nmarkers, ncol=3,byrow=TRUE)
  }
  shapex
}

#' extract the first PCA, select start and max, return initial and max frames
#'
#' @param sm The motion data array
#'
#' @return list containing the first PC score trajectory, the frame number where the maximum is obtained,
#' the initial frame, the the frame at the max, the geodesic initial frame and the geodesic max frame
#' @export
#'
first_PCA_traj = function(sm){
  dimsm = dim(sm)
  n <- dimsm[3]
  k <- dimsm[1]
  m <- dimsm[2]

  spca <- try(procGPA(sm, scale=TRUE, distances=FALSE))
  score1 = spca$scores[,1]
  if(score1[1] > 0) score1 = -score1
  framemax = which.max(score1)

  list(score1=score1,
       framemax = framemax,
       rawstart = sm[,,1],
       rawext = sm[,,framemax],
       geostart = spca$mshape+matrix(spca$scores[1,1]*spca$pcar[,1],nrow=k,ncol=m),
       geoext = spca$mshape+matrix(spca$scores[framemax,1]*spca$pcar[,1],nrow=k,ncol=m))
}

#' Plot a face or a pair of faces
#'
#' @param fm A face matrix
#' @param mfconfig motion capture config file
#' @param fm2 Another face matrix (optional)
#' @param view direction of the 2D view
#' @param num Marker numbers, rather than points, to be displayed
#' @param main option title
#' @param dirm when two faces are plotted, what should connect the corresponding markers
#'
#' @return nothing - a plot is constructed
#' @export
#'
plotface <- function(fm,mfconfig,fm2=NULL,view=c("front","side","top"), num=TRUE, main="",
                     dirm=c("arrows","points","segments")){
  n <- nrow(fm)
  dirm <- match.arg(dirm)
  view = match.arg(view)
  orthproj = switch(view, "front" = mfconfig$coordir[c(1,2)],
                          "side" = mfconfig$coordir[c(3,2)],
                          "top" = mfconfig$coordir[c(1,3)])

  latdir = mfconfig$coordir[1]
  if(latdir %in% orthproj){ # is lateral direction in selected view
    selpts = 1:n
  }else{ # which markers are on the right
    rightside = which(fm[,latdir] < mean(fm[mfconfig$midline,latdir]))
    rightside = union(rightside, mfconfig$midline)
    selpts = rightside
  }

  afm <- rbind(fm[selpts,],fm2[selpts,])
  xr <- range(afm[,orthproj[1]])
  yr <- range(afm[,orthproj[2]])

  maxsp <- max(diff(xr),diff(yr))
  nxlim <- mean(xr)+c(-1,1)*0.5*maxsp
  nylim <- mean(yr)+c(-1,1)*0.5*maxsp

  if(is.null(fm2)){
    if(num){
      plot(fm[selpts,orthproj[1]],fm[selpts,orthproj[2]],type="n",xlab="",ylab="",main=main,xlim=nxlim,ylim=nylim,asp=1)
      text(fm[selpts,orthproj[1]],fm[selpts,orthproj[2]],as.character(selpts))
    }else{
      plot(fm[selpts,orthproj[1]],fm[selpts,orthproj[2]],xlab="",ylab="",main=main,xlim=nxlim,ylim=nylim,asp=1)
    }
  }else{
    plot(fm[selpts,orthproj[1]],fm[selpts,orthproj[2]],xlab="",ylab="",main=main,xlim=nxlim,ylim=nylim,asp=1)
    if(dirm == "arrows"){
      arrows(fm[selpts,orthproj[1]],fm[selpts,orthproj[2]],
             fm2[selpts,orthproj[1]],fm2[selpts,orthproj[2]],length=0.1)
    }
    if(dirm == "points"){
      points(fm2[selpts,orthproj[1]],fm2[selpts,orthproj[2]],pch=20)
    }
    if(dirm == "segments"){
      segments(fm[selpts,orthproj[1]],fm[selpts,orthproj[2]],fm2[selpts,orthproj[1]],fm2[selpts,orthproj[2]])
    }
  }
}

#' Construct a single frame suitable for movie construction
#'
#' @param fc1 A face
#' @param mfconfig face configuration
#' @param fc2 Another face (optional)
#' @param box use boxes (default is FALSE)
#' @param header optional title for plot
#'
#' @return nothing - a plot is created
#' @export
faceframe <- function(fc1,mfconfig,fc2=NULL,rang=NULL,box=FALSE,header=FALSE){

  n <- mfconfig$nmarkers
  latdir = mfconfig$coordir[1]
  vertdir = mfconfig$coordir[2]
  frontdir = mfconfig$coordir[3]

  rightside = which(fc1[,latdir] < mean(fc1[mfconfig$midline,latdir]))
  leftside = setdiff(1:n, rightside)
  rightside = union(rightside, mfconfig$midline)
  leftside = union(leftside, mfconfig$midline)

  if(is.null(rang)) rang = apply(fc1,2,range) # should be set globally

  par(mar=rep(0.1,4))
  layout(matrix(c(1,2,3),ncol=3))

  if(!missing(header)) par(oma=c(0,0,2,0))

  orthoproj <- c(frontdir,vertdir)
  selpts <- rightside

  xr <- rang[,orthoproj[1]]
  yr <- rang[,orthoproj[2]]

  plot(fc1[selpts,orthoproj[1]],fc1[selpts,orthoproj[2]],ann=box,axes=FALSE,xlim=c(xr[2],xr[1]),ylim=yr,pch=1,frame.plot=box)
  if(!is.null(fc2)) points(fc2[selpts,orthoproj[1]],fc2[selpts,orthoproj[2]],pch=16,col=2)


  orthoproj <- c(latdir,vertdir)
  selpts <- 1:n

  xr <- rang[,orthoproj[1]]
  yr <- rang[,orthoproj[2]]

  plot(fc1[selpts,orthoproj[1]],fc1[selpts,orthoproj[2]],
       ann=box,axes=FALSE,xlim=xr,ylim=yr,pch=1, frame.plot=box)
  if(!is.null(fc2)) points(fc2[selpts,orthoproj[1]],fc2[selpts,orthoproj[2]],pch=16,col=2)

  if(!missing(header)){
    mtext(header,outer=TRUE,cex=1.5)
  }


  orthoproj <- c(frontdir,vertdir)
  selpts <- leftside

  xr <- rang[,orthoproj[1]]
  yr <- rang[,orthoproj[2]]

  plot(fc1[selpts,orthoproj[1]],fc1[selpts,orthoproj[2]],
       ann=box,axes=FALSE,xlim=xr,ylim=yr,pch=1,frame.plot=box)
  if(!is.null(fc2)) points(fc2[selpts,orthoproj[1]],fc2[selpts,orthoproj[2]],pch=16,col=2)
}

#' Extract information about the motion capture file
#'
#' @param fname A name of a typical motion capture file
#' @param nskip Number of lines at top containing the info
#'
#' @return list containing number of frames, number of markers, names of the
#'   markers, the markers which lie on the midline of the face and the
#'   coordinate directions as (lateral, vertical, frontal)
#' @export
extmfconfig = function(fname, nskip=5){
  hs = readLines(fname,n=nskip)
  l3 = strsplit(hs[3],"\t")[[1]]
  nframes = as.numeric(l3[3])
  nmarkers = as.numeric(l3[4])
  marklabs = strsplit(hs[4],"\t")[[1]][-c(1,2)]
  if(length(marklabs)/3 != nmarkers) warning("number of marker labels not equal to number of markers")
  marklabs = marklabs[3*(1:nmarkers)-2]
  fl = substring(marklabs,1,1)
  midline = which(!(fl %in% c("r","l","R","L"))) # guess which are midline markers but needs checking
  coordir = 1:3 # needs to be  checked by user
  list(nframes=nframes, nmarkers=nmarkers, marklabs=marklabs,midline=midline,coordir=coordir)
}

#' Make a facial motion movie
#'
#' @param shapem1 A shape motion
#' @param mfconfig motion config
#' @param shapem2 another shape motion (optional)
#' @param everynth sample every nth frame (default is 6)
#' @param adjust adjust second motion to start from same place as first motion
#' @param moviename name of the movie (default is "movie.mov")
#' @param header title displayed on movie
#'
#' @return nothing - find movie file in "animations" directory
#' @export
facemovie <- function(shapem1,mfconfig,shapem2=NULL,everynth=6,adjust=TRUE,moviename="movie.mov",header=NULL){

  nframes = mfconfig$nframes
  if(!missing(shapem2)){ # if two motions, attempt to coordinate them
    combf <- procGPA(abind(shapem1,shapem2))
    shapem1 <- combf$rotated[,,1:nframes]
    shapem2 <- combf$rotated[,,(nframes+1):(2*nframes)]
  }

  k <- mfconfig$nmarkers
  m <- 3
  nframes <- mfconfig$nframes

  isel <- seq(1,nframes,by=everynth) # winnow frames
  animframes <- length(isel)
  shapem1 <- shapem1[,,isel]
  if(!is.null(shapem2)) shapem2 <- shapem2[,,isel]

  if(adjust & !is.null(shapem2)){ # adjust to start from same spot
    adjm <- (shapem2[,,1] - shapem1[,,1])
    for(i in 1:dim(shapem1)[3]) shapem2[,,i] <- shapem2[,,i] - adjm
  }

  if(!is.null(shapem2)){
    rangeinfo <- apply(abind(shapem1,shapem2),2,range)
  }else{
    rangeinfo = apply(shapem1,2,range)
  }

  maxrang <- max(rangeinfo[2,]-rangeinfo[1,])
  for(i in 1:3){
    adj <- (rangeinfo[2,i]-rangeinfo[1,i] - maxrang)/2
    rangeinfo[1,i] <- rangeinfo[1,i]+adj
    rangeinfo[2,i] <- rangeinfo[2,i]-adj
  }

  system("rm -f movie/*")

  if(missing(header)){
    pnght <- 240
  }else{
    pnght <- 270
  }


  for(i in 1:animframes){
    png(filename=paste("movie/face",i+100,".png",sep=""),width=480,height=pnght)
    if(missing(shapem2)){
      faceframe(shapem1[,,i],mfconfig,rang=rangeinfo,header=header)
    }else{
      faceframe(shapem1[,,i],mfconfig,shapem2[,,i],rang=rangeinfo,header=header)
    }
    dev.off()
  }

  fps = animframes/(nframes/60)

  system(paste("rm -f animations/",moviename,sep=""))
  system(paste("~/bin/crtimgseq.py animations/",moviename," 1 ",fps," movie/*",sep=""))
}

#' Average faces and rotate so head is upright
#'
#' @param shapem An array of facial shapes
#' @param mconfig configuration info
#'
#' @return A facial shape matrix rotated to have minimum variance in frontal and lateral directions
#' @export
averecface = function(shapem, mconfig){
  latdir = mfconfig$coordir[1]
  vertdir = mfconfig$coordir[2]
  frontdir = mfconfig$coordir[3]
  mshape = procGPA(shapem, pcaoutput = FALSE, distances = FALSE)$mshape
  # rotate in the saggital plane
  aseq = seq(-pi/4, pi/4, length=1000)
  sdv = numeric(length(aseq))
  for(j in seq_along(aseq)){
    rshape = rotface(mshape, aseq[j], c(vertdir, frontdir))
    sdv[j] = sd(rshape[,frontdir])
  }
  bestang = aseq[which.min(sdv)]
  mshape = rotface(mshape, bestang, c(vertdir, frontdir))

  # rotate in the horizontal plane
  aseq = seq(-pi/4, pi/4, length=1000)
  sdv = numeric(length(aseq))
  for(j in seq_along(aseq)){
    rshape = rotface(mshape, aseq[j], c(latdir, frontdir))
    sdv[j] = sd(rshape[,frontdir])
  }
  bestang = aseq[which.min(sdv)]
  mshape = rotface(mshape, bestang, c(latdir, frontdir))

  # rotate in the frontal plane
  aseq = seq(-pi/4, pi/4, length=1000)
  for(j in seq_along(aseq)){
    rshape = rotface(mshape, aseq[j], c(latdir, vertdir))
    sdv[j] = sd(rshape[,latdir])+sd(rshape[,vertdir])
  }
  bestang = aseq[which.min(sdv)]
  rotface(mshape, bestang, c(latdir, frontdir))
}

#' Rotate a face matrix in SD
#'
#' @param sm The face matrix
#' @param theta Angle of rotation
#' @param dirs Coordinate axes of plane of rotation - first dir is theta=0
#'
#' @return a face matrix
#' @export
rotface = function(sm,theta,dirs){
  mr = diag(3)
  mr[dirs[1],dirs[1]] = cos(theta)
  mr[dirs[2],dirs[2]] = cos(theta)
  mr[dirs[1],dirs[2]] = -sin(theta)
  mr[dirs[2],dirs[1]] = sin(theta)
  sm %*% mr
}

#' Fill in missing values using linear interpolation
#'
#' @param x 1D motion trace
#'
#' @return 1D motion trace
#' @export
misfixlin <- function(x){
  mv <- is.na(x)
  if(sum(mv) == 0) return(x) # if nothing missing, do nothing
  if(sum(mv) > length(x)-2) return(x) # if everything missing, return missing
  n <- length(x)
  approx(1:n,x,rule=2,xout=1:n)$y
}

#' Process a facial motion to fill-in missing values and report diagnostics
#'
#' @param sm a facial motion
#' @param MISSMAX maximum allowable number of missing cases in a marker (default=50)
#'
#' @return a list containing a facial motion with missing values filled in, the number of
#' missing cases in each marker, the max jump in a marker between frames
#' @export
cleansm = function(sm,MISSMAX=50){
  nmiss = apply(sm,1,function(x) sum(is.na(x)))/3
  nobs = dim(sm)[3]
  missmark = which(nmiss == nobs)
  if(any(nmiss > MISSMAX)){
    warning(paste("More than",MISSMAX,"missing values in marker(s):", which(nmiss > MISSMAX)))
  }
  fm = aperm(apply(sm,c(1,2),misfixlin),c(2,3,1))
  dm = abs(fm[,,-nobs] - fm[,,-1])
  sc = max(apply(dm, 1, sum))/3
  list(sm=fm,nmiss=nmiss,jump=sc,allmiss=missmark)
}
