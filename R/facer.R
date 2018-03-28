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
  sm <-  data.matrix(sm)
  nfr <- nrow(sm)
  shapex <- array(NA, dim = c(mfconfig$nmarkers, 3, nfr))

  for(i in 1:nfr){
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
#' @param markerlabs optional marker labels
#' @param dirm when two faces are plotted, what should connect the corresponding markers
#' @param subset of markers to be shown
#'
#' @return nothing - a plot is constructed
#' @export
#'
plotface <- function(fm,mfconfig,fm2=NULL,view=c("front","side","top"), num=TRUE, main="",markerlabs=NULL,
                     dirm=c("arrows","points","segments"),subset=1:nrow(fm)){
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
  selpts = intersect(selpts,subset)
  if(is.null(markerlabs)) markerlabs = as.character(1:n)

  afm <- rbind(fm[selpts,],fm2[selpts,])
  xr <- range(afm[,orthproj[1]],na.rm=TRUE)
  yr <- range(afm[,orthproj[2]],na.rm=TRUE)

  maxsp <- max(diff(xr),diff(yr))
  nxlim <- mean(xr)+c(-1,1)*0.5*maxsp
  nylim <- mean(yr)+c(-1,1)*0.5*maxsp

  if(is.null(fm2)){
    if(num){
      plot(fm[selpts,orthproj[1]],fm[selpts,orthproj[2]],type="n",xlab="",ylab="",main=main,xlim=nxlim,ylim=nylim,asp=1)
      text(fm[selpts,orthproj[1]],fm[selpts,orthproj[2]],markerlabs[selpts])
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

  if(is.null(mfconfig$lrmat)){
    rightside = which(fc1[,latdir] < mean(fc1[mfconfig$midline,latdir]))
    leftside = setdiff(1:n, rightside)
  }else{
    rightside = mfconfig$lrmat[,1]
    leftside = mfconfig$lrmat[,2]
  }
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
#' @param missmat 0/1 nframes by nmarkers matrix with 1's for missing values
#'
#' @return nothing - find movie file in "animations" directory. Second sequence is in red
#' @export
facemovie <- function(shapem1,mfconfig,shapem2=NULL,everynth=6,adjust=TRUE,moviename="movie.mov",header=NULL,missmat=NULL){

  k <- mfconfig$nmarkers
  m <- 3
  nframes <- mfconfig$nframes

  if(dim(shapem1)[1] != k) stop("Number of markers inconsistent")
  if(dim(shapem1)[3] != nframes) stop("Number of frames inconsistent")

  if(adjust & !missing(shapem2)){ # if two motions, attempt to coordinate them
    combf <- faceGPA(abind(shapem1,shapem2))
    shapem1 <- combf$rotated[,,1:nframes]
    shapem2 <- combf$rotated[,,(nframes+1):(2*nframes)]
  }


  if(!dir.exists("movpngs")) dir.create("movpngs")
  if(!dir.exists("animations")) dir.create("animations")

  isel <- seq(1,nframes,by=everynth) # winnow frames
  animframes <- length(isel)
  shapem1 <- shapem1[,,isel]
  if(!is.null(shapem2)) shapem2 <- shapem2[,,isel]
  if(!is.null(missmat)) missmat = missmat[isel,]

  if(adjust & !is.null(shapem2)){ # adjust to start from same spot
    adjm <- (shapem2[,,1] - shapem1[,,1])
    for(i in 1:dim(shapem1)[3]) shapem2[,,i] <- shapem2[,,i] - adjm
  }

  if(!is.null(shapem2)){
    rangeinfo <- apply(abind(shapem1,shapem2),2,range,na.rm=TRUE)
  }else{
    rangeinfo = apply(shapem1,2,range,na.rm=TRUE)
  }

  maxrang <- max(rangeinfo[2,]-rangeinfo[1,])
  for(i in 1:3){
    adj <- (rangeinfo[2,i]-rangeinfo[1,i] - maxrang)/2
    rangeinfo[1,i] <- rangeinfo[1,i]+adj
    rangeinfo[2,i] <- rangeinfo[2,i]-adj
  }

  system("rm -f movpngs/*")

  if(missing(header)){
    pnght <- 240
  }else{
    pnght <- 270
  }


  for(i in 1:animframes){
    png(filename=paste("movpngs/face",i+100,".png",sep=""),width=480,height=pnght)
    if(!is.null(missmat)){
      shapem1[missmat[i,]==1,,i] = NA
    }
    if(missing(shapem2)){
      faceframe(shapem1[,,i],mfconfig,rang=rangeinfo,header=header)
    }else{
      faceframe(shapem1[,,i],mfconfig,shapem2[,,i],rang=rangeinfo,header=header)
    }
    dev.off()
  }

  fps = animframes/(nframes/60)

  system(paste("rm -f animations/",moviename,sep=""))
  system(paste("~/bin/crtimgseq.py animations/",moviename," 1 ",fps," movpngs/*",sep=""))
  system("rm -f movpngs/*")

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
#' missing cases in each marker, the max jump in a marker between frames, which cases
#' are all missing and a matrix of where the missing values occur
#' @export
cleansm = function(sm,MISSMAX=50){
  missmat = t(apply(sm,c(1,3),function(x) sum(is.na(x)))/3)
  nmiss = colSums(missmat)
  nobs = dim(sm)[3]
  missmark = which(nmiss == nobs)
  if(any(nmiss > MISSMAX)){
    warning(paste("More than",MISSMAX,"missing values in marker(s):", which(nmiss > MISSMAX)))
  }
  fm = aperm(apply(sm,c(1,2),misfixlin),c(2,3,1))
  dm = abs(fm[,,-nobs] - fm[,,-1])
  sc = max(apply(dm, 1, sum))/3
  list(sm=fm,nmiss=nmiss,jump=sc,allmiss=missmark,missmat=missmat)
}

#' Read in an LMK file
#'
#' @param filename - name of the file
#'
#' @return list containing the short labels, long labels and matrix of landmarks
#' @export
readlmk <- function(filename){
  ff = readLines(filename,warn=FALSE)
  n = length(ff)
  ff = ff[-c(1,2,n-1,n)]
  n=n-4

  shortlab = as.character(1:n)
  longlab = as.character(1:n)
  cx = numeric(n)
  cy = numeric(n)
  cz = numeric(n)

  Sys.setlocale('LC_ALL','C')

  for(i in 1:n){
    ll = strsplit(ff[[i]],"\"")[[1]]
    nn = as.numeric(strsplit(ll[6]," ")[[1]])
    shortlab[i] = ll[2]
    longlab[i] = ll[4]
    cx[i] = nn[1];cy[i]=nn[2];cz[i]=nn[3]
  }
  acp = !is.na(cx)
  list(shortlab[acp],longlab[acp],facem=cbind(cx[acp],cy[acp],cz[acp]))

}

#' Array of faces constructed from geodesic and profile information
#'
#' @param initshap shape matrix with starting posture
#' @param farshap shape matrix with extreme posture
#' @param profile vector of value in [0,1] for convex combination of postures
#'
#' @return array where each frame is a shape matrix. Number of frames is the length of the profile vector.
#' @export
constructfaceseq <- function(initshap,farshap,prof){
  ds <- dim(initshap)
  fseq <- array(NA,c(ds[1],ds[2],length(prof)))
  for(i in 1:length(prof)){
    fseq[,,i] <- prof[i]*farshap+(1-prof[i])*initshap
  }
  fseq
}

#' Array of faces constructed from geodesic and score information
#'
#' @param initshap initial shape
#' @param farshap extremal shape
#' @param score score profile
#'
#' @return facial motion sequence with same number frames as the score
#' @export
constructscorefaceseq <- function(initshap,farshap,score){
  ds <- dim(initshap)
  fseq <- array(NA,c(ds[1],ds[2],length(score)))
  ss <- score-min(score)
  ss <- ss/max(ss)
  for(i in seq_along(score)){
    fseq[,,i] <- ss[i]*farshap+(1-ss[i])*initshap
  }
  fseq
}


#' Plot an array of shapes with a different color for each marker
#'
#' @param fm array markers by directions by cases
#' @param mfconfig configuration information
#' @param view front, side or top
#' @param main optional title
#'
#' @return nothing
#' @export
multifaceplot <- function(fm, mfconfig, view=c("front","side","top"), main=NULL, subset=NULL){
  mdim = dim(fm)
 # pcolors = sample(colors(distinct = TRUE),mdim[1])
  pcolors = sample(rainbow(mdim[1]))
  view = match.arg(view)
  orthproj = switch(view, "front" = mfconfig$coordir[c(1,2)],
                    "side" = mfconfig$coordir[c(3,2)],
                    "top" = mfconfig$coordir[c(1,3)])
  alabs = switch(view, "front" = c("lateral","vertical"),
                    "side" = c("forward","vertical"),
                    "top" = c("lateral","forward"))
  if(is.null(main)) main <- paste(view, "view")
  if(is.null(subset)) subset = 1:mdim[1]

  latdir = mfconfig$coordir[1]
  if(latdir %in% orthproj){ # is lateral direction in selected view
    selpts = 1:mdim[1]
  }else{ # which markers are on the right
    rightside = which(fm[,latdir,1] < mean(fm[mfconfig$midline,latdir,1]))
    rightside = union(rightside, mfconfig$midline)
    selpts = rightside
  }
  selpts = intersect(subset, selpts)

  frame()
  plot.window(range(fm[selpts,orthproj[1],],na.rm=TRUE), range(fm[selpts,orthproj[2],],na.rm=TRUE), asp=1)
  title(main=main,xlab=alabs[1],ylab=alabs[2])
  axis(1)
  axis(2)
  box()
  for(i in selpts){
    points(fm[i,orthproj[1],], fm[i,orthproj[2],],col=pcolors[i],pch=20)
  }

}

#' Reflect face
#'
#' Reflect face by multiplying lateral direction by -1 and swapping the tags
#' of paired markers
#'
#' This function will work best when the face is centered with zero on
#' the midplane (lateral direction is tangent to this). But this is not
#' essential as reflection will happen anyway although the face may need to
#' be translated and/or rotated to match it up with the original face.
#'
#' @param sm face matrix
#' @param mfconfig face config info including left right pairs
#'
#' @return a face matrix
#' @export
reflectface <- function(sm, mfconfig){
  if(is.null(mfconfig$lrmat)) stop("Left Right pair matrix not defined")
  lrmat = mfconfig$lrmat
  nm = nrow(lrmat)
  latdir = mfconfig$coordir[1]
  midmn = abs(mean(sm[mfconfig$midline,latdir]))
  dd = sm[mfconfig$lrmat[,1],latdir] - sm[mfconfig$lrmat[,2],latdir]
  if(midmn > 0.1 * mean(abs(dd))) warning("Face may not be centered")
  for(i in 1:nm){
    x = sm[lrmat[i,1], latdir]
    sm[lrmat[i,1],latdir] = -sm[lrmat[i,2],latdir]
    sm[lrmat[i,2],latdir] = -x
    x = sm[lrmat[i,1], -latdir]
    sm[lrmat[i,1],-latdir] = sm[lrmat[i,2],-latdir]
    sm[lrmat[i,2],-latdir] = x
  }
  sm[mfconfig$midline,latdir] = -sm[mfconfig$midline,latdir]
  sm
}

#' Asymmetry scores
#'
#' Measures of asymmetry in object containing some matched pairs
#'
#' Face is reflected and the ordinary Procrustes applied to match up
#' with distances calculated between markers
#'
#' @param sm a face matrix
#' @param mfconfig face config info including left right pairs
#'
#' @return scores for each marker
#' @export
asymmetry <- function(sm, mfconfig){
  if(is.null(mfconfig$lrmat)) stop("Left Right pair matrix not defined")
  rsm = reflectface(sm, mfconfig)
  opr = faceOPA(sm, rsm, scale=FALSE)
  apply(opr$Ahat - opr$Bhat,1,function(x) sqrt(sum(x^2)))
}

#' Make a face symmetrical
#'
#' Make a face symmetrical by averaging with reflection
#'
#' Face is reflected and averaged with reflection
#'
#' @param sm a face matrix
#' @param mfconfig face config info including left right pairs
#'
#' @return scores for each marker
#' @export
makesym <- function(sm, mfconfig){
  if(is.null(mfconfig$lrmat)) stop("Left Right pair matrix not defined")
  rsm = reflectface(sm, mfconfig)
  (rsm + sm)/2.0
}

#' Location shift onto template
#'
#' Face is location-adjusted to make selected marker coincide with the template
#'
#' Usually do this with some (relatively) fixed landmark like the bridge of the nose
#'
#' @param templface the template face
#' @param face the face to be adjusted
#' @param marker the marker on the template to be used as the anchor
#'
#' @return the adjusted face
#' @export
bridgefix <- function(templface,face,marker=4){
  nd <- face[marker,]-templface[marker,]
  sweep(face,2,nd)
}

#' Compute size of face
#'
#' Average distance of markers from the origin
#'
#' Assumes the face has been centered
#'
#' @param f the face
#' @param subset of markers to be considerd
#'
#' @return the size of the face
#' @export
facenorm <- function(f,subset=1:nrow(f)){
  mean(sqrt(apply(f[subset,]^2,1,sum)))
}



#' Array mean
#'
#' computes array mean across the third dimension of an array
#'
#' surprisingly faster than using apply because colMeans is fast
#'
#' @param A a 3D array
#'
#' @return a matrix with dimensions equal to the first two dimensions of A
#' @export
arraymean = function(A){
  colMeans(aperm(A,c(3,1,2)))
}

#' Partial Procrustes rotation
#'
#' rotates B onto A
#'
#'
#' rotates B onto A
#'
#' Assumes that A and B have both been centered (and preferably scaled)
#'
#' @param A a face matrix
#' @param B a face matrix
#'
#' @return the rotated B
#' @export
parproc = function(A,B){
  k = ncol(A)
  ABCP <- crossprod(A,B)
  sv <- La.svd(ABCP,k,k)
  sv$u[,k] <- sv$u[,k]*determinant(ABCP)$sign
  B%*%tcrossprod(t(sv$vt),sv$u)
}

#' Ordinary Procrustes Analysis
#'
#' Fit one face onto another
#'
#' Face A is centered, Face B is centered then rotated to obtain the best Procrustes fit
#'
#' @param A a face matrix (no missing values allowed)
#' @param B another face matrix with the same dimension (missing values allowed)
#'
#' @return a list with the following components:
#' \item{Ahat}{centered A}
#' \item{Bhat}{B centered and rotated onto A}
#' \item{OSS}{sum of squares for the difference}
#' \item{rmsd}{RMS difference = sqrt(OSS/(no. of markers))}
#' \item{mss}{square distance per marker}
#' @export
faceOPA = function(A,B){
  if(any(is.na(A))) stop("Missing values in first argument")
  missrow = apply(B,1,function(x) any(is.na(x)))
  adim = dim(A)
  if(any(missrow)){
    Aorig = A
    A = A[!missrow,]
    B = B[!missrow,]
  }
  csA = center.scale(A)
  csB = center.scale(B)
  Ahat = csA$coords * csA$CS
  Bhat = parproc(csA$coords, csB$coords) * csB$CS
  if(any(missrow)){
    Afull = matrix(NA,adim[1],adim[2])
    Afull[!missrow] = Ahat
    Adiff = Afull - Aorig
    Ahat = Aorig
    Bfull = matrix(NA,adim[1],adim[2])
    Bfull[!missrow] = Bhat
    Bhat = Bfull - Adiff
  }
  mss = apply(Ahat - Bhat, 1,function(x) sum(x^2))
  OSS = sum(mss,na.rm=TRUE)
  rmsd = sqrt(OSS/sum(!missrow))

  list(Ahat=Ahat,Bhat=Bhat,OSS=OSS,rmsd=rmsd,mss=mss)
}

#' Generalized Procrustes Analysis
#'
#' @param A a 3D array arranged as markers x dimensions x faces
#' @param itmax maximum number of iterations (5 by default)
#'
#' @return a list containing the rotated array of faces, the mean face and the centroid sizes of the faces
#' @export
faceGPA = function(A, itmax=5, pcaoutput=FALSE){
  n = dim(A)[3]
  size = rep(0,n)
  prms = Inf
  for(i in 1:n){
    A[,,i] = center(A[,,i])
    size[i] = csize(A[,,i])
  }
  for(j in 1:itmax){
    for(i in 1:n){
      mshape = arraymean(A[,,-i])
      A[,,i] = parproc(mshape,A[,,i])
    }
    mshape = arraymean(A)
    rmsd = sqrt(mean(apply(A,3,function(x) centroid.size(x-mshape)))/n)
    if(prms - rmsd < 0.001) break
    prms = rmsd
  }
  if(pcaoutput){
    flatrot = t(apply(A,3,as.numeric))
    pcx = prcomp(flatrot)
  }else{
    pcx = list(x=NA,sdev=NA,rotation=NA)
  }
  list(rotated=A, mshape=mshape, size=size, pcasd=pcx$sdev, scores=pcx$x, pcar=pcx$rotation)
}

#' Conform a trajectory to a template
#'
#' @param sm The motion data array (may contain missing values)
#' @param mfi template shape (no missing values allowed)
#'
#' @return a list with the following components:
#' \item{sm}{conformed trajectory}
#' \item{dm}{matrix of squared distances to the template by frame and by marker}
#' \item{rawstart}{the initial frame}
#' \item{rawext}{the extreme frame}
#' @export
trajtemplate = function(sm, mfi){
  if(any(is.na(mfi))) stop("Missing values in the template shape")

  dm = matrix(NA, dim(sm)[3], dim(sm)[1])
  for(i in 1:dim(sm)[3]){
    ft = faceOPA(mfi, sm[,,i])
    sm[,,i] = ft$Bhat
    dm[i,] = ft$mss
  }
  dmf = apply(dm,2,medianfill)
  distem = apply(dmf,1,function(x) sqrt(sum(x,na.rm=TRUE)))
  maxi = which.max(distem)
  list(sm=sm,dm=dm,distem=distem,rawstart=sm[,,1],rawext=sm[,,maxi])
}

# median fill

medianfill = function(x){
  md = median(x,na.rm=TRUE)
  x[is.na(x)] = md
  x
}


# center
# centers a matrix faster than scale()
# used in other functions for gpagen; digitsurface
center <- function(x){
  if(is.vector(x)) x- mean(x) else {
    x <- as.matrix(x)
    x - rep(colMeans(x), rep.int(nrow(x), ncol(x)))
  }
}

# helper functions

# csize
# calculates centroid size
# digitsurface
csize <- function(x) sqrt(sum(center(as.matrix(x))^2))

# cs.scale
# divide matrices by centroid size
# used in other functions for gpagen
cs.scale <- function(x) x/csize(x)

# center.scale
# center and divide matrices by centroid size; faster than scale()
# used in other functions for gpagen
center.scale <- function(x) {
  x <- center(x)
  cs <- sqrt(sum(x^2))
  y <- x/cs
  list(coords=y, CS=cs)
}

