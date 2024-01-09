ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#libraries
packages <- c("VBTree", "tensorA", "Hmisc", "MVNBayesian", "stats", "graphics")
ipak(packages)


data <- as.data.frame(iris3)
names(data) <- mapply(gsub, colnames(data), pattern=" |\\.{2}", replacement="-")
dim(data)
cnames <- colnames(data)
cnames

dl <- chrvec2dl(cnames)
dl

vbt <- dl2vbt(dl)
vbt

ts <- vbt2ts(vbt)
# Generation from double list is alternative, using dl2ts();
# For conversion to array, call vbt2arr() or dl2arr().
ts

subvbt <- vbtsub(vbt, c(-1,-1,2))
subvbt
subdata <- datavisit(data, subvbt)
subdata

# for (i in 1:parts){
#   for (j in 1:species) {
#     plot body
#   }
# }

versicolorvbt <- vbtsub(vbt, c(-1,-1,2))
versicolorvbt

versicolor <- data[, vbt2arr(versicolorvbt)]
head(versicolor)

i <- 1
for (i in 1:3) {
  temp <- as.vector(vbt2arr(vbtsub(vbt, c(-1,-1, i))))
  print(temp)
}

autoscatter <- function(data, pltargs, ...){
  # note: pltargs is d3 vector, c(fixlayer, xlevel, ylevel)
  vbt <- dl2vbt(chrvec2dl(colnames(data)))
  dims <- vbt[[2]]
  rpt <- prod(dims)/dims[pltargs[1]]

  # xdata, ydata as chrvec
  xptr <- rep(-1, length(dims))
  yptr <- xptr
  xptr[pltargs[1]] <- pltargs[2]
  yptr[pltargs[1]] <- pltargs[3]

  xdata <- as.vector(vbt2arr(vbtsub(vbt, xptr)))
  ydata <- as.vector(vbt2arr(vbtsub(vbt, yptr)))

  dl <- vbt2dl(vbt)
  x_lab <- dl[[pltargs[1]]][pltargs[2]]
  y_lab <- dl[[pltargs[1]]][pltargs[3]]

  par(...) # plots arrangement

  i <- 1
  for (i in 1:rpt) {
    chrvectemp <- strsplit(xdata[i], split = "-")[[1]]
    title <- paste(chrvectemp[-pltargs[1]], collapse = ".")
    plot(x=data[, xdata[i]], y=data[, ydata[i]], pch = 20, col="red",
         xlab = x_lab, ylab = y_lab, main = title)
  }
}

autoscatter(data, c(2,2,1), mfcol=c(2,3))

require(tensorA)
groupscatter <- function(data, pltargs, grpby=NULL, ...){
  if(is.null(grpby)){
    autoscatter(data=data, pltargs = pltargs, ... = ...) # no group assignment
  }
  # note1: pltargs is d3 vector, c(fixlayer, xlevel, ylevel)
  # note2: grpby is group layer vector, with unfixed length
  vbt <- dl2vbt(chrvec2dl(colnames(data)))
  dims <- vbt[[2]]
  fixlayers <- append(grpby, pltargs[1])
  rpt <- prod(dims[-fixlayers]) # number of plots

  # build mapping
  overlaps <- prod(dims[grpby])
  trvsdim <- rep(1, length(dims))

  dl <- vbt2dl(vbt)
  x_lab <- dl[[pltargs[1]]][pltargs[2]]
  y_lab <- dl[[pltargs[1]]][pltargs[3]]

  # automatically generate query vector for vbtsub(), for x and y
  trvsdim[grpby] <- dims[grpby]
  xtrvs <- pos.tensor(trvsdim)
  # fix the bug of pos.tensor()
  if (length(xtrvs[,1])!=prod(trvsdim)){
    xtrvs <- matrix(1, prod(trvsdim), length(trvsdim))
    xtrvs[,length(trvsdim)] <- c(1:prod(trvsdim))
  }

  xtrvs[,-fixlayers] <- xtrvs[,-fixlayers]*(-1)
  ytrvs <- xtrvs
  xtrvs[, pltargs[1]] <- rep(pltargs[2], overlaps)
  ytrvs[, pltargs[1]] <- rep(pltargs[3], overlaps)

  # save character vectors in group
  emptlist <- list()
  xgroup <- list(emptlist)[rep(1,overlaps)]
  ygroup <- xgroup
  j <- 1
  for (j in 1:overlaps) {
    xgroup[[j]] <- as.vector(vbt2arr(vbtsub(vbt, as.vector(xtrvs[j,]))))
    ygroup[[j]] <- as.vector(vbt2arr(vbtsub(vbt, as.vector(ytrvs[j,]))))
  }

  par(...) # arrangement control

  # note3: rpt = number of plots
  # note4: overlaps = number of groups
  i <- 1
  for (i in 1:rpt) {
    j <- 1
    xtemp <- c()
    ytemp <- xtemp
    for (j in 1:overlaps) {
      xtemp <- append(xtemp, xgroup[[j]][i])
      ytemp <- append(ytemp, ygroup[[j]][i])
    }
    # axis scale for x and y
    xscale <- c(min(data[, unlist(xtemp)]), max(data[, unlist(xtemp)]))
    yscale <- c(min(data[, unlist(ytemp)]), max(data[, unlist(ytemp)]))

    j <- 1
    ptlegend <- c()
    for (j in 1:overlaps) {

      # making title and legend
      labtable <- unlist(chrvec2dl(xgroup[[j]][i]))
      notation <- paste(labtable[grpby], collapse = ".")
      ptlegend <- append(ptlegend, notation)
      mainlab <- paste(labtable[-fixlayers], collapse = ".")

      # plot in groups
      clrs <- palette()[j+1]
      plot(x=data[, xgroup[[j]][i]], y=data[, ygroup[[j]][i]], pch=20, col=clrs,
           xlab = x_lab, ylab = y_lab, xlim = xscale, ylim = yscale, main=mainlab)
      if(j==overlaps){
        legend("topleft", ptlegend, fill = palette()[2:(overlaps+1)],
               cex = 0.7, bg = "transparent", box.lty = 0)
        par(new=FALSE) # to next plot
      } else {
        par(new=TRUE) # overlap scatter points in next group
      }
    }
  }
}

groupscatter(data, c(2,2,1), 3, mfrow=c(1,2))
groupscatter(data, c(2,1,2), 1, mfrow=c(1,3))
groupscatter(data, c(2,2,1), c(1,3), mfrow=c(1,1))
groupscatter(data, c(2,2,1), mfcol=c(2,3))

set.seed(9)
data1 <- data
names(data1) <- as.vector(mapply(paste, cnames, "Canada", sep="-"))
data2 <- data + rnorm(prod(dim(data))) + 0.5
names(data2) <- as.vector(mapply(paste, cnames, "US", sep="-"))
newdata <- cbind(data1, data2)

require(Hmisc)
require(tensorA)
require(MVNBayesian)
groupmodeling <- function(data, ftlys, grpby){
  # note1: ftlys is an integer vector for layer(s) of features
  # note2: grpby is an integer vector for group conditions
  vbt <- dl2vbt(chrvec2dl(colnames(data)))
  dims <- vbt[[2]]

  rpt <- prod(dims[grpby])
  sumrpt <- prod(dims[-c(grpby, ftlys)])
  fts <- prod(dims[ftlys])

  initinq <- rep(1, length(dims))

  ftnameinq <- initinq
  ftnameinq[ftlys] <- ftnameinq[ftlys]*(-1)
  ftnamedl <- vbtinq(vbt, ftnameinq)
  i <- 1
  for (i in 1:length(ftnamedl)) {
    if(i %in% ftlys){
      ftnamedl[[i]] <- ftnamedl[[i]]
    } else {
      ftnamedl[[i]] <- "*"
    }
  }
  namevec <- as.vector(dl2arr(ftnamedl))
  namevec <- mapply(gsub, namevec, pattern="[*]|", replacement="")
  namevec <- mapply(gsub, namevec, pattern="^[-]{0,}|[-]{0,}$", replacement="")
  namevec <- as.vector(mapply(gsub, namevec, pattern="[-]{1,}", replacement="-"))

  grparr <- array(NA, c(rpt, sumrpt, fts))

  grpinq <- initinq
  trvsidx <- initinq
  trvsidx[-ftlys] <- dims[-ftlys]
  trvstb <- pos.tensor(trvsidx)
  trvstb[,ftlys] <- trvstb[,ftlys]*(-1)

  grpptr <- pos.tensor(grpby)

  i <- 1
  for (i in 1:rpt) {
    idx <- as.numeric(find.matches(t(grpptr[i,]),trvstb[,grpby])[[1]])
    j <- 1
    for (j in 1:sumrpt) {
      grparr[i,j,] <- as.vector(vbt2arr(vbtsub(vbt, trvstb[idx[j],])))
    }
  }

  emptlist <- list()
  result <- list(emptlist)[rep(1, rpt)]

  i <- 1
  for (i in 1:rpt) {
    j <- 1
    sumtb <- matrix(NA, 1, fts)
    colnames(sumtb) <- namevec
    for (j in 1:sumrpt) {
      clns <- grparr[i,j,] # sum for data
      temp <- data[, clns]
      colnames(temp) <- namevec
      sumtb <- rbind(sumtb, temp)
      grpdata <- sumtb[-1,]

      grpname <- unlist(chrvec2dl(grparr[i,j,1]))[grpby] # group name
      if (length(grpname)!=1){
        grpname <- paste(grpname, collapse = "-")
      }

      grpdim <- dim(grpdata) # Bayesian Posteriori using least informative
      Lambda <- rWishart(1, grpdim[2], (var(grpdata)/grpdim[2]))[,,1]
      BP <- MVN_BayesianPosteriori(grpdata, colMeans(grpdata), solve(Lambda*grpdim[1]))
    }
    result[[i]] <- list(group=grpname, mixpar=grpdim[1], mean=BP[[1]], var=BP[[2]])
  }

  # modify mixture parameters
  counts <- 0
  i <- 1
  for (i in 1:rpt) {
    counts <- counts + result[[i]][[2]]
  }
  i <- 1
  for (i in 1:rpt) {
    result[[i]][[2]] <- result[[i]][[2]]/counts
  }
  return(result)
}

groupmodeling(data, c(1,2), 3)
groupmodeling(newdata, c(1,2), 3)
groupmodeling(newdata, c(1,2), 4)
groupmodeling(newdata, c(1,2), c(3,4))

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
i <- 1
for (i in 1:3) {
  temp <- rbind(iris3[,1:2,i], iris3[,3:4,i])
  colnames(temp) <- gsub("Sepal |Petal |\\.", "", colnames(iris3[,,i]))[1:2]
  temp <- cbind(parts=rep(gsub("W|L|\\.", "", colnames(iris3[,,i])),each=25),
                species=rep(names(iris3[1,1,][i]), 100), temp)
  temp <- as.data.frame(temp)
  assign(paste("element", i, sep = ""), temp)
}
trddataframe <- rbind(element1, element2, element3)
trddataframe[,3:4] <- sapply(trddataframe[,3:4], as.numeric.factor)
head(trddataframe)
dim(trddataframe)

arr <- vbt2arr(vbt)
object.size(trddataframe)
object.size(data)
object.size(cnames)
object.size(ts)
object.size(arr)
object.size(vbt)
object.size(dl)
