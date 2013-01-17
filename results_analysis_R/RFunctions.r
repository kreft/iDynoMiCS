#library(lattice)
#library(MASS)
install.packages('stringr',repo="http://www.stats.bris.ac.uk/R/")
library(vegan)
library(rJava)       # these two are needed to read in pictures into R
library(ReadImages)
library(stringr) #needed to read the fikes in in the right order!
#require(utils)

#function to read all the agent states in. One arguement, path to folder containing the files
#path defaults to the working directory
readStates <- function(path="."){
  Numbers <- sort(as.numeric(str_extract(dir(path=".", pattern="agent_State.*.csv"), "[0-9]+")))
  Files <- sapply(Numbers, function(x) paste("agent_State(",x,").csv", sep=""))
  agent_States <- sapply(Files, function(x) read.csv(x,header= TRUE, dec=".", na.strings="NA", skip = 2, fill = TRUE));
  agent_States <- t(agent_States);
  return(agent_States);
}
readEnv <- function(path="."){
  Numbers <- sort(as.numeric(str_extract(dir(path=".", pattern="env_State.*.csv"), "[0-9]+")))
  Files <- sapply(Numbers, function(x) paste("env_State(",x,").csv", sep=""))
  env_States <- sapply(Files, function(x) read.csv(x,header= TRUE, dec=".", na.strings="NA", skip = 3, fill = TRUE));
  env_States <- t(env_States);
  return(env_States);
}

################################################################################
############# plot agents on x and y    corresponds to plotAgents in Matlab  ###
################################################################################
# the following function automates plotting a list as would be the result from reading in as above.
# the blobs at the respective agent sites are coloured according to the species, and a legend is provided
plotAgents <- function(filex) {
  if(is.numeric(filex)){
    name <- paste("agent_State(",filex,").csv", sep="")
    filex <- read.csv(name,header= TRUE, dec=".", na.strings="NA", skip = 2, fill = TRUE)
  }
  x_lim <- filex$resolution[1]*filex$grid_ni[1]
  y_lim <- filex$resolution[1]*filex$grid_nj[1]
  no.species <- length(unique(filex$species))
  plot(filex$locationY, filex$locationX, las= 1, xlab = "y [µm]", ylab = "x [µm]", xlim = c(0,x_lim), ylim=c(0,y_lim), cex=.6, col = as.factor(filex$species_id), pch=19)
  lege <- c()
  for(i in 1:no.species){
    lege[i] <- paste("sp.", unique(filex$species)[i])
  }
  legend("topright",  legend = lege,  text.col = c(1:no.species))
}

################################################################################
############# time course  
################################################################################

plotTimeCourseAgents <- function(y=NULL) {
  if(is.null(y)){
    y <- readStates()
  }
  spab <- matrix();
  mini <- 1e66;
  maxi <- 0;
  for(i in 1:dim(y)[1]){
    mini <- min(min(y[i][[1]]),mini)
    maxi <- max(max(y[i][[1]]),maxi)
  }
  for(i in 1:dim(y)[1]){
    temp <- c()
    for(j in mini:maxi){
      temp <- c(temp,length(which(y[i,]$species_id == j)))
    }
    if(i==1){
      spab <- temp
    }
    else{
      spab <- rbind(spab, temp)
    }
  }
  spab <- t(spab)
  plot(0, type = "n", main = "Individual abundances per species", las = 1, xlab = "Time course", ylab = "species abundance", xlim = c(0,dim(y)[1]), ylim = c(0,max(spab)))
  lapply(1:nrow(spab), function(p){points(spab[p, ], col = p, pch = 1) ; lines(spab[p, ], col = p)})
  lege <- c()
  for(i in 1:nrow(spab)){
    lege[i] <- paste("sp.", i-1)
  }
  legend("topleft",  legend = lege,  text.col = c(1:nrow(spab)))
  return(spab)
}

################################################################################
############# diversity over time -  corresponds to Matlab plotTime(type)  #####
################################################################################
simpsonIndex <-  function(y) {
  div_time <- NULL
  for (i in 1:dim(y)[2]) {
    # get the diversity of species
      div_time[i] <- diversity(y[,i], index = "simpson", MARGIN = 2)
      maxc <- max(div_time, na.rm = TRUE)
    }
    plot(div_time,  xlab = "Time course", ylab = "Species diversity (Simpson)", ylim= c(0,maxc), main = "Species diversity", las = 1, type = "b")
    return(div_time)
}

################################################################################
############# abundance over time  corresponds to Matlab plotTime(type)  #######
################################################################################

plotTimeCourseAbund <- function(y=NULL) {
  if(is.null(y)){
    y <- readStates()
  }
  abTimes <- NULL
  for (i in 1:dim(y)[1]){
      x <- y[i,]
      abTimes[i] <- sum(table(x[1]))
    }
    maxc <- max(abTimes, na.rm=TRUE)
    xmax <- length(y)
    plot(abTimes, col = 1, xlab = "time course", ylab = "Total abundance",  ylim= c(0,maxc), main = "Total abundances", las = 1, type = "n")
    points(abTimes, col = 1, pch = 1) ; lines(abTimes, col = 1)
}

################################################################################
############# contours    corresponds to Matlab plotContour(nFile,fieldName) ###
##########################R######################################################

#example_biomass <-  read.table("bac_000_x1_100000.txt")
plotContour <- function(iteration, solute){
  if(is.numeric(iteration)){
    name <- paste("env_State(",iteration,").csv", sep="")
    iteration <- read.csv(name,header= TRUE, dec=".", na.strings="NA", skip = 3, fill = TRUE)
  }
  nums <- which(iteration[[1]]==solute)
  minix <- min(iteration[[2]][nums])
  maxix <- max(iteration[[2]][nums])
  x <- minix:maxix 
  miniy <- min(iteration[[3]][nums])
  maxiy <- max(iteration[[3]][nums])
  y <- miniy:maxiy
  temp <- matrix(nrow=length(y), ncol=length(x))
  for(i in x){
    for(j in y){
	  if(minix==0){
	    iv <- i+1
	  }
	  else{
	    iv <- i
	  }
	  if(miniy==0){
	    jv <- j+1
	  }
	  else{
	    jv <- j
	  }
	  temp[iv,jv] <- iteration[[5]][nums][which(iteration[[2]][nums]==i)][which(iteration[[3]][nums][which(iteration[[2]][nums]==i)]==j)]
	}
  }
  filled.contour(1:dim(temp)[2],1:dim(temp)[1],t(temp))
  return(temp)
}


################################################################################
############# read in an image and overlay the contours ########################
################################################################################
Overlay <- function(picture, contours){
  #THIS DOES NOT WORK AT PRESENT!!!!
  #pic <- read.jpeg(picture)  # that s the picture copied from the website
  #plot(pic) #plot image
  #contour(1:dim(temp)[2]*(dim(pic)[2]/dim(temp)[2]),1:dim(temp)[1]*(dim(pic)[1]/dim(temp)[1]),t(temp), add=TRUE)
  print("This is non-functional atm.")
}

