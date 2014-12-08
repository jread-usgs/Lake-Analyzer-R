#'@title Calculates the center of buoyancy
#'@description Calculate the center of buoyancy using buoyancy frequency with 
#'a center of mass analysis. 
#'@details Brunt-Vaisala frequency is used for a temperature 
#'profile. Negative values for N2 are set to 0 (as they represent transient 
#'instabilities or sensor calibration issues) for this calculation.
#'@param wtr
#'@param depths
#'return Returns a value for the center of buoyancy.
#'@author Jordan S. Read
#'@seealso \code{\link{buoyancy.freq}}, \code{\link{ts.buoyancy.freq}}, 
#'\code{\link{ts.center.buoyancy}}
#'@examples 
#'# A vector of water temperatures
#'wtr = c(22.51, 22.42, 22.4, 22.4, 22.4, 22.36, 22.3, 22.21, 22.11, 21.23, 16.42, 
#'        15.15, 14.24, 13.35, 10.94, 10.43, 10.36, 9.94, 9.45, 9.1, 8.91, 8.58, 8.43)
#'
#'#A vector defining the depths
#'depths = c(0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
#'           17, 18, 19, 20)
#'
#'c.b = center.buoyancy(wtr, depths)
#'@export
center.buoyancy <- function(wtr, depths){
  
  if (depths[2] - depths[1] < 0 ){stop('depths must be in descending order')}
  N2 <- buoyancy.freq(wtr, depths)
  num.slices <- length(N2)
  areas <- vector('numeric',length = num.slices)
  cent.depths <- vector('numeric',length = num.slices)
  
  for (i in 1:num.slices){
    dz <- depths[i+1] - depths[i]
    areas[i] <- dz * N2[i] # assumes depths are in descending order
    cent.depths[i] <- mean(depths[i:(i+1)])
  }
  
  areas[areas < 0] <- 0
  cent.buoyancy <- sum(cent.depths*areas)/sum(areas)
  cent.buoyancy[cent.buoyancy == Inf] <- NA # division by zero
  return(cent.buoyancy)
}

#'@title Calculates the center of buoyancy for multiple temperature profiles.
#'@description Function for simplifying the calculation of the center of buoyancy. 
#'Can usually be called directly on data loaded directly using \code{\link{load.ts}} 
#'and \code{\link{load.bathy}}.
#'@param wtr a data frame of water temperatures (in Celsius). Loaded using \code{\link{load.ts}}
#'@return Returns a data frame with the timeseries of the center of buoyancy frequency. Includes a ‘datetime’ column.
#'@author Jordan S. Read
#'@seealso \code{\link{center.buoyancy}}, \code{\link{load.bathy}}, \code{\link{load.ts}}
#'@examples
#' #Get the path for the package example file included
#'wtr.path <- system.file('extdata', 'Sparkling.daily.wtr', package="rLakeAnalyzer")
#'
#'#Load data for example lake, Sparkilng Lake, Wisconsin.
#'sp.wtr = load.ts(wtr.path)
#'
#'#calculate and plot the thermocline depth
#'t.d = ts.thermo.depth(sp.wtr)
#'
#'center.N2 = ts.center.buoyancy(sp.wtr)
#'
#'plot(center.N2, type='l', ylab='Depth (m)', xlab='Date', ylim=c(19,0), lwd = 1.5)
#'lines(t.d, type='l', col='red', lwd = 1.5)
#'legend(x = t.d[3,1], y = .25,
#'       c('center of buoyancy','thermocline depth'),
#'       lty=c(1,1),
#'       lwd=c(1.5,1.5),col=c("black","red"), bty = "n")
#'@export
ts.center.buoyancy <- function(wtr){
  
  depths = get.offsets(wtr)
  
  n = nrow(wtr)
  
  wtr.mat = as.matrix(wtr[,-1])
  cent.n2 = rep(NA, n)
    
  for(i in 1:n){
    cent.n2[i] <- center.buoyancy(wtr.mat[i, ], depths)
  }
  
  cent.buoyancy = data.frame(wtr[,'datetime', drop=F], cent.n2)
  
  return(cent.buoyancy)
}