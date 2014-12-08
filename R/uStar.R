#'@title Calculates the water friction velocity, uStar
#'@description uStar is the water friction velocity due to wind stress at the lake surface, 
#'it is calculated following the methods of Imberger (1985) as a function of the shear stress 
#'of air (Fischer et al., 1979), drag coefficient for momentum (Hicks, 1972), and a 
#'dimensionless constant (von Karman constant) that decribes the logarithmic velocity 
#'profile at the air-water interface
#'@param wnd a numeric vector of wind speed in m s-1
#'@param wnd.height Height of the anemometer above the lake surface in meters
#'@param averageEpiDense a numeric vector of epilimnion density in kg m-3
#'@return a numeric vector of uStar
#'@author R. Iestyn. Woolway
#'@references Hicks, B.B., 1972. A procedure for the formulation of bulk transfer coefficients 
#'over water bodies of different sizes. Boundary-Layer Meterology 3: 201-213.
#'
#'Amorocho, J., DeVries, J.J., 1980. A new evaluation of the wind stress coefficient over water 
#'surfaces. Journal of Geophysical Research 85: 433-442.
#'
#'Fischer, H.B., List, E.J., Koh, R.C.Y., Imberger, J., Brooks, N.H., 1979. Mixing in inland and 
#'coastal waters. Academic Press.
#'
#'Imberger, J., 1985. The diurnal mixed layer. Limnology and Oceanography 30: 737-770.
#'
#'@seealso \code{\link{ts.uStar}}, \code{\link{layer.density}}
#'@examples
#'  wndSpeed  <- c(5.1,6.3,6.3,5.2,7,7.2)
#'  wndHeight  <-	2
#'  averageEpiDense	<- c(14,15,14.2,13,12,12)
#'  
#'  cat('uStar for input vector is: ')
#'  cat(uStar(wndSpeed,wndHeight,averageEpiDense))
#'@export
uStar <- function(wnd,wndHeight,averageEpiDense){
  
  # define constants
  rhoAir <- 1.2 # density of air
  vonK <- 0.4 # von Karman constant

  # -- calculate drag coefficient (from Hicks, 1972)
  Cd <- rep(0.0015, length(wnd))
  Cd[wnd < 5] <- 0.001

  # -- correct for wind measurement height if < 10 m (Amorocho and DeVries, 1980)
  if (wndHeight != 10){
    wnd <- wnd/(1-sqrt(Cd)/vonK*log(10/wndHeight))
  }
    
  # -- calculate shear stress of air (Fischer et al., 1979)
  tau <- Cd*rhoAir*wnd^2
  
  # -- calculate uStar following Imberger (1985)
  uStar <- sqrt(tau/averageEpiDense)
  return(uStar)
}

# -- References
#?? Hicks, B.B., 1972. A procedure for the formulation of bulk transfer
#?? coefficients over water bodies of different sizes. Boundary-Layer
#?? Meterology 3: 201-213

#?? Amorocho, J., DeVries, J.J., 1980. A new evaluation of the wind ??
#?? stress coefficient over water surfaces. Journal of Geophysical  ??
#?? Research 85: 433-442.

#?? Fischer, H.B., List, E.J., Koh, R.C.Y., Imberger, J., Brooks, N.H.,
#?? 1979. Mixing in inland and coastal waters. Academic Press.

#?? Imberger, J., 1985. The diurnal mixed layer. Limnology and Oceanography
#?? 30: 737-770.