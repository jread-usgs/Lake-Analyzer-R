#'@title Calculate the Schmidt stability
#'@description Schmidt stability, or the resistance to mechanical mixing 
#'due to the potential energy inherent in the stratification of the water column.
#'
#'@param wtr a numeric vector of water temperature in degrees C
#'@param depths a numeric vector corresponding to the depths (in m) of the wtr measurements
#'@param bthA a numeric vector of cross sectional areas (m^2) corresponding to bthD depths
#'@param bthD a numeric vector of depths (m) which correspond to areal measures in bthA
#'@param sal a numeric vector of salinity in Practical Salinity Scale units
#'
#'@details Schmidt stability was first defined by Schmidt (1928) and later modified by 
#'Hutchinson (1957). This stability index was formalized by Idso (1973) to reduce the 
#'effects of lake volume on the calculation (resulting in a mixing energy requirement 
#'per unit area).
#'
#'@return a numeric vector of Schmidt stability (J/m^2)
#'@author Luke Winslow, Jennifer Brentrup
#'@seealso \link{ts.schmidt.stability}, \link{lake.number}, \link{wedderburn.number}
#'@references Schmidt, W., 1928. Ueber Temperatur and Stabilitaetsverhaltnisse von Seen. 
#'Geo- graphiska Annaler 10, 145-177.
#'
#'Hutchinson, G.E., 1957. A Treatise on Limnology, vol. 1. John Wiley & Sons, Inc., New York.
#'
#'Idso, S.B., 1973. On the concept of lake stability. Limnology and Oceanography 18, 681-683.
#'@examples
#'bthA  <-	c(1000,900,864,820,200,10)
#'bthD	<-	c(0,2.3,2.5,4.2,5.8,7)
#'
#'wtr	<-	c(28,27,26.4,26,25.4,24,23.3)
#'depths	<-	c(0,1,2,3,4,5,6)
#'
#'cat('Schmidt stability for input is: ')
#'cat(schmidt.stability(wtr, depths, bthA, bthD))
#'@export
schmidt.stability = function(wtr, depths, bthA, bthD, sal = wtr*0){

if(length(wtr) != length(depths)){
	stop('water temperature array must be the same length as the depth array')
}

#Constants
g = 9.81
dz = 0.1

# Here is just some madeup data. This should 
# seem valid to the Schmidt Stability algorithm. Valid enough at least
#wtr = c(24,24,24,20,17,12,11,10,10)
#depths = 1:9
#sal = wtr*0
#bthD = 1:9
#bthA = seq(8,0,by=-1)

# if bathymetry has negative values, drop and interpolate to 0
if(min(bthD) < 0){
	useI = bthD >= 0
	
	if(any(bthD == 0)){
		depT = bthD[useI]
	}else{
		depT = c(0, bthD[useI])
	}
	
	bthA = approx(bthD, bthA, depT)$y
	bthD = depT
}

numD = length(wtr)
if(max(bthD) > depths[numD]){
	wtr[numD+1] = wtr[numD]
	sal[numD+1] = sal[numD]
	depths[numD+1] = max(bthD)
}else if(max(bthD) < depths[numD]) {
	bthD = c(bthD, depths[numD])
	bthA = c(bthA, 0)
}

if(min(bthD) < depths[1]) {
	wtr = c(wtr[1], wtr)
	sal = c(sal[1], sal)
	depths = c(min(bthD), depths)
}

Zo = min(depths)
Io = which.min(depths)
Ao = bthA[Io]

if(Ao == 0){
	stop('Surface area cannot be zero, check *.bth file')
}

#Calculate water density 
rhoL = water.density(wtr, sal)

#The approx (interp1 in matlab) just does linear interpolation
layerD = seq(min(depths), max(depths), by=dz)
layerP = approx(depths, rhoL, layerD)$y
layerA = approx(bthD, bthA, layerD)$y

Zv = layerD * layerA * dz
Zcv = sum(Zv)/sum(layerA)/dz

numInt = length(layerA)
st = layerA * NaN
for (i in 1:numInt){
	z = layerD[i]
	A = layerA[i]
	st[i] = -(Zcv-z)*layerP[i]*A*dz
}
St = g/Ao*sum(st)

return(St)

}

