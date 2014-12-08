## Helper functions for Lake Analyzer R
#'@title Gets depths from data frame containing profile info.
#'@description Extracts the depth information from a data frame 
#'containing multi-depth observation data. Relies on the format 
#'of the header to get information and may fail if your file 
#'format is incorrect.
#'@param data Data frame returned from \link{load.ts}.
#'@return A numeric vector of depth values. Should be the ncol(data) - 
#'1 in length as the first column contains date/time data.
#'@author Luke Winslow
#'@seealso \link{load.ts}
#'@examples
#'#Get the path for the package example file included
#'exampleFilePath <- system.file('extdata', 'Sparkling.wtr', package="rLakeAnalyzer")
#'
#'#Load
#'sparkling.temp = load.ts(exampleFilePath)
#'
#'#get the lake depths associated with each column
#'depths = get.offsets(sparkling.temp)
#'
#'print(depths)
#'@export
get.offsets <- function(data){
  
  header = names(data)
  
  #check for existence of datetime header and drop if there
  dt_indx = grep(pattern= "datetime", x= header, ignore.case= TRUE)
  if(length(dt_indx) > 0){
  	header = header[-dt_indx] #Drop datetime
  }
  
  matches = regexpr("(\\d+\\.?\\d*)" ,header)
  
  lengths = attr(matches,'match.length')
  offsets = vector(mode="numeric", length=length(matches))
  
  for(i in 1:length(matches)){
    offsets[i] = as.numeric(substr(header[i], matches[i], matches[i] + lengths[i]))
  }
  
  return(offsets)
}

# -- private function --
get.drho_dz <- function(wtr, depths){
	numDepths = length(wtr)
	
	rhoVar = water.density(wtr)
	
	drho_dz = vector(mode="double", length=numDepths-1);
	
	#Calculate the first derivative of density
	for(i in 1:numDepths-1){
		drho_dz[i] = ( rhoVar[i+1]-rhoVar[i] )/( depths[i+1] - depths[i] );
	}
	drho_dz
}
