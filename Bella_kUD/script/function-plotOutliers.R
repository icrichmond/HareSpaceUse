# Code was taken from Amanda Droghini - cite in paper
# Droghini, A. 2020. Southwest Alaska moose. Git Repository. Available: https://github.com/accs-uaa/southwest-alaska-moose

plotOutliers <- function(data,minIndex,maxIndex,output = NULL) {
  require(ggmap)
  
  temp <- data[minIndex:maxIndex,]
  print(unique(temp$tag.local.identifier))
  studyArea<-matrix(c(min(temp$location.long)-0.002,min(temp$location.lat)-0.002,
                      max(temp$location.long)+0.001,max(temp$location.lat))+0.001, nrow = 2)
  
  mapData <- get_map(studyArea, zoom=9, source="google", maptype="terrain")
  
  plotSubset <- ggmap(mapData)+
    geom_path(data=temp, aes(x=location.long, y=location.lat))
  
  print(plotSubset)
  if(!is.null(output)) 
    return(temp)
}
