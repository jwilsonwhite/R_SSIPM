#' Smurfdata
#'
#' Read in PISCO SMURF data and format into dataframe
#'
#' @return Dataframe
#'
#' @export
#
Smurfdata <- function(Species,Names=FALSE){
  D <- read.csv('data/PISCO_UCSB_subtidal_recruitment_fish_data.1.2.csv',header=TRUE)

  # Population with zeros (although it looks like this is not an issue for the SMURF data, at least not for all dates):
  D$Temp<-0
  D$Temp[D$Species==Species]<- D$Number.fish[D$Species==Species]

  # Summarize just for species of interest:
  D2<-aggregate(D$Temp,by=list(Island=D$Island,Site=D$SITE,Month=D$Month,Day=D$Day,Year=D$Year,NumDays=D$NumDays,N.Smurfs=D$N.Smurfs,
                               Lat=D$Lat,Lon=D$Long),FUN=sum)
  D2$fsd = D2$x / D2$NumDays / D2$N.Smurfs

  # Now sum up by year:
  D3.x<-aggregate(D2$x,by=list(Island=D2$Island,Site=D2$Site,Year=D2$Year,
                               Lat=D2$Lat,Lon=D2$Lon),FUN=sum)
  D3.fsd<-aggregate(D2$fsd,by=list(Island=D2$Island,Site=D2$Site,Year=D2$Year,
                                   Lat=D2$Lat,Lon=D2$Lon),FUN=sum)

  # Fix the names & combine frames
  D3<-D3.x
  D3$fsd <- D3.fsd$x

  # Make column names match species name:
  if (Names){
    colnames(D3)[colnames(D3)=='x']<-Species
    colnames(D3)[colnames(D3)=='fsd']<-paste(Species,'.fsd',sep="")
  } else {
    colnames(D3)[colnames(D3)=='x']<-'fish'
    colnames(D3)[colnames(D3)=='fsd']<-'fish.fsd'}

  #Combine some sites that were renamed and slightly relocated
  D3$Site[D3$Site=="HAZ-WEST"]<-"HAZ"
  D3$Site[D3$Site=='PEL-WEST']<-'PEL'
  D3$Site[D3$Site=='WIL']<-'VALLEY'
  D3$Site[D3$Site=='MORSE']<-'GULL'

  return(D3)
}
