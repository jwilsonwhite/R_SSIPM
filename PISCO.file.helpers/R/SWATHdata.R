#' SWATHdata
#'
#' Function to read in PISCO kelp swath data and format into dataframe
#'
#' @return Dataframe
#'
#' @export
# Read in SWATH data

SWATHdata <- function(){

  # Species = multi-digit name [THIS CODE BASICALLY ASSUMES WE WANT MACROCYSTIS STIPES]
  Species <- 'MACPYRAD'

  D <- read.csv('data/PISCO_kelpforest_swath.1.2.csv',header=TRUE)

  # Only take SBTL_SWATH
  D <- D[D$method=='SBTL_SWATH',]

  # Population with zeros:
  D$Temp<-0
  OKrows = D$classcode==Species
  OKrows[is.na(OKrows)] = FALSE
  D$Temp[OKrows]<- D$macstipes[OKrows] # Would need to change this line if you want something besides macrocystis

  # Tally up number of stipes per transect. Combine across sides in a given site.
  D2<-aggregate(D$Temp,by=list(Site=D$site,Year=D$year,Side=D$side,
                               Zone=D$zone,Transect=D$transect),FUN=sum)

  # Get total count, plus sample size.
  D2total <- aggregate(D$Temp,by=list(Site=D$site,Side=D$side,Month=D$month,Day=D$day,Year=D$year),FUN=sum)
  D2ntrans <- aggregate(D$Temp,by=list(Site=D$site,Side=D$side,Month=D$month,Day=D$day,Year=D$year),FUN=length)

  # Count per transect
  D2total$stipespt <- D2total$x/D2ntrans$x

  # Now sum up by year:
  D3.x<-aggregate(D2total$x,by=list(Site=D2total$Site,Side=D2total$Side,Year=D2total$Year),FUN=sum)
  D3.stipespt<-aggregate(D2total$stipespt,by=list(Site=D2total$Site,Side=D2total$Side,Year=D2total$Year),FUN=sum)

  # Fix the names & combine frames
  D3<-D3.x
  D3$stipespt <- D3.stipespt$x

  # Make column names match species name:
  colnames(D3)[colnames(D3)=='x']<-Species
  colnames(D3)[colnames(D3)=='fpt']<-paste(Species,'.stipespt',sep="")

  return(D3)
}


