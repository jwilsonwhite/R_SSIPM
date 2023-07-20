#' SBTLdata
#'
#' Read in PISCO subtidal visual transect data and format into dataframe
#'
#' @return Dataframe
#'
#' @export
#
# Read in SBTL data
SBTLdata <- function(Species,MinLen,pool.sides=TRUE,Names=FALSE){

  # Species = 4-digit name
  # MinLen = minimum length to be counted (e.g., to exclude PLCA YOY, choose 9)
  # Names = if TRUE, the resulting data frame will have data columns with the species name. Otherwise, 'fish'

  D <- read.csv('data/PISCO_kelpforest_fish.1.3.csv',header=TRUE)

  # Only take SBTL_FISH from UCSB
  D <- D[D$method=='SBTL_FISH',]

  # Population with zeros:
  D$Temp<-0
  OKrows = D$classcode==Species & D$fish_tl >= MinLen
  OKrows[is.na(OKrows)] = FALSE # there was one stray NA at SBI_SOUTHEAST_SEA_LION
  D$Temp[OKrows]<- D$count[OKrows]


  # Tally up number of fish per transect. Replicate across sides in a given site.
  # D2<-aggregate(D$Temp,by=list(Site=D$site,Year=D$year,Side=D$side,
  #                              Zone=D$zone,Level=D$level,Transect=D$transect),FUN=sum)

  # Get total count, plus sample size.
  D2total <- aggregate(D$Temp,by=list(Site=D$site,Side=D$side,Month=D$month,Day=D$day,Year=D$year),FUN=sum)
  D2ntrans <- aggregate(D$Temp,by=list(Site=D$site,Side=D$side,Month=D$month,Day=D$day,Year=D$year),FUN=length)

  # Count per transect
  D2total$fpt <- D2total$x/D2ntrans$x

  # Now sum up by year & side:
  if (!pool.sides){
    D3.x<-aggregate(D2total$x,by=list(Site=D2total$Site,Side=D2total$Side,Year=D2total$Year),FUN=sum)
    D3.fpt<-aggregate(D2total$fpt,by=list(Site=D2total$Site,Side=D2total$Side,Year=D2total$Year),FUN=sum)
  } else { # pool sides together
    D3.x<-aggregate(D2total$x,by=list(Site=D2total$Site,Year=D2total$Year),FUN=sum)
    D3.fpt<-aggregate(D2total$fpt,by=list(Site=D2total$Site,Year=D2total$Year),FUN=sum)
  }

  # Fix the names & combine frames
  D3<-D3.x
  D3$fpt <- D3.fpt$x

  # Make column names match species name:
  if (Names){
    colnames(D3)[colnames(D3)=='x']<-Species
    colnames(D3)[colnames(D3)=='fpt']<-paste(Species,'.fpt',sep="")
  } else {
    colnames(D3)[colnames(D3)=='x']<-'fish'
    colnames(D3)[colnames(D3)=='fpt']<-'fish.fpt'}

  return(D3)
}
