#' linkdata
#'
#' Function for combining SMURF, subtidal fish, and kelp swath dataframes into a single dataframe
#'
#' @param dSmurf SMURF dataframe
#' @param dFish  subtidal fish dataframe
#' @param Dswath kelp swath dataframe
#' @param Meta Metadata associated with site IDs
#' @param M Natural mortality rate (for calculating discounted recruitment index)
#' @param pool.sides Pool both 'sides' of the study site (default TRUE)
#' @param interp Interpolate missing data years (default FALSE)
#'
#' @return Combined dataframe. Also returns digits counting through the loop over sites and years, to show progress
#'
#' @export

linkdata <- function(Dsmurf,Dfish,Dswath,Meta,M,pool.sides=TRUE,interp=FALSE){

  # Now loop over SBTL sites/years and add in appropriate data
  OKSites = (Meta$sitename_short[Meta$include])
  OKSites = levels(factor(OKSites))

  Dfish$MPA <- as.factor('reference')
  levels(Dfish$MPA)<-c('reference','SMR','SMCA')
  Dfish$MACPYRAD <- NA
  Dfish$stipespt <- NA
  Dfish$smurf.year <- NA
  Dfish$smurf.wgt <- NA
  Dfish$smurf.wgt5 <- NA
  Dfish$smurf.year.kelp <- NA
  Dfish$smurf.wgt.kelp <- NA

  for (s in 1:length(OKSites)){
    OK <- Dfish$Site==OKSites[s] # subset data by site
    Dfish$MPA[OK] <- Meta$MPA_STATUS[Meta$sitename_short==OKSites[s]][1] # Assign MPA status
    Years <- unique(Dfish$Year[OK]) # number of years
    Smurfsite <- as.character(Meta$smurf_site[Meta$sitename_short==OKSites[s]][1]) # which smurfing site to use
    OKsmurf <- Dsmurf$Site == Smurfsite # Index the rows of Dsmurf to use
    OKswath <- Dswath$Site == OKSites[s] # Index the rows of Dswath to use
    SmurfSub <- Dsmurf[OKsmurf,]

    for (y in 1:length(Years)){
      OKy <- Dfish$Year == Years[y] & OK
      # add swath data
      swathyear <- Dswath$Year == Years[y] & OKswath
      #print(s)
      Sides = as.character(unique(Dfish$Side[OKy]))
      if (any(swathyear)){

        if (!pool.sides){
          for (d in 1:length(Sides)){
            if (any(  swathyear & Dswath$Side==Sides[d] )){ # if there is one for the corresponding side
              Dfish$MACPYRAD[OKy & Dfish$Side==Sides[d]] <- Dswath$MACPYRAD[swathyear & Dswath$Side==Sides[d]]
              Dfish$stipespt[OKy & Dfish$Side==Sides[d]] <- Dswath$stipespt[swathyear & Dswath$Side==Sides[d]]

              # interpolate if it is NA
              if (interp){
                if (is.na(Dfish$MACPYRAD[OKy & Dfish$Side==Sides[d]])) {
                  Dfish$MACPYRAD[OKy & Dfish$Side==Sides[d]] <- median(Dswath$MACPYRAD[OKswath & Dswath$Side==Sides[d]],na.rm=TRUE)
                  Dfish$stipespt[OKy & Dfish$Side==Sides[d]] <- median(Dswath$stipespt[OKswath & Dswath$Side==Sides[d]] ,na.rm=TRUE)
                }}

            }else{ # otherwise just take a site mean
              Dfish$MACPYRAD[OKy & Dfish$Side==Sides[d]] <- mean(Dswath$MACPYRAD[swathyear])
              Dfish$stipespt[OKy & Dfish$Side==Sides[d]] <- mean(Dswath$stipespt[swathyear])}

            # interpolate if it is NA
            if (interp){
              if (is.na(Dfish$MACPYRAD[OKy & Dfish$Side==Sides[d]])) {
                Dfish$MACPYRAD[OKy & Dfish$Side==Sides[d]] <- median(Dswath$MACPYRAD[OKswath],na.rm=TRUE)
                Dfish$stipespt[OKy & Dfish$Side==Sides[d]] <- median(Dswath$stipespt[OKswath] ,na.rm=TRUE)
              }}
          } # end loop over sides
        } else { # if pooling sides
          Dfish$MACPYRAD[OKy] <- mean(Dswath$MACPYRAD[swathyear],na.rm=TRUE)
          Dfish$stipespt[OKy] <- mean(Dswath$stipespt[swathyear],na.rm=TRUE)
        } # end if pool.sides

      } else{# end if any
        if (interp){
          for (d in 1:length(Sides)){
            Dfish$MACPYRAD[OKy & Dfish$Side==Sides[d]] <- median(Dswath$MACPYRAD[OKswath],na.rm=TRUE)
            Dfish$stipespt[OKy & Dfish$Side==Sides[d]] <- median(Dswath$stipespt[OKswath] ,na.rm=TRUE)}}}


      # Add smurf data (year-to-year)
      Smurf_y <- SmurfSub[SmurfSub$Year == Years[y],] # set of smurf data preceding the current SBTL year
      if (dim(Smurf_y)[1]>0){
        Dfish$smurf.year[OKy] <- Smurf_y$PCLA.fsd
        Dfish$smurf.year.kelp[OKy] <- Smurf_y$PCLA.fsd * log10(Dfish$stipespt[OKy]+1)
      } else {
        if (interp){
          Dfish$smurf.year[OKy] <- median(SmurfSub$PCLA.fsd,na.rm=T) # use median of all years
          Dfish$smurf.year.kelp[OKy] <- median(SmurfSub$PCLA.fsd,na.rm=T) * log10(Dfish$stipespt[OKy]+1)
        }}

      # Add smurf data (discounted)
      #Smurf_y <- SmurfSub[SmurfSub$Year < Years[y],] # set of smurf data preceding the current SBTL year
      Tmp.PCLA.fsd <- Dfish$smurf.year[OK & Dfish$Year < Years[y]]
      Tmp.PCLA.fsd.kelp <- Dfish$smurf.year.kelp[OK & Dfish$Year < Years[y]]
      Tmp.Year <- Dfish$Year[OK & Dfish$Year < Years[y]]
      Tmp.Year.kelp <- Tmp.Year
      # trim NAs
      Tmp.Year <- Tmp.Year[!is.na(Tmp.PCLA.fsd)]
      Tmp.PCLA.fsd <- Tmp.PCLA.fsd[!is.na(Tmp.PCLA.fsd)]
      Tmp.Year.kelp <- Tmp.Year.kelp[!is.na(Tmp.PCLA.fsd.kelp)]
      Tmp.PCLA.fsd.kelp <- Tmp.PCLA.fsd[!is.na(Tmp.PCLA.fsd.kelp)]

      if (length(Tmp.Year)>0){ # if there are any smurf data preceding this year
        Smurf_penal <- sum(Tmp.PCLA.fsd * exp(-M * (Years[y]- Tmp.Year) ),na.rm=T)
        Smurf_penal.kelp <- sum(Tmp.PCLA.fsd.kelp * exp(-M * (Years[y]- Tmp.Year.kelp) ),na.rm=T)
        Dfish$smurf.wgt[OKy] <- Smurf_penal
        Dfish$smurf.wgt.kelp[OKy] <- Smurf_penal.kelp
      } # end if there are any smurf data

      if (length(Tmp.Year)>=5){ # if there are at least 5 years of smurf data (for older fish)
        # zero out more recent years
        Tlen = length(Tmp.Year)
        Tmp.Year[(Tlen-3):Tlen] = 0
        Smurf_penal5 <- sum(Tmp.PCLA.fsd * exp(-M * (Years[y]- Tmp.Year) ),na.rm=T)
        Dfish$smurf.wgt5[OKy] <- Smurf_penal5
      }


    } # end loop over years

  } # end loop over sites
  DD <- Dfish
  return(DD)
}
