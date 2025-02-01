#' Posterier Predictive check
#'
#' @importFrom tidyverse
#'
#' @return List
#'
#' @export
#'
#'

post.pred.check <- function(Data, fix.param, Post.fit, ipm.burnin = FALSE, ipm.i = TRUE, ipm.r = TRUE){
  
# sampled posterior
  samp.dist <- data.frame(matrix(ncol = 1000, nrow = nrow(Data))) #set up data frame
  sample.posterior <- Post.fit$Values[sample(nrow(Post.fit$Values), 1000),] #sample 1000 from post.fit 

for (i in c(1:1000)) { #sample 1000 posterior
  
sample.cand <- sample.posterior[i,]
#take 
data.test <- Data
data.test[,c(1:ncol(data.test))] <- NA #use datastructure but no data fitting
samp.ipm <- run.IPM(fix.param = fix.param, cand.param = sample.cand, Data = data.test, burnin = ipm.burnin, ipm.i = ipm.i, ipm.r = ipm.r  )

#standardize and place in data frame
last.ipm <- samp.ipm$N[,ncol(samp.ipm$N)]*fix.param$correction[[ncol(Data)]] # abundance on the last year

#Poission distribution
# samp.dist[,i] <- rpois(n = 1:nrow(Data), lambda = last.ipm) #simulated count data 

# Negative binomial
samp.dist[,i] <- rnbinom(n = 1:nrow(Data), size = sample.cand[["nb.k"]], mu = last.ipm) #negative binomial
}
  samp.dist <- data.frame(samp.dist)
  
#  #split size at 20
#   samp.dist.20 <- samp.dist[c(20:nrow(samp.dist)),]
#   sam.dist.20.avg <- t(data.frame(colSums(samp.dist.20)))
#   
#   samp.dist.0 <- samp.dist[c(0:20),]
#   samp.dist.0.avg <- t(data.frame(colSums(samp.dist.0)))
#  
#   
#   # get quartile of each size 
#   samp.dist.quart.20 <- as.data.frame(matrix(ncol = 5,nrow = nrow(sam.dist.20.avg)))
#   names(samp.dist.quart.20) <- c('mean','lower_95','lower_50','upper_50','upper_95')
#   rownames(samp.dist.quart.20) <- "large"
#    samp.dist.quart.20['mean'] <- mean(as.numeric(sam.dist.20.avg[1,]))
#    samp.dist.quart.20['lower_95'] <- quantile(as.numeric(sam.dist.20.avg[1,]), 0.025)
#    samp.dist.quart.20['lower_50'] <- quantile(as.numeric(sam.dist.20.avg[1,]), 0.25)
#    samp.dist.quart.20['upper_50'] <- quantile(as.numeric(sam.dist.20.avg[1,]), 0.75)
#    samp.dist.quart.20['upper_95'] <- quantile(as.numeric(sam.dist.20.avg[1,]), 0.975)
# 
#   
#   samp.dist.quart.0 <- as.data.frame(matrix(ncol = 5,nrow = nrow(samp.dist.0.avg)))
#   names(samp.dist.quart.0) <- c('mean','lower_95','lower_50','upper_50','upper_95')
#   row.names(samp.dist.quart.0) <- "small"
#   samp.dist.quart.0['mean'] <- mean(as.numeric(samp.dist.0.avg[1,]))
#   samp.dist.quart.0['lower_95'] <- quantile(as.numeric(samp.dist.0.avg[1,]), 0.025)
#   samp.dist.quart.0['lower_50'] <- quantile(as.numeric(samp.dist.0.avg[1,]), 0.25)
#   samp.dist.quart.0['upper_50'] <- quantile(as.numeric(samp.dist.0.avg[1,]), 0.75)
#   samp.dist.quart.0['upper_95'] <- quantile(as.numeric(samp.dist.0.avg[1,]), 0.975)
#   
#   #join rows
# samp.quart <- data.frame(t(bind_rows(data.frame(samp.dist.quart.20),data.frame(samp.dist.quart.0)))) %>% 
#   mutate(stat = rownames(.)) %>% pivot_longer(cols = c(1:2), names_to = "size", values_to = "value") %>% 
#   pivot_wider(names_from = stat, values_from = value)
  
  #bin size classes to 10 bins
  if(!is.null( fix.param$ogive)){
  samp.dist.ogive <- samp.dist[c(fix.param$ogive:nrow(samp.dist)),]
  samp.bin <- samp.dist.ogive %>% mutate(size.bin = cut(c(fix.param$ogive :nrow(samp.dist)),breaks = 10)) #categorize bins
  samp.bin.name <- unique(samp.bin$size.bin)
  } else {
    samp.bin <- samp.dist %>% mutate(size.bin = cut(c(1:nrow(samp.dist)),breaks = 10)) #categorize bins
    samp.bin.name <- unique(samp.bin$size.bin)
  }
#integrate to bin
    samp.int.bin <- samp.bin %>%
      pivot_longer(cols = c(1:1000),names_to = "col",values_to = "dist") %>%
      summarise(sum.dist = sum(dist),.by = c(size.bin, col)) %>%  #summarize across a size bins
      mutate(sum.dist = sum.dist/fix.param$correction[[ncol(Data)]]) %>% #make abundance to density
      pivot_wider(names_from = col,values_from = sum.dist) %>%
      dplyr::select(!size.bin) #remove categroical bins

# get quartile of each size
  samp.dist.quart <- as.data.frame(matrix(ncol = 6,nrow = nrow(samp.int.bin)))
  names(samp.dist.quart) <- c('size','mean','lower_95','lower_50','upper_50','upper_95')
  samp.dist.quart$size <- samp.bin.name

  for (i in 1:nrow(samp.int.bin)) {
    samp.dist.quart[i,'mean'] <- mean(as.numeric(samp.int.bin[i,]))
    samp.dist.quart[i,'lower_95'] <- quantile(as.numeric(samp.int.bin[i,]), 0.025)
    samp.dist.quart[i,'lower_50'] <- quantile(as.numeric(samp.int.bin[i,]), 0.25)
    samp.dist.quart[i,'upper_50'] <- quantile(as.numeric(samp.int.bin[i,]), 0.75)
    samp.dist.quart[i,'upper_95'] <- quantile(as.numeric(samp.int.bin[i,]), 0.975)
  }
  
  
# density of last year of data by size
  # data.ipm <- Data %>% select(ncol(Data)) %>% #select the last column (last year of data)
  #   rename(count = everything()) %>%
  #   mutate(length = row.names(Data),
  #                   length = as.numeric(length)) %>% 
  # mutate(size = case_when(length<20 ~ "small" ,
  #                  length>=20 ~ "large")) %>%
  #   summarise(count = sum(count), .by = size) %>%
  #   mutate(density = count/ fix.param$correction[ncol(Data)])
  
# get quartile of each size 
  if(!is.null(fix.param$ogive)){
    data.ipm <- Data[fix.param$ogive:nrow(Data),]

data.ipm <- data.ipm %>% select(ncol(Data)) %>% #select the last column (last year of data)
  rename(count = everything()) %>% #rename so universial naming scheme
  mutate(length = fix.param$ogive:nrow(Data),
         length = as.numeric(fix.param$ogive:nrow(Data))) %>%
  mutate(size.bin = cut(c(fix.param$ogive:nrow(Data)),breaks = 10)) %>%  #set into bins
  summarise(sum.bins = sum(count), .by = size.bin) %>%
  mutate(length = as.numeric(c(1:10)),
         density = sum.bins/fix.param$correction[[ncol(Data)]]) # rename to have bins be numerical
  }else{
    data.ipm <- Data %>% select(ncol(Data)) %>% #select the last column (last year of data)
      rename(count = everything()) %>% #rename so universial naming scheme
      mutate(length = row.names(Data),
             length = as.numeric(length)) %>%
      mutate(size.bin = cut(c(1:nrow(.)),breaks = 10)) %>%  #set into bins
      summarise(sum.bins = sum(count), .by = size.bin) %>%
      mutate(length = as.numeric(c(1:10)),
             density = sum.bins/fix.param$correction[[ncol(Data)]]) # rename to have bins be numerical
    
  }
  
  
  # uncount(sum.bins) # setup to be used in denstiy

# den.dist <- density(data.ipm$length) #density distribution at that year
# af <- approxfun(den.dist$x,den.dist$y) #linear approximation of the values from density plot
# 
# mesh.len <- c(1:nrow(samp.int.bin)) #length interval
# approx.den <- af(mesh.len) # approximate density based on size bins
# approx.den[ is.na(approx.den)  ] = 1e-323 #give rest a minimum value
# std.den <- approx.den/sum(approx.den*fix.param$dx) # unity
#
# data.den <- data.frame(std.den) #make into dataframe for graphing

#plot poter simualtion + data
plot <- ggplot()+
  # geom_point(data = samp.dist.quart, mapping = aes(x = size, y =mean),color = "lightblue", size = 2)+
  geom_pointrange(data = samp.dist.quart, mapping = aes(x = size, y =mean,ymin = lower_95, ymax = upper_95), color = "skyblue2", linewidth = 1.5)+ #95% confidence Interval
  geom_pointrange(data = samp.dist.quart, mapping = aes(x = size, y =mean,ymin = lower_50, ymax = upper_50), color = "skyblue4", linewidth = 2, alpha = 0.5)+ #50% confidence Interval
  # geom_line(data = samp.dist %>% mutate(size = c(1:nrow(samp.dist))) %>%  pivot_longer(cols = c(1:1000),
  #                                             names_to = "iteration",
  #                                             values_to = "density"),
  #           aes(x = size, y = density, group = iteration), color = "grey66")+
  # geom_point(data = data.ipm, aes(x = size, y = density),color = "#9A32CD")+
  geom_point(data.ipm %>% mutate(size = size.bin),
  mapping = aes(x = size, y = density),color = "#9A32CD", size = 2)+
  # labs(x = "Size (20cm)")+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1,size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank())

return(plot)

}
