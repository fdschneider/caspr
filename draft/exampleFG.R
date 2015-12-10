# 
# This is an example file that illustrates the way to load, extract and 
#   process the output from the caspr simulations
# 
#####charge multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# We will need the following packages
library(devtools) 
library(ggplot2)


# Install the spatialwarnings package if it is not already installed (this 
# sometimes require restarting R).
if ( ! "spatialwarnings" %in% row.names(installed.packages())) { 
  install_github('fdschneider/spatial_warnings') # mind the "_"
}
library(spatialwarnings)

# Here we'll take the forestgap model as example
# 

# Let's merge the two branches together in one table
merge_two_branches <- function(up, low) { 
  data.frame(branch = c(rep("upper", nrow(up$DatBif)), 
                        rep("lower", nrow(low$DatBif))), 
             rbind(up$DatBif, low$DatBif))
} 
merged_results <- merge_two_branches(upper_branch, lower_branch)

#######FIG 1.i
# Do a graph "à la flo" (hard to read for now, needs more tweaking)
ggplot() + 
  geom_raster(aes(x = delta, y = d), fill = "grey70",
              data = subset(lower_branch$DatBif, mean_cover_. < .05)) + 
  geom_contour(aes(x = delta, y = d, z = mean_cover_.), 
               data = subset(upper_branch$DatBif), 
               color = "black", linesize = .1) 


#######
# Let's extract a 1D bifurcation diagram over delta (d in x-axis)
uny<-unique(merged_results$delta)
a1=list()
for (i in 1:length(uny)){
  data_subset <- subset(merged_results, delta==uny[i])
  a1[[i]]<-ggplot(data_subset) + 
    geom_point(aes(d, mean_cover_., color = branch)) +
    ggtitle(paste("delta =",uny[i]))
  print(a1[[i]])
  
}

####################################3duplicate
#Do it over d (delta in x-axis)
uny<-unique(merged_results$d)
a2=list()
for (i in 1:length(uny)){
  data_subset <- subset(merged_results, d==uny[i])
  a2[[i]]<-ggplot(data_subset) + 
    geom_point(aes(delta, mean_cover_., color = branch)) +
    ggtitle(paste("d =",uny[i]))
  print(a2[[i]])
  
}

#####
#So the transit will be
deltaSEL<-c(uny[2],uny[10])
#Selected transits (value of d for which we will subset the matrix)
dSEL<-c(uny[2],uny[9])


#subset transits of interest
var_dat <- subset(upper_branch$DatBif, d==dSEL)
idx<-which(upper_branch$DatBif$d==dSEL)
landscapes<-'['(upper_branch[[2]],idx)
#######################################################
# Computing common early warning: Skewness and Variance
#######################################################

#compute indicators. Output is length(idx)->nreplicates->result from generic_spews
genind<- generic_spews(landscapes,subsize = 8)

#extract summary. Output is length(idx)->summary tables
sumIND <- lapply(genind, spatialwarnings:::summary.generic_spews_list, 
       null_replicates = 0)
#extract variance. Output is length(idx)->summary tables, only variance
var<- lapply(sumIND, subset, indicator == 'Variance')
#transform into matrix and label the matrix
varDB<-t(sapply(1:length(var), function(x,y) rbind(y[[x]]$value),y=var))

##cleaner: varDB<-t(sapply(1:length(idx), function(x,y) rbind(y[[x]]$value),y=lapply(sumIND, subset, indicator == 'variance')))
colnames(varDB)<-row.names(var[[1]])

#extract Skewness Output is length(idx)->summary tables, only Skewness
Skew<- lapply(sumIND, subset, indicator == 'Skewness')
#transform into matrix and label the matrix
SkewDB<-t(sapply(1:length(Skew), function(x,y) rbind(y[[x]]$value),y=Skew))

##cleaner: SkewDB<-t(sapply(1:length(idx), function(x,y) rbind(y[[x]]$value),y=lapply(sumIND, subset, indicator == 'skewness')))
colnames(SkewDB)<-row.names(Skew[[1]])

#extract variance. Output is length(idx)->summary tables, only variance
MI<- lapply(sumIND, subset, indicator == "Moran's I")
#transform into matrix and label the matrix
MIDB<-t(sapply(1:length(MI), function(x,y) rbind(y[[x]]$value),y=MI))

##cleaner: varDB<-t(sapply(1:length(idx), function(x,y) rbind(y[[x]]$value),y=lapply(sumIND, subset, indicator == 'variance')))
colnames(MIDB)<-row.names(MI[[1]])





#Database format
#var_dat  <- cbind(var_dat,varDB,SkewDB)

##calculating avg +std
varDBmc<-cbind(apply(varDB,1,mean),apply(varDB,1,function(x) mean(x)-t.test(x)$conf.int[1]))
colnames(varDBmc)<-c("var_mean","var_CI")
SkewDBmc<-cbind(apply(SkewDB,1,mean),apply(SkewDB,1,function(x) {if(sum(is.na(x))<8){mean(x)-t.test(x)$conf.int[1]}else{NA}}))
colnames(SkewDBmc)<-c("Skew_mean","Skew_CI")
MIDBmc<-cbind(apply(MIDB,1,mean),apply(MIDB,1,function(x) {if(sum(is.na(x))<8){mean(x)-t.test(x)$conf.int[1]}else{NA}}))
colnames(MIDBmc)<-c("MI_mean","MI_CI")

#Database format
var_dat  <- cbind(var_dat,varDBmc,SkewDBmc,MIDBmc)

#Plotting
DBplot1<-subset(var_dat,d==dSEL[1])
DBplot2<-subset(var_dat,d==dSEL[2])

p1<-ggplot( DBplot1 , aes(delta,var_mean))+
  geom_line()+
  geom_ribbon(data=DBplot1,aes(ymin= var_mean - var_CI , ymax= var_mean + var_CI), alpha=0.3) +
  ggtitle(paste("Variance in transit 1. d =",dSEL[1]))

p2<-ggplot (DBplot2 ,aes(delta,var_mean))+
  geom_line()+
  geom_ribbon(data=DBplot2,aes(ymin=var_mean-var_CI,ymax=var_mean+var_CI),alpha=0.3) +
  ggtitle(paste("Variance in transit 2. d =",dSEL[2]))


p3<-ggplot(DBplot1 , aes(delta,Skew_mean))+
  geom_line()+
  geom_ribbon(data=DBplot1,aes(ymin=Skew_mean-Skew_CI,ymax=Skew_mean+Skew_CI),alpha=0.3) +
  ggtitle(paste("Skewness in transit 1. d =",dSEL[1]))


p4<-ggplot(DBplot2 , aes(delta,Skew_mean))+
  geom_line()+
  geom_ribbon(data=DBplot2,aes(ymin=Skew_mean-Skew_CI,ymax=Skew_mean+Skew_CI),alpha=0.3) +
  ggtitle(paste("Skewness in transit 2. d =",dSEL[2]))

p5<-ggplot(DBplot1 , aes(delta,MI_mean))+
  geom_line()+
  geom_ribbon(data=DBplot1,aes(ymin=MI_mean-MI_CI,ymax=MI_mean+MI_CI),alpha=0.3) +
  ggtitle(paste("Moran's I in transit 1. d =",dSEL[1]))

p6<-ggplot(DBplot2 , aes(delta,MI_mean))+
  geom_line()+
  geom_ribbon(data=DBplot2,aes(ymin=MI_mean-MI_CI,ymax=MI_mean+MI_CI),alpha=0.3) +
  ggtitle(paste("Moran's I in transit 2. d =",dSEL[2]))


multiplot(p1,p2,p3,p4,p5,p6,cols=3)

#######################################################
# Computing PSD
#######################################################
#Merging patch sizes
lcpsd<-list()
for(i in 1:length(landscapes)){
  try (lcpsd[[i]]<-indicator_cumpsd(landscapes[[i]]))
  
}
