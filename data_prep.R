# 1) Exploration of data to see how current benthic habitat clustering groups transects
# 2) Demonstration of R function to perform analyses similar to Kay et al(2018)

# load libraries
library(dplyr)
library(ggplot2)
library(cooccur)
library(igraph)
library(labdsv)
library(vegan)
library(ade4)
library(e1071)

####~~~~ Explore data ~~~~####
#read survey data
dat<-read.csv('C:/wallacea/Files for Mark and Maria/Data & scripts/Visual survey data/Wallacea_combined data_ALL.csv')
dat$code_trans<-paste(dat$Code, dat$Transect, sep='_')

sites<-read.csv('C:/wallacea/Files for Mark and Maria/Data & scripts/Wallacea_sites.csv')

# attrib sites with cluster and write for GIS
sites<-left_join(sites, dat%>%group_by(Code)%>%mutate(clusts=paste(unique(Cluster), collapse=' '))%>%
                   summarise_all(first)%>%select(Code, clusts, Impact_Score, Impact_category),by=c('SiteCode'='Code'))

#write.csv(sites, 'C:/wallacea/Files for Mark and Maria/Data & scripts/Wallacea_sites.csv', quote=F, row.names=F)

# explore
dat%>%group_by(Location, Site, Transect, id, Kingdom)%>%summarise(n())%>%View()
dat%>%group_by(Location, Site, Kingdom)%>%summarise(n())%>%View()

#figure out clusters (do they span locations/sites?)
dat%>%group_by(Cluster, Location, Site)%>%summarise(n())
# They actually span transects within sites!
dat%>%group_by(Code)%>%summarise(n_clust=length(unique(Cluster)))

# The plan:
# Run on Pres-Abs data
# run per species not FG
# run per transect
# run per taxa at the moment - could run bipartite for inter-taxa networks?
# outputting network metrics per taxa will allow us to test effect of Location/Site/habitat (Cluster) and Impacts/Env
# try seperate sites into good/degraded habitat per Location, then compare networks
# Figure out functional implications post-network creation e.g. what FG are hub species

####~~~~ * ~~~~~#####

####~~~~ Cluster sites into good/degraded habitat per site ~~~~####
# focus on CPCe data
# check variance within site, can we average across transects within site?
# check appropriate distance measure, use hclust and PCoA to visualise

####~~~~ * ~~~~~#####

####~~~~ Demonstration of R function to perform Kay et al(2018) analyses ~~~~####

# assign good/poor categories to individual sites
# exmaple below is just random e.g. Spermonde sites 1-6
# are assigned good and sites 7-12 poor. Also you would want to 
# assign individual transects to good/poor categories rather than 
# all transects within a site

dat$Loc_class<-'Lucipara_pristine'
dat[dat$Code %in% paste0('SPE', c(1,2,3,4,5,6)),]$Loc_class<-'Spermonde_good'
dat[dat$Code %in% paste0('SPE', c(7,8,9,10,11,12)),]$Loc_class<-'Spermonde_poor'
dat[dat$Code %in% paste0('AMB', c(1,2,3,4,5,6)),]$Loc_class<-'Ambon_good'
dat[dat$Code %in% paste0('AMB', c(7,8,9,10,11,12)),]$Loc_class<-'Ambon_poor'
dat[dat$Code %in% paste0('HAL', c(1,2,3,4,5,6)),]$Loc_class<-'Halmahera_good'
dat[dat$Code %in% paste0('HAL', c(7,8,9,10,11,12)),]$Loc_class<-'Halmahera_poor'

# Run fish data first

fish<-dat[dat$Kingdom=='Fish',]

#sum invidual species observations at transect level
fish_sum<-fish%>%group_by(Location, Loc_class, code_trans, Species)%>%summarise(sumAbun=sum(Count))%>%as.data.frame()
# source function
source('C:/wallacea/Wallacea_reefs/main_analyses_function.R')
# run function per region
# function requires data (surv_dat) be formatted identically to fish_sum.
# Column1: 'Location'= name of region/location
# Column2: 'Loc_class'= good/poor category assigned to transect
# Column3: 'code_trans'= transect name
# Column3: 'Species'= species name
# Column3: 'sumAbun'= species abundance per transect

ambon_compare<-compare_comms(surv_dat=filter(fish_sum, Location=='Ambon'),
                             good_sites='Ambon_good', poor_sites='Ambon_poor') # run function!
# can ignore warnings
# returns a list of outputs
ambon_compare[[1]] # plot data for igraph network plot
ambon_compare[[2]] # node metrics
ambon_compare[[3]] # network metrics
ambon_compare[[4]] # species accumulation curve plot
ambon_compare[[5]] # Multivariate PCoA plot

# Make igraph plots for both threat treatments with little function
plot_network<-function(regional_compare_output=ambon_compare[[1]],site_type='Ambon_good')
{
ig_ready<-graph_from_data_frame(regional_compare_output[regional_compare_output$run==site_type,c(2:4)], directed=F,
                                 vertices=unique(c(regional_compare_output$spe1,regional_compare_output$spe2)))
  
  la <- layout_in_circle(ig_ready)
  par(mar=c(8,6,6,6))
  plot(ig_ready, layout=la, edge.width=abs(E(ig_ready)$effects),
       edge.color='black', vertex.size=2,
       vertex.label="", edge.curved=0.2)
  
  ## Apply labels manually
  #Specify x and y coordinates of labels, adjust outward as desired
  x = la[,1]*1.3
  y = la[,2]*1.3
  vat<-vertex.attributes(ig_ready)$`name`
  threat_lev_sp<-unique(c(regional_compare_output[regional_compare_output$run==site_type,]$spe1,
                          regional_compare_output[regional_compare_output$run==site_type,]$spe2))
  
  #create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
  angle = ifelse(atan(-(la[,1]/la[,2]))*(180/pi) < 0,  90 + atan(-(la[,1]/la[,2]))*(180/pi), 270 + atan(-la[,1]/la[,2])*(180/pi))
  
  #Apply the text labels with a loop with angle as srt
  for (i in 1:length(x)) {
    text(x=x[i], y=y[i], labels=vat[i],
         adj=NULL, pos=NULL, cex=.7, col=ifelse(vat[i]%in%threat_lev_sp, 'black', 'grey'), srt=angle[i], xpd=T)}
}

#make Ambon plots with plot_network function
plot_network(regional_compare_output=ambon_compare[[1]],site_type='Ambon_good')
plot_network(regional_compare_output=ambon_compare[[1]],site_type='Ambon_poor')

# run for Halmahera
halmahera_compare<-compare_comms(surv_dat=filter(fish_sum, Location=='Halmahera'),
                             good_sites='Halmahera_good', poor_sites='Halmahera_poor')
halmahera_compare[[3]]
halmahera_compare[[4]]
halmahera_compare[[5]]
plot_network(regional_compare_output=halmahera_compare[[1]],site_type='Halmahera_good')
plot_network(regional_compare_output=halmahera_compare[[1]],site_type='Halmahera_poor')

# run for Spermonde
spermonde_compare<-compare_comms(surv_dat=filter(fish_sum, Location=='Spermonde'),
                                 good_sites='Spermonde_good', poor_sites='Spermonde_poor')
spermonde_compare[[3]]
spermonde_compare[[4]]
spermonde_compare[[5]]
plot_network(regional_compare_output=spermonde_compare[[1]],site_type='Spermonde_good')
plot_network(regional_compare_output=spermonde_compare[[1]],site_type='Spermonde_poor')

### OLD  stuff ###

#heatmap of positive significant interactions only
heatmap(get.adjacency(graph_from_data_frame(fish_out[fish_out$sig>0 &
                                                       fish_out$effects>0,12:15], directed=F), attr="effects", sparse=F))


fish<-dat[dat$Kingdom=='Fish',]
#sum invidual observations
fish2<-fish%>%group_by(code_trans, Species)%>%summarise(sumAbun=sum(Count))%>%as.data.frame()
fish_mat<-matrify(fish2)

# cluster transects/sites
fish_mat_clust<-fish_mat
fish_mat_clust[fish_mat_clust > 0] = 1 
fish_mat_clust
bdist<-dist.binary(fish_mat_clust, method=1) # jaccard dist
#shouldn't use ward centroid or median methods for jaccard dist
plot(hclust(bdist, 'single'))  #


fish_mat<-t(fish_mat)#transpose so columns become rows for cooccur
fish_mat<-fish_mat[-which(rowSums(fish_mat)<2),] #remove rare species with 1 observation across all transects
fish_mat[fish_mat > 0] = 1 # convert to P/A

fish$site_gp<-'SP1'
fish[fish$Code %in% c("SPE5",  "SPE6",  "SPE7"  ,"SPE8"),]$site_gp<-'SP2'
fish[fish$Code %in% c("SPE9",  "SPE10",  "SPE11"  ,"SPE12"),]$site_gp<-'SP3'
fish[fish$Code %in% c("AMB1",  "AMB2",  "AMB3"  ,"AMB4"),]$site_gp<-'AM1'
fish[fish$Code %in% c("AMB5",  "AMB6",  "AMB7"  ,"AMB8"),]$site_gp<-'AM2'
fish[fish$Code %in% c("AMB9",  "AMB10",  "AMB11"  ,"AMB12"),]$site_gp<-'AM3'
fish[fish$Code %in% c("HAL1",  "HAL2",  "HAL3"  ,"HAL4"),]$site_gp<-'HL1'
fish[fish$Code %in% c("HAL5",  "HAL6",  "HAL7"  ,"HAL8"),]$site_gp<-'HL2'
fish[fish$Code %in% c("HAL9",  "HAL10",  "HAL11"  ,"HAL12"),]$site_gp<-'HL3'
fish[fish$Code %in% c("LUC1",  "LUC2",  "LUC3"  ,"LUC4"),]$site_gp<-'LC1'

ok2<-NULL
for(i in unique(fish$site_gp))
{
  trans<-unique(fish[fish$site_gp==i,]$code_trans)

fish_co = cooccur(mat=fish_mat[,which(dimnames(fish_mat)[[2]]%in%trans)], 
                  type="spp_site", 
                  thresh=TRUE, 
                  spp_names=TRUE, 
                  prob="comb")

signif.cooc.C2 = subset(prob.table(fish_co) , p_gt <= 0.049 | 
                         p_lt <= 0.049 ) 

fish_co2 = cooccur(mat=fish_mat[,which(dimnames(fish_mat)[[2]]%in%trans)], 
                   type="spp_site", 
                   thresh=TRUE, 
                   spp_names=TRUE, 
                          only_effects=TRUE, 
                          eff_standard=TRUE)

ok1<-NULL
for(j in 1:nrow(signif.cooc.C2))
{
ok1<-rbind(ok1, fish_co2[fish_co2$sp1==as.character(signif.cooc.C2[j,]$sp1_name) &
           fish_co2$sp2==as.character(signif.cooc.C2[j,]$sp2_name) ,])
}

ok2<-rbind(ok2, data.frame(sp_gp=i, mean_effect=mean(ok1$effects, na.rm=TRUE), 
                           sd_effect=sd(ok1$effects, na.rm=TRUE)))
print(i)
}


f2<-fish%>%group_by(site_gp)%>%summarise(Cluster=median(Cluster), Impact_Score=mean(Impact_Score),
                                         abun=sum(Count))%>%as.data.frame()
f2<-cbind(f2, ok2%>%arrange(sp_gp))

ggplot(data=f2, aes(x=Impact_Score, y=mean_effect))+geom_point()+geom_smooth(method='lm')+
  xlab('Impact score')+ylab('Mean co-occurence network strength')

ggplot(data=f2, aes(x=log(abun), y=mean_effect))+geom_point()+geom_smooth(method='lm')+
  xlab('log Abundance')+ylab('Mean co-occurence network strength')

ggplot(data=f2, aes(x=as.factor(Cluster), y=mean_effect))+geom_boxplot()+geom_point()+
  xlab('Benthic habitat cluster')+ylab('Mean co-occurence network strength')


# plotting
ig1<-graph_from_data_frame(ok1, directed=F)
plot(ig1, curved=T)

plot(ig1, curved=T, layout=layout.circle)


net.sp.C2 <- delete.edges(net.C2, which(E(net.C2)$weight <0.005))
