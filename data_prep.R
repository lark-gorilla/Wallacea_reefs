# Data preparation and exploration 
library(dplyr)
library(cooccur)
library(igraph)
library(adonis) 
library(labdsv)
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

# OK so for a trial lets:
# Run on Pres-Abs data
# run per species not FG
# run per transect
# run per taxa at the moment
# outputting network metrics per taxa will allow us to test effect of Location/Site/habitat (Cluster) and Impacts/Env
# Figure out functional implications post-network creation e.g. what FG are hub species
# use Maarten's cluster code as template

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
