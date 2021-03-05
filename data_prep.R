# Data preparation and exploration 
# 1) cluster Sites within Locations based on threat/habitat into good/degraded
# 2) Develop code to splot networks and extract network metrics

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
# run per taxa at the moment - run bipartite for inter-taxa networks?
# outputting network metrics per taxa will allow us to test effect of Location/Site/habitat (Cluster) and Impacts/Env
# try seperate sites into good/degraded habitat per Location, then compare networks
# Figure out functional implications post-network creation e.g. what FG are hub species
# use Maarten's cluster code as template

####~~~~ Cluster sites into good/degraded habitat per site ~~~~####
# focus on CPCe data
# check variance within site, can we average across transects within site?
# check appropriate distance measure, use hclust and PCoA to visualise

####~~~~ * ~~~~~#####

# assign clusters to species data

dat$Loc_class<-'Lucipara_pristine'
dat[dat$Code %in% paste0('SPE', c(1,2,3,4,5,6)),]$Loc_class<-'Spermonde_good'
dat[dat$Code %in% paste0('SPE', c(7,8,9,10,11,12)),]$Loc_class<-'Spermonde_poor'
dat[dat$Code %in% paste0('AMB', c(1,2,3,4,5,6)),]$Loc_class<-'Ambon_good'
dat[dat$Code %in% paste0('AMB', c(7,8,9,10,11,12)),]$Loc_class<-'Ambon_poor'
dat[dat$Code %in% paste0('HAL', c(1,2,3,4,5,6)),]$Loc_class<-'Halmahera_good'
dat[dat$Code %in% paste0('HAL', c(7,8,9,10,11,12)),]$Loc_class<-'Halmahera_poor'

# Run fish data first

fish<-dat[dat$Kingdom=='Fish',]
#sum invidual observations
fish_sum<-fish%>%group_by(Loc_class, code_trans, Species)%>%summarise(sumAbun=sum(Count))%>%as.data.frame()
# summarise fish per Location
fish_loc_sum<-fish%>%group_by(Location, Species)%>%summarise(sumAbun=sum(Count))%>%as.data.frame()


#make networks per treatment per site 

fish_mat<-matrify(fish_sum[fish_sum$Loc_class == 'Halmahera_good', 2:4])
fish_mat<-t(fish_mat)#transpose so columns become rows for cooccur
fish_mat<-fish_mat[-which(rowSums(fish_mat)<2),] #remove rare species with 1 observation across all transects
fish_mat[fish_mat > 0] = 1 # convert to P/A

fish_co = cooccur(mat=fish_mat, 
                  type="spp_site", 
                  thresh=TRUE, 
                  spp_names=TRUE, 
                  prob="comb")

fish_eff<-effect.sizes(fish_co, standardized = TRUE, matrix = FALSE)
names(fish_eff)[1:2]<-c('spe1', 'spe2')

fish_out<-cbind(prob.table(fish_co), fish_eff)

#add columns for igraph
fish_out$sig<-ifelse(fish_out$p_gt<0.05 | fish_out$p_lt<0.05, 1, 0)

#format for igraph 

ig<-graph_from_data_frame(fish_out[,12:15], directed=F)

#if omitting non observed interactions then do clustering without weights
# if omitting non sig interactions then do cluster with weights
# if omitting non sig and negative effect interactions then do cluster with weights

#ig<-graph_from_data_frame(fish_out[fish_out$obs_cooccur>0,12:15], directed=F)
#ig<-graph_from_data_frame(fish_out[fish_out$sig==1,12:15], directed=F)
#ig<-graph_from_data_frame(fish_out[fish_out$sig==1 & fish_out$effects>0,12:15], directed=F)

# Select only edges that have significant positive cooccurence,
# but include all ns species as unconnected nodes
ig<-graph_from_data_frame(fish_out[fish_out$sig==1 & fish_out$effects>0,12:15], directed=F,
                          vertices=unique(fish_out$sp1_name))

length(unique(fish_out$sp1_name)) # n all species
# n species with at least 1 sig positive coocurrence
length(unique(fish_out[fish_out$sig==1 & fish_out$effects>0,12:15]$spe1))

#calc degree (n connections) and strength (sum edge weight 'effects') and
# eigen centrality (connected to other connected nodes aka 'hubs')
node_metrics<-data.frame(species=V(ig)$name, degree=as.vector(degree(ig)), 
           strength=as.vector(strength(ig, weights = E(ig)$effects)),
           eig_cent=as.vector(eigen_centrality(ig, weights = E(ig)$effects,
                                               scale=TRUE)[[1]]))


#calculate modularity springlass without weights
springL_noweight <- cluster_spinglass(ig, weights = NA, spins = 100, parupdate = FALSE,
                                           start.temp = 1, stop.temp = 0.01, cool.fact = 0.99, update.rule = "config", gamma = 1, implementation = "orig", gamma.minus = 1)
# Spinglass with negative weights
springL_weight <- cluster_spinglass(ig, weights = E(ig)$effects, spins = 100, parupdate = FALSE,
                                           start.temp = 1, stop.temp = 0.01, cool.fact = 0.99, update.rule = "config",
                                    gamma = 1, implementation = "neg", gamma.minus = 1)

# Spinglass with positive-only weights
springL_weight <- cluster_spinglass(ig, weights = E(ig)$effects, spins = 100, parupdate = FALSE,
                                    start.temp = 1, stop.temp = 0.01, cool.fact = 0.99, update.rule = "config",
                                    gamma = 1, implementation = "orig", gamma.minus = 1)


NewmanGir_weight <-edge.betweenness.community(ig, weights = E(ig)$effects, directed=F) 


la <- layout_in_circle(ig, order=order(membership(NewmanGir_weight)))

par(mar=c(8,6,6,6))


plot(ig, layout=la, edge.width=abs(E(ig)$effects),
     edge.color=adjustcolor(ifelse(E(ig)$effects<0, 'grey', 'black'),0.5), vertex.size=2,
     vertex.label="", edge.curved=0.2)

## Apply labels manually
#Specify x and y coordinates of labels, adjust outward as desired
x = la[,1]*1.3
y = la[,2]*1.3

#create vector of angles for text based on number of nodes (flipping the orientation of the words half way around so none appear upside down)
angle = ifelse(atan(-(la[,1]/la[,2]))*(180/pi) < 0,  90 + atan(-(la[,1]/la[,2]))*(180/pi), 270 + atan(-la[,1]/la[,2])*(180/pi))


pal1<-sample(c("#564138","#2e86ab","#f6f5ae","#f5f749","#f24236","#f18805","#f7996e","#285238","#138a36","#5f00ba",
  "#b0d0d3","#c08497","#f7af9d","#f7e3af","#f3eec3","#393e41","#d6ff79","#5f00ba","#3f3244","#2f2235",
  "#ddd78d","#dcbf85","#8b635c","#60594d","#93a29b","#050517","#197bbd","#125e8a","#204b57","#610f7f",
  "#3a405a","#aec5eb","#f9dec9","#e9afa3","#685044","#40c9a2","#e55381","#a31621","#db222a","#767522",
  "#628395","#96897b","#dbad6a","#cf995f","#d0ce7c","#f5efff","#48e5c2","#f991cc","#092327","#a8201a"),
  39)


pal2<-sample(rainbow(39), 39)
#Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
  text(x=x[i], y=y[i], labels=attr(membership(NewmanGir_weight)[i],'names'),
  adj=NULL, pos=NULL, cex=.7, col=pal1[membership(NewmanGir_weight)[i]], srt=angle[i], xpd=T)
}


#Apply the text labels with a loop with angle as srt
for (i in 1:length(x)) {
  text(x=x[i], y=y[i], labels=V(ig)[order(membership(NewmanGir_weight))]$name[i], adj=NULL, pos=NULL, cex=.7, col="black", srt=angle[i], xpd=T)
}

#https://gist.github.com/ajhmohr/5337a5c99b504e4a243fad96203fa74f


# igraph plotting

ig1<-graph_from_data_frame(fish_out[,12:15], directed=F,
                           vertices=filter(fish_loc_sum, Location=='Halmahera')%>%select(Species))


ig2<-delete_edges(ig1, which(E(ig1)$sig ==0))

la <- layout_in_circle(ig2)

plot(ig2, layout=la, edge.width=abs(E(ig2)$effects),
     edge.color=ifelse(E(ig2)$effects>0, 'dark green', 'dark red'), vertex.size=1.5,
     vertex.label="")


#heatmap of positive significant interactions only
heatmap(get.adjacency(graph_from_data_frame(fish_out[fish_out$sig>0 &
                                                       fish_out$effects>0,12:15], directed=F), attr="effects", sparse=F))

### OLD ###


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
