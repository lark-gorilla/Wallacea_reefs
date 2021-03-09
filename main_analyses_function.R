# function to carry out paper analysis per Location
# 1) split Location into good and degraded sites = 'treatments'
# 2) run species accumulation curve to get rarefied species richness per treatment
# 3) run jaccard dist/PCoA on Location, colour PCoA by treatment 
# 4) run cooccur and igraph functions per treatment
# 5) per treatment calculate: node metrics (degree, strength, eig centrality)
# 6) summarise per treatment network:
# Modularity, skewness of degree and strength (Kay & Mont-Cent)
# connectance (sum of sp degree)/(n nodes^2) (Kay, sensu Gilbert)
# OR/AND links per sp (sum of sp degree)/(n species) (Kay)
# OR/AND Network density is the portion of potential connections
# in a network that are realized (Mont-Cent sensu Wasserman & Faust 1994)
# rarefied species richness, species prevelence
# 7) combine treatment network results, pairwise sp connections
# Change in pairwise species co-occurrence connections (degree) (Kay)
# Stable (or no link) VS  Links restructured + breakdown:
# Links lost (because species lost from modified landscape)
# Links gained (because species gain in modified landscape)
# Links lost (species present across both landscapes)
# Links gained (species present across both landscapes)
# 8) PLot change in sp prevelence against change in sp degree (Kay)
# colour by quad-sector AND by funtional group (+ kernel)

compare_comms<-function(surv_dat=sp_abun_dat, good_sites='my_good_sites', poor_sites='my_poor_sites',
                         run_spaccu=TRUE, run_multiv=TRUE, run_netwrk=TRUE)
{  
  # Setup Location level sp occurrence matrix
  
  fish_mat_all<-matrify(surv_dat[c('code_trans', 'Species', 'sumAbun')])
  fish_mat_all[fish_mat_all > 0] = 1 # convert to P/A
  fish_mat_all<-fish_mat_all[,-which(colSums(fish_mat_all )<2)] 
  #remove rare species with 1 observation across all transects (~1% occ Kay, Tulloch)

  # do PCoA
  #all_nmds <- metaMDS(fish_mat_all, distance="jaccard", binary=T,
  #                              trace=FALSE, trymax=100)
  
  #pc2<-dudi.pco(d =vegdist(fish_mat_all, method='jaccard', binary=T),
  #              scannf = FALSE, nf = 4)
  
  pc2<-dudi.pco(d =dist.binary(fish_mat_all, method=1),
                scannf = FALSE, nf = 4)
  
  pc2.dfs <- data.frame(pc2$li, data.frame(code_trans=dimnames(fish_mat_all)[[1]])) # combine PCoA axes with trait/site data
  pc2.dfs<-left_join(pc2.dfs, surv_dat%>%group_by(code_trans)%>%summarise(Loc_class=first(Loc_class)), 
                     by='code_trans')
  
  #make global hull
  glob_hull<-pc2.dfs[chull(pc2.dfs$A1, pc2.dfs$A2),]
  
  # code to setup ggplot enviroment
  ppp <- ggplot() + coord_fixed() +
    labs(x="PCoA1, Axis1", y="PCoA2, Axis2") +
    geom_hline(yintercept=0, col="darkgrey") +
    geom_vline(xintercept=0, col="darkgrey") 
  
  pco_plot<-ppp+
    geom_polygon(data=glob_hull,aes(x=A1,y=A2),fill=NA,colour="grey70")+
    geom_point(data=pc2.dfs, aes(x=A1, y=A2, colour=Loc_class))+
    geom_polygon(data=pc2.dfs%>%group_by(Loc_class)%>%slice(chull(A1, A2)),
                 aes(x=A1,y=A2, fill=Loc_class),alpha=0.08,colour="black")+
   theme_bw()
  
  
  # do species accumulation curves
  spacc_good<-specaccum(fish_mat_all, method='exact', subset=
                   row.names(fish_mat_all)%in% 
                   surv_dat[surv_dat$Loc_class==good_sites,]$code_trans)
  spacc_poor<-specaccum(fish_mat_all, method='exact', subset=
                          row.names(fish_mat_all)%in% 
                          surv_dat[surv_dat$Loc_class==poor_sites,]$code_trans)
  
  spacc_pdat=rbind(data.frame(Sites=spacc_good$sites, Sp_richness=spacc_good$richness, 
             sd=spacc_good$sd, Loc_class=good_sites), 
             data.frame(Sites=spacc_poor$sites, Sp_richness=spacc_poor$richness, 
                        sd=spacc_poor$sd, Loc_class=poor_sites))
  
  spacc_plot<-ggplot(data=spacc_pdat, aes(x=Sites, y=Sp_richness))+
           geom_ribbon(aes(ymin=Sp_richness-sd, ymax=Sp_richness+sd,fill=Loc_class), alpha=0.5)+ 
             geom_line(aes(linetype=Loc_class),size=1)+theme_bw()+
        ylab('Species richness')
  
  
  for(i in c(good_sites, poor_sites))
  
  fish_mat<-fish_mat_all[row.names(fish_mat_all)%in% 
    surv_dat[surv_dat$Loc_class==i,]$code_trans,] # select gd/poor sites
  
  fish_mat<-t(fish_mat)#transpose so columns become rows for cooccur

  
  fish_co = cooccur(mat=fish_mat, 
                    type="spp_site", 
                    thresh=TRUE, 
                    spp_names=TRUE, 
                    prob="comb")
  
  fish_eff<-effect.sizes(fish_co, standardized = TRUE, matrix = FALSE)
  names(fish_eff)[1:2]<-c('spe1', 'spe2')
  
  fish_out<-cbind(prob.table(fish_co), fish_eff)
  
  #add columns for igraph, only positive significance taken forward
  # interested in co-occurence not exclusion
  fish_out$sig<-ifelse(fish_out$p_gt<0.05 , 1, 0)
  
  #format for igraph 
  
  # Select only edges that have significant positive cooccurence,
  # but include all ns species as unconnected nodes (from threat level community)
  ig<-graph_from_data_frame(fish_out[fish_out$sig==1 ,12:15], directed=F,
                            vertices=rownames(fish_mat))

  #calc Node metrics degree (n connections) and strength (sum edge weight 'effects') and
  # eigen centrality (connected to other connected nodes aka 'hubs')
  node_metrics<-data.frame(species=V(ig)$name, degree=as.vector(degree(ig)), 
                           strength=as.vector(strength(ig, weights = E(ig)$effects)),
                           eig_cent=as.vector(eigen_centrality(ig, weights = E(ig)$effects,
                                                               scale=TRUE)[[1]]))
  
  # calc Network metrics
  mod_newm<-modularity(edge.betweenness.community(ig, weights = E(ig)$effects, directed=F)) 
  mod_clwk<-modularity(cluster_walktrap(ig, weights = E(ig)$effects))
  skewness(degree) skewness(strength)
  # connectance (sum of sp degree)/(n nodes^2) (Kay, sensu Gilbert)
  # OR/AND links per sp (sum of sp degree)/(n species) (Kay)
  # OR/AND Network density is the portion of potential connections
  # in a network that are realized (Mont-Cent sensu Wasserman & Faust 1994)
