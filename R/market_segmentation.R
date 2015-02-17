#############################################
### Market Segmentation of Quantcast data ###
#############################################

#### The purpose of this study is cluster DMAs or 'Designated Market Areas'
#### using Quantcast data for a variety of sites.
#### A handful of sites with Quantcast data were chosen due to data availability,
#### but more sites can be added.

#### Hierarchical and K-means clusterings are tested, but there are other
#### algorithms that warrant investigation (i.e. E-M algorithms in 'mclust')

setwd("~/GitHub/Market Segmentation")

options(stringsAsFactors = F)

#### list of sites that were used in the analysis
#### looking at the Gawker, Vox, and Onion Media Networks
#### Unfortunately vox.com, my personal favorite within these networks, does not have msa data available
#### Also included are a handful of NBC websites, and 3 sites chosen for the affinity
#### to Gawker's cink.hu, Vox's vox.com, and The Onion's gameological.com, all of which had insufficient data
sites <- list("deadspin", "gawker", "gizmodo", "io9", "jalopnik", "jezebel", "kinja", "kotaku", "lifehacker",
              "sbnation", "theverge", "polygon", "curbed", "eater", "racked",
              "theonion", "avclub", "clickhole",
              "bringatrailer", "fool", "kickstarter",
              "nbcnews", "nbcsports", "eonline", "nbc", "msnbc", "cnbc")

library(data.table)
library(plyr)
library(dplyr)
library(magrittr)
library(cluster)
library(fastcluster)
library(ggplot2)
library(fpc)
library(DT)
library(rmarkdown)

#### each dataset is unique by DMA or metropolitan area
#### need to fix column names
dmas <- lapply(sites, function(site) {
	data <- fread(paste0("DMAs/", site, "_dma.csv"))
	setnames(data,
		   c("dma",
		     paste0("uniques_", site),
         paste0("pct_uniques_", site), paste0("index_uniques_", site) ,
		     paste0("pct_impressions_", site) , paste0("index_impressions_", site)))
	setkey(data, dma)
	select(data, dma, contains("index"))
})

dmas

#### Merge all datasets into one data.table
all_dmas <- Reduce(merge, dmas)
### should be 210 media markets
tbl_dt(all_dmas)

#### currently only using euclidean, but other distances can be used
distances <- list("euclidean")

#### We'll start with hierarchical clustering
#### all the hierarchical clustering methods
#### note, this script uses 'fastclusters' 'hclust' function,
#### which returns the same results as base R 'hclust', but with less computation time
hc_methods <- list("ward.D", "ward.D2", "single", "average", "mcquitty", "centroid")
cluster_numbers <- 2:30

#### test a variety of distances (just euclidean for now), clustering methods,
#### with various numbers of clusters.
#### test statistic is the silhouette coefficient, which has some issues,
#### (https://web.njit.edu/~yl473/papers/ICDM10CLU.pdf)
#### but will suffice for the time being.
#### average width of 1 indicates 'perfect' clusters, -1 completely 'imperfect'
hc_clusters <- lapply(distances, function(d) {
  dis <- dist(select(all_dmas, -dma), method = d,)
  lapply(hc_methods, function(hc) {
    lapply(cluster_numbers, function(n) {
      clus <- hclust(dis, method = hc)
      cut <- cutree(clus, n)
      sil <- summary(silhouette(cut, dis))$avg.width
      data.table(di = d, clus_method = hc, num_clus = n, silhou = sil)
    })
  })
})

#### turn the nested list into a data.table using magrittr and dplyr
hc_clusters %<>% lapply(., function(d) {
  lapply(d, function(hc) {
    rbindlist(hc)
  })
}) %>% lapply(function(d) {
  rbindlist(d)
}) %>% rbindlist()

arrange(hc_clusters, silhou)

#### boxplot of silhouette coefficients by hierarchicial clustering algorithm
ggplot(data = hc_clusters, aes(factor(clus_method), silhou)) + 
  geom_boxplot() + geom_jitter()
#### hard to say which hierarchical clustering method has the most potential,
#### but 'average' seems to have a slight edge

#### the clusGap function requires argumetns of k and cluster
new_hc_func <- function(x, k, hc, d) {
  list(cluster = cutree(hclust(dist(x, method = d), method = hc), k),
       k = k)
}

#### The Cluster Gap Statistic, developed by Robert Tibsharini
#### (http://web.stanford.edu/~hastie/Papers/gap.pdf) is more robust
#### than the silhouette coefficient, but more computationally expensive.
#### The silhouette coefficients were used to determine the best algorithm
cg_hc <- clusGap(select(all_dmas, -dma), FUN = new_hc_func, K.max = 30, B = 70,
                 hc = "average", d = "euclidean")
print(cg_hc, method = "globalSEmax")
#### Returns 30 clusters fairly reliably, but occasionally we see 14
plot(cg_hc)

#### now we'll try K-means clustering
km_methods <- list("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")

#### test a variety of distances (just euclidean for now), clustering methods,
#### with various numbers of clusters.
#### test statistic is the silhouette coefficient, which has some issues,
#### but will suffice for the time being.
#### average width of 1 indicates 'perfect' clusters, -1 completely 'imperfect'
km_clusters <- lapply(distances, function(d) {
  dis <- dist(select(all_dmas, -dma), method = d)
  lapply(km_methods, function(km) {
    lapply(cluster_numbers, function(n) {
      clus <- kmeans(select(all_dmas, -dma), algorithm = km,
                     centers = n, nstart = n, iter.max = 100)
      cut <- clus$cluster
      sil <- summary(silhouette(cut, dis))$avg.width
      data.table(di = d, clus_method = km, num_clus = n, silhou = sil)
    })
  })
})

km_clusters %<>% lapply(., function(d) {
  lapply(d, function(km) {
    rbindlist(km)
  })
}) %>% lapply(function(d) {
  rbindlist(d)
}) %>% rbindlist()

arrange(km_clusters, silhou)

#### boxplot of silhouette coefficients by k-means algorithm
ggplot(data = km_clusters, aes(factor(clus_method), silhou)) + 
  geom_boxplot() + geom_jitter()
#### hard to tell which k-means algorithm is best,
#### but 'Hartigan-Wong' seems to have a slight edge

#### Gap statistic again
cg_km <- clusGap(select(all_dmas, -dma), FUN = kmeans, K.max = 30, B = 70,
              algorithm = "Hartigan-Wong", nstart = 30, iter.max = 100)
print(cg_km, method = "globalSEmax")
#### Returns 27-29 clusters pretty reliably

#### Now let's look at the Calinski-Harabasz criterion of the best
#### hierarchical and kmeans results according to gap statistic
#### The Calinski-Harabasz criterion is still imperfect,
#### (https://web.njit.edu/~yl473/papers/ICDM10CLU.pdf)
#### but it is a strong internal evaluation measure
best_hc <- hclust(dist(select(all_dmas, -dma), method = "euclidean"),
                  method = "average")
hc_cut <- cutree(best_hc, 30)
best_km <- kmeans(select(all_dmas, -dma), algorithm = "Hartigan-Wong",
                  nstart = 29, centers = 29, iter.max = 100)
km_cut <- best_km$cluster

#### hc is ~23
calinhara(select(all_dmas, -dma), hc_cut)
#### km is ~37
calinhara(select(all_dmas, -dma), km_cut)

#### The algorithm that has the highest Calinski-Harabasz criterion
#### is the strongest clustering algorithm.
#### Thus, we choose K-Means Clustering (Hartigan-Wong method), with 27 clusters

#### Now let's see what DMAs are clustered together

all_dmas[, cluster := as.character(km_cut)]

display <- select(all_dmas, dma, cluster)
setnames(display, toupper(names(display)))

#### Let's see what our results are
DT::datatable(display)

#### Some fascinating findings:
#### a) San Francisco grouped with Los Angeles (NorCal-SoCal divide not as strong
#### as one might think) and Washington DC (West Coast meets East Coast)
#### b) Financial powerhouses Chicago and New York are found together
#### c) Southeastern powerhouses Atlanta and Miami are in the same cluster

#### Of course, these findings should be taken with a heavy grain of salt,
#### There is strong selection bias in the websites used, and we are not
#### examining the entire web. But what we can take away is that there are
#### not-so-obvious differences/similarities across geographic media markets,
#### even when it comes to a place that has no borders (in the US at least),
#### the Internet.