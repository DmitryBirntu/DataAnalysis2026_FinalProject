library(vegan)
#--------------------start-------------------------------
# Get current working directory
getwd()
#----------------read dataset----------------------------
data<-read.table("data.txt",header=TRUE,sep="\t",check.names = FALSE)
summary (data)
data
rownames(data)<-data[,1]
data<-data[,-1]
data
# Non-metric multidimensional scaling (NMDS) using Euclidean distances
# Perform NMDS ordination based on Euclidean distances.
# Euclidean distance is suitable for quantitative data but may be less appropriate
# for ecological community data due to sensitivity to double zeros.
# The metaMDS function automatically runs multiple random starts and scales the solution.
ord <- metaMDS(data, distance = "euclidean")

#Visualisation of NMDS results
plot(ord, type = "n")                     # Create empty plot (axes only)
points(ord, disp = "sites", pch = 21, 
       cex = 2.5, lwd = 2.5, col = "red") # Add site points (samples)
text(ord, display = "site", cex = 0.7, 
     col = "red", pos = 3)                         # Label sites with their names


# Cluster analysis with Euclidean distances
# Compute Euclidean distance matrix between sites
d<-vegdist(data,method="euclidean")
# Perform hierarchical agglomerative clustering using average linkage (UPGMA)
fit<-hclust(d, method="average")

# Visualise the dendrogram with labels aligned at the baseline
plot(fit, hang =-1)
# Simple dendrogram plot (default parameters)
plot(fit)


#/////////////////////////////////////////////////////////////
#Non-metric multidimensional scaling (NMDS) with Bray–Curtis dissimilarities
# Bray–Curtis is a standard ecological distance measure that ignores double zeros
# and is robust for abundance data.
# The metaMDS function automatically runs multiple random starts and scales the final solution.
ord <- metaMDS(data, distance = "bray")

#Visualisation of NMDS results
plot(ord, type = "n")
points(ord, disp="sites", pch=21, cex=2.5, lwd=2.5, col = "red")
text(ord, display = "site", cex=0.7, col="red", pos = 3)

# Cluster analysis with Bray–Curtis distances
# Compute Bray–Curtis dissimilarity matrix (standard for ecological community data)
d<-vegdist(data,method="bray")
# Perform hierarchical agglomerative clustering using average linkage (UPGMA)
fit<-hclust(d, method="average")

# Visualise dendrogram with labels aligned at the same horizontal level
plot(fit, hang =-1)
# Simple dendrogram plot (default R style)
plot(fit)


#/////////////////////////////////////////////////////////////
# Non-metric multidimensional scaling (NMDS) with Jaccard distance
ord <- metaMDS(data, distance = "jaccard")


plot(ord, type = "n")
points(ord, disp="sites", pch=21, cex=2.5, lwd=2.5, col = "red")
text(ord, display = "site", cex=0.7, col="red", pos = 3)

# Cluster analysis with Jaccard distance
# Compute Jaccard dissimilarity matrix
d<-vegdist(data,method="jaccard")
# Hierarchical clustering using average linkage (UPGMA)
fit<-hclust(d, method="average")
# Visualise dendrogram with labels aligned
plot(fit, hang =-1)
# Default dendrogram plot
plot(fit)

#---------------------HOMEWORK--------------------------

ord_bray  <- metaMDS(data, distance = "bray", trymax = 50)

set.seed(123)
fit_sp <- envfit(ord_bray, data, permutations = 999)
sig_sp <- fit_sp$vectors$p.val <= 0.05
cat("значимые виды ( p <= 0.05):\n")
print(fit_sp$vectors$arrows[sig_sp, , drop = FALSE])
print(fit_sp$vectors$p.val[sig_sp])

d_bray <- vegdist(data, method = "bray")
hc_bray <- hclust(d_bray, method = "average")

k_clusters <- 2
clusters <- cutree(hc_bray,k = k_clusters)
clusters <- factor(clusters, labels = c("cluster 1","cluster 2"))

cols <- c("blue","red")[as.numeric(clusters)]
plot(ord_bray,type = "n",
     main = paste0("NMDS (Bray-curtis), stress =", round(ord_bray$stress, 4)),
     sub = "кластеры выделены по UPGMA (k=2)")
points(ord_bray, col = cols, pch = 16, cex = 2.5)
ordiellipse(ord_bray, groups = clusters, col = c("blue","red"),
            kind = "sd", conf = 0.95, lwd = 2)
orditorp(ord_bray, display = "sites", cex = 0.8, col = "black",air = 0.5)
plot(fit_sp, p.max = 0.05, col ="darkgreen", cex = 0.9, add = TRUE)
legend("topright", legend = levels(clusters), col = c("blue", "red"),
       pch = 16, title = "cluster", bty = "n")


set.seed(456)
adonis_result <- adonis2(data ~ clusters, method = "bray", permutations = 999)

print(adonis_result)

R2 <- adonis_result$R2[1]
p_val <- adonis_result$'Pr(>F)'[1]

cat("\n=============================\n")
cat("Permanova results :\n")
cat("R squared =  ",round(R2,4), "\n")
cat("p-value = ", p_val, "\n")
if(p_val < 0.05){
  cat("Вывод: различия между кластерами статистически значимы")
} else {
  cat("Вывод: различия между кластерами не значимы")
}
cat("=======================\n")

cat("\nСостав кластеров:\n")
for (i in 1:k_clusters) {
  cat(names(clusters[clusters == levels(clusters)[i] ]),"\n")
}

