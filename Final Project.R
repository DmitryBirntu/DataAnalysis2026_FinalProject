#пакеты
install.packages("skimr", dependencies = TRUE)
install.packages("DataExplorer", dependencies = TRUE)
install.packages("Hmisc", dependencies = TRUE)
install.packages("vegan", dependencies = TRUE)
install.packages("factoextra", dependencies = TRUE)
install.packages("knitr", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("tidyr", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("cluster", dependencies = TRUE)

library(skimr)
library(DataExplorer)
library(Hmisc)
library(vegan)
library(factoextra)
library(knitr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cluster)
#----------ПОДГОТОВКА ДАННЫХ--------------
#читаем датасет
omul <- read.table("data_omul.txt",
                   header = TRUE,
                   sep = "\t",
                   dec = ".",
                   stringsAsFactors = FALSE)


str(omul)
omul$sex <- as.factor(omul$sex)
summary(omul)
skim(omul)

#количественные переменные
num_vars <- c("con_28", "con_52", "con_101",
              "con_118_d", "con_153", "con_138", "con_180")


#проверка на пропуски
sum(is.na(omul))
colSums(is.na(omul))


#EDA отчёт
create_report(
  data = omul,
  output_file = "EDA_omul.html",
  output_dir = getwd(),
  report_title = "EDA Report: Omul Dataset"
)
#---------------------Аналитические методы-------------------------
#описательный анализ
#гистограммы
par(mfrow = c(3, 3))
for (v in num_vars) {
  hist(omul[[v]],
       main = paste("Histogram of", v),
       xlab = v,
       col = "lightblue")
}
par(mfrow = c(1, 1))

#боксплоты
par(mfrow = c(3, 3))
for (v in num_vars) {
  boxplot(omul[[v]],
          main = paste("Boxplot of", v),
          col = "tomato")
}
par(mfrow = c(1, 1))

#выбросы
outliers_data <- omul %>%
  select(where(is.numeric)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

ggplot(outliers_data, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(title = "Outlier Detection",
       x = "Variables",
       y = "Value") +
  theme_minimal()

ggplot(outliers_data, aes(y = value)) +
  geom_boxplot(fill = "tomato", alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  labs(title = "Boxplots for all numeric variables") +
  theme_minimal()
#----------------Проверка гипотезы---------------------------
#проверка на нормальность распределения
normality_results <- data.frame(
  variable = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (v in num_vars) {
  shapiro_res <- shapiro.test(omul[[v]])
  normality_results <- rbind(
    normality_results,
    data.frame(variable = v, p_value = shapiro_res$p.value)
  )
}

normality_results
kable(normality_results, digits = 4,
      caption = "Shapiro-Wilk normality test results")

# Если p < 0.05, распределение не нормальное.


#Выбросы по Z-score
omul_num <- omul[, num_vars]
omul_z <- scale(omul_num)

# наблюдения с минимум 1 |z| > 3
outlier_idx <- apply(abs(omul_z) > 3, 1, any)
table(outlier_idx)

#новая таблица с записями о выбросах
omul_clean <- omul

#стандартизация
omul_conc <- omul_clean[, num_vars]

omul_scaled <- decostand(omul_conc,
                         method = "range",
                         MARGIN = 2)

summary(omul_scaled)

#анализ корреляции
#корреляция Спирмана
data_matrix <- as.matrix(omul_scaled)
rcorr_result <- rcorr(data_matrix, type = "spearman")

DD <- rcorr_result$r
DP <- rcorr_result$P

DD_sig <- DD
DD_sig[DP > 0.05] <- 0
diag(DD_sig) <- 1

#таблицы
kable(round(DD, 2),
      caption = "Spearman correlation coefficients")

kable(round(DP, 4),
      caption = "P-values for Spearman correlations")

kable(round(DD_sig, 2),
      caption = "Significant Spearman correlations (p < 0.05)")

#хитмапа для корреляции
heatmap(DD,
        Rowv = NA,
        Colv = NA,
        scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(50),
        main = "Spearman correlation matrix")

#анализ кластеров, взял из 5 задания

#евклидово расстояние
d_eucl <- vegdist(omul_scaled, method = "euclidean")
hc_eucl <- hclust(d_eucl, method = "average")

plot(hc_eucl,
     hang = -1,
     main = "Hierarchical clustering (Euclidean, average linkage)",
     xlab = "",
     sub = "")

plot(hc_eucl)

#расстояние Брэй-Кёртиса
d_bray <- vegdist(omul_scaled, method = "bray")
hc_bray <- hclust(d_bray, method = "average")

plot(hc_bray,
     hang = -1,
     main = "Hierarchical clustering (Bray-Curtis, average linkage)",
     xlab = "",
     sub = "")

plot(hc_bray)

#дендрограмма
hc_clusters_3 <- cutree(hc_bray, k = 3)
table(hc_clusters_3)


omul_clustered <- omul_clean
omul_clustered$hc_cluster3 <- factor(hc_clusters_3)

#k-means кластеризация и обоснвание модели
#метод локтя
wss <- numeric()

for (k in 1:10) {
  set.seed(123)
  km_tmp <- kmeans(omul_scaled, centers = k, nstart = 25)
  wss[k] <- km_tmp$tot.withinss
}

plot(1:10, wss, type = "b",
     xlab = "Number of clusters k",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow method for choosing k")

#метод силуэта
fviz_nbclust(omul_scaled, kmeans, method = "silhouette")

#итоговая k-means модель
set.seed(123)
k <- 3
km_fit <- kmeans(omul_scaled, centers = k, nstart = 25)

km_fit$size
km_fit$centers

omul_clustered$kmeans_cluster <- factor(km_fit$cluster)

#Сравниваем методы кластеризации
table(HAC = omul_clustered$hc_cluster3,
      KMEANS = omul_clustered$kmeans_cluster)

# k-means
cluster_means_kmeans <- aggregate(omul_clustered[, num_vars],
                                  by = list(cluster = omul_clustered$kmeans_cluster),
                                  FUN = mean)

cluster_means_kmeans_rounded <- cluster_means_kmeans
cluster_means_kmeans_rounded[ , num_vars] <-
  round(cluster_means_kmeans_rounded[ , num_vars], 2)

cluster_means_kmeans_rounded
kable(cluster_means_kmeans_rounded,
      caption = "Mean concentrations by k-means cluster")

#агломеративная кластеризация
cluster_means_hac <- aggregate(omul_clustered[, num_vars],
                               by = list(cluster = omul_clustered$hc_cluster3),
                               FUN = mean)

cluster_means_hac_rounded <- cluster_means_hac
cluster_means_hac_rounded[ , num_vars] <-
  round(cluster_means_hac_rounded[ , num_vars], 2)

cluster_means_hac_rounded
kable(cluster_means_hac_rounded,
      caption = "Mean concentrations by HAC cluster")

table(Kmeans_Cluster = omul_clustered$kmeans_cluster,
      Sex = omul_clustered$sex)

table(HAC_Cluster = omul_clustered$hc_cluster3,
      Sex = omul_clustered$sex)

#Визуализация кластеров
# k-means кластеры в пространстве переменных
plot(omul_clustered$con_28, omul_clustered$con_180,
     col = omul_clustered$kmeans_cluster,
     pch = ifelse(omul_clustered$sex == "m", 16, 17),
     xlab = "con_28",
     ylab = "con_180",
     main = "k-means clusters in con_28 - con_180 space")

legend("topright",
       legend = levels(omul_clustered$kmeans_cluster),
       col = 1:length(levels(omul_clustered$kmeans_cluster)),
       pch = 16,
       title = "k-means cluster")

# визуализация k-means
fviz_cluster(km_fit, data = omul_scaled)

#PCA
pca_fit <- prcomp(omul_scaled, center = FALSE, scale. = FALSE)

summary(pca_fit)
pca_fit$rotation

eig_values <- get_eigenvalue(pca_fit)
eig_values
kable(round(eig_values, 2),
      caption = "Eigenvalues and explained variance")

plot(pca_fit, type = "l", main = "Scree plot for PCA")

# PCA по k-means кластерам
fviz_pca_biplot(pca_fit,
                habillage = omul_clustered$kmeans_cluster,
                addEllipses = TRUE,
                ellipse.level = 0.95,
                title = "PCA biplot by k-means cluster")

# PCA по агломеративным кластерам
fviz_pca_biplot(pca_fit,
                habillage = omul_clustered$hc_cluster3,
                addEllipses = TRUE,
                ellipse.level = 0.95,
                title = "PCA biplot by HAC cluster")

# PCA по полу омуля
fviz_pca_biplot(pca_fit,
                habillage = omul_clustered$sex,
                addEllipses = TRUE,
                ellipse.level = 0.95,
                title = "PCA biplot by sex")
#-----------------------------Проверка на робастность----------------------------
#подтверждение
dist_mat <- dist(omul_scaled)
sil <- silhouette(km_fit$cluster, dist_mat)

summary(sil)

plot(sil,
     border = NA,
     main = "Silhouette plot for k-means clusters")

#сохранение результатов
write.csv(normality_results, "normality_results_omul.csv", row.names = FALSE)
write.csv(as.data.frame(DD), "spearman_correlations_omul.csv", row.names = TRUE)
write.csv(as.data.frame(DP), "spearman_pvalues_omul.csv", row.names = TRUE)
write.csv(cluster_means_kmeans, "cluster_means_kmeans_omul.csv", row.names = FALSE)
write.csv(cluster_means_hac, "cluster_means_hac_omul.csv", row.names = FALSE)
write.csv(omul_clustered, "omul_clustered.csv", row.names = FALSE)
