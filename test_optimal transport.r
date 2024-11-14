
#####latest exec code####
library(transport)
library(cluster)
library(pheatmap)
library(SummarizedExperiment)
library(approxOT)
library(diffusionMap)
library(dplyr);library(tidyr);library(tidyverse)
library(NMF)

set.seed(8)
n_genes <- 100
n_blood_patients <- 400
n_liver_patients <- 550
set.seed(8)
# Generate synthetic data with integer counts
blood_data <- matrix(rpois(n_genes * n_blood_patients, lambda = exp(rnorm(n_genes * n_blood_patients) + 5)), 
                     nrow = n_genes, ncol = n_blood_patients)
colnames(blood_data) <- paste0("blood_Patient_", 1:n_blood_patients)
rownames(blood_data) <- paste0("Gene_", 1:n_genes)

liver_data <- matrix(rpois(n_genes * n_liver_patients, lambda = exp(rnorm(n_genes * n_liver_patients) + 5)), 
                     nrow = n_genes, ncol = n_liver_patients)
colnames(liver_data) <- paste0("liver_Patient_", 1:n_liver_patients)
rownames(liver_data) <- paste0("Gene_", 1:n_genes)

clusters_blood <- sample(c("healthy", "disease"), n_blood_patients, replace = TRUE)
clusters_liver <- sample(c("healthy", "disease"), n_liver_patients, replace = TRUE)

#checking
head(blood_data[, 1:5])
head(liver_data[, 1:5])

####create synthetic SE
se_blood <- SummarizedExperiment(
  assays = list(counts = blood_data),
  colData = data.frame(
    sampleID = colnames(blood_data),
    clusters = clusters_blood
  )
)

se_liver <- SummarizedExperiment(
  assays = list(counts = liver_data),
  colData = data.frame(
    sampleID = colnames(liver_data),
    clusters = clusters_liver
  )
)

list_SE <- list('se_blood'=se_blood,'se_liver'=se_liver)
mat_list<-list(as.data.frame(blood_data),as.data.frame(liver_data))
nmf_res=run_nmf(mat_list, k = 2:4,nrun=2,seed=8)
nmf_res$nmf_object
#check optimal k for NMF
plot(nmf_res)
#higher  Cophenetic Correlation means better
#lower dispersion of intra cluster
#higher evar explained
#lower residuals and RSS
#higher sparseness more interpretability

nmf_selected<-selecting_nmf_model(nmf_res,k=2)
names(nmf_selected)
nmf_selected$nmf.clusters



H_matrix<-nmf_res$H
kmeans_dat=kmeans(H_matrix,4,nstart=25)
nmf_res$Hkmeans(H_matrix)
require(NMF)



??NMF()




elbow_kmeans(H_matrix,20)
H_matrix
testcost<-cost_matrix_centroid_new(list_SE, kmeans_dat,H_matrix, metric = "cosine")
testcost

p<-plot(1:max_clusters, wcss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters",
     ylab = "Total Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Plot for Optimal Number of Clusters")
p

dist_profiles <- Cluster_Representations_bulk(list_SE, kmeans_dat)

dist_profiles
#The rows of the cost matrix in the transport_plan_given_C function should correspond to the elements in mass_x,
#and the columns should correspond to the elements in mass_y

