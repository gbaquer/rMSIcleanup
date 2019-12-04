tsne_test <- function()
{
  library(caret)
  library(Rtsne)

  ######################################################################
  ## The WHOLE post is in: https://github.com/pablo14/post_cluster_tsne
  ######################################################################

  ## Download data from: https://github.com/pablo14/post_cluster_tsne/blob/master/data_1.txt (url path inside the gitrepo.)
  data_tsne=read.delim("data_1.txt", header = T, stringsAsFactors = F, sep = "\t")

  ## Rtsne function may take some minutes to complete...
  set.seed(9)
  tsne_model_1 = Rtsne(as.matrix(data_tsne), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)

  ## getting the two dimension matrix
  d_tsne_1 = as.data.frame(tsne_model_1$Y)

  ## plotting the results without clustering
  ggplot(d_tsne_1, aes(x=V1, y=V2)) +
    geom_point(size=0.25) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("t-SNE") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_brewer(palette = "Set2")

  ## keeping original data
  d_tsne_1_original=d_tsne_1

  ## Creating k-means clustering model, and assigning the result to the data used to create the tsne
  fit_cluster_kmeans=kmeans(scale(d_tsne_1), 3)
  d_tsne_1_original$cl_kmeans = factor(fit_cluster_kmeans$cluster)

  ## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
  fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

  ## setting 3 clusters as output
  d_tsne_1_original$cl_hierarchical = factor(cutree(fit_cluster_hierarchical, k=3))

  plot_cluster=function(data, var_cluster, palette)
  {
    ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
      geom_point(size=0.25) +
      guides(colour=guide_legend(override.aes=list(size=6))) +
      xlab("") + ylab("") +
      ggtitle("") +
      theme_light(base_size=20) +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            legend.direction = "horizontal",
            legend.position = "bottom",
            legend.box = "horizontal") +
      scale_colour_brewer(palette = palette)
  }


  plot_k=plot_cluster(d_tsne_1_original, "cl_kmeans", "Accent")
  plot_h=plot_cluster(d_tsne_1_original, "cl_hierarchical", "Set1")

  ## and finally: putting the plots side by side with gridExtra lib...
  library(gridExtra)
  grid.arrange(plot_k, plot_h,  ncol=2)
}
