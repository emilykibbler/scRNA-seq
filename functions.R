# This function takes a data set (Seurat object) and a list of desired features and returns a ggarrange plot object
# Features must be something that can be retrieved by Seurat's FetchData function
# Heavy lifting performed by Seurat's VlnPlot function; formatting and tweaking by E. Kibbler

myVlnPlot <- function(dat, feats) {
  # Structure as a list even if only one feature is given
  if (length(feats) == 1) {
    feats <- c(feats)
  }
  # Initiate empty list to put the plot(s) in
  plots <- c()
  # Plot each feature separately
  for (i in 1:length(feats)) {
    # Use suppressWarnings because the VlnPlot() has a default argument which uses a deprecated argument
    temp <- suppressWarnings(VlnPlot(dat,
                                     features = feats[i])) +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none")
    plots <- c(plots, list(temp))
  }
  # Generate and return a plot of panel(s)
  # Bottom text will display what the "project" is (defined in the Seurat object)
  sup <- paste("Data from", dat@meta.data$orig.ident[1], "project")
  res <- ggarrange(plotlist = plots, nrow = 1) %>%
    annotate_figure(bottom = text_grob(sup, 
                                       size = 12, 
                                       hjust = 0, 
                                       x = 0.1))
  return(res)
}