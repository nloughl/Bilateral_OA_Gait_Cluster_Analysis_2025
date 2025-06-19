# Single Component Reconstruction
library(matrixStats)
library(pracma) # needed for repmat
library(ggplot2)

scr_plot <- function(input, PCnum, PlotTitle) {
  loading <- as.matrix(input$rotation)
  # results in rows as Loading Vectors for each PC
  loading <-t(loading) 
  upperBand <- colQuantiles(input$x, probs = c(.95))
  lowerBand <- colQuantiles(input$x, probs = c(.05))
  # results in columns as 95th Percentile vectors for each PC.
  XU <- repmat(input$center, dim(input$x)[2], 1) + (loading * upperBand)
  # results in columns as 5th Percentile vectors for each PC.
  XL <- repmat(input$center, dim(input$x)[2], 1) + (loading * lowerBand)
  t_rownames <- rownames(XU)
  XU <- t(XU)
  XL <- t(XL)
  colnames(XU) <- t_rownames
  colnames(XL) <- t_rownames
  XU <- as.data.frame(XU) %>% rowid_to_column(var = "x")
  XL <- as.data.frame(XL) %>% rowid_to_column(var = "x")
  
  
  meanx <- as.data.frame(input$center)
  colnames(meanx) <- 'xBar'
  meanx <- meanx %>% rowid_to_column(var = "x")
  
  p <- ggplot() +
    #geom_line(data = as.data.frame(XU), aes(x = x, y = PC1), color = 'green') +
    #geom_line(data = as.data.frame(XL), aes(x = x, y = PC1), color = 'red') +
    geom_line(data = as.data.frame(XU), aes_string(x = "x", y = PCnum), linewidth = 1, color = 'green') +
    geom_line(data = as.data.frame(XL), aes_string(x = "x", y = PCnum), linewidth = 1, color = 'red') +
    geom_line(data = as.data.frame(meanx), aes(x = x, y = xBar), linewidth = 1, color = 'black') +
    geom_hline(yintercept = 0,
               color = "black",
               linewidth = 1) +
    theme_minimal() +
    theme(
      axis.line = element_line(linewidth = 1, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = 'none'
    ) +
    scale_x_continuous(expand = c(0, 0))+
    ggtitle(PlotTitle)
  
  return(p)
  
}