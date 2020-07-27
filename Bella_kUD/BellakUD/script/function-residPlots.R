# create a function to make the residual plots later on in the script 
# code originally obtained and edited from github.com/aosmith and https://rstudio-pubs-static.s3.amazonaws.com/43190_f9bd428cf5b94562a1b46e00b454e05a.html
resid_plots <- function(model, modelname){
  # augment is from broom package, it adds columns of residuals, predictions, etc.
  output <- broom.mixed::augment(model)
  # first plot is a residual vs fitted plot
  p1 <- ggplot(output, aes(.fitted, .resid))+geom_point()
  p1 <- p1+stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")
  p1 <- p1+xlab("Fitted values")+ylab("Residuals")
  p1 <- p1+ggtitle("Residual vs Fitted Plot")+theme_bw()+theme(plot.title = element_text(size=9),axis.title = element_text(size=8))
  # second plot is a qqplot to test normality
  p2 <- ggplot(output, aes(sample=.resid))+stat_qq()
  p2 <- p2+geom_qq_line(col='red')+xlab("Theoretical Quantiles")+ylab("Stdized Residuals")
  p2 <- p2+ggtitle("Normal Q-Q")+theme_bw()+theme(plot.title = element_text(size=9),axis.title = element_text(size=8))
  # third plot is a Cook's distance plot to assess outliers 
  p3 <- ggplot(output, aes(seq_along(.cooksd), .cooksd))+geom_bar(stat="identity", position="identity")
  p3 <- p3+xlab("Obs. Number")+ylab("Cook's distance")
  p3 <- p3+ggtitle("Cook's distance")+theme_bw()+theme(plot.title = element_text(size=9),axis.title = element_text(size=8))
  # last plot is a check of influential data points, similar to Cook's distance
  p4 <- ggplot(output, aes(.hat, .resid))+geom_point(aes(size=.cooksd), na.rm=TRUE)
  p4 <- p4+stat_smooth(method="loess", na.rm=TRUE)
  p4 <- p4+xlab("Leverage")+ylab("Stdized Residuals")
  p4 <- p4+ggtitle("Residual vs Leverage Plot")
  p4 <- p4+scale_size_continuous("Cook's Distance", range=c(1,5))
  p4 <- p4+theme_bw()+theme(plot.title = element_text(size=9),axis.title = element_text(size=8), legend.position="bottom")
  # plot with patchwork 
  (p1 | p2) /
    (p3 | p4) +
    patchwork::plot_annotation(title = paste("Diagnostic plots", modelname))
}
