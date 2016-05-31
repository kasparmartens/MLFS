create_plot_label_switching = function(mat, switched_indexes){
  mat.m = melt(mat)
  p1 = ggplot(mat.m, aes(rev(Var1), (Var2), fill=value)) + 
    geom_tile() + 
    scale_x_continuous(breaks=NULL) + scale_y_continuous(breaks=NULL) + 
    xlab("") + ylab("") + 
    theme(panel.background = element_blank()) + 
    scale_fill_gradient(low = "white", high = "steelblue")
  
  offdiag = mat[upper.tri(mat) | lower.tri(mat)]
  switched_estimates = mat[cbind(switched_indexes, rev(switched_indexes))]
  p2 = ggplot(data.frame(x = offdiag)) + 
    geom_density(aes(x = x)) + 
    geom_vline(xintercept = switched_estimates, col="red") +
    xlim(0, 1) + theme_bw() + xlab("Off-diagonal elements")
  g = grid.arrange(p1, p2, nrow=1)
  return(g)
}
