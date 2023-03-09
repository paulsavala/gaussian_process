library(ggplot2)
library(dplyr)

# Lattice for predictions -------------------------------------------------------
form_grid = function(df, min_radius_to_data, min_points_in_radius=5, num.x=25, num.y=25, plot.grid=FALSE) {
  x.lattice = seq(min(df$x), max(df$x), length.out=num.x)
  y.lattice = seq(min(df$y), max(df$y), length.out=num.y)
  
  grid.df = make.surface.grid(list(x=x.lattice, y=y.lattice)) %>% as.data.frame
  
  # Remove lattice points that have no training points to the left/right or above/below
  rows_to_delete = c()
  for (i in 1:nrow(grid.df)) {
    # Get lattice point coordinates
    r = grid.df[i, ]
    
    # Find distance to all training points
    D = rdist(r, df.train)
    
    # Find how many points it's within half the radius from
    r.idx = which(D < min_radius_to_data / 2)
    
    if (length(r.idx) < min_points_in_radius) {
      rows_to_delete = c(rows_to_delete, i)
    }
  }
  
  grid.df = grid.df[-rows_to_delete, ]
  
  if (plot.grid) {
    ggplot(grid.df) +
      geom_point(aes(x=x, y=y)) +
      geom_point(data=df.train, aes(x=x, y=y), color='red')
  }
  
  return(grid.df)
}


plot_preds_on_lattice = function(grid.preds, grid.df) {
  # Combine into a data frame
  grid.df.preds = cbind(grid.df, grid.preds)
  
  # Visualize predictions on lattice
  ggplot(grid.df.preds) +
    geom_tile(aes(x=x, y=y, fill=median_pred)) +
    geom_point(data=df.train, aes(x=x, y=y, fill=pm2_5), size=3, pch=21, color='black') +
    scale_fill_gradient2(low='blue', high='red', midpoint=0) +
    ggtitle('Median prediction (tile) vs actual (point) PM 2.5')
  
  ggplot(grid.df.preds) +
    geom_tile(aes(x=x, y=y, fill=se_pred)) +
    geom_point(data=df.train, aes(x=x, y=y), size=3, pch=21, color='black') +
    scale_fill_gradient(low='white', high='red') +
    ggtitle('SE prediction of PM 2.5') 
}
  
  
  
  
  
  
  