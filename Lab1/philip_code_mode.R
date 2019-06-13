dx <- density(G)
dn <- cumsum(dx$y)/sum(dx$y)
li <- which(dn>=0.05)[1]
ui <- which(dn>=0.95)[1]
dx$x[c(li,ui)]



# Method 2:
# idea 2 - h-line
density_G = density(G)
density_y = density_G$y
density_x = density_G$x
density_df = data.frame(nr = 1:length(density_x),
                        density_x = density_x,
                        density_y = density_y)
density_df_ordered = density_df[order(density_y, decreasing = TRUE),]
density_df_ordered$cumsum_y = cumsum(density_df_ordered$density_y)
density_df_ordered$cumsum_y_proportional_percent = density_df_ordered$cumsum_y/sum(density_y)*100


density_df_ordered$in_ci = density_df_ordered$cumsum_y_proportional_percent <= 95

#data frame with just true values 
density_df_ordered_trueci = density_df_ordered[(density_df_ordered$in_ci == TRUE),]
HPD = density_df_ordered_trueci[nrow(density_df_ordered_trueci),]
HPD
