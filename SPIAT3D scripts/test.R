number_of_cells_within_radius <- function(data, reference_celltype, 
                                          target_celltype, radius = 20, 
                                          feature_colname) 
{
  
  all.df <- list()
  for (i in reference_celltype) {
    reference_cells <- data[which(data[,feature_colname] == i), ]
    reference_cell_cords <- reference_cells[, c("Cell.ID", "Cell.X.Position","Cell.Y.Position")]
    dataframe <- tibble::remove_rownames(reference_cell_cords)
    dataframe <- dataframe %>% tibble::column_to_rownames("Cell.ID")
    reference_cell_cords <- reference_cells[, c( "Cell.X.Position","Cell.Y.Position")]
    for (j in target_celltype) {
      target_cells <- data[which(data[,feature_colname] == j),     ]
      target_cell_cords <- target_cells[, c("Cell.ID", "Cell.X.Position","Cell.Y.Position")]
      target_cell_cords <- tibble::remove_rownames(target_cell_cords)
      target_cell_cords <- target_cell_cords %>% tibble::column_to_rownames("Cell.ID")
      reference_target_result <- dbscan::frNN(target_cell_cords, eps = radius,
                                              query = reference_cell_cords, sort =FALSE)
      n_targets <- rapply(reference_target_result$id, length)
      dataframe[,j] <- n_targets
    }
    all.df[[i]] <- dataframe
  }
  return(all.df)
}
