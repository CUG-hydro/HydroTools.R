tree_depth <- function(self) {
  if (isLeaf(self)) {
    1
  } else {
    ans <- sapply(self$children, tree_depth) %>% max()
    ans + 1
  }
}

find_children <- function(df, root = -1) {
  children <- df[iddown == root, id] %>% unique()
  n <- length(children)
  # print(root)
  if (n > 0) {
    l <- list()
    for (i in 1:n) {
      l[[i]] <- find_children(df, root = children[i])
    }
    l <- set_names(l, children)
    return(l)
  } else {
    return(list())
  }
}
# df <- fread(file_basinId)
# lst <- find_children(root = -1, df)

#' plot_tree
#' 
#' @param info A data.table with the column of `["id", "iddown", "name"]`
#' @param title The title of the tree
#' 
#' @export 
plot_StreamNet <- function(info, title = "", fout = "tree.pdf", show = FALSE, root = -1) {
  # df <- fread(file_basinId)
  lst <- find_children(info, root = root)

  tree <- FromListSimple(lst)
  tree$Do(function(node) node$depth <- tree_depth(node))

  ids <- tree$Get("name")
  depths <- tree$Get("depth")
  names <- info$name[match(ids, info$id)]
  names[ids == "Root"] = title
  # depths
  tree$Set(site = names)
  tree$Set(label = sprintf("%s\r\n%s", ids, names))

  SetNodeStyle(tree, label = \(node) node$label, style = "filled,rounded")
  Do(tree$leaves, \(node) SetNodeStyle(node, 
    inherit = FALSE, shape = "box", fillcolor = "GreenYellow"))
  
  DiagrammeR::export_graph(tree %>% ToDiagrammeRGraph(dir = "descend"), file_name = fout)
  if (show) file.show(fout)
  tree
}

# 修复嵌套的流域边界
#' correct_basin_net_shapefile
#' 
#' @param ID the ID of basins
#' @param shp basin shapefile, sf polygon object
#' @param tree returned by [plot_StreamNet()]
#' @param outdir outdir of shapefile
#' @param prefix prefix of shapefile
#' 
#' @export 
correct_basin_net_shapefile <- function(ID, shp, tree, 
  outdir = "basins/subbasin", prefix = "shed_珠江", overwrite = FALSE) 
{
  Ipaper::check_dir(outdir)
  if ("grid" %in% names(shp)) shp = dplyr::rename(shp, id = grid)  
  
  x <- FindNode(tree, ID)
  name <- x$Get("site")[1]
  outfile <- sprintf("%s/%s_[%03d]_%s.shp", outdir, prefix, ID, name)

  if (!file.exists(outfile) || overwrite) {
    cat(outfile, "\n")
    # graph$Get('name', filterFun = isLeaf)
    ids <- x$Get("name") # ids = ID
    print(ids)

    shp_sub <- subset(shp, id %in% ids)
    shp_sub$id <- ID

    sf2::st_dissolve(shp_sub) %>% write_shp(outfile)
  }
}

# #' @import data.tree
# fix_basinId <- function(file_basinId, file_pour = NULL, plot = FALSE, overwrite = TRUE) {
#   df <- fread(file_basinId)
#   lst <- find_children(root = -1, df)
#   tree <- FromListSimple(lst)
#   tree$Do(function(node) node$depth <- tree_depth(node))

#   file_basinId_new <- gsub(".txt", "_fixed.txt", file_basinId)
#   if (!file.exists(file_basinId_new) || overwrite) {
#     # id <- tree$Get("name")
#     depth <- tree$Get("depth")
#     df_depth <- depth[-1] %>% {
#         data.table(id = names(.) %>% as.numeric(), depth = .)
#       }
#     df_new <- merge(df, df_depth, all.x = TRUE, sort = FALSE)

#     if (!is.null(file_pour)) {
#       pours <- read_sf(file_pour) %>% as.data.table()
#       df_new <- merge(df_new, pours, sort = FALSE)
#     }
#     fwrite(df_new, file_basinId_new)
#   } else {
#     fprintf("[ok]")
#   }
#   # df_new
# }
