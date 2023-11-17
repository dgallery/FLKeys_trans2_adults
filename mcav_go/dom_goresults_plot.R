##### PLOTS #######
library(ape)
library(ggplot2)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(RColorBrewer)
setwd("mcav_go/")

# What to change:
# Classifications in classify_labels()
# text_pad in go_plot() (if labels are cut off)
# Note that three text labels are cut off because their pvals are not <= 0.1
# 

results_rdata = "mcav_go_results.rdata"
out_dir = "outputs"

png_width = 6 # Width of the output files

## Functions ####
classify_labels = function(graph){
  
  labels = c('mitochondria and ribosomes', 'cyto? maybe?', 'metabolism & shit', 'Something Else')
  
  graph |> mutate(class = case_when(
    # Feel free to update these labels to something correct
    str_detect(label, 'translat|ribosom|rRNA|ribonu|mito') ~ labels[1], # Feel free to update these to other co
    str_detect(label, 'actin|bind|chromo|condens|cytoplas') ~ labels[2], 
    str_detect(label, 'metab|proteo|peptidase|secretory|catab') ~ labels[3], 
    TRUE ~ labels[4]
  ) |> factor(levels = labels))
}
tidy_clusters = function(data_list) {
  # data_list should be a list, where each element is list(data.frame, hclust)
  # First, pull out the metadata & reformat
  node_metadata = data_list |> map(1) |> 
    map(as_tibble, rownames = 'label')
  # Then tidy the hclusts
  tidy_hclust = data_list |> map(2) |>
    map(as_tbl_graph) |> # convert to tidygraph format
    # Join the metadata to each element
    map2(node_metadata, left_join, by = 'label')   |> 
    # Add the module name as a column
    imap(\(x, name) mutate(x, module = str_remove(name, 'results'))) |> 
    # Classify the labels
    map(classify_labels) 
  tidy_hclust
}
join_clusters = function(tidy_clusts) {
  do.call(bind_graphs, args = unname(tidy_clusts))
}

bin_pval = \(x, thresholds =  c(0.1, 0.05, 0.01)) {
  sort_thresh = sort(thresholds)
  thresh = c(sort_thresh, -1)
  cut(x, thresh, right = FALSE, include.lowest = TRUE, labels = c(paste('p ≤', sort_thresh))) |> 
    fct_rev()
}

go_plot = function(tidy_graph, text_pad = 20, 
                   pval_thresh = c(0.2,0.1, 0.05, 0.01),
                   pval_sizes = c(1, 1.25, 1.5, 1.75)*2,
                   pval_faces = c('plain','italic', 'plain', 'bold')) {
  # layout = create_layout(tidy_graph, 'dendrogram', direction = 'out')
  data = tidy_graph |> 
    mutate(pval_class = bin_pval(pval, pval_thresh))
  # browser()
  # key_pval = \(data = list(), ...) {
  #   data$label = c(levels(data$pval_class), '')
  #   data$fontface =  c(pval_faces, 'plain')
  #   data$sizes = c(pval_sizes, 0)
  #   draw_key_text(data, ...)
  # }
  data |> 
    ggraph(layout = 'dendrogram') + 
    geom_node_text(aes(color = pval_class, label = label, 
                       fontface = pval_class, size = pval_class), hjust = 0,
                   key_glyph = draw_key_rect) + 
    # geom_node_point(aes(color = class)) + 
    geom_edge_elbow() + 
    # Align graph correctly & make space for labels
    coord_flip(ylim = c(NA, -text_pad), clip = FALSE) + scale_y_reverse() + 
    facet_graph(module~., row_type = 'node', scales = 'free', space = 'free', switch = 'both') +
    theme(panel.background = element_blank(), strip.background = element_blank(), 
          strip.text = element_text(size = 12, face = 'bold'),
          legend.position = 'bottom', 
          legend.text = element_text(size = 8)
    )  + 
    #scale_color_brewer(type = 'seq', palette = 6 , drop = FALSE) + 
    scale_color_manual(values = colorRampPalette(brewer.pal(9, "Greys"))(12)[9:12]) +
    scale_size_manual(values = pval_sizes, na.value = 0, guide = 'none') + 
    scale_discrete_manual('fontface', values = pval_faces,
                          guide = 'none',
                          na.value = 'plain') 
  # ylim(c(10, -10))
  #+ facet_wrap(~module, strip.position = 'left')
}


### Input and clean data ####

input_names = load(results_rdata) 
# tibble::lst() auto-names the list elements
GO_terms_list <- lst(lightcyanMFresults,lightcyanBPresults,lightcyanCCresults, lightgreenMFresults,lightgreenBPresults,lightgreenCCresults, darkgreyBPresults, 
                     darkgreyCCresults, greenyellowMFresults, greenyellowBPresults, greenyellowCCresults)
tidy_GO_terms_list = tidy_clusters(GO_terms_list)
tidy_GO_terms = join_clusters(tidy_GO_terms_list) # don't end up using this, but I left it here incase it's of use

# make the plots
go_plot_lst = tidy_GO_terms_list |>
  map(go_plot,text_pad = 30)

# Create your output directory
dir.create(out_dir, showWarnings = FALSE)
# Save your plots
iwalk(go_plot_lst, \(plot, name) {
  n_terms = plot$data |> filter(leaf) |> nrow()
  out_file = file.path(out_dir, paste0(name, '.png'))
  ggsave(out_file, plot, dpi = 300, width = png_width, height = n_terms/5+1)
})

