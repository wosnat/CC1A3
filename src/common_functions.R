library(ComplexHeatmap)
library(tidyverse)

#####################################################################################
#####################################################################################
#####################################################################################



category.cols <- c(
  "Photosynthesis/Carbon fixation" = "#66a61e",  # Olive green
  "high light inducible" = "#e6ab02",       # Mustard yellow
  "Energy/Carbohydrate/Glycan" = "#00bfc4", # Bright Cyan
  "AA/Nucleotide;Energy/Carbohydrate/Glycan" = "#9e9ac8", # Lavender Blue
  "AA/Nucleotide" = "#7570b3",              # Purple
  "Nitrogen metabolism" = "#1b9e77",         # Teal green
  "Metabolism" = "#a6761d",                 # Brown
  "Membrane transport" = "#e7298a",         # Pink
  "Motility" = "#f781bf",                   # Light pink
  "Env. Info/Cellular Process" = "#e34a33",  # Warm Red
  "Genetic Info" = "#1f78b4",               # Deep blue
  "Other" = "#b0b0b0",                      # Gray
  "Uncharacterized" = "#d95f02"            # Orange
)

category.order <- c(
		"Photosynthesis/Carbon fixation", 
		"high light inducible",        
		"Energy/Carbohydrate/Glycan", 
		"AA/Nucleotide;Energy/Carbohydrate/Glycan", 
		"AA/Nucleotide", 
		"Nitrogen metabolism", 
		"Metabolism", 
		"Membrane transport", 
		"Motility",
		"Env. Info/Cellular Process", 
		"Genetic Info", 
		"Other",
		"Uncharacterized"
	)

#####################################################################################
#####################################################################################
#####################################################################################

contrast_map_to_label_alt_continues = c(
    "A2vsA1_prot" = 'AX\nvs 1st', 
    "A2vsA1_rna" = 'AX\nvs 1st', 
    "A3vsA1_prot" = 'AX\nvs 1st', 
    "A3vsA1_rna" = 'AX\nvs 1st', 
    "A3vsA2_prot" = 'AX\nvs prev', 
    "A3vsA2_rna" = 'AX\nvs prev', 
    "A5vsA1_prot" = 'AX\nvs 1st', 
    "A5vsA3_prot" = 'AX\nvs prev', 
    "C1vsA1_rna" = 'other', 
    "C2vsC1_prot" = 'CC\nvs 1st', 
    "C2vsC1_rna" = 'CC\nvs 1st', 
    "C3vsC1_prot" = 'CC\nvs 1st', 
    "C3vsC1_rna" = 'CC\nvs 1st', 
    "C3vsC2_prot" = 'CC\nvs prev', 
    "C3vsC2_rna" = 'CC\nvs prev', 
    "C4vsC1_prot" = 'CC\nvs 1st', 
    "C4vsC1_rna" = 'CC\nvs 1st', 
    "C4vsC3_prot" = 'CC\nvs prev', 
    "C4vsC3_rna" = 'CC\nvs prev', 
    "C5vsC1_prot" = 'CC\nvs 1st', 
    "C5vsC1_rna" = 'CC\nvs 1st', 
    "C5vsC4_prot" = 'CC\nvs prev', 
    "C5vsC4_rna" = 'CC\nvs prev', 
    "LATEvsA1_rna" = 'other', 
    "LATEvsA2_rna" = 'other', 
    "LATEvsA3_rna" = 'other', 
    "LATEvsC1_prot" = 'other', 
    "LATEvsC1_rna" = 'other', 
    "LATEvsC2_prot" = 'other', 
    "LATEvsC2_rna" = 'other', 
    "LATEvsC3_prot" = 'other', 
    "LATEvsC3_rna" = 'other'
    
    
)
 
contrast_map_to_label_pro_continues = c(
    "C1vsP1_rna" = 'other', 
    "C2vsC1_prot" = 'CC\nvs 1st', 
    "C2vsC1_rna" = 'CC\nvs 1st', 
    "C3vsC1_prot" = 'CC\nvs 1st', 
    "C3vsC1_rna" = 'CC\nvs 1st', 
    "C3vsC2_prot" = 'CC\nvs prev', 
    "C3vsC2_rna" = 'CC\nvs prev', 
    "C4vsC1_prot" = 'CC\nvs 1st', 
    "C4vsC1_rna" = 'CC\nvs 1st', 
    "C4vsC2_prot" = 'other', 
    "C4vsC3_prot" = 'CC\nvs prev', 
    "C4vsC3_rna" = 'CC\nvs prev', 
    "C5vsC1_prot" = 'CC\nvs 1st', 
    "C5vsC1_rna" = 'CC\nvs 1st', 
    "C5vsC2_prot" = 'other', 
    "C5vsC3_prot" = 'other', 
    "C5vsC3_rna" = 'Survival', 
    "C5vsC4_prot" = 'CC\nvs prev', 
    "C5vsC4_rna" = 'CC\nvs prev', 
    "LATEvsC1_prot" = 'Survival\nvs growth', 
    "LATEvsC1_rna" = 'Survival\nvs growth', 
    "LATEvsC2_prot" = 'other', 
    "LATEvsC2_rna" = 'other', 
    "LATEvsC3_prot" = 'other', 
    "LATEvsC3_rna" = 'other', 
    "LATEvsP1_prot" = 'Decompsing\nvs growth', 
    "LATEvsP1_rna" = 'other', 
    "LATEvsP2_prot" = 'Decompsing\nvs death', 
    "LATEvsP2_rna" = 'Survival\nvs death', 
    "P2vsP1_prot" = 'AX\nvs 1st', 
    "P2vsP1_rna" = 'AX\nvs 1st', 
    "P3vsP1_prot" = 'AX\nvs 1st', 
    "P5vsP1_prot" = 'AX\nvs 1st'
)

contrast_list_pro_continues = c(
    "P2vsP1_prot", 
    "P3vsP1_prot", 
    "P5vsP1_prot",
    "P2vsP1_rna", 

    "C2vsC1_prot", 
    "C3vsC1_prot", 
    "C4vsC1_prot", 
    "C5vsC1_prot",
	
    "C2vsC1_rna", 
    "C3vsC1_rna", 
    "C4vsC1_rna", 
    "C5vsC1_rna", 
	
    "C3vsC2_prot", 
    "C4vsC3_prot", 
    "C5vsC4_prot", 
	
    "C3vsC2_rna", 
    "C4vsC3_rna", 
    "C5vsC4_rna"
	
)
  
contrast_map_to_label_pro_survival_death = c(
    "C1vsP1_rna" = 'other', 
    "C2vsC1_prot" = 'Limitation', 
    "C2vsC1_rna" = 'Limitation', 
    "C3vsC1_prot" = 'Day 30', 
    "C3vsC1_rna" = 'Day 30', 
    "C3vsC2_prot" = 'other', 
    "C3vsC2_rna" = 'other', 
    "C4vsC1_prot" = 'Day 60', 
    "C4vsC1_rna" = 'Day 60', 
    "C4vsC2_prot" = 'other', 
    "C4vsC3_prot" = 'other', 
    "C4vsC3_rna" = 'other', 
    "C5vsC1_prot" = 'Day 90', 
    "C5vsC1_rna" = 'Day 90', 
    "C5vsC2_prot" = 'other', 
    "C5vsC3_prot" = 'other', 
    "C5vsC3_rna" = 'Survival', 
    "C5vsC4_prot" = 'other', 
    "C5vsC4_rna" = 'other', 
    "LATEvsC1_prot" = 'Survival\nvs growth', 
    "LATEvsC1_rna" = 'Survival\nvs growth', 
    "LATEvsC2_prot" = 'other', 
    "LATEvsC2_rna" = 'other', 
    "LATEvsC3_prot" = 'other', 
    "LATEvsC3_rna" = 'other', 
    "LATEvsP1_prot" = 'Decompsing\nvs growth', 
    "LATEvsP1_rna" = 'other', 
    "LATEvsP2_prot" = 'Decompsing\nvs death', 
    "LATEvsP2_rna" = 'Survival\nvs death', 
    "P2vsP1_prot" = 'Death\nvs growth', 
    "P2vsP1_rna" = 'Death\nvs growth', 
    "P3vsP1_prot" = 'other', 
    "P5vsP1_prot" = 'other'
)



#####################################################################################
#####################################################################################
#####################################################################################


generate_goseq_heatmap_merged <- function(goseq_res, contrast_map_to_label, contrast_list, enrichment_type, stage_order=NULL) {
    toppaths = goseq_res %>% 
        filter(contrast %in% contrast_list) %>%
        filter(
            enrich %in% enrichment_type, 
            padj < 0.05
        ) %>% 
        distinct(PATH)

    num_paths = length(toppaths$PATH)
    print(num_paths)
    heatmap_height = max(8,num_paths/2)
    goseq_filter_df = goseq_res %>% 
        filter(contrast %in% contrast_list) %>%
        filter(PATH %in% toppaths$PATH) %>% 
        replace_na(list(padj=1)) %>%
        mutate(score = -log10(padj)) %>%
        mutate(score = if_else(score > 6, 6, score)) %>%
        mutate(score = if_else(type == 'down', -score, score)) %>%
        mutate(star = if_else(padj < 0.05, '*', '')) %>%
        group_by(PATH, contrast) %>%
        arrange(over_represented_pvalue) %>%
           filter(row_number()==1) %>%
        ungroup 

    

    goseq_filter_meta = goseq_filter_df %>% 
        distinct(contrast) %>% 
        separate_wider_delim(contrast, "_", names = c(NA, "Assay"), cols_remove = FALSE) %>%
        mutate(Group = case_when(
            str_detect(contrast, "C")~ "Coculture", 
            str_detect(contrast, "P3|P5|LATEvsP")~ "Axenic late", 
            TRUE~ "Axenic"
        )) %>% 
        mutate(Stage=contrast_map_to_label[contrast]) 

    if (is.null(stage_order)) {
        stage_order = sort(unique(goseq_filter_meta$Stage))
    }
    goseq_filter_meta = goseq_filter_meta %>%
        mutate(
            Group = factor(Group, levels=c("Coculture", "Axenic", "Axenic late")),
            Assay = factor(Assay, levels=c("rna", "prot")),
            Stage = factor(Stage, levels=stage_order),
            contrast = factor(contrast, levels=contrast_list),
              ) %>%
        arrange(contrast) %>%
        column_to_rownames('contrast') 

    print(dput(unique(goseq_filter_meta$Stage)))
    # row_ha = rowAnnotation(foo2 = runif(10), ))
    
    module_mat = goseq_filter_df %>% 
        distinct(PATH, Category) %>%
        mutate(Category = factor(Category, levels=category.order)) %>%
        arrange(Category,PATH) %>%
        column_to_rownames('PATH') 

    goseq_filter_df = goseq_filter_df %>% 
        mutate(PATH = factor(PATH, levels=rownames(module_mat))) %>%
        mutate(Category = factor(Category, levels=category.order)) %>%
        arrange(Category, PATH)
    
    
    # Heatmap(small_mat, name = "mat", col = col_fun, 
    #     layer_fun = function(j, i, x, y, width, height, fill) {
    #         v = pindex(small_mat, i, j)
    #         l = v > 0
    #         grid.text(sprintf("%.1f", v[l]), x[l], y[l], gp = gpar(fontsize = 10))
    # })
    
    stars = goseq_filter_df %>% #filter(type=='up') %>%
        pivot_wider(names_from = contrast, values_from = star, id_cols =PATH) %>% 
        arrange(PATH) %>%
        relocate(rownames(goseq_filter_meta)) %>%
        column_to_rownames('PATH') 
    
    up_mat = goseq_filter_df %>% #filter(type=='up') %>%
        pivot_wider(names_from = contrast, values_from = score, id_cols =PATH) %>% 
        arrange(PATH) %>%
        relocate(rownames(goseq_filter_meta)) %>%
        column_to_rownames('PATH') %>% as.matrix() %>%
    Heatmap(
        name = 'GOSEQ score',
        col=goseq_pal_updown_fun,
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_columns  = FALSE,
        cluster_rows  = FALSE,
        show_row_dend = FALSE,
        #row_order=module_mat$PATH,
        column_order = rownames(goseq_filter_meta),
        row_split = module_mat$Category,
        column_split = goseq_filter_meta$Stage,
        #column_title_gp = gpar(col = group_pal),
        #column_names_gp = gpar(col = group_pal),
        column_title_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),
        row_names_gp = gpar(fontsize=8),
        border=TRUE,
        width=unit(8, 'cm'),
        height=unit(heatmap_height, 'cm'),
        heatmap_legend_param = list(direction = "horizontal"),
        row_title=NULL,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%s", stars[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        
        #column_title = 'HOT1A3: GOSEQ Score - upregulated',
    ) #, annotation_col = goseq_filter_meta) #, annotation_row = path_meta)
    
    
    ann_mat = Heatmap(
        as.matrix(module_mat), 
        #row_order=module_mat$PATH,
        name='Category',
        width=unit(0.5, 'cm'),
        heatmap_legend_param = list(direction = "horizontal", ncol=3),
        #col = structure(category.cols, names = category.order),
        col = category.cols,
        row_title=FALSE,
    )
    ht_list = up_mat + ann_mat
    return (draw(ht_list, heatmap_legend_side = "bottom", merge_legend = TRUE))
}




#####################################################################################
#####################################################################################
#####################################################################################


generate_goseq_heatmap <- function(contrast_list, enrichment_type, testcolname) {
    toppaths = goseq_res %>% 
        filter(contrast %in% contrast_list) %>%
        filter(
            enrich %in% enrichment_type, 
            padj < 0.05
        ) %>% 
        distinct(PATH)

    goseq_filter_df = goseq_res %>% 
        filter(contrast %in% contrast_list) %>%
        filter(PATH %in% toppaths$PATH) %>% 
        replace_na(list(padj=1)) %>%
        mutate(score = -log10(padj)) %>%
        mutate(score = if_else(score > 6, 6, score)) %>%
        mutate(star = if_else(padj < 0.05, '*', '')) 
    

    goseq_filter_meta = goseq_filter_df %>% 
        distinct(contrast) %>% 
        separate_wider_delim(contrast, "_", names = c(NA, "Assay"), cols_remove = FALSE) %>%
        mutate(Group = case_when(
            str_detect(contrast, "C")~ "Coculture", 
            str_detect(contrast, "P3|P5|LATEvsP")~ "Axenic late", 
            TRUE~ "Axenic"
        )) %>% 
        mutate(Stage=contrast_map_to_label[contrast]) %>%
        mutate(
            Group = factor(Group, levels=c("Coculture", "Axenic", "Axenic late")),
            Assay = factor(Assay, levels=c("rna", "prot"))
            
              ) %>%
        arrange(Stage, Group, Assay, contrast) %>%
        column_to_rownames('contrast') 
    
    # row_ha = rowAnnotation(foo2 = runif(10), ))
    
    module_mat = goseq_filter_df %>% 
        distinct(PATH, Category) %>%
        mutate(Category = factor(Category, levels=category.order)) %>%
        arrange(PATH) %>%
        column_to_rownames('PATH') 
    
    # Heatmap(small_mat, name = "mat", col = col_fun, 
    #     layer_fun = function(j, i, x, y, width, height, fill) {
    #         v = pindex(small_mat, i, j)
    #         l = v > 0
    #         grid.text(sprintf("%.1f", v[l]), x[l], y[l], gp = gpar(fontsize = 10))
    # })
    
    up_stars = goseq_filter_df %>% filter(type=='up') %>%
        pivot_wider(names_from = contrast, values_from = star, id_cols =PATH) %>% 
        arrange(PATH) %>%
        relocate(rownames(goseq_filter_meta)) %>%
        column_to_rownames('PATH') 
    down_stars = goseq_filter_df %>% filter(type=='down') %>%
        pivot_wider(names_from = contrast, values_from = star, id_cols =PATH) %>% 
        arrange(PATH) %>%
        relocate(rownames(goseq_filter_meta)) %>%
        column_to_rownames('PATH') 
    
    up_mat = goseq_filter_df %>% filter(type=='up') %>%
        pivot_wider(names_from = contrast, values_from = score, id_cols =PATH) %>% 
        arrange(PATH) %>%
        relocate(rownames(goseq_filter_meta)) %>%
        column_to_rownames('PATH') %>% as.matrix() %>%
    Heatmap(
        name = 'Up score',
        col=goseq_pal_up_fun,
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_columns  = FALSE,
        show_row_dend = FALSE,
        #column_order = rownames(goseq_filter_meta),
        row_split = module_mat$Category,
        column_split = goseq_filter_meta$Stage,
        #column_title_gp = gpar(col = group_pal),
        #column_names_gp = gpar(col = group_pal),
        column_title_gp = gpar(fontsize=6),
        column_names_gp = gpar(fontsize=6),
        row_names_gp = gpar(fontsize=6),
        border=TRUE,
        width=unit(6, 'cm'),
        height=unit(12, 'cm'),
        heatmap_legend_param = list(direction = "horizontal"),
        row_title=NULL,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%s", up_stars[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        
        #column_title = 'HOT1A3: GOSEQ Score - upregulated',
    ) #, annotation_col = goseq_filter_meta) #, annotation_row = path_meta)
    
    down_mat = goseq_filter_df %>% filter(type=='down') %>% 
        pivot_wider(names_from = contrast, values_from = score, id_cols =PATH) %>% 
        relocate(rownames(goseq_filter_meta)) %>%
        arrange(PATH) %>%
        column_to_rownames('PATH') %>% as.matrix() %>%
    Heatmap(
        name = 'Down Score',
        col=goseq_pal_down_fun,
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_columns  = FALSE,
        show_row_dend = FALSE,
        #column_order = rownames(goseq_filter_meta),
        column_split = goseq_filter_meta$Stage,
        # column_title_gp = gpar(col = group_pal),
        # column_names_gp = gpar(col = group_pal),
        column_title_gp = gpar(fontsize=6),
        column_names_gp = gpar(fontsize=6),
        row_names_gp = gpar(fontsize=6),
        border=TRUE,
        width=unit(6, 'cm'),
        height=unit(15, 'cm'),
        heatmap_legend_param = list(direction = "horizontal"),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%s", down_stars[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        #column_title = 'MED4: GOSEQ Score - upregulated',
    ) #, annotation_col = goseq_filter_meta) #, annotation_row = path_meta)
    
    ann_mat = Heatmap(
        as.matrix(module_mat), 
        name='Category',
        width=unit(0.5, 'cm'),
        heatmap_legend_param = list(direction = "horizontal", ncol=3),
        col = structure(brewer.pal(length(unique(module_mat$Category)), "Paired"), names = unique(module_mat$Category)),
        row_title=FALSE,
    )
    
    ht_list = up_mat + down_mat + ann_mat
    return (draw(ht_list, heatmap_legend_side = "bottom", merge_legend = TRUE))
}