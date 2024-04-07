df$Cell.Type

?match

df$Color <- color_order[match(df$Cell.Type, cell_type_order)]
