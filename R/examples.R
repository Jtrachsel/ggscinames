require(grid)
require(ggplot2)

# Examples can be run after loading allll of the functions
#
# df <- data.frame(V1 = c(1, 1.25, 1.75, 2), V2 = c(2,3,3,4),
#                  V3 = c("in vivo/.* ", "in vivo", "another", "another"),
#                  V4 = c("point", "strain", "point 7887", "strain CAG:34 ?/"),
#                  V5 = c(TRUE, TRUE, FALSE, FALSE))
#
#
# ggplot(data = df, aes(x = V1, y = V2)) +
#   geom_point() +
#   geom_text_sciname(aes(sci = V3, nonsci = V4, important = V5), hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
#
# ggplot(data = df, aes(x = V1, y = V2)) +
#   geom_point() +
#   geom_text_sciname(aes(sci = V3, important = V5), hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
#
# ggplot(data = df, aes(x = V1, y = V2)) +
#   geom_point() +
#   geom_text_sciname(aes(nonsci = V4, important = V5), hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
#
# ggplot(data = df, aes(x = V1, y = V2)) +
#   geom_point() +
#   geom_text_sciname(aes(sci = V3, nonsci = V4, important = V5), hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
#
# ggplot(data = df, aes(x = V1, y = V2)) +
#   geom_point() +
#   geom_text_sciname(hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
#
#
# ########## four different scenarios#######
# #
# # # non-sci and important present
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_label_sciname(aes(sci = V3, nonsci = V4, important = V5), hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
# #
# #
# # # neither present
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_label_sciname(aes(sci = V3), hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
# #
# # # sci and important only
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_label_sciname(aes(sci = V3, important = V5), hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
# #
# # # sci and nonsci, no important
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_label_sciname(aes(sci = V3, nonsci = V4), hjust =1, vjust = 2) + xlim(.75,2) + ylim(1.5,4)
# #
# #
# #
# #
# #
# #
# #
# # # make sure you have the github version!!!
# # # devtools::install_github("slowkow/ggrepel")
# # library(ggrepel)
# #
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_label_sciname_repel(aes(sci = V3, nonsci = V4, important = V5))
# #
# #
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_label_sciname_repel(aes(sci = V3))
# #
# #
# #
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_label_sciname_repel(aes(sci = V3, nonsci = V4))
# #
# #
# #
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_label_sciname_repel(aes(sci = V3, important = V5))
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_text_sciname_repel(aes(sci = V3, nonsci = V4, important = V5))
# #
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_text_sciname_repel(aes(sci = V3))
# #
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_text_sciname_repel(aes(sci = V3, nonsci = V4))
# #
# # ggplot(data = df, aes(x = V1, y = V2)) +
# #   geom_point() +
# #   geom_text_sciname_repel(aes(sci = V3, important = V5))
# #
# #
# #
# # # ggplot(data = df, aes(x = V1, y = V2)) +
# # #   geom_point() +
# # #   geom_label_repel(aes(label = V3))
# #
# #
# # # aes(genus=gen, species=spec, strain=strn, bold=important_ones)
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
