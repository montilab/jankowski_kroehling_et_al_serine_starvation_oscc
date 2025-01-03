library(ggplot2)

ggstyle <- function(font = "Helvetica", scale = 1) {
  fs <- function(x) x * scale # Dynamic font scaling
  ggplot2::theme(
    plot.title = ggplot2::element_text(family = font, size = fs(26), face = "bold", color = "#222222"),
    plot.subtitle = ggplot2::element_text(family = font, size = fs(18), margin = ggplot2::margin(0, 0, 5, 0)),
    plot.caption = ggplot2::element_blank(),
    legend.position = "right",
    legend.text.align = 0,
    legend.background = ggplot2::element_blank(),
    #legend.title = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family = font, size = fs(18), color = "#222222"),
    axis.title = ggplot2::element_text(family = font, size = fs(18), color = "#222222"),
    axis.text = ggplot2::element_text(family = font, size = fs(18), color = "#222222"),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10)),
    # axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color = "#222222"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "white"),
    strip.text = ggplot2::element_text(size = fs(22), hjust = 0)
  )
}

setwd("/restricted/projectnb/montilab-p/projects/oralcancer/serine_starvation/results/metabolomics/metaboanalystonline_out")
kegg <- read.csv('./msea_ora_resultKEGG.csv')
kegg$ratio <- kegg$hits/kegg$total

kegg1 <- kegg[order(kegg$Raw.p),]
kegg2 <- kegg1[1:25,]
kegg2$X <- factor(kegg2$X, levels = rev(kegg2$X))
mid<-mean(kegg2$Raw.p)

pdf('./keggmetab.pdf', height = 8, width = 10)
ggplot(kegg2, aes(x=X, y=-log10(Raw.p), label=Raw.p)) + 
  geom_point(stat='identity', aes(col=Raw.p), size = kegg2$ratio*6)  +
  scale_size_continuous(range = c(5, 10)) +
  scale_color_gradient(low="red", high="yellow")+
  labs(title="KEGG", col = "p-value") + 
  coord_flip() +
  ylab("-log10 (p-value)") +
  ggstyle()
dev.off()

smdb <- read.csv('./SMDB.csv')
smdb$ratio <- smdb$hits/smdb$total

smdb1 <- smdb[order(smdb$Raw.p),]
smdb2 <- smdb1[1:25,]
smdb2$X <- factor(smdb2$X, levels = rev(smdb2$X))

pdf('./smdbmetab.pdf', height = 8, width = 13)
ggplot(smdb2, aes(x=X, y=-log10(Raw.p), label=Raw.p)) + 
  geom_point(stat='identity', aes(col=Raw.p), size = smdb2$ratio*6)  +
  scale_size_continuous(range = c(5, 10)) +
  scale_color_gradient(low="red", high="yellow")+
  labs(title="SMDB", col = "p-value") + 
  coord_flip() +
  ylab("-log10 (p-value)") +
  ggstyle()
dev.off()
