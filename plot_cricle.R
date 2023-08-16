############1、UP
edges <- read.csv("E:/workspace/R/case4_micro/data/up_plot_group.csv")
# 在叶子(个体)之间创建一个相互连接的数据框
connect <- read.csv("E:/workspace/R/case4_micro/data/up_plot.csv", row.names = 1)
connect <- connect[,c(1,2)]
colnames(connect) <- c("from", "to")
# 子叶之间相互连接的程度
connect$value <- runif(nrow(connect))
# 创建顶点数据框，层次结构中的每个对象一行：
vertices  <-  data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to))) , 
  value = runif(123)
) 
# 用每个名称的组添加一个列，这对后面的色点是有用的  
vertices$group  <-  edges$from[match(vertices$name, edges$to)]
# 添加关于我们将要添加的标签的信息:角度，水平调整和潜在翻转计算标签的角度  
vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, edges$from) ))
nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)
vertices$angle <- 100 - 360 * vertices$id / nleaves

# 计算标签的对齐方式:左对齐还是右对齐  
# 如果在图的左侧，标签当前的角度为< -90  
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)

# 翻转角度BY使它们可读
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)


# 创建graph对象
mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )

# 连接对象必须引用叶节点的id:  
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)



p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, width=0.9, aes(colour=..index..)) +
  scale_edge_colour_distiller(palette = "RdPu") +
  geom_node_text(aes(x = x*1.25, y=y*1.25,label=name, angle=angle), size=3.5, alpha=1) +
  geom_node_point(aes(filter = leaf, x = x*1.1, y=y*1.1, colour=group, size=1, alpha=0.2)) +
  scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
  scale_size_continuous( range = c(0.1,10) ) +
  
  theme_void() +
  theme(
    legend.position="top",
    plot.margin=unit(c(0,0,0,0),"cm"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) +
  scale_color_discrete(name="Group", labels=c("Target-gene", "GO", "KEGG", "miRNA"),) 
  


pdf("./picture/plot_up.pdf", height = 16, width = 15)
p
dev.off()







############2、DOWN
edges <- read.csv("E:/workspace/R/case4_micro/data/down_plot_group.csv")
# 在叶子(个体)之间创建一个相互连接的数据框
connect <- read.csv("E:/workspace/R/case4_micro/data/down_plot.csv", row.names = 1)
connect <- connect[,c(1,2)]
colnames(connect) <- c("from", "to")
# 子叶之间相互连接的程度
connect$value <- runif(nrow(connect))
# 创建顶点数据框，层次结构中的每个对象一行：
vertices  <-  data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to))) , 
  value = runif(77)
) 
# 用每个名称的组添加一个列，这对后面的色点是有用的  
vertices$group  <-  edges$from[match(vertices$name, edges$to)]
# 添加关于我们将要添加的标签的信息:角度，水平调整和潜在翻转计算标签的角度  
vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, edges$from) ))
nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)
vertices$angle <- 100 - 360 * vertices$id / nleaves

# 计算标签的对齐方式:左对齐还是右对齐  
# 如果在图的左侧，标签当前的角度为< -90  
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)

# 翻转角度BY使它们可读
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)


# 创建graph对象
mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )

# 连接对象必须引用叶节点的id:  
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)



p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, width=0.9, aes(colour=..index..)) +
  scale_edge_colour_distiller(palette = "RdPu") +
  geom_node_text(aes(x = x*1.25, y=y*1.25,label=name, angle=angle), size=3.5, alpha=1) +
  geom_node_point(aes(filter = leaf, x = x*1.1, y=y*1.1, colour=group, size=1, alpha=0.2)) +
  scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
  scale_size_continuous( range = c(0.1,10) ) +
  
  theme_void() +
  theme(
    legend.position="top",
    plot.margin=unit(c(0,0,0,0),"cm"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12)
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3)) +
  scale_color_discrete(name="Group", labels=c("Target-gene", "GO", "KEGG", "miRNA"),) 



pdf("./picture/plot_down.pdf", height = 16, width = 15)
p
dev.off()
