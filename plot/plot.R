table1 <- read.table("example5.exp", header=T)
table2 <- read.table("example5_1_1.gdat", header=T)
table3 <- read.table("example5_2_1.gdat", header=T)

table1$source <- rep("exp1", nrow(table1))
table2$source <- rep("model1 (EGF=1.0nM)", nrow(table2))
table3$source <- rep("model2 (EGF=5.0nM)", nrow(table3))
table1$size <- rep(1, nrow(table1))
table2$size <- rep(1, nrow(table2))
table3$size <- rep(1, nrow(table2))


table4 <- rbind(table1, table2[,c(1,6,7,8,9)], table3[,c(1,6,7,8,9)])

require("ggplot2")

ggplot(data=table4, aes(x=time, y=pR, group = source, colour = source)) +
  #geom_line(aes(linetype=source, color=source)) +
  geom_point(aes(color=source, shape=source, fill=source), size=2) +
  geom_smooth(aes(group=source), method = "loess", se=FALSE)
  #geom_point( size=1, shape=21, fill="white")


#test later
#scale_linetype_manual(values=c("solid","twodash", "dotted"))+