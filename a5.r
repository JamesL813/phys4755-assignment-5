library(ggplot2)

data1 <- read.csv("5-1.csv", header = TRUE, sep = ",")

ggplot(data1) +
    geom_line(aes(v0, i1), color = "lightblue") +
    geom_line(aes(v0, i2), color = "pink") +
    geom_line(aes(v0, i3), color = "lightgreen") +
    theme_bw()
ggsave(file = "5-1.pdf")    
