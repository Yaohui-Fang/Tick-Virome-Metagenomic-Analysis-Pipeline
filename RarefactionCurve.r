library(vegan)
library(ggplot2)

# 读取并整理数据
otu <- read.delim("vOTU.txt", row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu1 <- as.matrix(t(otu))

# 获取 rarecurve 数据
rare_data <- rarecurve(otu1, step = 1, label = TRUE, tidy = TRUE)

# 作图（不指定颜色，避免映射问题）
ggplot(rare_data, aes(x = Sample, y = Species, color = Site)) +
  geom_smooth(se = FALSE, method = "lm", formula = y ~ log(x + 1), linewidth = 1.2) +
  labs(x = "Sequencing depth", y = "Viral Species") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linetype = 2),
    panel.background = element_rect(fill = NA, color = "black"),
    legend.key = element_rect(fill = "transparent")
  )
