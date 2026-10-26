library(ggplot2)
library(dplyr)
library(ggsidekick); theme_set(theme_sleek())

home <- here::here()

B0 <- readRDS(paste0(home, "/output/b0.rds"))

ggplot(B0, aes(year, B0)) +
  geom_line(color = "grey30") +
  facet_wrap(~sid, scales = "free_y", ncol = 6) +
  labs(x = "Year", y = expression(B[0])) +
  scale_x_continuous(breaks = c(1980, 2020)) +
  theme(
    strip.text.x.top = element_text(size = 7),
    axis.text.y = element_blank(),
    #axis.title.y = element_blank(),
    axis.ticks.y = element_blank()
    )

ggsave(paste0(home, "/figures/supp/b0_by_stock.pdf"), width = 21, height = 30, units = "cm")
