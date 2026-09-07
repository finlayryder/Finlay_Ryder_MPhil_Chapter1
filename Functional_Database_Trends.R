library(here)
library(ggplot2)
library(tidyr)


data <- read.csv("Functional_Database_Trends.csv")

data_long <- pivot_longer(
  data,
  cols = c(Publications, Datasets),
  names_to = "Type",
  values_to = "Count"
)

ggplot(data_long, aes(Year, Count, colour = Type)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = seq(2000, 2024, by = 2)
  ) +
  scale_colour_manual(
    values = c("Publications" = "#2C7BB6",
               "Datasets" = "#D7191C")
  ) +
  labs(
    x = "Year",
    y = "Number",
    colour = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

