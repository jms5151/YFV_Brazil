# Initialize a list to store results for each yfv_params_idx
num_params <- length(yfv_params_list)
grouped_results <- vector("list", num_params)

i=9
test <- resultsNew[[1]]

# Extract the third element from each list within grouped_results
third_elements <- lapply(test, function(x) x[[3]])

# Combine all third elements into a single data frame using do.call and rbind
combined_df <- do.call(rbind, third_elements)

library(tidyverse)

test2 <- combined_df %>%
  na.omit() %>%
  group_by(time) %>%
  summarise(
    I_h_median = quantile(I_h, 0.5)
    , I_h_25 = quantile(I_h, 0.25)
    , I_h_75 = quantile(I_h, 0.75)
  )

plot.ts(test2$I_h_75, lty = 2, col = 'red')
lines(test2$I_h_median)
lines(test2$I_h_25, lty = 2, col = 'red')
