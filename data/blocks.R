library(dplyr)

blocks <- read.csv("blocks.csv")

head(blocks)

blocks %>% group_by(main_chain) %>% summarise(n())
## several blocks not in the main chain - what does that mean?

blocks %>% group_by(ver) %>% summarise(n())
## version numbers are strange (is 'ver' version numbers?)

blocks %>% summarise(mean(size), sd(size), mean(fees), sd(fees),
                     mean(transactions_count), sd(transactions_count))

## cleaned up version of blocks
## removes unnecessary variables
## orders by block height
## changes variable names and creates useful variables
blocks_clean <- blocks %>%
  select(-hash, -prev_block, -mrkl_root, -nonce, -bits) %>%
  arrange(height) %>%
  rename(transactions = transactions_count) %>%
  mutate(main_chain = (main_chain == "true"),
         size_mb = size / 1e6,
         time_min = (time - min(time))/60)

write.csv(blocks_clean, file = "blocks_clean.csv", row.names = FALSE)
