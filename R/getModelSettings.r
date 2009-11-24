#
# model settings
# mean: E, D
# between: E, D 
# within: E, D
# cov: 0, D        # 0 - zero diagonal, D - same as variance

getModelStructure <- function(mean = "D", between = "D", within = "D", cov = "D") {
  list(mean = mean, between = between, within = within, cov = cov)
}
