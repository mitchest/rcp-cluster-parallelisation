# purity & entropy functions, from:
# Manning, C. D., Raghavan, P., & Schütze, H. (2008). Introduction to information retrieval (Vol. 1, No. 1, p. 496). Cambridge: Cambridge university press.


cluster_purity <- function(conf_mat) {
  print(round(apply(conf_mat, 2, max) / apply(conf_mat, 2, sum), 2))
  print(round(apply(conf_mat, 1, max) / apply(conf_mat, 1, sum), 2))
  invisible()
}

max_plus_2ndmax <- function(x) {
  sum(sort(x, decreasing = T)[1:2])
}

cluster_purity_2nd <- function(conf_mat) {
  print(round(apply(conf_mat, 2, max_plus_2ndmax) / apply(conf_mat, 2, sum), 2))
  print(round(apply(conf_mat, 1, max_plus_2ndmax) / apply(conf_mat, 1, sum), 2))
  invisible()
}

log2_zero <- function(x) {
  zero_idx <- x == 0
  logs <- log2(x)
  logs[zero_idx] <- 0
  logs
}

mle_prob <- function(x) {
  prob <- x / sum(x)
  sum(abs(prob * log2_zero(prob)))
}

cluster_entropy <- function(conf_mat) {
  print(round(apply(conf_mat, 2, mle_prob),2))
  print(round(apply(conf_mat, 1, mle_prob),2))
}


# load confusion matrices
nsw_formations <- read.csv("nsw_confusion.csv", header = F)
cwl_formations <- read.csv("cwl_confusion_form.csv", header = F)
cwl_dsf <- read.csv("cwl_confusion_dsf.csv", header = F)

cluster_purity(nsw_formations)
cluster_purity(cwl_formations)
cluster_purity(cwl_dsf)

cluster_purity_2nd(nsw_formations)
cluster_purity_2nd(cwl_formations)
cluster_purity_2nd(cwl_dsf)

cluster_entropy(nsw_formations)
cluster_entropy(cwl_formations)
cluster_entropy(cwl_dsf)