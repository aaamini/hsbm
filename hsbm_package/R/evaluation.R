#' Testing in-group versus random pairs
#'
#' @param grps  grps a list of indices within each group
#' @param A a list of (sparse) adj. matrices -- should be a list; assuming all
#'   networks have the same number of nodes
#' @export
pair_hamming_test = function(A, grps, npairs = 5000, type = 1) {
  pri = sapply(grps, function(x) choose(length(x),2))
  pri = pri/sum(pri)
  npri = round(pri*npairs)

  ingrp_pairs = do.call(rbind, lapply(seq_along(grps), function(j) {
    if (length(npri) > 0){
      matrix(sample(grps[[j]], npri[j]*2, T), ncol=2)
    }
  }))
  npairs = nrow(ingrp_pairs) # could be different than the original target due to rounding errors
  rand_pairs = matrix(sample(nrow(A[[1]]), npairs*2, T), ncol=2)

  df = rbind(
    data.frame(type = "Within", nham = get_multi_nhamming(A, ingrp_pairs-1, type = type)),
    data.frame(type = "Random", nham = get_multi_nhamming(A, rand_pairs-1, type = type))
  )
  # df$type = factor(df$type, c("Within","Random"))
  p = df %>% ggplot(aes(x = nham, color = type, fill = type)) +
    geom_density(alpha = .25) +
    ylab("Density") +
    xlab("Average pariwise N-Hamming distance") +
    theme_minimal() +
    ggplot2::theme(
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position = c(0.8, 0.9),
      legend.spacing.y = unit(3.0, 'cm')
      # legend.text = ggplot2::element_text(size=18),
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 1, keyheight = 1))

  list(plot = p, df = df)
}

#' Compactify mulilayer labels
#'
#' The function takes as input a list of lists of labels and collapses it to a
#' matrix. The outer list iterates over chain iterations, and the inner list
#' iterates over the layers. The output is a matrix of dimension
#' number_of_node_layer_pairs x number_of_iterations
#'
#' @param zh_list a list of lists of labels
#' @export
compactify_chain_multi_labels = function(zh_list) {
  # sapply(zh_list, function(zh) as.vector(do.call(rbind, zh)))
  sapply(zh_list, function(zh) do.call(c, zh))
}


#' @export
seq_nmi_plot = function(zh_list) {
  zh_compact = compactify_chain_multi_labels(zh_list)
  niter = ncol(zh_compact)
  df = data.frame(iter = 2:niter,
             nmi =  sapply(2:niter, function(t) nett::compute_mutual_info(zh_compact[,t], zh_compact[,t-1])))
    ggplot2::ggplot(df, ggplot2::aes(iter, nmi)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(trans="log10") +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab("Squential Agg. NMI")
}
