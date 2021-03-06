
# measure the strength of the arcs in a directed graph.
arc.strength = function(x, data, criterion = NULL, ..., debug = FALSE) {

  # check x's class.
  check.bn(x)
  # arc strength is undefined in partially directed graphs.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # check the data are there.
  data.info = check.data(data)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)
  # check criterion.
  if (is.null(criterion)) {

    # if no criterion is specified use either the default one or the
    # one used by the learning algorithm.
    if (x$learning$test == "none")
      criterion = check.test(criterion, data)
    else
      criterion = x$learning$test

  }#THEN
  else  {

    criterion = check.criterion(criterion, data)

  }#ELSE

  # set the test/score counter.
  reset.test.counter()

  # expand and sanitize score-specific arguments and the alpha threshold.
  if (criterion %in% available.tests) {

    # sanitize the alpha threshold.
    alpha = check.alpha(list(...)$alpha, network = x)

    # sanitize B (the number of bootstrap/permutation samples).
    B = check.B(list(...)$B, criterion)

    # warn about unused arguments.
    check.unused.args(list(...), c("alpha", "B"))

    res = arc.strength.test(network = x, data = data, alpha = alpha,
            test = criterion, B = B, debug = debug,
            complete = data.info$complete.nodes)

    # add extra information for strength.plot().
    res = structure(res, method = "test", threshold = alpha)

  }#THEN
  else if (criterion %in% available.scores) {

    # expand and sanitize score-specific arguments.
    extra.args = check.score.args(score = criterion, network = x,
                   data = data, extra.args = list(...), learning = FALSE)

    res = arc.strength.score(network = x, data = data, score = criterion,
            extra = extra.args, debug = debug)

    # add extra information for strength.plot().
    res = structure(res, method = "score", threshold = 0)

  }#THEN

  # set the class of the return value.
  res = structure(res, class = c("bn.strength", class(res)))

  # reset the row names if there are rows.
  if (nrow(res) > 0)
    res  = structure(res, row.names = seq(nrow(res)))

  return(res)

}#ARC.STRENGTH

# compute an approximation of arc and direction strength from the Bayes factors
# that can be computed from a single MAP network.
bf.strength = function(x, data, score, ..., debug = FALSE) {

  # check whether Rmpfr is loaded.
  check.and.load.package("Rmpfr")
  # check x's class.
  check.bn(x)
  # arc strength is undefined in partially directed graphs.
  if (is.pdag(x$arcs, names(x$nodes)))
    stop("the graph is only partially directed.")
  # check the data are there.
  data.info = check.data(data)
  # check the network against the data.
  check.bn.vs.data(x, data)
  # check debug.
  check.logical(debug)

  # make sure the score function is suitable for computing a Bayes factor.
  if (missing(score)) {

    if (data.info$type %in% discrete.data.types)
      score = "bde"
    else if (data.info$type %in% continuous.data.types)
      score = "bge"
    else if (data.info$type %in% mixed.data.types)
      score = "bic-cg"

  }#THEN
  else {

    score = check.score(score, data,
              allowed = c(available.discrete.bayesian.scores,
                          available.continuous.bayesian.scores,
                          grep("bic", available.scores, value = TRUE)))

  }#ELSE

  # expand and sanitize score-specific arguments.
  extra.args = check.score.args(score = score, network = x,
                 data = data, extra.args = list(...), learning = FALSE)

  bf.strength.backend(x = x, data = data, score = score, debug = debug,
    extra.args = extra.args)

}#BF.STRENGTH

# compute the strength of all possible arcs from a list of network
# structures/arc sets.
custom.strength = function(networks, nodes, weights = NULL, cpdag = TRUE,
    debug = FALSE) {

  # check debug.
  check.logical(debug)
  # check cpdag.
  check.logical(cpdag)
  # check networks.
  if (is(networks, c("bn.kcv", "bn.kcv.list"))) {

    extract = function(x) bn.net(x$fitted)

    if (is(networks, "bn.kcv"))
      networks = lapply(networks, extract)
    else if (is(networks, "bn.kcv.list")) {

      networks = lapply(networks, function(x) lapply(x, extract))
      networks = do.call("c", networks)

    }#THEN

    # extract the node labels from the networks, and ignore the argument.
    if (!missing(nodes))
      warning("the labels from the 'nodes' argument will be ignored.")

    nodes = names(networks[[1]]$nodes)

  }#THEN
  else {

    # check the node labels.
    check.nodes(nodes)
    # check the networks.
    check.customlist(networks, nodes = nodes)

  }#ELSE
  # check the weights.
  weights = check.weights(weights, length(networks))

  res = arc.strength.custom(custom.list = networks, nodes, cpdag = cpdag,
          arcs = NULL, weights = weights, debug = debug)

  # add extra information for strength.plot().
  res = structure(res, method = "bootstrap", threshold = threshold(res),
          class = c("bn.strength", class(res)))

  return(res)

}#CUSTOM.STRENGTH

# build the averaged network structure using arc strengths and a
# significance threshold.
averaged.network = function(strength, nodes, threshold) {

  # check the strength threshold.
  threshold = check.threshold(threshold, strength)
  # check nodes.
  if (missing(nodes)) {

    # check the strength parameter.
    check.bn.strength(strength, valid = c("bootstrap", "bayes-factor"))
    # use the bn.strength object to get a node set.
    nodes = unique(c(strength[, "from"], strength[, "to"]))

  }#THEN
  else {

    # sanitize the node set.
    check.nodes(nodes = nodes)
    # check the strength object and whether it agrees with the node set.
    check.bn.strength(strength, nodes = nodes,
      valid = c("bootstrap", "bayes-factor"))

  }#ELSE

  avg = averaged.network.backend(strength = strength, nodes = nodes,
          threshold = threshold)

  # add the metadata for the print() method.
  avg$learning$algo = "averaged"
  avg$learning$args = list(threshold = threshold)

  return(avg)

}#AVERAGED.NETWORK

# average multiple bn.strength objects.
mean.bn.strength = function(x, ..., weights = NULL) {

  # check the bn.strength objects.
  strength = c(list(x), list(...))

  method = character(length(strength))
  nodes = unique(c(x[, "from"], x[, "to"]))

  for (s in seq_along(strength)) {

    # check the nodes are the same, and that the object has the right structure.
    check.bn.strength(strength[[s]], nodes = nodes)
    # check the objects have the same number of arcs.
    if (length(strength[[s]][, "from"]) != length(x[, "from"]))
      stop("the bn.strength objects have different numbers of arcs.")
    # check that the objects have the same arcs.
    if (!all(which.listed(as.matrix(strength[[s]][, c("from", "to")]),
                         as.matrix(x[, c("from", "to")]))))
      stop("the bn.strength objects have different arcs.")

    method[s] = attr(strength[[s]], "method")

  }#FOR

  # check all strengths are either score differences, p-values or
  # strength/direction probabilities.
  if (any(method != method[1]))
    stop("all the bn.strength objects must be computed from either a score, ",
      "a conditional independence test, or bootstrap resampling.")

  # check the weights.
  weights = check.weights(weights, length(strength))
  # average the objects.
  res = mean.strength(strength, nodes, weights)
  # set the attributes of the return value.
  attributes(res) = attributes(x)
  if (attr(res, "method") == "bootstrap")
    attr(res, "threshold") = threshold(res)

  return(res)

}#MEAN.BN.STRENGTH

