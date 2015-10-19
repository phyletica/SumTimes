.. role:: bolditalic
.. role:: hlight 

.. _background:

**********
Background
**********

In phylogenetics, the relationships among divergence times within or across
trees are often of interest. For example, you might want to know whether
Species X and Y divergenced before, after, or around the same time as Species A
and B. Or, you might want to know whether both diverged within a particular
window of time, like a geological period, for example.

|sumtimes|_ is a simple tool that allows you to approximate the posterior
probability of such "divergence-time scenarios" by summarizing one or more
posterior samples of time-relative or time-absolute ultrametric phylogenies.
These posterior samples of trees are collected with some other phylogenetic
software package, like |revbayes|_ or |beast|_.

One strength of |sumtimes|_ is flexibility. You can approximate the posterior
probability of (almost) any divergence-time scenario that entails an arbitrary
number of nodes across an arbitrary number of posterior tree samples.


.. _how-it-works:

How |sumtimes|_ works
=====================

The main inputs to |sumtimes|_ are:

*   One or more posterior samples of trees
*   Nodes of interest in those trees
*   Expressions of divergence-time scenarios of interest

You specify a node of interest by listing a subset of the tips found on each
tree in the posterior sample.
From each tree in the posterior sample, |sumtimes|_ will extract the age of the
node representing the most recent common ancestor (MRCA) of those tips (or the
parent node of the MRCA if you specify the node is "stem-based").


.. _within-tree-comparisons:

Comparing nodes within a tree
-----------------------------

The simplest possible use case for |sumtimes|_ is if you are interested in
comparing two divergence times across a "single" posterior sample of trees.
By "single" I mean a sample collected from a Bayesian phylogenetic analysis of
one dataset (it could actually consist of tens of thousands of trees found in
several different files on your computer).

In such a case, we would specify the two MRCA nodes (or their parent nodes)
that we want to compare. Let's call these odes A and B, and assume we want to
determine the probability that Node A diverged before (further back in the
past) than Node B.
In this case, |sumtimes|_ will check every tree in the posterior sample
(barring any burn-in samples ignored) to see if the divergence time of Node A is
greater than Node B. It will report the proportion of trees for which that is
true, which is an approximation of the posterior probability that Node A
diverged before Node B.

Note, when comparing node ages from the same posterior sample, |sumtrees|_ will
always do this *jointly*. I.e., it always compares the node-age values from the
same tree to account for the potential correlation between the ages of the
nodes.


Comparing nodes across trees
----------------------------

With |sumtimes|_ you can also compare nodes across "multiple" posterior samples
of trees.
By "multiple," I mean samples collected from multiple Bayesian phylogenetic
analyses, potentially of different datasets with completely different taxa.
In other words, you will be comparing the age of nodes in *different*
trees.
How |sumtimes|_ does this is very similar to how it :ref:`compares nodes within
a tree<within-tree-comparisons>`.
The main difference is that it can *not* do the comparisons *jointly*, because
the posterior samples from the different analyses are independent of one
another.

How does |sumtimes|_ decide which trees to compare between posterior samples?
If both posteriors have the same number of sampled trees (after burn-in is
removed), then it will simply compare them in the order they appear in the tree
files (i.e., it will use all the trees in both posterior samples).
If one posterior has a larger number of samples than the other, then
|sumtimes|_ will use the ages in all the trees from the posterior with the
smaller sample size, and randomly subsample (with replacement) the same number
of trees from the larger posterior sample.
Thus, which trees end up being compared will be random; the results will likely
differ each time you run the program, unless you specify a seed for the random
number generator.
