# Development notes — 29 August 2026 (v0.22 → v0.23)

Working record of a session spent on the annotated-tree readers and two
neighbouring functions. The user-facing summary is in the version history at the
end of `README.md`; this file keeps the reasoning, the evidence, and what was
deliberately left alone.

Commits, oldest first: `a01d7b7`, `13e82fc`, `1a9895a`, `94fa884`, `2d06ea9`,
`dbeb304`.

## How this started

Summarising BEAST 2.7 posterior tree files with TreeAnnotator into MCC trees,
then reading those trees in R. The reading step is what exposed the bugs below.
The MCC trees and the script that made them live in another project
(`~/Dropbox/projects_WORKING/infinite_sites_simona/GBE_R2`); one of them was
copied here as `example_data/h1n1_beast2_mcc.tree` (34 influenza A/H1N1
sequences from 2009, TreeAnnotator `-burnin 10 -height keep -topology MCC`).

## What was wrong, and how it was found

### 1. `read.annotated.nexus` returned the wrong tip labels

Translated labels were assigned positionally, `tree$tip.label <- TRANS[, 2]`,
which assumes the tips appear in the Newick string in TRANSLATE-block order.
BEAST writes them in topological order, so the tree came back with the right
topology and branch lengths and permuted tip labels.

Found by rebuilding the tree independently — `ape::read.tree` on the raw Newick,
then mapping the numbers through the TRANSLATE table by lookup — and comparing:

| | vs. independently built tree |
|---|---|
| `ape::read.nexus` | `all.equal.phylo` TRUE, RF = 0 |
| old `read.annotated.nexus` | FALSE, **RF = 64** |

This affected the package's own `example_data/hiv_A_env.tree`, whose tips run
`19 13 16 12 ...`. Fixed by re-deriving the tip order from the Newick string and
mapping through the TRANSLATE table.

### 2. Annotations were attached to nodes by counting, not by identity

`.annotated.tree.build()` collected the `[&...]` blocks in file order and handed
them out one per node visited, which assumes every node is annotated. A tree
with only some nodes annotated failed with `subscript out of bounds`, or, when
the counts happened to line up, **silently attached annotations to the wrong
nodes**. A tree with no annotations at all also failed.

`.strip.annotations()` already leaves a numbered `[i]` placeholder where it took
each annotation out, so each node's annotation is now read from the placeholder
in that node's own token. Identifying which tokens are nodes needs care:

- a token is a node's own token when the delimiter closing it is `,`, `)` or `;`;
- the empty tokens produced by runs of `(` are not nodes;
- the root's token is **empty** when the root has neither label nor branch
  length, so filtering tokens on `nzchar` drops it.

The first attempt at this fix aligned placeholders with the non-empty tokens and
so silently dropped the **root's** annotation on `hiv_A_env.tree`. Caught by
diffing every annotation against the previous version — worth repeating for any
future change here.

### 3. `find.sister` only worked for a single tip

Given an internal node it returned that node's own descendants rather than its
sister group. Separately, the final `else` branch referenced `all_descendants`,
a variable never defined there, so any set of two or more tips whose MRCA is not
a polytomy failed with `object 'all_descendants' not found`. On a strictly
bifurcating tree the function therefore worked for one tip and nothing else; it
appeared to work only because the example tree has a polytomy, which routes
through a different branch.

Now computed one way in every case: the tips descending from the parent of the
node, minus the tips descending from the node itself.

### 4. `get.oldest.branch.length` — not a bug

Flagged as a possible naming slip because it returns `min()`, then checked and
withdrawn. On a time-scaled tree a child of the root is reached sooner along a
shorter branch, so the shortest root-descending branch subtends the *oldest* of
the root's children. `git log` confirms it has been `min()` since 2018-06-06.
Documentation only: the reasoning is now in the roxygen block, together with the
caveat that on a phylogram the name should not be read literally.

## Other changes

- `$annotations` is a **named** list: `names(x$annotations)[i]` is
  `as.character(x$edge[i, 2])`, with the root's annotations last under
  `Ntip + 1`. Retrieval works positionally (by edge) or by node number.
- Multi-tree files return a `multiPhylo` **without** annotations and with a
  warning, rather than a mislabelled annotated object. Fast (~3 s for 1094 trees
  of 328 taxa) because the annotations are stripped and ape's parser does the
  building.
- A TREES block with no closing `End;` — as left by a run stopped part way
  through writing, which is how the files that prompted this session were — is
  read up to its last complete tree. Both the old code and `ape::read.nexus`
  fail on those with `Error in start:end : NA/NaN argument`.
- Tree names are no longer mangled for tab-indented NEXUS files.
- New arguments `simplify` (interval annotations such as `height_95%_HPD` as
  length-two vectors) and `check` (compare against `ape::read.nexus` and warn).
- `read.annotated.tree` names its annotations the same way and takes `simplify`.
- README: tutorial section 10 on BEAST 2 summary trees, and the version history
  filled in so that every name in `NAMESPACE` is mentioned somewhere.

## Verifying changes here

```r
Rscript run_tests.R      # castor 94, annotated-tree-reader 47, find-sister 19, branch-length-stats 6
R CMD INSTALL .          # then restart R: library() holds the old namespace otherwise
```

Two techniques that earned their keep and are worth reusing:

- **A/B against the installed package.** `sys.source` the working tree into an
  environment and compare its output element by element with the installed
  version on the same file. This is what caught the dropped root annotation.
- **Self-identifying test trees.** A tree whose tips appear in an order that is
  *not* the TRANSLATE order, with a unique `[&id=]` on every node, pins down both
  the label mapping and the annotation-to-node mapping at once. See
  `tagged_tree_file()` in `tests/testthat/test-annotated-tree-reader.R`.

## Environment

- `geiger` 2.0.11 installed from CRAN (pulling in `subplex` and `ncbit`).
  `run_tests.R` calls `library(geiger)` and aborted before this. One test in the
  castor suite had been skipping for the same reason and now runs.
- NELSI 0.23 installed from this working tree.

## Left alone, deliberately

- **`roxygen2` is not installed**, so every `man/*.Rd` touched here was edited by
  hand to match its roxygen block. Re-run `roxygen2::roxygenise()` when it is
  available and check the two agree.
- **The root's token in `.annotated.tree.build`.** Node labels and branch lengths
  are still read through the `nzchar`-filtered token vector, so a root with
  neither label nor branch length gives `NA`. Pre-existing, and papered over by
  `if (is.na(node.label[1])) node.label[1] <- ""` and by omitting `root.edge`.
  Only the annotation mapping was changed; the label and length handling was
  left as it was to avoid altering long-standing behaviour.
- **`.annotated.clado.build`** is still a `stop()` — trees without branch lengths
  are not supported by either reader.
- **`get.rtt.dist` and `node.to.tip.dist` are byte-for-byte identical.** The
  documentation calls the second an alias. A candidate for consolidation.
- **No chronogram check** in `get.oldest.branch.length`: an ultrametricity test
  would wrongly reject heterochronous trees, and a chronogram cannot be told
  from a phylogram by inspection.
- **`find.sister` now returns sorted tip indices**, where it previously returned
  them in traversal order. A behaviour change for anything indexing the result
  positionally.
