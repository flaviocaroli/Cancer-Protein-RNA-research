# Gene identifier mapping

For a network-independent Collection A build, place a tab-separated file here:

`ensembl_to_symbol.tsv`

Accepted Ensembl column names: `ensembl_gene_id`, `ensembl`, `query`, `gene_id`.
Accepted symbol column names: `gene_symbol`, `symbol`, `external_gene_name`, `hgnc_symbol`.

Example:

```text
ensembl_gene_id	gene_symbol
ENSG00000141510	TP53
```

If no local mapping is present, notebook 01 can query MyGene once and cache the result as
`ensembl_to_symbol_cache.tsv`. A local versioned mapping is preferable for reproducibility.
