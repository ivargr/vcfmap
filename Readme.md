## Tools for working with vcf haplotype information on genome graphs

### Usage
```python
vcfmap make -g graph.nobg -l linear_ref.interval -v variants.vcf.gz -o myvcfmap
```

This stores a `VcfMap` that allows for efficient lookup of which haplotypes are present at edges in the graph:

```python
from vcfmap import VcfMap
m = VcfMap.from_file("myvcfmap.npz")
haplotypes = m.get_haplotypes_on_edge(1, 2)
haplotypes = m.get_haplotypes_on_edge(1, 3)
```

