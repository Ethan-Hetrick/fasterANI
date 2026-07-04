# Benchmarks

# 100 v 100 prokaryote

2 sets of 100 random prok genomes from GTDB r207, all performed with 20 threads in query-list vs reference-list mode.

## Minmer

- Tested minmer count in testing-guided intervals between 1-1000
- No minmer count testedresulted in more emitted pairs than baseline
- Candidates began to increase at minmer count 260
- Emitted pairs began to decrease at minmer count minmer count 225
- No minmer count significantly improved execution time without severe drop in recall
- Pairs recovered with low minmer counts were exclusively divergent (<80% ANI)

## Confidence

- The default mash confidence is 0.9, which emitted 6003 pairs ranging from 73.228-88.591 ANI, and 0-0.619
- Increasing the default increased emitted pairs
- Confidence of 9x10^-10 recovered 9999/1000 possible pairs
- Confidence of <=0.2 recovered 1332/1000 possible pairs
- Confidence of 9x10^-10 recovered pairs ranging from 69.865-81.978 ANI, and AF 0.001-0.837 (might be value in this extreme)
- Confidence of 0.00001 recovered pairs ranging from 80.008-92.609 ANI, and up to 0.494 AF
- Decreasing confidence significantly lowers candidantes and improves execution time, at the cost of fewer emitted pairs
- Lower confidence values produce pairs in a lower range
- Decreasing confidence 


