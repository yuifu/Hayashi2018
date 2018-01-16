# Custom codes for Hayashi et al. (2018)

This repository includes custom R and Julia codes for generating and processing of the data described in:

- Tetsutaro Hayashi†, Haruka Ozaki†, Yohei Sasagawa, Mana Umeda, Hiroki Danno and Itoshi Nikaido, Single-cell full-length total RNA sequencing uncovers dynamics of recursive splicing and enhancer RNAs, accepted.

#### Citation

Tetsutaro Hayashi, Haruka Ozaki, Yohei Sasagawa, Mana Umeda, Hiroki Danno and Itoshi Nikaido, Single-cell full-length total RNA sequencing uncovers dynamics of recursive splicing and enhancer RNAs, accepted.

#### License
Copyright (c) 2017 Haruka Ozaki Released under the MIT license

----

## Brief description of custom codes
### Analysis of recursive splicing (Fig. 4, Supplementary Fig. 16 and 17)

The directory `recursive_splicing` includes custom codes for preparing of coverage data, fitting sawtooth patterns, replotting fitted data, and calculating RS detection probability.


### Analysis of enahncer RNA (Fig. 5, Supplementary Fig. 18 and 19)

The directory `enhancer_rna` includes custom codes for generating aggregation plots and heatmaps from output files of [deepTools](https://deeptools.readthedocs.io/en/latest/).



### Millefy (Fig. 3e and 4bc, Supplementary Fig. 10)

`millefy`, Genome browser-like visualization of single-cell RNA-seq dataset, is available from https://github.com/yuifu/millefy.
