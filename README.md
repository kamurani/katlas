# katlas

Python implementation of the kinase motif atlas for substrate specificity analysis.

## TODO:

- support `asterisk` and `central` positional formats. 
- apply distribution statistics to find percentiles across kinome. 
- include Ser/Thr favourability into motif score. 
- clustering of profiles / sequence based clustering supported. 
- able to plot heatmaps of profiles / generate sequence logos. 

Somehow show importance of residues across kinases? (enrichment analysis? entropy?) 
- and visualise the 'direction' or position of each residue in feature space (e.g. negatively charged)

Able to work backwards and produce 'optimum' motif for a given kinase? (This may be what we use for 3D analysis) 


## Reference

An atlas of substrate specificities for the human serine/threonine kinome

Jared L. Johnson, Tomer M. Yaron, Emily M. Huntsman, Alexander Kerelsky, Junho Song, Amit Regev, Ting-Yu Lin, Katarina Liberatore, Daniel M. Cizin, Benjamin M. Cohen, Neil Vasan, Yilun Ma, Konstantin Krismer, Jaylissa Torres Robles, Bert van de Kooij, Anne E. van Vlimmeren, Nicole Andrée-Busch, Norbert F. Käufer, Maxim V. Dorovkov, Alexey G. Ryazanov, Yuichiro Takagi, Edward R. Kastenhuber, Marcus D. Goncalves, Benjamin D. Hopkins, Olivier Elemento, Dylan J. Taatjes, Alexandre Maucuer, Akio Yamashita, Alexei Degterev, Mohamed Uduman, Jingyi Lu, Sean D. Landry, Bin Zhang, Ian Cossentino, Rune Linding, John Blenis, Benjamin E. Turk, Michael B. Yaffe, and Lewis C. Cantley

Nature 613, 759–766 (2023). DOI: 10.1038/s41586-022-05575-3

## Error message 

### Single-letter amino acid code
Valid sequences in the Kinase Library can contain single-letter codes for the 20 canonical amino acids, 'X' to mask a position, '_' to mark truncation, and lower-case s/t/y to denote phosphorylated residues (if the phospho-priming setting is enabled).

### 20 canonical amino acids
A: Alanine
R: Arginine
N: Asparagine
D: Aspartic acid
C: Cysteine
E: Glutamic acid
Q: Glutamine
G: Glycine
H: Histidine
I: Isoleucine
L: Leucine
K: Lysine
M: Methionine
F: Phenylalanine
P: Proline
S: Serine
T: Threonine
A: Tryptophan
Y: Tyrosine
V: Valine
### Phosphorylated residues
s: Phosphoserine
t: Phosphothreonine
y: Phosphotyrosine

