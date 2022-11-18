# ASR_CC
ASR_CC.py ranked the sequences based on the correlation coefficients calculated by analyzing “the number of occurrences of 20 amino acids for the template sequence and the sequence” and “Mean difference for number of AAs per thousands in mesophiles and thermophiles (Gregory A.C. Singer and Donal A. Hickey, Gene 317, 2003, 39-47)”.

[Required calculation environment]
Linux (CentOS 6, 7)
Python >3.7, BioPython, Numpy, Scikit-learn

[Required input data]
・Sequence data in Fasta format (Multiple number of sequences are available)
・Template protein sequence data (only one sequence, fasta format)
