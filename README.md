# MINT
# Language: R
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 1.1, R 4.0.0
# Dependency: knitr_1.28, mixOmics_6.12.1


Plugin to integrate multiple independent studies using Multivariate INTegrative analysis (Rohart et al, 2017).

The plugin accepts a TXT file of keyword-value pairs, tab-delimited:
data: CSV file representing the data collected (rows correspond to samples, columns to quantities measured)
outcome: CSV file containing two columns; the first is the sample, the second is an outcome of the study (i.e. healthy or diseased)
study: CSV file containing two columns; the first is the sample, the second is a unique identifier for which independent study it represents.

The plugin will output several statistics as CSV file commencing with the user-specified output prefix, including principle components, importance of variables, and error rates.  For more information, see this site:
http://mixomics.org/mixmint/

All CSV data in example/ was obtained using datasets from mixOmics, which also includes MINT.
