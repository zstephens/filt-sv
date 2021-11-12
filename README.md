# sv-annotate
A workflow for filtering and annotating SVs from PacBio CLR long reads

**NOTES/TODO-LIST:**
- Only PBSV VCFs are presently supported
- No additional dependencies / executables are presently required
- Final merged VCF
- Better documentation (+ annotated figure)

## example usage

`python3 sv-annotate.py \ `  
`    -i input.vcf \ `  
`    -o output_dir/ \ `  

## --control-freq

The `--control-freq` input option determines the minimum frequency required for SVs from the control set to be used as a substractive filter. Examples:

`--control-freq 0.0`: Filter SVs that were found in *any* of the control samples.  
`--control-freq 0.5`: Filter SVs that were found in at least half of the control samples.  
`--control-freq 1.0`: Only filter SVs that were found in every single control sample.  

Note that the highest SV frequency in the control set is ~0.6, so any value above this effectively disables this filter.
