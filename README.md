# Imputation Pipeline

This pipeline is designed for imputing and filtering genotyping data using Plink, Beagle, and BCFtools. It includes steps for filtering panels, checking duplicated markers, and imputing missing data.


## Dependencies

The following software tools are required to run the pipeline:
- [Plink](https://www.cog-genomics.org/plink/) - Version  v1.90b6.26 64-bit (2 Apr 2022)
- [Beagle](https://faculty.washington.edu/browning/beagle/) - (beagle.r1399.jar)
- [BCFtools](http://samtools.github.io/bcftools/)
- [Java](https://www.oracle.com/java/technologies/javase-downloads.html)

#Usage
./imputation_pipeline.sh -p /path/to/plink -r reference_panel -o /path/to/output -i /path/to/input -t target_panel -B /path/to/beagle.jar -b /path/to/bcftools
