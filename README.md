# metabolomics-tools
A compilation of code helpful in filtering, assessing data quality, and conducting first pass metabolomics analysis.

Divided into:

* [*First pass filtering*](https://github.com/CooperstoneLab/metabolomics-tools/blob/main/Metabolomics-filtering-template_github.md): takes data directly from MZmine, filters based on retention time, cleans up file names, removes unwanted samples, filtering based on presence in the QCs (both by number and by CV), removing features present in the process blanks.
* [*First pass data analysis*](https://github.com/CooperstoneLab/metabolomics-tools/blob/main/Metabolomics-data-analysis-template_github.md): takes data directly from first pass filtering and surveys and handles missing data, imputes missing values, assesses data quality, clusters features coming from the same metabolite using [`notame`](https://github.com/antonvsdata/notame), transforms data, runs and visualizes PCAs, and conducts basic univarite analysis.
