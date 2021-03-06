SpatialCPie
===========
SpatialCPie is an R package designed to facilitate cluster evaluation for
spatial transcriptomics data by providing intuitive visualizations that
display the relationship between clusters in order to guide the user during
cluster identification, selection and further downstream applications.

Installation
------------
Using `devtools`, execute the following from the R console:
```r
devtools::install_github(
    "jbergenstrahle/SpatialCPie"
   ,build_opts=c("--no-resave-data", "--no-manual")
)
```

Bioconductor version
--------------------
SpatialCPie is currently on the developmental version of Bioconductor (3.9):
https://bioconductor.org/packages/devel/bioc/html/SpatialCPie.html

Usage
-----
See the vignette:
```r
vignette("SpatialCPie", package="SpatialCPie")
```
