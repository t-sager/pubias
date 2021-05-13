# pubias
Code for a R Package which implements identification &amp; correction for Publication Bias using the approach by Andrews &amp; Kasy (2019).


# Installing the package

To install the package, use the following command:

- `devtools::install_github("t-sager/pubias", build_vignettes = TRUE)`

In case you run into problems, try updating R to the most recent version or run `Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")` before trying to install the package. 

To get more information about the procedures implemented in the package and 
the available functions, please check the vignette with `utils::vignette("pubias")`.
