
# 2020-05-08-variability

<!-- badges: start -->
<!-- badges: end -->

This project contains the code and results for the analysis of the methylation variability analysis.

This is a [drake](https://github.com/ropensci/drake) project which is a framework for data analysis. For more information about drake see [documentation](https://docs.ropensci.org/drake/) and [manual](https://books.ropensci.org/drake/).

The main parts of this project is the the `packages.R` file, `R/` folder, `doc/` folder and `data-raw/` folder.

`packages.R` is a file which contains all the `library()` calls that will be used in the project.

`data-raw/` is where we keep the data we are reading from. This data should not be modified, if you want to clean these file, do so and save the results to a folder called `data/` to avoid confusion.

the `R/` folder contains a variety of R scripts. The most important one is `plan.R` which outlines the [drake plan](https://books.ropensci.org/drake/plans.html) which contains all the code and structure of the analysis. This will contain code to generate the intermediate objects. Some of the functions that are used are placed in seperate R scrips also located in the `R/` folder.

The `doc/` folder contains the final .Rmd file that contains and illustrates the results. This .Rmd file is knitted along with the drake project to keep up to date with the intermediate objects.
