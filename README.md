# ggscinames
ggplot2 geoms for scientific names in figures (italics)
  
  This package provides altered versions of ggplot2 and ggrepel geoms with the intent of making the labeling of scientific 
  names in R figures easier.
  
  #### Install and Dependencies:  
  `install.packages('devtools')`  
  `library(devtools)`  
  `install.packages('grid')`  
  `install_github('jtrachsel/ggscinames')`  
  `library(grid)`  
  `library(ggscinames)`  
  
  If you want to use the repel geoms you will need to install `ggrepel` as well  
  `install.packages('ggrepel')`
  
  #### Usage
  ggscinames provides 4 new geoms for ggplot2 plots, these 'new' geoms are really just altered versions of the existing 
  text and label geoms.  These new geoms take the same arguments as their analogous geoms with a few exceptions.
  Within the aesthetic `aes()` call instead of the usual `label = ...`,  
  These geoms require at least one of the following:  
  + `sci`  
    - This variable will be displayed in italics  
  + `nonsci`  
    - This variable will be displayed in regular text  
    
  If both of these parameters are specified then the resulting labels on the plot will be a concatenation of
  "`sci` `nonsci`".  
  Additionally there is another parameter:
  + `important`  
    - Should be a logical vector of the same length as either sci or nonsci.  
    - Either `sci` or `nonsci` must be provided with `important`
    - values in `sci` and/or `nonsci` corresponding to `TRUE` values in this vector will be **bolded**
  
  |ggscinames geom | based on|
  |---|---|
  |`geom_text_sciname` | `ggplot2::geom_text`|
  |`geom_label_sciname` | `ggplot2::geom_label`|
  |`geom_text_sciname_repel`| `ggrepel::geom_text_repel`|
  |`geom_label_sciname_repel` | `ggrepel::geom_label_repel`|
  
  
  
