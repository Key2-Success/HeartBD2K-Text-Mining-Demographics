# HeartBD2K Text Mining Demographics
This was the main project I tackled Summer of 2017 at my internship at HeartBD2K, which resulted in an R Shiny App. This project was essentially a culmination of my summer internship. This site allows a clinician, for instance, to compare two different demographics and find which drugs (RN terms) and MeSH terms (which include phrases such as diseases, demographics, and symptoms) are most statistically distinct between each demographic. The purpose of this project was to centralize all of the back-end analyzing, instead of manually running scrips for each demographic of interest.

## Content
The resulting website includes the following tabs
- *Home* to input the two demographics as well as more information about the project
- *Heat Maps* to show which MesH and RN terms are statistically distinct between each demographic and by how much
- *Tree Maps* to represent the top MeSH terms within each demographic by size and color (color represents the distinctness to that search)
- *PCA Plots* to spatially visualize how similar each MeSH term is to each other within each demographic
- *Correlation Plots* to demonstrate how MeSH terms are correlated within each demographic

I have included screenshots of the website of each tab in the `screenshots` folder.

## Future Work
This project is unfortunately not yet complete, however still functional. This is, essentially, my to-do list:
- find hclust or other way of measuring distance for heatmap
- learn if you can figure out how long a code will run before it's done executing
- compare multiple 10k (for RN) and 1k (for MH) runs to see if terms are consistent/sub is representative
- if mtx code too slow, filter out to cluster that matches demographics and manually calculate distance from each male and female 
to each MH
- fix named_RN parsed where all ) gets cut off
- error handling
- improve UI
- optimize speed
- change middle-aged to middle aged
- re-read in correct RN terms (those with ) at end are all getting cut off)
- fix "can't take in sample greater than population when replace = F" error
- add gif spinner when plots are loading
- get inspired by the ucla pubgraph project
- implement tabs for RN features too
- add confidence intervals on the jaccard index based on sample size
- extend left panel grey box to the bottom
- treemap on one page
- change scale of sliderInput for number of case reports
- make pca plot interactive so it can be hovered upon (plotly?)
- make sure the color of treemap is correct logically
- make sure the correlation plot is actually outputting top MeSH terms
- rewrite prevalence factors as "odds ratio"
- try plotly pca zoom in
- reduce redundancy via tree nodes
- use only single patients
- determine case report size needed for statistical significance based on number of entries in search
- topic modeling on MH to cluster/reduce
- remove searched MH terms in outputs (heart failure vs stroke)
- i think pca plot left side top slider controls both plots
- add in underscore for correlation plot to read in phrases
