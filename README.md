# Lab8uIowaMicrobiology

<h3> View the App </h3>
The shiny app can be accessed from this <a href="https://ecelias.shinyapps.io/Lab8uIowaMicrobiology/" target="_blank">link.</a> <br>
Please use the files linked below if you would like to try this Shiny app:<br>
<a href="https://github.com/ecelias/Lab8uIowaMicrobiology/blob/main/level5_example.csv" target="_blank"> Level 5 File Example </a><br>
<a href="https://github.com/ecelias/Lab8uIowaMicrobiology/blob/main/metadata_example.csv" target="_blank"> Metadata File Example </a><br>
**Please note, If viewing the app from the above link, it was deployed on a free shiny.io plan so if usage exceeds 25 hours per month, the app will not be available. 

<h3>Project Overview</h3>

<p>This RShiny app was reconfigured from source code written by Carolyne Huang from Emory University.</p>

<p>While the original app was more than capable of serving its purpose, high traffic made it difficult for students at the University of Iowa to use properly. This app will be used by the University of Iowa specifically to reduce site traffic and offer new features such as the ability to save any graphs as a PNG and view multiple taxonomic levels side-by-side rather than one at a time. Additionally, this version contains improvement to figure format so that they are in a more publishable format. </p>

<p>Unfortunately, the code base provided by Emory University did not provide the version of R that was used and contained multiple errors which rendered the app non-functional. Although the version of the app previously used by the University of Iowa functioned as expected, a large portion of the source code had to be altered or updated with non-depracted functions to work correctly.  </p>

<p>As the University of Iowa's Microbiology & Immunology department values giving students a wide variety of opportunities to explore careers in different aspects of the field, this project also sought to provide detailed comments and robust documentation within this repository so any future students interested in bioinformatics are able to play with this code should they desire a starting point. </p>

<h3>R Version Information and Project Dependencies</h3>
<strong>R Version: 4.4.2 "Pile of leaves"</strong> <br><br>
Platform: aarch64-apple-darwin20 <br><br>
Project Dependencies:
<ul>
  <li>BiocManager, version=3.20</li>
  <li>shiny, version=1.10.0</li>
  <li>tidyverse, version=2.0.0</li>
  <li>ggplot2, version=3.5.1</li>
  <li>vegan, version=2.6-10</li>
  <li>data.table, version=1.16.4</li>
  <li>readr, version=2.1.5</li>
  <li>phyloseq, version=1.50.0</li>
  <li>phyloseqCompanion, version=1.1</li>
  <li>broom, version=1.0.7</li>
  <li>shinycssloaders, version=1.1.0</li>
  <li>bslib, version=0.8.0</li>
  <li>rsconnect, version=1.3.4</li>
  <li>igraph, version=2.1.2</li>
  <li>purr, version=1.0.2</li>
</ul>

<h3>App Deployment:</h3>
This app was deployed to a server using <a href="https://www.shinyapps.io" target="_blank"> shinyapps.io</a> <br>
<a href="https://www.youtube.com/watch?v=1g7IAUWD7P0&t=2s" target="_blank"> Tutorial</a> for running the app locally and deploying to shinyapps.io server 

<h3>Contact:</h3>
For questions regarding this software, please contact <strong>Regina McGrane</strong> (University of Iowa, Microbiology & Immunology) at regina-mcgrane@uiowa.edu
