# CL_Create_Project
![ADBio][adbio-logo]

Command line program that generates the files for an Adbio project repository
## Want to Contribute?
We welcome contibutions to this repository. Please refer to the [Contibuting Guide](CONTRIBUTING.md).
## Using this program
###Requirements
* Java 8
* R

After setting up all requirements download ![CLCreateProject.jar][CL-file]. Refer to [How to prepare for my data](#prepare) for instructions on how to... well prepare your data. When your data is prepared run the CLCreateProject.jar file by using
```cmd
java -jar CLCreateProject.jar <Oganism abrivation e.g. hsa> <Destination path> <RData file path> <metadata.tsv path>
```
This will create all the required file in the destination folder. After this you need to create a repository on git and push this folder to it. Then follow the adbio guide to view you data.(link coming soon)

##<a name="prepare"></a> How to prepare for my data
In order to create a new ADBio project using your own data, you need to prepare these two files.

* RData file: it contains the data matrix and clustering information
* Tab-delimited text file: it contains the meta data for each experiment (or each sample)  (See '[example_metadata.tsv](https://github.com/ActiveDataBio/adbio_tutorial/blob/master/example_metadata.tsv)')

### How to generate RData from my data file.
You need to install [R](https://cran.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/download/) first. If you have R, please clone or download this repository in your local computer. Then, please open '[generate_rdata_hc.R](https://github.com/ActiveDataBio/adbio_tutorial/blob/master/generate_rdata_hc.R)' and run it with '[example_dataset.csv](https://github.com/ActiveDataBio/adbio_tutorial/blob/master/example_dataset.csv)'. Or please go to [here](https://github.com/ActiveDataBio/adbio_tutorial/blob/master/tutorial_1_generate_rdata.ipynb) and follow the guidelines.

### How to use a meta data file.
please go to [here](https://github.com/ActiveDataBio/adbio_tutorial/blob/master/tutorial_2_metadata.ipynb) and refer to the tutorial.

<!--
## How to upload my data
with some pictures
* Go to 'myproject' page
* Click the 'Create' tab
How to fill in the forms
How to assign tests
Bug reports
-->

[adbio-logo]:https://adbio.pnnl.gov/bioviz/images/activeData-biglogo.png
[CL-file]:https://github.com/ActiveDataBio/CL_Create_Project/raw/master/CLCreateProject/src/target/CLCreateProject.jar "CLCreateProject.jar"
