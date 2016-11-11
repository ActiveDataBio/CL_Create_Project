![ADBio][adbio-logo]

# CL_Create_Project
Command line program that generates the files for an Adbio data project repository.

## Want to Contribute?
We welcome contibutions to this repository. Please refer to the [Contributing Guide](CONTRIBUTING.md).
## Using this program
###Requirements
* Java 8
* R

After setting up all requirements download [CLCreateProject.jar][CL-file]. Refer to [How to prepare for my data](#prepare) for instructions on how to... well prepare your data. When your data is prepared run the CLCreateProject.jar file by using
```cmd
java -jar CLCreateProject.jar <Oganism abrivation e.g. hsa> <Destination path> <RData file path> <metadata.tsv path>
```
This will create all the required file in the destination folder. You will need to create a README.md file in this directory, use this template so the information in it will be parse properly and displayed correctly.
```markdown
<!--adbio-title-->
<provide a title for your project>
<!--adbio-description-->
<provide a description for your project>
<!--adbio-funding-->
<!--adbio-publication-->
<!--adbio-organism-->
<provide an organism with this format>
Organism [<name to display in readme file e.g. Homo sapiens (human)>](<link to information about organism e.g. http://www.genome.jp/kegg-bin/show_organism?org=hsa>)
<!------------------------------------------------------------------------------>
<!--you can add any other information here-->

```
Next you need to create a repository on git and push this folder to it. You can also add a png file called icon.png to this repository, this will be displayed in the information for project. Then refer the [Adbio repository](https://github.com/ActiveDataBio/Adbio) to install the visual analytic tool to view your data.

##<a name="prepare"></a> How to prepare for my data
In order to create a new ADBio project using your own data, you need to prepare these two files.

* RData file: it contains the data matrix and clustering information (See '[example data.RData](CLCreateProject/src/target/data.RData)') 
* Tab-delimited text file: it contains the meta data for each experiment (or each sample)  (See '[example_metadata.tsv](CLCreateProject/src/target/metadata.tsv)')

### How to generate RData from my data file.
Please refer to [this repository](https://github.com/ActiveDataBio/adbio_tutorial/) for instructions and tutorials to convert your data to a RData.

### How to use a meta data file.
please go to [here](https://github.com/ActiveDataBio/adbio_tutorial/blob/master/tutorial_2_metadata.ipynb) and refer to the tutorial.

[adbio-logo]:https://adbio.pnnl.gov/bioviz/images/activeData-biglogo.png
[CL-file]:https://github.com/ActiveDataBio/CL_Create_Project/raw/master/CLCreateProject/src/target/CLCreateProject.jar "CLCreateProject.jar"
