# Genome Wide Cardiovascular Association

This project is an attempt to categorize gene variations into cardiovascular disease bins using MATLAB 9.7 to build a neural network (NN).

## The Database

Gene Expression Omnibus (GEO) is a searchable database, where keyword searches can query data sets of the human genome. The data used in this code was found by quereing the terms **cardiovascular** and **heart disease** One platform included data sets of samples taken from the left ventricular wall of human males. 5 patients of normal human heart physiology, 12 with ischemic cardiomyopathy, and 12 with dilated cardiomyopathy.

Affymetrix collected the data sets used. Their protocol for transcription profiling is done so with micro arrays. They run distinct sequences in parallel which makes the technique resilient to issues detecting and measuring low abundance transcripts, or rare alternative splicing events.

The database is hosted [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1053922)

## Genome-Wide Association Studies (GWAS)

GWAS is a method of associating specefic genetic variations with particular diseases. This association technique involves using many individuals' entire genome to look for genetic markers that can be used to predict the presence of certain diseases.

#### The ```GWAS.m``` file 

This file uses the GWAS method for association. It builds and trains the NN, then displays the probablilty, time, and % error. The ouput is going to be 7 txt files that include an array of output values. These values can then be graphed to see relational data more clearly.

---
**Note:** Lines 11-43 can be commented out after the first run as these lines are just collecting the data from the database.


## Sequential Forward Selection (SFS)

SFS is a method of associating specefic genetic variations with particular diseases. This association technique differs from GWAS because instead of using an individual's entire genome it focuses on a seubset of the genome that has high entropuy. It uses this subset to look for genetic markers that can be used to predict the presence of certain diseases.

#### The ```SFS.m``` file 

This file uses the SFS method for association. It builds and trains the NN, then displays the probablilty, time, and % error. The ouput is going to be 6 txt files that include an array of output values. These values can then be graphed to see relational data more clearly.

---
**Note:** Lines 10-42 can be commented out after the first run as these lines are just collecting the data from the database.

## Downloading the Code

1. Open Console or Terminal
2. ```cd``` into working directory 
3. Run ```git clone``` + the clone url

By default the ```master``` branch will be checked out. A different branch can be checked out using ```git checkout``` + the branch name. 

## Data From the First Run

In the first run GWAS had a slightly lower mean squared error as well as being considerably faster than SFS. This means that based on our research on genome based machine learning, there is a significant advantage to using the GWAS method over SFS.

![alt text](https://static.wixstatic.com/media/7d0401_1907af8ddb83451b90d81248626ec280~mv2.png/v1/fill/w_703,h_430,al_c,lg_1,q_80/GWAS_SFS.webp)

## Going Forward

Ideally this project would be expanded to include more cardiovascular disease datasets and additional species datasets. We would also like to compare the outputs to proven genetic indicators to verify our results. **Any updates to the code base are welcome in the form of commit or merge requests.**
