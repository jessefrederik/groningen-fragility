Metadata
Quantifying the probability of light damage to masonry structures 
An exploration of crack initiation and progression due to seismic vibrations on masonry buildings with existing damage
Author: Paul Korswagen

See the dissertation for additional details.
---------------------------------------------------

In this repository, there are three types of files:

1) ExperimentalData.mat
Containing measurements from the wall experiments for all masonry wall tests.

2) Zip files with DIC processed data
Compressed .csv files with the horizontal (X) and vertical (Z) displacement fields as determined with DIC. Refer to ExperimentalData.mat for folder and file names associated with the experimental data.

3) ComputationalData.mat
With the results of the extrapolation models. Tabular Diana output files and processed Ψ values are collected.

---------------------------------------------------
Details .mat files

1) WallData_Share is a Matlab struct class table with rows for each wall test identified with its Name such as "Comp45". Each wall is associated with various tests within the "Test" field.
Each wall will have multiple tests. This substructure presents a table with the .Date and .Type of the test. The associated DIC_Ponter is the name of the .zip file with the DIC data. The experimental measurements (Force, Displacement, Drift, etc) are inside the .Data table. The rows of the table linked to the Ψ from DIC are gathered in another table .PsiData.

3) FEM_ResultTable is a Matlab table class. Uncompressed, it is 1.7GB. The table has several variables detailed in the dissertation. Additionally, .DianaTable includes the output from the Diana FEA model. These tables have six variables:
Elmnr. Element number or integration point number.
Eknn. Cracking strain normal direction.
Gknt. Cracking strain in tangential direction.
X0-Z0. Coordinates of the integration point.