# sholl_processing
 A script to automatically process and plot sholl data, comparing genotypes.

## BEFORE running the script: 
1) Perform sholl analysis, e.g. in Fiji. After running it, store all tracing profiles to be analysed in the same directory.

2) Naming of the tracing profiles. the name of each tracing profile should start like this: 
                "cellx_"
    where "x" is a number that is unique in the whole dataset and doesn't occur twice for cells (this is the ID for the cell).

3) create a .csv file (delimiter ;) that contains information on each cell (age, genotype, etc.). include a column with the assigned number x for each cell, and call it "cellnr". 
    3a) when specifying the genotype, call them just "Ctrl" and "KO". 

4) note down the cell numbers you want to exclude from your analysis. The script will ask you to give these numbers later on. 

## AFTER running the script: 
The script has produced: 
1) a table that combines all tracing profiles with information on each cell. 
2) a graph that plots the average number of intersections for each sholl radius and genotype, +- standard error. 
3) a graph of p values, comparing control and KO with each other. 

Note: you can change the comparisons in the intersections plot by switching out the "genotype" variable with e.g. "age". 


