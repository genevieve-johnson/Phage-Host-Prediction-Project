# Phage-Host-Prediction-Project
A pipeline script for the prediction of phage host ranges.

This pipeline is not yet perfected, but the basic concept layout of the pipeline is within the code. 

This pipeline takes in one or multiple fasta files of bacteriophage genomes. It then uses BLAST+ to remotely compare the bacteriophage sequences to a database of CRISPR spacers from bacterial genomes. The top result of the CRISPR homology is output as a possible host for the bacteriophage. 

Then the pipeline moves to look at genetic homology by running BLAST locally against known bacteriophages in the nr/nt Virsues[ORGN]. The pipeline then takes the accession numbers of the closest related bacteriophages and finds the known bacterial host of each of the related phage and outputs the bacterial hosts as possible hosts for the query bacteriophage entered into the pipeline. 

Finally, the pipeline uses k-mer frequency and correlation of the input bacteriophage sequence and a set of all bacterial sequences. This section of the pipeline then outputs the bacterial species that have the highest k-mer correlation to the input bacteriophage as the possible hosts for the bacteriophage.

The CRISPR and Genetic Homology sections were completed by Genevieve Johnson (me) and the k-mer section was completed by Taylor Miller-Ensminger


