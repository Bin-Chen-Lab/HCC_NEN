This pipeline is used to identify drug hits for hepatocellular carcinoma ([I'm an inline-style link with title](www.gastrojournal.org/article/S0016-5085(17)30264-0/abstract "Chen B. Gastroenterology, 2017")). It can be applied to other cancers as well. 
![alt text][overview]
[overview]: https://github.com/Bin-Chen-Lab/HCC_NEN/blob/master/workflow_general.png "workflow overview"


The pipeline includes four components:
##tumor_cell_line
Correlate tumor samples and cell lines

Select tumor samples that are not correlated to cell lines

Select cell lines that are correlated to tumor samples

##disease_sig
Create disease gene expression signatures

Validate signatures using external sets

#drug_prediction
Predict drugs using CMap and LINCS

Analyze drug hits including niclosamide and NEN

#reverse_genes
Examine correlation between disease gene signatures and drug gene signatures derived from in vitro/in vivo studies

Identify genes reversed by the drugs

#To run the scripts.
1) dowload data from https://www.synapse.org/#!Synapse:syn6173892/files/

2) set up workspace in main.R

3) run workflow.R under each component. 

Contact Bin Chen (bin.chen@ucsf.edu) for any questions.