#download raw data from GDAC

tumors = c('LIHC') #, 'PAAD', 'THCA', 'UCS', 'UCEC','OV','GBM','LGG' LIHC', 'BRCA','LUAD','LUSC','PRAD', 'COAD', 'LAML','SKCM','BLCA','STAD','KIRC') #STAD no rnaseq2
for (tumor in tumors){
  url = paste('http://gdac.broadinstitute.org/runs/stddata__2014_05_18/data/', tumor, '/20140518/gdac.broadinstitute.org_', tumor, '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2014051800.0.0.tar.gz', sep="")
  if (!file.exists("raw/gdac/rnaseq")){
    dir.create("gdac")
  }
  system(paste('curl ', url, ' -o raw/gdac/rnaseq/temp.raw.tar.gz '))
  system(paste('tar -zxvf', 'raw/gdac/rnaseq/temp.raw.tar.gz -C raw/gdac/rnaseq/'))
}