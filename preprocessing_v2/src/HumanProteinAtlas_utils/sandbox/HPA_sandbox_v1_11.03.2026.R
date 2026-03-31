if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("HPAanalyze")
library(HPAanalyze)
# pobranie wybranych datasetów HPA
downloadedData <-  HPAanalyze::hpaDownload(downloadList='histology', version='example')
summary(downloadedData)


EGFR <- HPAanalyze::hpaXml(inputXml='ENSG00000146648')
names(EGFR)

FKBP5 <- HPAanalyze::hpaXml(inputXml = "ENSG00000096060")


EGFRxml <- HPAanalyze::hpaXmlGet('ENSG00000146648')

hpaXmlProtClass(EGFRxml)

hpaXmlTissueExprSum(EGFRxml)


tissueExpression <- hpaXmlTissueExpr(EGFRxml)
summary(tissueExpression)


tissueExpression[[1]]
