library(qs2)
library(edgeR)
# Read Final clinical data and Raw expression data
gex <- qs_read("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCAN_B_Raw_Expression_Final.qs2")
clin_complete <- qs_read("C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCAN_B_Final_Clinical_Data.qs2")

# TMM normalization for CPM and LCPM data

dgelist <- DGEList(as.matrix(gex))
keep.genes <- filterByExpr(dgelist)
dgelist <- dgelist[keep.genes, , keep=FALSE]
dgelist <- calcNormFactors(dgelist, method = "TMM")
cpm <- cpm(dgelist, log = F)
lcpm <- cpm(dgelist, log = T, prior.count = 1)
qs_save(cpm, "C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCANB_Final_CPM.qs2")
qs_save(lcpm, "C:/Users/empor8949/UiT Office 365/O365-Thesis-Emma - General/Saikat/SCANB_Data/SCANB_Final_LCPM.qs2")
