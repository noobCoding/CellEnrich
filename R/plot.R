myf <- function(text){
  return(strsplit(text, '\n')[[1]])
}


{
scM01_B <- myf('Allograft rejection
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Cell adhesion molecules (CAMs)
Epstein-Barr virus infection
Graft-versus-host disease
Hematopoietic cell lineage
Human T-cell leukemia virus 1 infection
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Kaposi sarcoma-associated herpesvirus infection
Leishmaniasis
Leukocyte transendothelial migration
Phagosome
Rheumatoid arthritis
Ribosome
Spliceosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Th1 and Th2 cell differentiation
Th17 cell differentiation
Toxoplasmosis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')}

{
scM01_CD14Mono <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Bacterial invasion of epithelial cells
Cardiac muscle contraction
Cell adhesion molecules (CAMs)
Drug metabolism
Endocytosis
Epstein-Barr virus infection
Fc gamma R-mediated phagocytosis
Fluid shear stress and atherosclerosis
Graft-versus-host disease
Huntington disease
Influenza A
Intestinal immune network for IgA production
Kaposi sarcoma-associated herpesvirus infection
Leishmaniasis
Leukocyte transendothelial migration
Necroptosis
Non-alcoholic fatty liver disease (NAFLD)
Osteoclast differentiation
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Pertussis
Phagosome
Proteasome
Retrograde endocannabinoid signaling
Rheumatoid arthritis
Ribosome
Salmonella infection
Shigellosis
Staphylococcus aureus infection
Systemic lupus erythematosus
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Vibrio cholerae infection
Viral myocarditis
')

}

scM01_CD8T <- myf('Natural killer cell mediated cytotoxicity
Human immunodeficiency virus 1 infection
Th1 and Th2 cell differentiation
Th17 cell differentiation
Leukocyte transendothelial migration
Antigen processing and presentation
Proteasome
Spliceosome
Shigellosis
Viral myocarditis
Parkinson disease
Graft-versus-host disease
Human T-cell leukemia virus 1 infection
Salmonella infection
Ribosome
')

scM01_DC <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Bacterial invasion of epithelial cells
Cardiac muscle contraction
Cell adhesion molecules (CAMs)
Drug metabolism
Endocytosis
Epstein-Barr virus infection
Fc gamma R-mediated phagocytosis
Fluid shear stress and atherosclerosis
Glycolysis / Gluconeogenesis
Graft-versus-host disease
Hematopoietic cell lineage
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
Huntington disease
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Kaposi sarcoma-associated herpesvirus infection
Leishmaniasis
Leukocyte transendothelial migration
Necroptosis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Regulation of actin cytoskeleton
Retrograde endocannabinoid signaling
Rheumatoid arthritis
Ribosome
RNA transport
Salmonella infection
Shigellosis
Spliceosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Vibrio cholerae infection
Viral myocarditis
')

scM01_fcgr <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Bacterial invasion of epithelial cells
Cardiac muscle contraction
Cell adhesion molecules (CAMs)
Drug metabolism
Endocytosis
Epstein-Barr virus infection
Fc gamma R-mediated phagocytosis
Fluid shear stress and atherosclerosis
Graft-versus-host disease
Human immunodeficiency virus 1 infection
Huntington disease
Influenza A
Kaposi sarcoma-associated herpesvirus infection
Leishmaniasis
Leukocyte transendothelial migration
Lysosome
Natural killer cell mediated cytotoxicity
Necroptosis
Non-alcoholic fatty liver disease (NAFLD)
Osteoclast differentiation
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Platelet activation
Proteasome
Retrograde endocannabinoid signaling
Rheumatoid arthritis
Ribosome
Salmonella infection
Shigellosis
Spliceosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Vibrio cholerae infection
Viral myocarditis
')

scM01_memory <- myf('Alzheimer disease
Autophagy
Human immunodeficiency virus 1 infection
Huntington disease
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Ribosome
Salmonella infection
Shigellosis
Spliceosome
Thermogenesis
')

scM01_naive <- myf('Ribosome
')

scM01_NK <- myf('Antigen processing and presentation
Autophagy
Bacterial invasion of epithelial cells
Endocytosis
Graft-versus-host disease
Human immunodeficiency virus 1 infection
Kaposi sarcoma-associated herpesvirus infection
Leukocyte transendothelial migration
Natural killer cell mediated cytotoxicity
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Proteasome
Protein processing in endoplasmic reticulum
Regulation of actin cytoskeleton
Ribosome
Salmonella infection
Shigellosis
Spliceosome
Viral myocarditis
')

scM01_platelet <- myf('Allograft rejection
Autoimmune thyroid disease
Bacterial invasion of epithelial cells
Dilated cardiomyopathy (DCM)
Endocytosis
Fc gamma R-mediated phagocytosis
Ferroptosis
Focal adhesion
Gap junction
Gastric acid secretion
Glutathione metabolism
Hypertrophic cardiomyopathy (HCM)
Leukocyte transendothelial migration
Necroptosis
Pancreatic secretion
Pathogenic Escherichia coli infection
Phagosome
Platelet activation
Rap1 signaling pathway
Regulation of actin cytoskeleton
Retrograde endocannabinoid signaling
Salmonella infection
Shigellosis
Tight junction
Vasopressin-regulated water reabsorption
Vibrio cholerae infection
Viral carcinogenesis
')


scM <- data.frame()
scM <- rbind(scM, data.frame(Pathway = scM01_B, Cell='B'))
scM <- rbind(scM, data.frame(Pathway = scM01_CD14Mono, Cell='CD14Mono'))
scM <- rbind(scM, data.frame(Pathway = scM01_CD8T, Cell='CD8T'))
scM <- rbind(scM, data.frame(Pathway = scM01_DC, Cell='DC'))
scM <- rbind(scM, data.frame(Pathway = scM01_fcgr, Cell='fcgr'))
scM <- rbind(scM, data.frame(Pathway = scM01_memory, Cell='memory'))
scM <- rbind(scM, data.frame(Pathway = scM01_naive, Cell='naive'))
scM <- rbind(scM, data.frame(Pathway = scM01_NK, Cell='NK'))
scM <- rbind(scM, data.frame(Pathway = scM01_platelet, Cell='platelet'))
{
fisher_B <- myf('Drug metabolism
Ferroptosis
Phagosome
Ribosome
Tuberculosis
')
fisher_CD14Mono <- myf('Allograft rejection
Antigen processing and presentation
Autoimmune thyroid disease
Epstein-Barr virus infection
Graft-versus-host disease
Hematopoietic cell lineage
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
Natural killer cell mediated cytotoxicity
Primary immunodeficiency
Ribosome
T cell receptor signaling pathway
Type I diabetes mellitus
Viral myocarditis
')
fisher_CD8T <- myf('Antigen processing and presentation
Ferroptosis
Ribosome
Tuberculosis
')

fisher_DC <- myf('Allograft rejection
Apoptosis
Graft-versus-host disease
Natural killer cell mediated cytotoxicity
Ribosome
Type I diabetes mellitus
Viral protein interaction with cytokine and cytokine receptor
')
fisher_fcgr <- myf('Hematopoietic cell lineage
PD-L1 expression and PD-1 checkpoint pathway in cancer
Primary immunodeficiency
Ribosome
T cell receptor signaling pathway
Th1 and Th2 cell differentiation
Th17 cell differentiation
')

fisher_memory <-myf('Allograft rejection
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Cell adhesion molecules (CAMs)
Epstein-Barr virus infection
Graft-versus-host disease
Hematopoietic cell lineage
Herpes simplex virus 1 infection
Human T-cell leukemia virus 1 infection
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
Phagosome
Rheumatoid arthritis
Ribosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Th1 and Th2 cell differentiation
Th17 cell differentiation
Toxoplasmosis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')

fisher_naive = myf('Epstein-Barr virus infection
Hematopoietic cell lineage
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Primary immunodeficiency
Ribosome
T cell receptor signaling pathway
Th1 and Th2 cell differentiation
Th17 cell differentiation
')

fisher_NK <- myf('Allograft rejection
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Cell adhesion molecules (CAMs)
Epstein-Barr virus infection
Graft-versus-host disease
Hematopoietic cell lineage
Herpes simplex virus 1 infection
Human T-cell leukemia virus 1 infection
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
Phagosome
Rheumatoid arthritis
Ribosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Th1 and Th2 cell differentiation
Th17 cell differentiation
Toxoplasmosis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')
fisher_platelet <- myf('Chemokine signaling pathway
ECM-receptor interaction
Hematopoietic cell lineage
Platelet activation
Viral protein interaction with cytokine and cytokine receptor
')
}

fisher <- data.frame()

fisher <- rbind(fisher, data.frame(Pathway = fisher_B, Cell='B'))
fisher <- rbind(fisher, data.frame(Pathway = fisher_CD14Mono, Cell='CD14Mono'))
fisher <- rbind(fisher, data.frame(Pathway = fisher_CD8T, Cell='CD8T'))
fisher <- rbind(fisher, data.frame(Pathway = fisher_DC, Cell='DC'))
fisher <- rbind(fisher, data.frame(Pathway = fisher_fcgr, Cell='fcgr'))
fisher <- rbind(fisher, data.frame(Pathway = fisher_memory, Cell='memory'))
fisher <- rbind(fisher, data.frame(Pathway = fisher_naive, Cell='naive'))
fisher <- rbind(fisher, data.frame(Pathway = fisher_NK, Cell='NK'))
fisher <- rbind(fisher, data.frame(Pathway = fisher_platelet, Cell='platelet'))
scM$Method = 'scM'
fisher$Method = 'fisher'

datas <- rbind(scM, fisher)
int <- data.frame()
for(i in unique(datas$Cell)){
  filtered <- datas %>% filter(Cell == i)
  Pathways <- filtered[filtered %>% select(Pathway)  %>% duplicated() %>% which(),] %>% select(Pathway)

  int <- rbind(int, data.frame(Pathways = Pathways, Cell = i, Method = "scM-Fisher" ))
}

datas <- rbind(datas, int)
{
pbmc_B <- myf('Allograft rejection
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
B cell receptor signaling pathway
Cell adhesion molecules (CAMs)
Epstein-Barr virus infection
Graft-versus-host disease
Hematopoietic cell lineage
Human T-cell leukemia virus 1 infection
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
Phagosome
Protein export
Rheumatoid arthritis
Spliceosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Th1 and Th2 cell differentiation
Th17 cell differentiation
Toxoplasmosis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')

pbmc_CD14Mono <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Bacterial invasion of epithelial cells
Cardiac muscle contraction
Drug metabolism
Endocytosis
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Glycolysis / Gluconeogenesis
Graft-versus-host disease
Hematopoietic cell lineage
Huntington disease
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
Leukocyte transendothelial migration
Lysosome
Necroptosis
Non-alcoholic fatty liver disease (NAFLD)
Osteoclast differentiation
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
Retrograde endocannabinoid signaling
Rheumatoid arthritis
Salmonella infection
Shigellosis
Staphylococcus aureus infection
Systemic lupus erythematosus
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Vibrio cholerae infection
Viral myocarditis
')

pbmc_CD8T <- myf('Antigen processing and presentation
Cell adhesion molecules (CAMs)
Endocytosis
Epstein-Barr virus infection
Graft-versus-host disease
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
Leukocyte transendothelial migration
Natural killer cell mediated cytotoxicity
Parkinson disease
Primary immunodeficiency
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Spliceosome
Th1 and Th2 cell differentiation
Th17 cell differentiation
')


pbmc_DC <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Bacterial invasion of epithelial cells
Cardiac muscle contraction
Cell adhesion molecules (CAMs)
Drug metabolism
Endocytosis
Epstein-Barr virus infection
Glycolysis / Gluconeogenesis
Graft-versus-host disease
Hematopoietic cell lineage
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
Huntington disease
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
Leukocyte transendothelial migration
Necroptosis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
Rheumatoid arthritis
Ribosome
RNA transport
Salmonella infection
Shigellosis
Spliceosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Thermogenesis
Toxoplasmosis
Tuberculosis
Type I diabetes mellitus
Vibrio cholerae infection
Viral myocarditis
')

pbmc_fcgr <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Bacterial invasion of epithelial cells
Cardiac muscle contraction
Cell adhesion molecules (CAMs)
Drug metabolism
Endocytosis
Epstein-Barr virus infection
Fc gamma R-mediated phagocytosis
Ferroptosis
Graft-versus-host disease
Human immunodeficiency virus 1 infection
Huntington disease
Influenza A
Intestinal immune network for IgA production
Kaposi sarcoma-associated herpesvirus infection
Leishmaniasis
Leukocyte transendothelial migration
Lysosome
Natural killer cell mediated cytotoxicity
Necroptosis
NOD-like receptor signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Osteoclast differentiation
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Platelet activation
Proteasome
Protein export
Regulation of actin cytoskeleton
Retrograde endocannabinoid signaling
Rheumatoid arthritis
Salmonella infection
Shigellosis
Spliceosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Vibrio cholerae infection
Viral myocarditis
')

pbmc_memory <- myf('Alzheimer disease
Autophagy
Human immunodeficiency virus 1 infection
Huntington disease
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Protein processing in endoplasmic reticulum
Ribosome
Spliceosome
Thermogenesis
')

pbmc_naive <- myf('Ribosome
')

pbmc_NK <- myf('Alzheimer disease
Antigen processing and presentation
Apoptosis
Autophagy
Bacterial invasion of epithelial cells
Endocytosis
Graft-versus-host disease
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Leukocyte transendothelial migration
Natural killer cell mediated cytotoxicity
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
Salmonella infection
Shigellosis
Spliceosome
Viral myocarditis
')

pbmc_platelet <- myf('Bacterial invasion of epithelial cells
Cardiac muscle contraction
Dilated cardiomyopathy (DCM)
Endocytosis
Fc gamma R-mediated phagocytosis
Ferroptosis
Focal adhesion
Gap junction
Glutathione metabolism
Human cytomegalovirus infection
Hypertrophic cardiomyopathy (HCM)
Leukocyte transendothelial migration
Necroptosis
Pancreatic secretion
Pathogenic Escherichia coli infection
Phagosome
Platelet activation
Rap1 signaling pathway
Regulation of actin cytoskeleton
Salmonella infection
Shigellosis
Synaptic vesicle cycle
Tight junction
Vascular smooth muscle contraction
Vibrio cholerae infection
Viral carcinogenesis
')

}

pbmc<- data.frame()
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_B, Cell='B'))
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_CD14Mono, Cell='CD14Mono'))
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_CD8T, Cell='CD8T'))
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_DC, Cell='DC'))
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_fcgr, Cell='fcgr'))
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_memory, Cell='memory'))
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_naive, Cell='naive'))
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_NK, Cell='NK'))
pbmc <- rbind(pbmc, data.frame(Pathway = pbmc_platelet, Cell='platelet'))
pbmc$Method = 'CellEnrich'
datas <- rbind(datas, pbmc)


int2 <- data.frame()

for(i in unique(datas$Cell)){
  filtered <- datas %>% filter(Cell == i, Method %in% c('fisher', 'CellEnrich'))

  Pathways <- filtered[filtered %>% select(Pathway)  %>% duplicated() %>% which(),] %>% select(Pathway)
  if(nrow(Pathways)>0){
    int2 <- rbind(int2, data.frame(Pathways = Pathways, Cell = i, Method = "CellEnrich-Fisher" ))
  }
}
datas <- rbind(datas, int2)

ggplot(datas %>% group_by(Method, Cell) %>% dplyr::summarise(n = n()), aes(x = Cell, y = n, fill = Method) ) +
  geom_bar(stat = 'identity', position = position_dodge())
