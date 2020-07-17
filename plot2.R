library(dplyr)
library(ggplot2)

myf <- function(text){
  return(strsplit(text, '\n')[[1]])
}



# plot2

{

c1b <- myf('Allograft rejection
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

c1cd14 <- myf('Allograft rejection
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
c1cd8 <- myf('Antigen processing and presentation
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

c1dc <- myf('Allograft rejection
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

c1fcgr <- myf('Allograft rejection
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

c1mem <- myf('Alzheimer disease
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

c1naive <- myf('Ribosome
')

c1nk <- myf('Alzheimer disease
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

c1plate <- myf('Bacterial invasion of epithelial cells
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


{

  c3b <- myf('Allograft rejection
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
Leishmaniasis
Phagosome
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


  c3cd14 <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Cardiac muscle contraction
Drug metabolism
Graft-versus-host disease
Huntington disease
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
Non-alcoholic fatty liver disease (NAFLD)
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
Viral myocarditis
')

  c3cd8 <- myf('Antigen processing and presentation
Epstein-Barr virus infection
Huntington disease
Parkinson disease
Proteasome
Spliceosome
')

  c3dc <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Cardiac muscle contraction
Cell adhesion molecules (CAMs)
Drug metabolism
Epstein-Barr virus infection
Glycolysis / Gluconeogenesis
Graft-versus-host disease
Hematopoietic cell lineage
Huntington disease
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
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
Salmonella infection
Shigellosis
Spliceosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')

  c3fcgr <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Bacterial invasion of epithelial cells
Cardiac muscle contraction
Drug metabolism
Endocytosis
Epstein-Barr virus infection
Fc gamma R-mediated phagocytosis
Graft-versus-host disease
Huntington disease
Influenza A
Leishmaniasis
Leukocyte transendothelial migration
Necroptosis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
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
Viral myocarditis
')

  c3mem <- myf('Alzheimer disease
Huntington disease
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Ribosome
Spliceosome
Thermogenesis
')

  c3naive <- myf('Ribosome
')

  c3nk <- myf('Alzheimer disease
Antigen processing and presentation
Autophagy
Graft-versus-host disease
Human immunodeficiency virus 1 infection
Natural killer cell mediated cytotoxicity
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Salmonella infection
Spliceosome
Viral myocarditis
')

  c3plate <- myf('Dilated cardiomyopathy (DCM)
Focal adhesion
Leukocyte transendothelial migration
Pathogenic Escherichia coli infection
Phagosome
Platelet activation
Regulation of actin cytoskeleton
Shigellosis
Tight junction
')
}


{
  c5b <- myf('Allograft rejection
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
Leishmaniasis
Phagosome
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

  c5cd14 <- myf('Alzheimer disease
Antigen processing and presentation
Asthma
Autophagy
Huntington disease
Leishmaniasis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Phagosome
Rheumatoid arthritis
Salmonella infection
Thermogenesis
Tuberculosis
Viral myocarditis
')

  c5cd8 <- myf('Parkinson disease
Spliceosome
')

  c5dc <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Cardiac muscle contraction
Drug metabolism
Epstein-Barr virus infection
Graft-versus-host disease
Huntington disease
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
Retrograde endocannabinoid signaling
Rheumatoid arthritis
Ribosome
Salmonella infection
Shigellosis
Spliceosome
Systemic lupus erythematosus
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')

  c5fcgr <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Cardiac muscle contraction
Epstein-Barr virus infection
Huntington disease
Influenza A
Leishmaniasis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
Retrograde endocannabinoid signaling
Salmonella infection
Shigellosis
Spliceosome
Staphylococcus aureus infection
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')

  c5mem <- myf('Autophagy
Huntington disease
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Ribosome
Spliceosome
')

  c5naive <- myf('Ribosome
')

  c5nk <- myf('Autophagy
Natural killer cell mediated cytotoxicity
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Spliceosome
')

  c5plate <- myf('Focal adhesion
Leukocyte transendothelial migration
Pathogenic Escherichia coli infection
Phagosome
Platelet activation
Regulation of actin cytoskeleton
Shigellosis
Tight junction
')


}

c1 <- data.frame()
c1 <- rbind(c1, data.frame(Pathway = c1b, Cell='B'))
c1 <- rbind(c1, data.frame(Pathway = c1cd14, Cell='CD14Mono'))
c1 <- rbind(c1, data.frame(Pathway = c1cd8, Cell='CD8T'))
c1 <- rbind(c1, data.frame(Pathway = c1dc, Cell='DC'))
c1 <- rbind(c1, data.frame(Pathway = c1fcgr, Cell='fcgr'))
c1 <- rbind(c1, data.frame(Pathway = c1mem, Cell='memory'))
c1 <- rbind(c1, data.frame(Pathway = c1naive, Cell='naive'))
c1 <- rbind(c1, data.frame(Pathway = c1nk, Cell='NK'))
c1 <- rbind(c1, data.frame(Pathway = c1plate, Cell='platelet'))
c1$Param <- '0.1'

c3 <- data.frame()
c3 <- rbind(c3, data.frame(Pathway = c3b, Cell='B'))
c3 <- rbind(c3, data.frame(Pathway = c3cd14, Cell='CD14Mono'))
c3 <- rbind(c3, data.frame(Pathway = c3cd8, Cell='CD8T'))
c3 <- rbind(c3, data.frame(Pathway = c3dc, Cell='DC'))
c3 <- rbind(c3, data.frame(Pathway = c3fcgr, Cell='fcgr'))
c3 <- rbind(c3, data.frame(Pathway = c3mem, Cell='memory'))
c3 <- rbind(c3, data.frame(Pathway = c3naive, Cell='naive'))
c3 <- rbind(c3, data.frame(Pathway = c3nk, Cell='NK'))
c3 <- rbind(c3, data.frame(Pathway = c3plate, Cell='platelet'))
c3$Param <- '0.3'

c5 <- data.frame()
c5 <- rbind(c5, data.frame(Pathway = c5b, Cell='B'))
c5 <- rbind(c5, data.frame(Pathway = c5cd14, Cell='CD14Mono'))
c5 <- rbind(c5, data.frame(Pathway = c5cd8, Cell='CD8T'))
c5 <- rbind(c5, data.frame(Pathway = c5dc, Cell='DC'))
c5 <- rbind(c5, data.frame(Pathway = c5fcgr, Cell='fcgr'))
c5 <- rbind(c5, data.frame(Pathway = c5mem, Cell='memory'))
c5 <- rbind(c5, data.frame(Pathway = c5naive, Cell='naive'))
c5 <- rbind(c5, data.frame(Pathway = c5nk, Cell='NK'))
c5 <- rbind(c5, data.frame(Pathway = c5plate, Cell='platelet'))
c5$Param <- '0.5'
datas <- rbind(c1,c3,c5)
datas2 <- datas %>% group_by(Param, Cell) %>% summarise(Count = n())
ggplot(datas2, aes(x = Cell, y = Count, fill = Param)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  geom_text(aes(label = Count), position = position_dodge(0.9), vjust = -1 )


{
m5b <- myf('Allograft rejection
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Cell adhesion molecules (CAMs)
Epstein-Barr virus infection
Graft-versus-host disease
Hematopoietic cell lineage
Human T-cell leukemia virus 1 infection
Inflammatory bowel disease (IBD)
Intestinal immune network for IgA production
Leishmaniasis
Phagosome
Rheumatoid arthritis
Ribosome
Staphylococcus aureus infection
Systemic lupus erythematosus
Th1 and Th2 cell differentiation
Type I diabetes mellitus
Viral myocarditis
')
m5cd14 <- myf('Allograft rejection
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Graft-versus-host disease
Huntington disease
Leishmaniasis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Rheumatoid arthritis
Ribosome
Salmonella infection
Shigellosis
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')
m5cd8 <- myf('Antigen processing and presentation
Parkinson disease
Ribosome
Salmonella infection
Viral myocarditis
')
m5dc <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Cardiac muscle contraction
Drug metabolism
Epstein-Barr virus infection
Graft-versus-host disease
Huntington disease
Influenza A
Intestinal immune network for IgA production
Leishmaniasis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
Retrograde endocannabinoid signaling
Rheumatoid arthritis
Ribosome
Salmonella infection
Shigellosis
Spliceosome
Systemic lupus erythematosus
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')
m5fcgr <- myf('Allograft rejection
Alzheimer disease
Antigen processing and presentation
Asthma
Autoimmune thyroid disease
Autophagy
Cardiac muscle contraction
Graft-versus-host disease
Huntington disease
Leishmaniasis
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
Rheumatoid arthritis
Ribosome
Salmonella infection
Shigellosis
Thermogenesis
Tuberculosis
Type I diabetes mellitus
Viral myocarditis
')
m5mem <- myf('Autophagy
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Salmonella infection
')
m5naive <- myf('Ribosome
')
m5nk <- myf('Antigen processing and presentation
Autophagy
Graft-versus-host disease
Natural killer cell mediated cytotoxicity
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Ribosome
Salmonella infection
Viral myocarditis
')
m5plate <- myf('Bacterial invasion of epithelial cells
Dilated cardiomyopathy (DCM)
Endocytosis
Ferroptosis
Focal adhesion
Leukocyte transendothelial migration
Pathogenic Escherichia coli infection
Phagosome
Platelet activation
Regulation of actin cytoskeleton
Salmonella infection
Shigellosis
Tight junction
')
}

m5 <- data.frame()
m5 <- rbind(m5, data.frame(Pathway = m5b, Cell='B'))
m5 <- rbind(m5, data.frame(Pathway = m5cd14, Cell='CD14Mono'))
m5 <- rbind(m5, data.frame(Pathway = m5cd8, Cell='CD8T'))
m5 <- rbind(m5, data.frame(Pathway = m5dc, Cell='DC'))
m5 <- rbind(m5, data.frame(Pathway = m5fcgr, Cell='fcgr'))
m5 <- rbind(m5, data.frame(Pathway = m5mem, Cell='memory'))
m5 <- rbind(m5, data.frame(Pathway = m5naive, Cell='naive'))
m5 <- rbind(m5, data.frame(Pathway = m5nk, Cell='NK'))
m5 <- rbind(m5, data.frame(Pathway = m5plate, Cell='platelet'))
m5$Param <- '0.5'
{
  fb <- myf('Drug metabolism
Ferroptosis
Phagosome
Ribosome
Tuberculosis
')

  fcd14 <- myf('Allograft rejection
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

  fcd8 <- myf('Antigen processing and presentation
Ferroptosis
Ribosome
Tuberculosis
')

  fdc <- myf('Allograft rejection
Apoptosis
Graft-versus-host disease
Natural killer cell mediated cytotoxicity
Ribosome
Type I diabetes mellitus
Viral protein interaction with cytokine and cytokine receptor
')

  ffcgr <- myf('Hematopoietic cell lineage
PD-L1 expression and PD-1 checkpoint pathway in cancer
Primary immunodeficiency
Ribosome
T cell receptor signaling pathway
Th1 and Th2 cell differentiation
Th17 cell differentiation
')

  fmem <- myf('Allograft rejection
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

  fnaive <- myf('Epstein-Barr virus infection
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

  fnk <- myf('Allograft rejection
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

  fplate <- myf('Chemokine signaling pathway
ECM-receptor interaction
Hematopoietic cell lineage
Platelet activation
Viral protein interaction with cytokine and cytokine receptor
')
}

f <- data.frame()
f <- rbind(f, data.frame(Pathway = fb, Cell='B'))
f <- rbind(f, data.frame(Pathway = fcd14, Cell='CD14Mono'))
f <- rbind(f, data.frame(Pathway = fcd8, Cell='CD8T'))
f <- rbind(f, data.frame(Pathway = fdc, Cell='DC'))
f <- rbind(f, data.frame(Pathway = ffcgr, Cell='fcgr'))
f <- rbind(f, data.frame(Pathway = fmem, Cell='memory'))
f <- rbind(f, data.frame(Pathway = fnaive, Cell='naive'))
f <- rbind(f, data.frame(Pathway = fnk, Cell='NK'))
f <- rbind(f, data.frame(Pathway = fplate, Cell='platelet'))
f$Param <- 'Fisher'

compare <- rbind(c1, c3, c5, f)
ggplot(
  compare %>%
    group_by(Pathway, Cell) %>%
    dplyr::summarise(Count = ifelse(n()==1, 'Unique', 'Duplicate')) %>%
    group_by(Cell, Count) %>% dplyr::summarise(Count2 = n()) ,
    aes(x = Cell, y = Count2, fill = Count)) +
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('') +
  ylab('Number of Pathways')

ggplot( compare , aes(x = Cell, y = Param, colour = Cell) ) +
  geom_jitter()

compare2 <- compare %>%
  inner_join(
    compare %>% group_by(Cell, Pathway) %>% summarise(Count = n())
  )

compare2 %>% filter(Count ==1)
compare2 %>% filter(Count==1) %>% group_by(Pathway, Param) %>% summarise(Count2 = n()) %>% filter(Count2 == 1)


ggplot(
  compare2 %>%
    filter(Count ==1) %>%
    group_by(Cell, Param) %>%
    summarise(Count2 = n()),
  aes(x = Cell, y = Count2, fill = Param)) +
  geom_bar(stat = 'identity', position = position_dodge())

compare2 %>% filter(Count == 1) %>% filter(Cell %in% c('Delta', 'Epsilon', 'Gamma'))



ggplot(
  compare %>%
    inner_join(
      compare %>% group_by(Cell, Pathway) %>% summarise(Count = n())
    )
  , aes(x = Cell, y = Param, colour = as.factor(Count))
  ) +
  geom_jitter(size = 3) +
  scale_colour_manual(values = c('#eb4d4b', '#f0932b','#f9ca24',  '#7ed6df'))


data0.5 <- datas %>% filter(Param=='0.5') %>% group_by(Pathway) %>% summarise(Count = n())
data0.3 <- datas %>% filter(Param=='0.3') %>% group_by(Pathway) %>% summarise(Count = n())
data0.1 <- datas %>% filter(Param=='0.1') %>% group_by(Pathway) %>% summarise(Count = n())

data0.5$Param <- '0.5'
data0.3$Param <- '0.3'
data0.1$Param <- '0.1'

datai <- rbind(data0.5, data0.3, data0.1)
datai$Param <- factor(datai$Param, levels = c('0.5', '0.3', '0.1'), ordered = TRUE)
ggplot(datai , aes(x = Param, y = Count, group = Param, fill = Param)) +
  geom_boxplot(width = 0.3)+
  scale_fill_manual(values = c('#74b9ff', '#a29bfe', '#fd79a8'))

m5 %>% inner_join()

## PANCREAS

{
  c1acinar <- myf('Adherens junction
Adipocytokine signaling pathway
Amino sugar and nucleotide sugar metabolism
Aminoacyl-tRNA biosynthesis
Apoptosis
Bacterial invasion of epithelial cells
Basal transcription factors
Base excision repair
Bladder cancer
Cell cycle
Cellular senescence
Choline metabolism in cancer
Chronic myeloid leukemia
Colorectal cancer
Cysteine and methionine metabolism
Endocytosis
Endometrial cancer
Epithelial cell signaling in Helicobacter pylori infection
Epstein-Barr virus infection
ErbB signaling pathway
Fatty acid degradation
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Fructose and mannose metabolism
Glioma
Glutathione metabolism
Glycolysis / Gluconeogenesis
Glyoxylate and dicarboxylate metabolism
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
HIF-1 signaling pathway
Hippo signaling pathway
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Influenza A
Insulin resistance
Insulin signaling pathway
Kaposi sarcoma-associated herpesvirus infection
Legionellosis
Longevity regulating pathway
Lysine degradation
MAPK signaling pathway
Measles
Mitophagy
Necroptosis
Neurotrophin signaling pathway
NF-kappa B signaling pathway
N-Glycan biosynthesis
NOD-like receptor signaling pathway
Non-small cell lung cancer
Notch signaling pathway
Nucleotide excision repair
Pancreatic cancer
Pancreatic secretion
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Pentose phosphate pathway
Peroxisome
Propanoate metabolism
Prostate cancer
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Pyruvate metabolism
Regulation of actin cytoskeleton
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
Steroid biosynthesis
Thyroid cancer
Tight junction
TNF signaling pathway
Toxoplasmosis
Valine, leucine and isoleucine degradation
Various types of N-glycan biosynthesis
Viral carcinogenesis
Yersinia infection
')

  c1activestell <- myf('Acute myeloid leukemia
Adherens junction
Adipocytokine signaling pathway
AGE-RAGE signaling pathway in diabetic complications
Amino sugar and nucleotide sugar metabolism
Aminoacyl-tRNA biosynthesis
Amoebiasis
AMPK signaling pathway
Apelin signaling pathway
Apoptosis
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Basal transcription factors
Base excision repair
Bladder cancer
Breast cancer
Cell cycle
Cellular senescence
Central carbon metabolism in cancer
Choline metabolism in cancer
Chronic myeloid leukemia
Circadian rhythm
Colorectal cancer
C-type lectin receptor signaling pathway
Cysteine and methionine metabolism
DNA replication
ECM-receptor interaction
Endocytosis
Endometrial cancer
Epstein-Barr virus infection
ErbB signaling pathway
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Gastric cancer
Glioma
Glutathione metabolism
Glycolysis / Gluconeogenesis
Glycosaminoglycan biosynthesis
Glycosaminoglycan biosynthesis
Hedgehog signaling pathway
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
HIF-1 signaling pathway
Hippo signaling pathway
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Influenza A
Insulin resistance
Kaposi sarcoma-associated herpesvirus infection
Longevity regulating pathway
MAPK signaling pathway
Measles
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Necroptosis
Neurotrophin signaling pathway
N-Glycan biosynthesis
NOD-like receptor signaling pathway
Non-small cell lung cancer
Notch signaling pathway
Nucleotide excision repair
p53 signaling pathway
Pancreatic cancer
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Phagosome
Phospholipase D signaling pathway
PI3K-Akt signaling pathway
Platelet activation
Prostate cancer
Proteasome
Protein digestion and absorption
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Regulation of actin cytoskeleton
Relaxin signaling pathway
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA polymerase
RNA transport
Salmonella infection
Shigellosis
Signaling pathways regulating pluripotency of stem cells
Small cell lung cancer
SNARE interactions in vesicular transport
Sphingolipid signaling pathway
Spliceosome
Terpenoid backbone biosynthesis
TGF-beta signaling pathway
Thyroid cancer
Thyroid hormone signaling pathway
TNF signaling pathway
Toxoplasmosis
Transcriptional misregulation in cancer
Ubiquitin mediated proteolysis
Vasopressin-regulated water reabsorption
Viral carcinogenesis
Wnt signaling pathway
Yersinia infection
')

  c1alpha <- myf('Alcoholism
Alzheimer disease
Cardiac muscle contraction
Citrate cycle (TCA cycle)
Dopaminergic synapse
Endocytosis
Glucagon signaling pathway
Glycosylphosphatidylinositol (GPI)-anchor biosynthesis
Growth hormone synthesis, secretion and action
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Huntington disease
Inositol phosphate metabolism
Insulin signaling pathway
Long-term potentiation
Lysosome
mRNA surveillance pathway
Non-alcoholic fatty liver disease (NAFLD)
Oocyte meiosis
Other glycan degradation
Oxidative phosphorylation
Parkinson disease
Phosphatidylinositol signaling system
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Purine metabolism
Retrograde endocannabinoid signaling
RNA degradation
RNA polymerase
SNARE interactions in vesicular transport
Sphingolipid signaling pathway
Spliceosome
Steroid biosynthesis
Synaptic vesicle cycle
Terpenoid backbone biosynthesis
Thermogenesis
Ubiquitin mediated proteolysis
Various types of N-glycan biosynthesis
Vasopressin-regulated water reabsorption
Vibrio cholerae infection
')

  c1beta <- myf('Adipocytokine signaling pathway
Alzheimer disease
AMPK signaling pathway
Autophagy
Autophagy
Basal transcription factors
Cardiac muscle contraction
Cell cycle
Circadian rhythm
Dopaminergic synapse
FoxO signaling pathway
Glucagon signaling pathway
Growth hormone synthesis, secretion and action
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Huntington disease
Inositol phosphate metabolism
Insulin resistance
Insulin secretion
Insulin signaling pathway
Longevity regulating pathway
Longevity regulating pathway
Lysine degradation
Maturity onset diabetes of the young
Mitophagy
mTOR signaling pathway
Nucleotide excision repair
Oxidative phosphorylation
Parkinson disease
Phosphatidylinositol signaling system
Prolactin signaling pathway
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA transport
Spinocerebellar ataxia
Synaptic vesicle cycle
Thermogenesis
Thyroid hormone signaling pathway
Type II diabetes mellitus
Ubiquitin mediated proteolysis
')

  c1delta <- myf('Circadian rhythm
Dopaminergic synapse
Growth hormone synthesis, secretion and action
Insulin secretion
Retrograde endocannabinoid signaling
Spinocerebellar ataxia
Synaptic vesicle cycle
Vasopressin-regulated water reabsorption
')

  c1ductal <- myf('Adherens junction
Adipocytokine signaling pathway
AGE-RAGE signaling pathway in diabetic complications
Amino sugar and nucleotide sugar metabolism
AMPK signaling pathway
Antigen processing and presentation
Apoptosis
Arginine and proline metabolism
Autophagy
Autophagy
Axon guidance
Bacterial invasion of epithelial cells
Base excision repair
Bladder cancer
Cell cycle
Cellular senescence
Central carbon metabolism in cancer
Choline metabolism in cancer
Chronic myeloid leukemia
Circadian rhythm
Citrate cycle (TCA cycle)
Colorectal cancer
C-type lectin receptor signaling pathway
Cysteine and methionine metabolism
Endocytosis
Endometrial cancer
Epithelial cell signaling in Helicobacter pylori infection
Epstein-Barr virus infection
ErbB signaling pathway
Fatty acid degradation
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Fructose and mannose metabolism
Galactose metabolism
Glioma
Glutathione metabolism
Glycolysis / Gluconeogenesis
Glycosylphosphatidylinositol (GPI)-anchor biosynthesis
Glyoxylate and dicarboxylate metabolism
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Hippo signaling pathway
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Influenza A
Kaposi sarcoma-associated herpesvirus infection
Legionellosis
Longevity regulating pathway
Longevity regulating pathway
Lysine degradation
Lysosome
MAPK signaling pathway
Measles
Mitophagy
Necroptosis
Neurotrophin signaling pathway
NF-kappa B signaling pathway
N-Glycan biosynthesis
NOD-like receptor signaling pathway
Non-small cell lung cancer
Notch signaling pathway
Nucleotide excision repair
Other glycan degradation
p53 signaling pathway
Pancreatic cancer
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Pentose phosphate pathway
Peroxisome
Phagosome
Propanoate metabolism
Prostate cancer
Proteoglycans in cancer
Pyruvate metabolism
Regulation of actin cytoskeleton
Renal cell carcinoma
Ribosome biogenesis in eukaryotes
Salmonella infection
Shigellosis
Small cell lung cancer
Sphingolipid signaling pathway
Steroid biosynthesis
Terpenoid backbone biosynthesis
Thyroid cancer
Thyroid hormone signaling pathway
Tight junction
TNF signaling pathway
Toxoplasmosis
Valine, leucine and isoleucine degradation
Viral carcinogenesis
Viral myocarditis
Yersinia infection
')

  c1endo <- myf('Acute myeloid leukemia
Adherens junction
Adipocytokine signaling pathway
AGE-RAGE signaling pathway in diabetic complications
Alzheimer disease
Apelin signaling pathway
Apoptosis
Autophagy
Autophagy
Axon guidance
B cell receptor signaling pathway
Bacterial invasion of epithelial cells
Basal transcription factors
Base excision repair
Bladder cancer
Breast cancer
Cell cycle
Cellular senescence
Chagas disease (American trypanosomiasis)
Choline metabolism in cancer
Chronic myeloid leukemia
Colorectal cancer
Endocytosis
Endometrial cancer
Epithelial cell signaling in Helicobacter pylori infection
Epstein-Barr virus infection
ErbB signaling pathway
Fc epsilon RI signaling pathway
Fc gamma R-mediated phagocytosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Gastric cancer
Glioma
Glycosaminoglycan biosynthesis
Glycosaminoglycan biosynthesis
Growth hormone synthesis, secretion and action
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Hippo signaling pathway
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Huntington disease
Influenza A
Inositol phosphate metabolism
Kaposi sarcoma-associated herpesvirus infection
Leukocyte transendothelial migration
Lysine degradation
MAPK signaling pathway
Measles
Mitophagy
mTOR signaling pathway
Necroptosis
Neurotrophin signaling pathway
NOD-like receptor signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Non-small cell lung cancer
Notch signaling pathway
Nucleotide excision repair
Osteoclast differentiation
Oxidative phosphorylation
Pancreatic cancer
Parathyroid hormone synthesis, secretion and action
Parkinson disease
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Phagosome
Phosphatidylinositol signaling system
Phospholipase D signaling pathway
PI3K-Akt signaling pathway
Platelet activation
Prolactin signaling pathway
Prostate cancer
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Rap1 signaling pathway
Ras signaling pathway
Regulation of actin cytoskeleton
Relaxin signaling pathway
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA polymerase
RNA transport
Salmonella infection
Shigellosis
Signaling pathways regulating pluripotency of stem cells
Small cell lung cancer
Sphingolipid signaling pathway
Spliceosome
T cell receptor signaling pathway
TGF-beta signaling pathway
Th17 cell differentiation
Thermogenesis
Thyroid cancer
Thyroid hormone signaling pathway
Tight junction
TNF signaling pathway
Toxoplasmosis
Transcriptional misregulation in cancer
Ubiquitin mediated proteolysis
VEGF signaling pathway
Vibrio cholerae infection
Viral carcinogenesis
Viral myocarditis
Yersinia infection
')
  c1epsilon <- myf('Alzheimer disease
Antigen processing and presentation
Bacterial invasion of epithelial cells
Basal transcription factors
Cardiac muscle contraction
Citrate cycle (TCA cycle)
Endocytosis
Epithelial cell signaling in Helicobacter pylori infection
ErbB signaling pathway
Fatty acid elongation
Ferroptosis
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Human T-cell leukemia virus 1 infection
Huntington disease
Inositol phosphate metabolism
Lysosome
mRNA surveillance pathway
mTOR signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Nucleotide excision repair
Oocyte meiosis
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Peroxisome
Phagosome
Phosphatidylinositol signaling system
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
RNA degradation
RNA polymerase
RNA transport
Salmonella infection
SNARE interactions in vesicular transport
Spliceosome
Terpenoid backbone biosynthesis
Thermogenesis
Ubiquitin mediated proteolysis
Valine, leucine and isoleucine degradation
Vasopressin-regulated water reabsorption
Vibrio cholerae infection
')

  c1gamma <- myf('Alcoholism
Alzheimer disease
Bacterial invasion of epithelial cells
Cardiac muscle contraction
Circadian rhythm
Citrate cycle (TCA cycle)
Dopaminergic synapse
Endocrine and other factor-regulated calcium reabsorption
Endocytosis
Epithelial cell signaling in Helicobacter pylori infection
Ferroptosis
Glucagon signaling pathway
Growth hormone synthesis, secretion and action
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Huntington disease
Inositol phosphate metabolism
Insulin secretion
Insulin signaling pathway
Long-term potentiation
Lysosome
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Neurotrophin signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Oocyte meiosis
Oxidative phosphorylation
Parkinson disease
Phosphatidylinositol signaling system
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Renal cell carcinoma
Retrograde endocannabinoid signaling
RNA degradation
RNA transport
Salmonella infection
Shigellosis
SNARE interactions in vesicular transport
Sphingolipid signaling pathway
Spinocerebellar ataxia
Spliceosome
Synaptic vesicle cycle
Thermogenesis
Thyroid hormone signaling pathway
Valine, leucine and isoleucine degradation
Various types of N-glycan biosynthesis
Vasopressin-regulated water reabsorption
Vibrio cholerae infection
Viral carcinogenesis
')

  c1macro <- myf('Acute myeloid leukemia
Adipocytokine signaling pathway
AGE-RAGE signaling pathway in diabetic complications
Allograft rejection
Amino sugar and nucleotide sugar metabolism
Amoebiasis
Amyotrophic lateral sclerosis (ALS)
Antigen processing and presentation
Apelin signaling pathway
Apoptosis
Asthma
Autoimmune thyroid disease
Autophagy
Autophagy
B cell receptor signaling pathway
Bacterial invasion of epithelial cells
Base excision repair
Cell cycle
Cellular senescence
Central carbon metabolism in cancer
Chagas disease (American trypanosomiasis)
Chemokine signaling pathway
Choline metabolism in cancer
Chronic myeloid leukemia
Collecting duct acid secretion
Colorectal cancer
C-type lectin receptor signaling pathway
Cysteine and methionine metabolism
DNA replication
Endometrial cancer
Epithelial cell signaling in Helicobacter pylori infection
Epstein-Barr virus infection
ErbB signaling pathway
Fc epsilon RI signaling pathway
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Galactose metabolism
Glioma
Glycosaminoglycan biosynthesis
Glycosaminoglycan biosynthesis
GnRH signaling pathway
Graft-versus-host disease
Growth hormone synthesis, secretion and action
Hematopoietic cell lineage
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Homologous recombination
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Inflammatory bowel disease (IBD)
Influenza A
Inositol phosphate metabolism
Insulin resistance
Intestinal immune network for IgA production
Kaposi sarcoma-associated herpesvirus infection
Legionellosis
Leishmaniasis
Leukocyte transendothelial migration
Longevity regulating pathway
Lysosome
Malaria
MAPK signaling pathway
Measles
Natural killer cell mediated cytotoxicity
Necroptosis
Neurotrophin signaling pathway
NF-kappa B signaling pathway
NOD-like receptor signaling pathway
Non-small cell lung cancer
Oocyte meiosis
Osteoclast differentiation
Other glycan degradation
p53 signaling pathway
Pancreatic cancer
Parathyroid hormone synthesis, secretion and action
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Pertussis
Phagosome
Phosphatidylinositol signaling system
Phospholipase D signaling pathway
Platelet activation
Prolactin signaling pathway
Propanoate metabolism
Prostate cancer
Proteoglycans in cancer
Pyrimidine metabolism
Regulation of actin cytoskeleton
Relaxin signaling pathway
Renal cell carcinoma
Rheumatoid arthritis
Salmonella infection
Shigellosis
Small cell lung cancer
SNARE interactions in vesicular transport
Sphingolipid signaling pathway
Spinocerebellar ataxia
Staphylococcus aureus infection
Systemic lupus erythematosus
T cell receptor signaling pathway
Th1 and Th2 cell differentiation
Th17 cell differentiation
Thyroid hormone signaling pathway
TNF signaling pathway
Toll-like receptor signaling pathway
Toxoplasmosis
Transcriptional misregulation in cancer
Tuberculosis
Type I diabetes mellitus
Type II diabetes mellitus
VEGF signaling pathway
Vibrio cholerae infection
Viral carcinogenesis
Viral myocarditis
Yersinia infection
')

  c1mast <- myf('Acute myeloid leukemia
AGE-RAGE signaling pathway in diabetic complications
Apoptosis
Autophagy
B cell receptor signaling pathway
Base excision repair
Cellular senescence
Central carbon metabolism in cancer
Choline metabolism in cancer
Chronic myeloid leukemia
C-type lectin receptor signaling pathway
Epstein-Barr virus infection
Fc epsilon RI signaling pathway
Fc gamma R-mediated phagocytosis
Glycolysis / Gluconeogenesis
Hepatitis B
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
Influenza A
Leishmaniasis
Long-term potentiation
MAPK signaling pathway
Measles
Natural killer cell mediated cytotoxicity
Neurotrophin signaling pathway
NF-kappa B signaling pathway
NOD-like receptor signaling pathway
Non-small cell lung cancer
Notch signaling pathway
Nucleotide excision repair
Osteoclast differentiation
Oxidative phosphorylation
Pancreatic cancer
PD-L1 expression and PD-1 checkpoint pathway in cancer
Phosphatidylinositol signaling system
Platelet activation
Regulation of actin cytoskeleton
Ribosome
Salmonella infection
Shigellosis
Sphingolipid signaling pathway
Spliceosome
T cell receptor signaling pathway
Th1 and Th2 cell differentiation
Th17 cell differentiation
TNF signaling pathway
Transcriptional misregulation in cancer
Yersinia infection
')

  c1quiestell <- myf('Acute myeloid leukemia
Adherens junction
Adipocytokine signaling pathway
AGE-RAGE signaling pathway in diabetic complications
Alzheimer disease
Aminoacyl-tRNA biosynthesis
Amoebiasis
Antigen processing and presentation
Apelin signaling pathway
Apoptosis
Autophagy
Autophagy
Axon guidance
B cell receptor signaling pathway
Bacterial invasion of epithelial cells
Basal transcription factors
Base excision repair
Bladder cancer
Breast cancer
Cardiac muscle contraction
Cell cycle
Cellular senescence
Central carbon metabolism in cancer
cGMP-PKG signaling pathway
Chagas disease (American trypanosomiasis)
Choline metabolism in cancer
Chronic myeloid leukemia
Citrate cycle (TCA cycle)
Colorectal cancer
C-type lectin receptor signaling pathway
Cushing syndrome
DNA replication
ECM-receptor interaction
Endocytosis
Endometrial cancer
Epstein-Barr virus infection
ErbB signaling pathway
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Gap junction
Gastric cancer
Glioma
Glutathione metabolism
Glycosaminoglycan biosynthesis
Glyoxylate and dicarboxylate metabolism
GnRH signaling pathway
Growth hormone synthesis, secretion and action
Hedgehog signaling pathway
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Hippo signaling pathway
Homologous recombination
Human cytomegalovirus infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Influenza A
Inositol phosphate metabolism
Insulin signaling pathway
Kaposi sarcoma-associated herpesvirus infection
Longevity regulating pathway
MAPK signaling pathway
Measles
Melanoma
Mismatch repair
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Neurotrophin signaling pathway
N-Glycan biosynthesis
Non-alcoholic fatty liver disease (NAFLD)
Non-small cell lung cancer
Notch signaling pathway
Nucleotide excision repair
Osteoclast differentiation
Oxytocin signaling pathway
p53 signaling pathway
Pancreatic cancer
Parathyroid hormone synthesis, secretion and action
Parkinson disease
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Phagosome
Phosphatidylinositol signaling system
Phospholipase D signaling pathway
PI3K-Akt signaling pathway
Platelet activation
Prolactin signaling pathway
Prostate cancer
Proteasome
Protein digestion and absorption
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Purine metabolism
Pyrimidine metabolism
Rap1 signaling pathway
Regulation of actin cytoskeleton
Relaxin signaling pathway
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA polymerase
RNA transport
Salmonella infection
Shigellosis
Signaling pathways regulating pluripotency of stem cells
Small cell lung cancer
Sphingolipid signaling pathway
Spinocerebellar ataxia
Spliceosome
T cell receptor signaling pathway
Terpenoid backbone biosynthesis
TGF-beta signaling pathway
Thermogenesis
Thyroid cancer
Thyroid hormone signaling pathway
TNF signaling pathway
Toxoplasmosis
Ubiquitin mediated proteolysis
Vascular smooth muscle contraction
Vasopressin-regulated water reabsorption
Vibrio cholerae infection
Viral carcinogenesis
Yersinia infection
')

  c1schwann <- myf('Acute myeloid leukemia
Adherens junction
AGE-RAGE signaling pathway in diabetic complications
Alzheimer disease
AMPK signaling pathway
Apelin signaling pathway
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Basal transcription factors
Breast cancer
Cell cycle
Cellular senescence
Central carbon metabolism in cancer
Chagas disease (American trypanosomiasis)
Choline metabolism in cancer
Chronic myeloid leukemia
Colorectal cancer
C-type lectin receptor signaling pathway
ECM-receptor interaction
Endocytosis
Endometrial cancer
Epstein-Barr virus infection
ErbB signaling pathway
Fatty acid elongation
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Gap junction
Gastric cancer
Glucagon signaling pathway
Growth hormone synthesis, secretion and action
Hepatitis B
Hepatocellular carcinoma
HIF-1 signaling pathway
Hippo signaling pathway
Human cytomegalovirus infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Huntington disease
Inositol phosphate metabolism
Insulin secretion
Kaposi sarcoma-associated herpesvirus infection
Longevity regulating pathway
Long-term depression
Lysine degradation
Lysosome
MAPK signaling pathway
Measles
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Neurotrophin signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Non-small cell lung cancer
Oxidative phosphorylation
Pancreatic cancer
Parathyroid hormone synthesis, secretion and action
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Phosphatidylinositol signaling system
Phospholipase D signaling pathway
PI3K-Akt signaling pathway
Platelet activation
Prolactin signaling pathway
Prostate cancer
Proteasome
Protein digestion and absorption
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Regulation of actin cytoskeleton
Relaxin signaling pathway
Renal cell carcinoma
Retrograde endocannabinoid signaling
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA transport
Salmonella infection
Shigellosis
Signaling pathways regulating pluripotency of stem cells
Small cell lung cancer
Sphingolipid signaling pathway
Spinocerebellar ataxia
Spliceosome
T cell receptor signaling pathway
Thermogenesis
Thyroid cancer
Thyroid hormone signaling pathway
Thyroid hormone synthesis
Tight junction
Toxoplasmosis
Ubiquitin mediated proteolysis
Vasopressin-regulated water reabsorption
VEGF signaling pathway
Viral carcinogenesis
Yersinia infection
')
}

{
  c3acinar <- myf('Adherens junction
Amino sugar and nucleotide sugar metabolism
Apoptosis
Bacterial invasion of epithelial cells
Cellular senescence
Chronic myeloid leukemia
Colorectal cancer
Cysteine and methionine metabolism
Endocytosis
Endometrial cancer
Epstein-Barr virus infection
ErbB signaling pathway
Ferroptosis
FoxO signaling pathway
Glutathione metabolism
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
Human T-cell leukemia virus 1 infection
Insulin resistance
Insulin signaling pathway
Kaposi sarcoma-associated herpesvirus infection
Longevity regulating pathway
Mitophagy
Neurotrophin signaling pathway
N-Glycan biosynthesis
Pancreatic cancer
Pancreatic secretion
Pathogenic Escherichia coli infection
Peroxisome
Propanoate metabolism
Prostate cancer
Protein export
Protein processing in endoplasmic reticulum
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
Tight junction
TNF signaling pathway
Valine, leucine and isoleucine degradation
Various types of N-glycan biosynthesis
Viral carcinogenesis
Yersinia infection
')
  c3activestell <- myf('Adherens junction
AGE-RAGE signaling pathway in diabetic complications
Amino sugar and nucleotide sugar metabolism
AMPK signaling pathway
Apoptosis
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Basal transcription factors
Cell cycle
Cellular senescence
Central carbon metabolism in cancer
Chronic myeloid leukemia
Colorectal cancer
Cysteine and methionine metabolism
Endocytosis
Endometrial cancer
Epstein-Barr virus infection
ErbB signaling pathway
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Glycosaminoglycan biosynthesis
Glycosaminoglycan biosynthesis
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
HIF-1 signaling pathway
Hippo signaling pathway
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Insulin resistance
Kaposi sarcoma-associated herpesvirus infection
Longevity regulating pathway
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Neurotrophin signaling pathway
N-Glycan biosynthesis
Non-small cell lung cancer
Notch signaling pathway
Nucleotide excision repair
Pancreatic cancer
Pathogenic Escherichia coli infection
Prostate cancer
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Regulation of actin cytoskeleton
Relaxin signaling pathway
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA polymerase
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
SNARE interactions in vesicular transport
Sphingolipid signaling pathway
Spliceosome
TGF-beta signaling pathway
Thyroid cancer
Thyroid hormone signaling pathway
TNF signaling pathway
Ubiquitin mediated proteolysis
Viral carcinogenesis
Yersinia infection
')
  c3alpha <- myf('Alzheimer disease
Citrate cycle (TCA cycle)
Dopaminergic synapse
Endocytosis
Human immunodeficiency virus 1 infection
Huntington disease
Insulin signaling pathway
Lysosome
mRNA surveillance pathway
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
RNA degradation
Spliceosome
Thermogenesis
Ubiquitin mediated proteolysis
Various types of N-glycan biosynthesis
Vibrio cholerae infection
')
  c3beta <- myf('Alzheimer disease
AMPK signaling pathway
Autophagy
Autophagy
Dopaminergic synapse
Huntington disease
Insulin resistance
Insulin signaling pathway
Longevity regulating pathway
Longevity regulating pathway
Mitophagy
mTOR signaling pathway
Oxidative phosphorylation
Parkinson disease
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA transport
Thermogenesis
Thyroid hormone signaling pathway
Ubiquitin mediated proteolysis
')
  c3delta <- myf('Retrograde endocannabinoid signaling
')
  c3ductal <- myf('Adherens junction
AMPK signaling pathway
Apoptosis
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Cellular senescence
Chronic myeloid leukemia
Citrate cycle (TCA cycle)
Colorectal cancer
Endocytosis
Endometrial cancer
Epithelial cell signaling in Helicobacter pylori infection
Epstein-Barr virus infection
ErbB signaling pathway
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
Glutathione metabolism
Glyoxylate and dicarboxylate metabolism
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
Human immunodeficiency virus 1 infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Influenza A
Kaposi sarcoma-associated herpesvirus infection
Longevity regulating pathway
Longevity regulating pathway
Lysosome
Mitophagy
Necroptosis
Neurotrophin signaling pathway
N-Glycan biosynthesis
Non-small cell lung cancer
Pancreatic cancer
Pathogenic Escherichia coli infection
Peroxisome
Prostate cancer
Proteoglycans in cancer
Pyruvate metabolism
Regulation of actin cytoskeleton
Renal cell carcinoma
Ribosome biogenesis in eukaryotes
Salmonella infection
Shigellosis
Small cell lung cancer
Steroid biosynthesis
Tight junction
TNF signaling pathway
Toxoplasmosis
Valine, leucine and isoleucine degradation
Viral carcinogenesis
Yersinia infection
')
  c3endo <- myf('Acute myeloid leukemia
Adherens junction
AGE-RAGE signaling pathway in diabetic complications
Alzheimer disease
Apelin signaling pathway
Apoptosis
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Cellular senescence
Chronic myeloid leukemia
Colorectal cancer
Endocytosis
Endometrial cancer
Epithelial cell signaling in Helicobacter pylori infection
Epstein-Barr virus infection
ErbB signaling pathway
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Glioma
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Huntington disease
Kaposi sarcoma-associated herpesvirus infection
Leukocyte transendothelial migration
MAPK signaling pathway
Mitophagy
mTOR signaling pathway
Neurotrophin signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Non-small cell lung cancer
Notch signaling pathway
Oxidative phosphorylation
Pancreatic cancer
Parkinson disease
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Phagosome
PI3K-Akt signaling pathway
Platelet activation
Prostate cancer
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Rap1 signaling pathway
Ras signaling pathway
Regulation of actin cytoskeleton
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
Sphingolipid signaling pathway
Spliceosome
TGF-beta signaling pathway
Thermogenesis
Thyroid hormone signaling pathway
TNF signaling pathway
Ubiquitin mediated proteolysis
Vibrio cholerae infection
Viral carcinogenesis
Yersinia infection
')
  c3epsilon <- myf('Alzheimer disease
Bacterial invasion of epithelial cells
Citrate cycle (TCA cycle)
Endocytosis
Ferroptosis
Huntington disease
Lysosome
mRNA surveillance pathway
mTOR signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Peroxisome
Phagosome
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
RNA degradation
RNA transport
Salmonella infection
Spliceosome
Thermogenesis
Ubiquitin mediated proteolysis
Valine, leucine and isoleucine degradation
Vibrio cholerae infection
')
  c3gamma <- myf('Alzheimer disease
Bacterial invasion of epithelial cells
Citrate cycle (TCA cycle)
Dopaminergic synapse
Endocytosis
Human immunodeficiency virus 1 infection
Huntington disease
Insulin signaling pathway
Lysosome
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Renal cell carcinoma
Retrograde endocannabinoid signaling
RNA degradation
RNA transport
Salmonella infection
Shigellosis
Spliceosome
Synaptic vesicle cycle
Thermogenesis
Thyroid hormone signaling pathway
Valine, leucine and isoleucine degradation
Various types of N-glycan biosynthesis
Vasopressin-regulated water reabsorption
Vibrio cholerae infection
Viral carcinogenesis
')
  c3macro <- myf('Acute myeloid leukemia
Adipocytokine signaling pathway
Allograft rejection
Amino sugar and nucleotide sugar metabolism
Antigen processing and presentation
Apoptosis
Asthma
Autophagy
Autophagy
B cell receptor signaling pathway
Bacterial invasion of epithelial cells
Cellular senescence
Chagas disease (American trypanosomiasis)
Chemokine signaling pathway
Choline metabolism in cancer
Chronic myeloid leukemia
Colorectal cancer
C-type lectin receptor signaling pathway
Epithelial cell signaling in Helicobacter pylori infection
Epstein-Barr virus infection
ErbB signaling pathway
Fc epsilon RI signaling pathway
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Glioma
Graft-versus-host disease
Hematopoietic cell lineage
Hepatitis B
Hepatitis C
Herpes simplex virus 1 infection
HIF-1 signaling pathway
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
Inflammatory bowel disease (IBD)
Influenza A
Insulin resistance
Intestinal immune network for IgA production
Kaposi sarcoma-associated herpesvirus infection
Legionellosis
Leishmaniasis
Longevity regulating pathway
Lysosome
Measles
Necroptosis
Neurotrophin signaling pathway
NF-kappa B signaling pathway
NOD-like receptor signaling pathway
Non-small cell lung cancer
Osteoclast differentiation
Other glycan degradation
Pancreatic cancer
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Pertussis
Phagosome
Phosphatidylinositol signaling system
Platelet activation
Regulation of actin cytoskeleton
Renal cell carcinoma
Rheumatoid arthritis
Salmonella infection
Shigellosis
SNARE interactions in vesicular transport
Sphingolipid signaling pathway
Staphylococcus aureus infection
Systemic lupus erythematosus
T cell receptor signaling pathway
Th1 and Th2 cell differentiation
Th17 cell differentiation
TNF signaling pathway
Toll-like receptor signaling pathway
Toxoplasmosis
Transcriptional misregulation in cancer
Tuberculosis
Type I diabetes mellitus
Vibrio cholerae infection
Viral carcinogenesis
Viral myocarditis
Yersinia infection
')
  c3mast <- myf('Autophagy
B cell receptor signaling pathway
Chronic myeloid leukemia
Epstein-Barr virus infection
Fc epsilon RI signaling pathway
Fc gamma R-mediated phagocytosis
Hepatitis B
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
Influenza A
Neurotrophin signaling pathway
Osteoclast differentiation
Oxidative phosphorylation
Pancreatic cancer
PD-L1 expression and PD-1 checkpoint pathway in cancer
Platelet activation
Ribosome
Salmonella infection
Shigellosis
Sphingolipid signaling pathway
Spliceosome
Yersinia infection
')
  c3quiestell <- myf('Acute myeloid leukemia
Adherens junction
AGE-RAGE signaling pathway in diabetic complications
Alzheimer disease
Apelin signaling pathway
Apoptosis
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Basal transcription factors
Cell cycle
Cellular senescence
Central carbon metabolism in cancer
Chronic myeloid leukemia
Citrate cycle (TCA cycle)
Colorectal cancer
Endocytosis
Endometrial cancer
Epstein-Barr virus infection
ErbB signaling pathway
Fc gamma R-mediated phagocytosis
Ferroptosis
Fluid shear stress and atherosclerosis
Focal adhesion
FoxO signaling pathway
Glioma
Hepatitis B
Hepatitis C
Hepatocellular carcinoma
HIF-1 signaling pathway
Hippo signaling pathway
Human cytomegalovirus infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Insulin signaling pathway
Kaposi sarcoma-associated herpesvirus infection
Longevity regulating pathway
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Neurotrophin signaling pathway
N-Glycan biosynthesis
Non-alcoholic fatty liver disease (NAFLD)
Non-small cell lung cancer
Notch signaling pathway
Pancreatic cancer
Parkinson disease
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Platelet activation
Prolactin signaling pathway
Prostate cancer
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Regulation of actin cytoskeleton
Relaxin signaling pathway
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA polymerase
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
Sphingolipid signaling pathway
Spliceosome
TGF-beta signaling pathway
Thermogenesis
Thyroid cancer
Thyroid hormone signaling pathway
TNF signaling pathway
Ubiquitin mediated proteolysis
Vasopressin-regulated water reabsorption
Vibrio cholerae infection
Viral carcinogenesis
Yersinia infection
')
  c3schwann <- myf('Adherens junction
Alzheimer disease
AMPK signaling pathway
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Cellular senescence
Chronic myeloid leukemia
Endocytosis
ErbB signaling pathway
Ferroptosis
Focal adhesion
HIF-1 signaling pathway
Human papillomavirus infection
Huntington disease
Longevity regulating pathway
Lysosome
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Neurotrophin signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Pancreatic cancer
Parkinson disease
Pathogenic Escherichia coli infection
Phagosome
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Regulation of actin cytoskeleton
Renal cell carcinoma
Retrograde endocannabinoid signaling
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
Spliceosome
Thermogenesis
Thyroid hormone signaling pathway
Tight junction
Ubiquitin mediated proteolysis
Viral carcinogenesis
')
}

{
  c5acinar <- myf('Apoptosis
Chronic myeloid leukemia
Endocytosis
Glutathione metabolism
Mitophagy
N-Glycan biosynthesis
Pathogenic Escherichia coli infection
Protein processing in endoplasmic reticulum
Ribosome
RNA transport
Salmonella infection
Shigellosis
Tight junction
Valine, leucine and isoleucine degradation
Viral carcinogenesis
')
  c5activestell <- myf('Adherens junction
AGE-RAGE signaling pathway in diabetic complications
AMPK signaling pathway
Apoptosis
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Cellular senescence
Chronic myeloid leukemia
Colorectal cancer
Endocytosis
Endometrial cancer
Epstein-Barr virus infection
Focal adhesion
Glycosaminoglycan biosynthesis
Glycosaminoglycan biosynthesis
Hepatitis B
Hepatocellular carcinoma
Human cytomegalovirus infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Longevity regulating pathway
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Notch signaling pathway
Pancreatic cancer
Pathogenic Escherichia coli infection
Prostate cancer
Proteasome
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Regulation of actin cytoskeleton
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
Spliceosome
TGF-beta signaling pathway
Ubiquitin mediated proteolysis
Viral carcinogenesis
Yersinia infection
')
  c5alpha <- myf('Alzheimer disease
Endocytosis
Huntington disease
Lysosome
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Protein processing in endoplasmic reticulum
Spliceosome
Thermogenesis
Ubiquitin mediated proteolysis
')
  c5beta <- myf('Autophagy
Autophagy
Protein processing in endoplasmic reticulum
RNA transport
Ubiquitin mediated proteolysis
')
  c5delta <- myf('')
  c5ductal <- myf('Apoptosis
Autophagy
Autophagy
Chronic myeloid leukemia
Endocytosis
Epstein-Barr virus infection
Focal adhesion
Lysosome
Pancreatic cancer
Pathogenic Escherichia coli infection
Salmonella infection
Shigellosis
Small cell lung cancer
Tight junction
Viral carcinogenesis
')
  c5endo <- myf('Alzheimer disease
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Cellular senescence
Chronic myeloid leukemia
Colorectal cancer
Endocytosis
Epstein-Barr virus infection
Fluid shear stress and atherosclerosis
Focal adhesion
HIF-1 signaling pathway
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Huntington disease
Kaposi sarcoma-associated herpesvirus infection
Neurotrophin signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Pathogenic Escherichia coli infection
Prostate cancer
Proteasome
Protein processing in endoplasmic reticulum
Renal cell carcinoma
Ribosome biogenesis in eukaryotes
RNA degradation
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
Spliceosome
Thermogenesis
Ubiquitin mediated proteolysis
Viral carcinogenesis
Yersinia infection
')
  c5epsilon <- myf('Alzheimer disease
Endocytosis
Huntington disease
Lysosome
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
RNA degradation
RNA transport
Salmonella infection
Spliceosome
Thermogenesis
Ubiquitin mediated proteolysis
Valine, leucine and isoleucine degradation
Vibrio cholerae infection
')
  c5gamma <- myf('Alzheimer disease
Dopaminergic synapse
Endocytosis
Huntington disease
Lysosome
Mitophagy
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Protein processing in endoplasmic reticulum
Retrograde endocannabinoid signaling
RNA transport
Salmonella infection
Shigellosis
Spliceosome
Thermogenesis
Vibrio cholerae infection
Viral carcinogenesis
')
  c5macro <- myf('Acute myeloid leukemia
Allograft rejection
Antigen processing and presentation
Apoptosis
Asthma
Autophagy
Autophagy
B cell receptor signaling pathway
Chagas disease (American trypanosomiasis)
Chemokine signaling pathway
Chronic myeloid leukemia
Colorectal cancer
C-type lectin receptor signaling pathway
Epstein-Barr virus infection
Fc epsilon RI signaling pathway
Fc gamma R-mediated phagocytosis
Ferroptosis
Graft-versus-host disease
Hepatitis B
Human cytomegalovirus infection
Human immunodeficiency virus 1 infection
Human T-cell leukemia virus 1 infection
Inflammatory bowel disease (IBD)
Influenza A
Intestinal immune network for IgA production
Kaposi sarcoma-associated herpesvirus infection
Leishmaniasis
Lysosome
Measles
Neurotrophin signaling pathway
NF-kappa B signaling pathway
NOD-like receptor signaling pathway
Non-small cell lung cancer
Osteoclast differentiation
Pancreatic cancer
Pathogenic Escherichia coli infection
PD-L1 expression and PD-1 checkpoint pathway in cancer
Pertussis
Phagosome
Platelet activation
Rheumatoid arthritis
Salmonella infection
Shigellosis
Sphingolipid signaling pathway
Systemic lupus erythematosus
T cell receptor signaling pathway
Th1 and Th2 cell differentiation
Th17 cell differentiation
TNF signaling pathway
Toll-like receptor signaling pathway
Toxoplasmosis
Tuberculosis
Type I diabetes mellitus
Viral carcinogenesis
Viral myocarditis
Yersinia infection
')
  c5mast <- myf('Autophagy
Fc epsilon RI signaling pathway
Oxidative phosphorylation
Platelet activation
Ribosome
Salmonella infection
Shigellosis
Spliceosome
')
  c5quiestell <- myf('Adherens junction
Alzheimer disease
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Cell cycle
Cellular senescence
Chronic myeloid leukemia
Colorectal cancer
Endocytosis
Epstein-Barr virus infection
Fluid shear stress and atherosclerosis
Focal adhesion
Hippo signaling pathway
Human cytomegalovirus infection
Human papillomavirus infection
Human T-cell leukemia virus 1 infection
Insulin signaling pathway
Longevity regulating pathway
Mitophagy
mRNA surveillance pathway
mTOR signaling pathway
Non-alcoholic fatty liver disease (NAFLD)
Pancreatic cancer
Parkinson disease
Pathogenic Escherichia coli infection
Prostate cancer
Proteasome
Protein processing in endoplasmic reticulum
Proteoglycans in cancer
Regulation of actin cytoskeleton
Renal cell carcinoma
Ribosome
Ribosome biogenesis in eukaryotes
RNA degradation
RNA transport
Salmonella infection
Shigellosis
Small cell lung cancer
Spliceosome
Thermogenesis
Thyroid hormone signaling pathway
Ubiquitin mediated proteolysis
Viral carcinogenesis
Yersinia infection
')
  c5schwann <- myf('Adherens junction
Alzheimer disease
AMPK signaling pathway
Autophagy
Autophagy
Bacterial invasion of epithelial cells
Chronic myeloid leukemia
Endocytosis
Focal adhesion
Huntington disease
Lysosome
Mitophagy
mRNA surveillance pathway
Non-alcoholic fatty liver disease (NAFLD)
Oxidative phosphorylation
Parkinson disease
Proteasome
Protein export
Protein processing in endoplasmic reticulum
Renal cell carcinoma
Retrograde endocannabinoid signaling
Ribosome
RNA degradation
RNA transport
Salmonella infection
Shigellosis
Spliceosome
Thermogenesis
Tight junction
Ubiquitin mediated proteolysis
Viral carcinogenesis
')
}

{
  facinar <- myf('AGE-RAGE signaling pathway in diabetic complications
Amoebiasis
ECM-receptor interaction
Focal adhesion
Human papillomavirus infection
PI3K-Akt signaling pathway
Protein digestion and absorption
Proteoglycans in cancer
Relaxin signaling pathway
Small cell lung cancer
')
  factivestell <- myf('Maturity onset diabetes of the young
')
  falpha <- myf('cAMP signaling pathway
Glucagon signaling pathway
Maturity onset diabetes of the young
Type I diabetes mellitus
')
  fbeta <- myf('Insulin secretion
Maturity onset diabetes of the young
Type I diabetes mellitus
')
  fdelta <- myf('AMPK signaling pathway
Insulin secretion
Maturity onset diabetes of the young
Type I diabetes mellitus
Type II diabetes mellitus
')
  fductal <- myf('Bile secretion
Cell adhesion molecules (CAMs)
ECM-receptor interaction
Pancreatic secretion
Proximal tubule bicarbonate reclamation
Staphylococcus aureus infection
')
  fendo <- myf('Focal adhesion
Leukocyte transendothelial migration
PI3K-Akt signaling pathway
')
  fepsilon <- myf('AGE-RAGE signaling pathway in diabetic complications
Amoebiasis
ECM-receptor interaction
Focal adhesion
Human papillomavirus infection
Pathways in cancer
PI3K-Akt signaling pathway
Platelet activation
Protein digestion and absorption
Proteoglycans in cancer
Relaxin signaling pathway
')
  fgamma <- myf('Chemical carcinogenesis
Drug metabolism
Drug metabolism
Fat digestion and absorption
Fluid shear stress and atherosclerosis
Glutathione metabolism
Metabolism of xenobiotics by cytochrome P450
Pancreatic secretion
Protein digestion and absorption
Renin-angiotensin system
')
  fmacro <- myf('Fc epsilon RI signaling pathway
')
  fmast <- myf('Cell adhesion molecules (CAMs)
Complement and coagulation cascades
Fc gamma R-mediated phagocytosis
Osteoclast differentiation
Pertussis
Prion diseases
Rheumatoid arthritis
Staphylococcus aureus infection
Systemic lupus erythematosus
Viral myocarditis
')
  fquiestell <- myf('Arrhythmogenic right ventricular cardiomyopathy (ARVC)
Cell adhesion molecules (CAMs)
Dilated cardiomyopathy (DCM)
ECM-receptor interaction
Focal adhesion
Human papillomavirus infection
Hypertrophic cardiomyopathy (HCM)
PI3K-Akt signaling pathway
')
  fschwann <- myf('Cholinergic synapse
Circadian entrainment
Insulin secretion
Maturity onset diabetes of the young
Signaling pathways regulating pluripotency of stem cells
Transcriptional misregulation in cancer
Type I diabetes mellitus
')
}

c1 <- data.frame()
c1 <- rbind(c1, data.frame(Pathway = c1acinar, Cell='Acinar'))
c1 <- rbind(c1, data.frame(Pathway = c1activestell, Cell='Activestell'))
c1 <- rbind(c1, data.frame(Pathway = c1alpha, Cell='Alpha'))
c1 <- rbind(c1, data.frame(Pathway = c1beta, Cell='Beta'))
c1 <- rbind(c1, data.frame(Pathway = c1delta, Cell='Delta'))
c1 <- rbind(c1, data.frame(Pathway = c1ductal, Cell='Ductal'))
c1 <- rbind(c1, data.frame(Pathway = c1endo, Cell='Endo'))
c1 <- rbind(c1, data.frame(Pathway = c1epsilon, Cell='Epsilon'))
c1 <- rbind(c1, data.frame(Pathway = c1gamma, Cell='Gamma'))
c1 <- rbind(c1, data.frame(Pathway = c1macro, Cell='Macro'))
c1 <- rbind(c1, data.frame(Pathway = c1mast, Cell='Mast'))
c1 <- rbind(c1, data.frame(Pathway = c1quiestell, Cell='Quiestell'))
c1 <- rbind(c1, data.frame(Pathway = c1schwann, Cell='Schwann'))
c1$Param <- '0.1'


c3 <- data.frame()
c3 <- rbind(c3, data.frame(Pathway = c3acinar, Cell='Acinar'))
c3 <- rbind(c3, data.frame(Pathway = c3activestell, Cell='Activestell'))
c3 <- rbind(c3, data.frame(Pathway = c3alpha, Cell='Alpha'))
c3 <- rbind(c3, data.frame(Pathway = c3beta, Cell='Beta'))
c3 <- rbind(c3, data.frame(Pathway = c3delta, Cell='Delta'))
c3 <- rbind(c3, data.frame(Pathway = c3ductal, Cell='Ductal'))
c3 <- rbind(c3, data.frame(Pathway = c3endo, Cell='Endo'))
c3 <- rbind(c3, data.frame(Pathway = c3epsilon, Cell='Epsilon'))
c3 <- rbind(c3, data.frame(Pathway = c3gamma, Cell='Gamma'))
c3 <- rbind(c3, data.frame(Pathway = c3macro, Cell='Macro'))
c3 <- rbind(c3, data.frame(Pathway = c3mast, Cell='Mast'))
c3 <- rbind(c3, data.frame(Pathway = c3quiestell, Cell='Quiestell'))
c3 <- rbind(c3, data.frame(Pathway = c3schwann, Cell='Schwann'))
c3$Param <- '0.3'


c5 <- data.frame()
c5 <- rbind(c5, data.frame(Pathway = c5acinar, Cell='Acinar'))
c5 <- rbind(c5, data.frame(Pathway = c5activestell, Cell='Activestell'))
c5 <- rbind(c5, data.frame(Pathway = c5alpha, Cell='Alpha'))
c5 <- rbind(c5, data.frame(Pathway = c5beta, Cell='Beta'))
# c5 <- rbind(c5, data.frame(Pathway = c5delta, Cell='Delta'))
c5 <- rbind(c5, data.frame(Pathway = c5ductal, Cell='Ductal'))
c5 <- rbind(c5, data.frame(Pathway = c5endo, Cell='Endo'))
c5 <- rbind(c5, data.frame(Pathway = c5epsilon, Cell='Epsilon'))
c5 <- rbind(c5, data.frame(Pathway = c5gamma, Cell='Gamma'))
c5 <- rbind(c5, data.frame(Pathway = c5macro, Cell='Macro'))
c5 <- rbind(c5, data.frame(Pathway = c5mast, Cell='Mast'))
c5 <- rbind(c5, data.frame(Pathway = c5quiestell, Cell='Quiestell'))
c5 <- rbind(c5, data.frame(Pathway = c5schwann, Cell='Schwann'))
c5$Param <- '0.5'

f <- data.frame()
f <- rbind(f, data.frame(Pathway = facinar, Cell='Acinar'))
f <- rbind(f, data.frame(Pathway = factivestell, Cell='Activestell'))
f <- rbind(f, data.frame(Pathway = falpha, Cell='Alpha'))
f <- rbind(f, data.frame(Pathway = fbeta, Cell='Beta'))
f <- rbind(f, data.frame(Pathway = fdelta, Cell='Delta'))
f <- rbind(f, data.frame(Pathway = fductal, Cell='Ductal'))
f <- rbind(f, data.frame(Pathway = fendo, Cell='Endo'))
f <- rbind(f, data.frame(Pathway = fepsilon, Cell='Epsilon'))
f <- rbind(f, data.frame(Pathway = fgamma, Cell='Gamma'))
f <- rbind(f, data.frame(Pathway = fmacro, Cell='Macro'))
f <- rbind(f, data.frame(Pathway = fmast, Cell='Mast'))
f <- rbind(f, data.frame(Pathway = fquiestell, Cell='Quiestell'))
f <- rbind(f, data.frame(Pathway = fschwann, Cell='Schwann'))
f$Param <- 'Fisher'

compare <- rbind(c1, c3, c5, f)

compare2 <- compare %>%
  inner_join(
    compare %>% group_by(Cell, Pathway) %>% summarise(Count = n())
  )


ggplot(
  compare %>%
    inner_join(
      compare %>% group_by(Cell, Pathway) %>% summarise(Count = n())
    )
  , aes(x = Cell, y = Param, colour = as.factor(Count))
) +
  geom_jitter(size = 3) +
  scale_colour_manual(values = c('#eb4d4b', '#f0932b','#f9ca24','#7ed6df', '#6c5ce7'))


compare2 %>% filter(Cell %in% c('Delta', 'Epsilon', 'Gamma')) %>% filter(Count == 1)
