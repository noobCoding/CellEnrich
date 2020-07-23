# plot3.R

# compare between methods


plot3 <- function(compare, Params){
  Params <- c(Params, 'Fisher')
  compare2 <- compare %>%
    inner_join(
      compare %>%
        filter(Param %in% Params) %>%
        group_by(Cell, Pathway) %>%
        summarise(Count = n()) %>%
        mutate(Count2 = ifelse(Count == 1, 'Unique', 'Duplicate')) %>%
        select(-Count)
    ) %>%
    filter(Param %in% Params) %>%
    mutate(Param = ifelse(Count2 == 'Unique', Param, 'Intersect')) %>%
    distinct() %>%
    select(-Count2) %>%
    group_by(Cell, Param) %>%
    summarise(Count = n()) %>%
    mutate(Count = ifelse(Param=='Fisher', -Count, Count)) %>%
    mutate(Param = ifelse(!Param %in% c('Fisher', 'Intersect'), 'CellEnrich', Param))

  maxV <- max(abs(compare2$Count))
  gobj <- ggplot(compare2, aes(x = Cell, y = Count, fill = Param)) +
    geom_bar(
      stat = 'identity',
      position = 'identity',
      width = 0.6,
      colour = '#2d3436'
    ) +
    ylim(-50, 50) +
    scale_fill_manual(values = c('#74b9ff', '#fdcb6e', '#00b894')) +
    labs(title = paste(Params, collapse = ' VS '))
  return(gobj)
}

# UNIQ FUNCTION

  compare %>% inner_join(
    compare %>%
      filter(Param %in% c('Fisher', '0.5')) %>%
      group_by(Pathway, Param) %>%
      summarise(Count = n())
  ) %>%
    filter(Count==1) %>%
    mutate(Param = ifelse(Param=='Fisher', Param, 'CellEnrich')) %>%
    select(-Count) %>%
    write.csv(quote = FALSE, row.names = FALSE)

plot4 <- function(compare, Params){
  Params <- c(Params, 'Fisher')

  compare2 <- compare %>%
    inner_join(
      compare %>%
        filter(Param %in% Params) %>%
        group_by(Param, Pathway) %>%
        summarise(Count = n()) %>%
        mutate(Count2 = ifelse(Count == 1, 'Unique', 'Duplicate')) %>%
        select(-Count)
    ) %>%
    filter(Param %in% Params) %>%
    mutate(Param = ifelse(grepl('Fisher', Param),Param, 'CellEnrich' )) %>%
    mutate(Param = ifelse(Count2 == 'Unique', paste0(Param, 'U'), Param)) %>%
    select(-Count2) %>%
    group_by(Cell, Param) %>%
    summarise(Count = n()) %>%
    mutate(Count = ifelse(grepl('Fisher',Param), -Count, Count))
  #
  #compare2 <- compare %>%
    #filter(Param %in% Params) %>%
    #group_by(Cell, Param) %>%
    #summarise(Count = n()) %>%
    #mutate(Param = ifelse(Param=='Fisher', Param, 'CellEnrich')) %>%
    #rbind(
      #compare %>%
        #inner_join(
          #compare %>%
            #filter(Param %in% Params) %>%
            #group_by(Pathway, Param) %>%
            #summarise(Count = n()) %>%
            #mutate(Unique = ifelse(Count==1, 'Unique', 'Not')) %>%
            #select(-Count)
        #) %>%
        #group_by(Cell, Param) %>%
        #summarise(Count = n()) %>%
        #mutate(Param = ifelse(Param=='Fisher', paste0(Param, 'U'), paste0('CellEnrich', 'U')))
    #) %>%
    #mutate(Count = ifelse(grepl('Fisher', Param), -Count, Count ))

  maxV <- max(abs(compare2$Count))
  gobj <- ggplot(compare2, aes(x = Cell, y = Count, fill = Param)) +
    geom_bar(
      stat = 'identity',
      position = 'identity',
      width = 0.6,
      colour = '#2d3436'
    ) +
    ylim(-50, 50) +
    scale_fill_manual(values = c('#ffeaa7', '#fdcb6e', '#74b9ff',  '#0984e3')) +
    labs(title = paste(Params, collapse = ' VS '))
  return(gobj)
}
plot3(compare, 0.5)
plot3(compare, 0.3)
plot3(compare, 0.1)


plot4(compare, 0.5)
plot4(compare, 0.3)
plot4(compare, 0.1)



# plot5 boxplot odd ratio?

compare <- read.csv('mixture.csv')
compare[which(compare$Cell=='NA?VE'),'Cell'] <- 'NAIVE'
plot5 <- function(value){
  Params <- paste0( c('CELLENRICH', 'SCMERGE'), value)
  compare2 <- compare %>%
    inner_join(
      compare %>%
        filter(Param %in% Params) %>%
        group_by(Cell, Pathway) %>%
        dplyr::summarise(Count = n()) %>%
        mutate(Count2 = ifelse(Count == 1, 'Unique', 'Duplicate')) %>%
        select(-Count)
    ) %>%
    filter(Param %in% Params) %>%
    mutate(Param = ifelse(Count2 == 'Unique', Param, 'Intersect')) %>%
    distinct() %>%
    select(-Count2) %>%
    group_by(Cell, Param) %>%
    dplyr::summarise(Count = n()) %>%
    mutate(Count = ifelse(Param == paste0('SCMERGE', value), -Count, Count)) %>%
    mutate(Param = ifelse(!Param %in% c(paste0('SCMERGE', value), 'Intersect'), 'CellEnrich', Param))


  maxV <- max(abs(compare2$Count))
  gobj <- ggplot(compare2, aes(x = Cell, y = Count, fill = Param)) +
    geom_bar(
      stat = 'identity',
      width = 0.6,
      colour = '#2d3436'
    ) +
    ylim(-50, 50) +
    scale_fill_manual(values = c('#74b9ff', '#fdcb6e', '#00b894')) +
    labs(title = paste(Params, collapse = ' VS '))
  return(gobj)
}
plot5(0.1)
plot5(0.3)
plot5(0.5)

Cells <- unique(compare$Cell)

res <- data.frame()
for(i in 1:length(Cells)){
  res <- rbind(res,
  compare %>%
    filter(grepl('0.1',Param)) %>%
    inner_join(
      compare %>%
        filter(Cell==Cells[i]) %>%
        filter(grepl('0.1',Param)) %>%
        arrange(Pathway) %>%
        group_by(Cell, Pathway) %>%
        dplyr::summarise(Count = n()) %>%
        filter(Count == 1) %>% select(-Count)
    )
  )
}

res %>%
  mutate(Param = ifelse(grepl('CELLENRICH', Param), 'CellEnrich', 'scMerge')) %>%
  arrange(Param) %>%
  write.csv(quote = FALSE, row.names = FALSE, file = 'tmp.csv')

compare %>% inner_join(
  compare %>%
    filter(grepl('0.3', Param)) %>%
    group_by(Pathway, Param) %>%
    dplyr::summarise(Count = n())
) %>%
  filter(Count==1) %>%
  #mutate(Param = ifelse(Param=='Fisher', Param, 'CellEnrich')) %>%
  select(-Count) %>%
  write.csv(quote = FALSE, row.names = FALSE, file = 'tmp.csv')

