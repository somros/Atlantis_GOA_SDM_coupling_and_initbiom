# validation tables for Alaska and Canada

library(data.table)
library(dplyr)

val_files_ak <- list.files('../../sdmTMB_Alaska_stages/output/validation_tables/', full.names = T)
val_files_bc <- list.files('../../sdmTMB_Canada_stages/output/validation_tables/', full.names = T)

val_ak <- rbindlist(lapply(val_files_ak, 
                        function(x) (read.csv(x) %>% mutate(Stage=substr(x,(nchar(x)-9), (nchar(x)-9))) %>% distinct() ))) %>%
  select(Group,Stage,Convergence:NRMSE...)

val_bc <- rbindlist(lapply(val_files_bc, 
                           function(x) (read.csv(x) %>% mutate(Stage=substr(x,(nchar(x)-8), (nchar(x)-8))) %>% distinct() ))) %>%
  select(Group,Stage,Convergence:NRMSE...)

### convergence messages

val_ak %>% group_by(Message) %>% tally()
val_bc %>% group_by(Message) %>% tally() # about 50:50 4 and 5 

### max gradients being <0.001

val_ak %>% filter(abs(Max.gradient)<0.001) %>% nrow() %>% `/` (nrow(val_ak)) # all good
val_bc %>% filter(abs(Max.gradient)<0.001) %>% nrow() %>% `/` (nrow(val_bc)) # all good

### range
val_ak %>% filter(Practical.range..km.>20) %>% nrow() %>% `/` (nrow(val_ak)) # 1 model has range smaller than 20 km
val_bc %>% filter(Practical.range..km.>20) %>% nrow() %>% `/` (nrow(val_bc)) # 8 models have range smaller than 20 km

### pearson's correlation

val_ak %>% filter(Pearson.s.correlation>0.7) %>% nrow() %>% `/` (nrow(val_ak))
val_bc %>% filter(Pearson.s.correlation>0.7) %>% nrow() %>% `/` (nrow(val_bc)) # only 35% >0.7

# BC models show in general lower R and not as good stats as Alaska models, mainly due to WCHG.
# WCHG helps with predictions at depth of slope-dwelling species

write.csv(val_ak, '../output/validation/val_ak.csv', row.names = F)
write.csv(val_bc, '../output/validation/val_bc.csv', row.names = F)
