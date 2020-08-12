require(survival)
require(data.table)
################################## generate csv for simulation  ###################################################
## prepare test data  
data(lung)
lung2=lung[,1:5]
lung2$inst[lung$inst>=2] = 2
lung2$inst[lung$inst>=12] = 3
lung2$inst[is.na(lung2$inst)] <- 3
lung2$inst <- c('site2', 'site3', 'site1')[lung2$inst]   # LETTERS[lung2$inst]
table(lung2$inst)
lung2$status <- ifelse(lung2$status == 2, 1, 0)
lung2$sex <- lung2$sex-1
lung2 <- data.table(lung2)
write.csv(lung2[inst=='site1', -'inst'], file='data/Lung_site1.csv', row.names = F)
write.csv(lung2[inst=='site2', -'inst'], file='data/Lung_site2.csv', row.names = F)
write.csv(lung2[inst=='site3', -'inst'], file='data/Lung_site3.csv', row.names = F)
write.csv(lung2[,-'inst'], file='data/Lung.csv')
