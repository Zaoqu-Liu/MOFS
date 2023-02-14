library(neuralnet)
load('fit.rda')
d <- read.csv('InputExample.csv')
if (ncol(d)>2&colnames(d)[1]!='Feature'&colnames(d)[2]!='Value') {
  message('The data you entered does not meet the format, please refer to the InputExample data!')
}
if (sum(d$Feature%in%feature$ID)<22) {
  message('The data you entered does not fit the model!')
}
d <- d[d$Feature%in%feature$ID,]
d2 <- as.data.frame(matrix(d$Value,nrow = 1))
colnames(d2) <- d$Feature
d2 <- d2[,feature$ID]
mlppre <- predict(mlpcla,d2)
mlpprelab <- apply(mlppre , 1, which.max)

message('Probability:\nMOFS1 = ',round(mlppre[1],6),'\nMOFS2 = ',round(mlppre[2],6),'\nMOFS3 = ',round(mlppre[3],6))
message('Hence, this sample was identified as MOFS',mlpprelab)











