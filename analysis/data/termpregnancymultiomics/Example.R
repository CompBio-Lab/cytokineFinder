##This script serves as a time-efficient example, demonstrating the overall LOOCV process applied to one omics dataset. The complete analysis was performed using a similar strategy but using parallel processing to optimize lambda and alpha, followed by a stack generalization layer as described in the article.
load('data/termpregnancymultiomics/Data.Rda')
x=InputData[[7]]
y=featureweeks

##exclude the postpartum samples
index=which(featureweeks<0)
x=x[-index,]
y=y[-index]

library(glmnet)
results=vector()
for (i in seq(length(y))){
    xtrain=x[-i,]
    ytrain=y[-i]
    cv=cv.glmnet(xtrain, ytrain)
    results[i]=predict(cv, t(x[i,]))
}

print(cor.test(results, y))
