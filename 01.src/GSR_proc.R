##Load raw data for resampling

#Data should have measurements and categorizations.
#For example, 
#  measurement category
#1      54.808       11
#2      45.968        8
#3      65.416       10
#4      54.808       10
#5      56.576       10
#6      43.316       11

##Set for resampling 
resam.size = 2
resam.time = 120



##Process resampling

#resam.result is for final resampling result. 
resam.result = c()

for(category.i in min(resam.raw$category):max(resam.raw$category)){
  #Extracting measurements of each category as a data frame format.
  category.n = resam.raw[which(resam.raw$category==category.i), ]
  
  #Iteration is necessary until successful verification is achieved.
  repeat{
    
    resam.iter = c()
    for(resam.i in 1:resam.time){
      #Selection of measurements for resampling according to resam.size is performed
      # through random selection among the rows of the data frame.
      selected.measurement = sample(nrow(category.n), size=resam.size, replace=T)
      
      #Generate resample value, empirically the value is mean value.
      resam.value = mean(category.n[selected.measurement, 'measurement'])
      
      #Gathering resample values.
      resam.iter = c(resam.iter, resam.value)
      
      #Clean spent things.
      rm(resam.i, selected.measurement, resam.value)
    }
    
    #Sort resample values in ascending order.
    resam.iter = resam.iter[order(resam.iter)]
    
    #Verification 1: whether resample values are satisfy a normality or not
    #Shapiro-Wilk test was used for normality verification.
    if(shapiro.test(resam.iter)$p.value<0.05) next
    
    #Verification 2: Whether resample values are satisfy a Gaussian distribution
    #package 'moments' is used for calculating kurtosis and skewness.
    #Gaussian distribution satisfies that kurtosis is 3 and skewness is 0.
    if(round(kurtosis(resam.iter))!=3 | round(skewness(resam.iter))!=0) next
    
    #Verification 3: 95% confidence interval(CI) of resample value encompass more than 95% of original measurements
    #95% CI is range from 2.5th percentile to 97.5th percentile.
    #Becuase number of resample value is 120, 2.5th percentile is 3rd value and 97.5th percentile is 118th value.
    if(length(which(category.n$measurement>=resam.iter[3] & category.n$measurement<=resam.iter[118]))<nrow(category.n)*0.95) next
    
    #If resam.iter successfully passes the three-stage verification, it is stored in resam.result.
    #Becuase resam.iter is character vector, this is stored with category.i in data frame format.
    resam.iter = data.frame(category=category.i, resam.value=resam.iter)
    resam.result= rbind(resam.result, resam.iter)
    
    #Finally, repeat must be stopped.
    #before break, Clean spent things.
    rm(resam.iter)
    break
  }
  #Clean spent things.
  rm(category.i, category.n)
}

##Generating optimal polynomial regression curve
#The previous "resam.result" must be necessary.

##Set for conducting regression analysis.
#Empirically, upper limit is 97.5th percentile, and lower limit is 2.5th percentile.
#At this time, upper limit is generated.
reg.target = 97.5

#For polynomial regression, k-fold cross-validation is used.
#Since the entire gestational period is divided into trimester, k is set to 3.
reg.k=3

#Empirically, degree is searched between 1 and 5.
reg.degree = 5

#Regression analysis is iterated thousand times.
reg.iter = 1000


##Extract values corresponding to reg.target for each category from resam.result.
target.values = c()
for(category.i in min(resam.result$category):max(resam.result$category)){
  category.n = resam.result[which(resam.result$category==category.i), ]
  
  category.target = category.n$resam.value[(nrow(category.n)+1)*reg.target/100]
  category.target = data.frame(
    target.values=category.target,
    category=category.i
  )
  
  target.values = rbind(target.values, category.target)
  
  #Clean spent things.
  rm(category.i, category.n, category.target)
}

##Mean square error(MSE) is used to define the best degree.
#Regression analysis is repeated in reg.iter
iter.degrees = c()
for(iter.i in 1:reg.iter){
  suffle.values = target.values[sample(nrow(target.values)), ]
  k.fold = cut(seq(1,nrow(suffle.values)), breaks=reg.k, labels=F)
  
  mse = matrix(data=NA, nrow=reg.k, ncol=reg.degree)
  for(i in 1:reg.k){
    index.test = which(k.fold==i, arr.ind=T)
    data.test = suffle.values[index.test, ]
    data.train = suffle.values[-index.test, ]
    
    for(j in 1:reg.degree){
      fit.train = lm(target.values ~ poly(category, j), data=data.train)
      fit.test = predict(fit.train, newdata=data.test)
      
      mse[i, j] = mean((fit.test-data.test$target.values)^2)
      
      #Clean spent things.
      rm(j, fit.train, fit.test)
    }
    #Clean spent things.
    rm(i, index.test, data.test, data.train)
  }
  
  ##Best degree is the degree which has the minimum MSE value
  best.degree = which(colMeans(mse)==min(colMeans(mse)))
  
  #Gather all best degrees
  iter.degrees = c(iter.degrees, best.degree)
  
  #Clean spent things.
  rm(iter.i, suffle.values, k.fold, mse, best.degree)
}

#Among thousand results, most frequent degree is defined as optimal degree.
optimal.degree = data.frame(table(iter.degrees))
optimal.degree = optimal.degree[order(optimal.degree$Freq, decreasing=T), ]
optimal.degree = as.numeric(optimal.degree$iter.degrees[1])

#Generate polynomial regression curve equation
equation.target = lm(data=target.values, formula=target.values~poly(category, optimal.degree, raw=T))
