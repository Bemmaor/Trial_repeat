#This code implements the simulated Models 1 (NBD), 2 (no covariates + static changepoint), 3 (dynamic changepoint + no covariates), 4 (no changepoint + covariates), 5 (covariates + static changepoint) and 6 (covariates + dynamic changepoint) of Fader/Hardie/Huang (2004). It relies on the parameter estimates that are reported in the article. As input, the user can use the Kiwi Bubbles marketing mix data file (kiwibubbles_mkt_mix.txt)that is available at Professor Bruce Hardie website. The file includes 104 observations (2 markets * 52 weeks).  
# http://www.brucehardie.com/datasets/
# In its current state, the output applies to Model 6. To obtain output for Model 1 for example, the user can comment all the commands that apply to the other models, including the output commands. 
# The code simulates the behavior of 100 samples of 2,799 households each over 52 weeks. There are five output files per run. These files are relatively self-explanatory. The manuscript provides information on the expected output. 
# In a first step, the user can select the first two lines of code to upload the input file. In a second step, he/she can select the rest of the code to run it.  
market<-
  read.table(file.choose(),sep = "",col.names = c("week","market","coupon_stock","advertising_stock","%ACV_any_promotion"))

t1<-matrix(0,nrow=52,ncol=1300)
t2<-matrix(0,nrow=52,ncol=1499)
prob<-matrix(0,nrow=52,ncol=52)
lambda<-matrix(0,nrow=52,ncol=1)
gammaj<-matrix(0,nrow=52,ncol=1)
u0<-matrix(0,nrow=52,ncol=1)
u1<-matrix(0,nrow=52,ncol=1)
transaction_per_week<-matrix(0,nrow=52,ncol=100)
result_simulation<-matrix(0,nrow=100,ncol=6)
result_transaction_times<-matrix(0,nrow=52,ncol=100)
result_cum_week<-matrix(0,nrow=52,ncol=600)
names_result<-c("total_transactions","non_buyers","one-time buyer","two-time buyer","three-time buyer","four-time buyer")
dimnames(result_cum_week)<-list(1:52,rep(names_result,100))
dimnames(result_simulation)<-list(1:100,c("total_transactions","non_buyers","one-time buyer","two-time buyer","three-time buyer","four-time buyer"))

#Model 6
r<-.061
alpha<-80.228/7
#Model 3
#r<-.047
#alpha<-24.057/7
#Model 4
#r<-.076
#alpha<-138.239/7
#Model 5
#r<-.066
#alpha<-97.661/7
#Model 6 - estimated over 12 weeks
#r<-.067
#alpha<-86.573/7
#Model 6 - estimated over 20 weeks
#r<-.071
#alpha<-97.682/7
#Model 2
#r<-.049
#alpha<-26.797/7
#Model 1
#r<-.0787
#alpha<-71.3745/7

#Model 6
AI<-exp(5.204*market$coupon_stock+.012*market[,5])
#Model 3
#AI<-exp(0*market$coupon_stock+0*market[,5])
#Model 4
#AI<-exp(5.182*market$coupon_stock+.014*market[,5])
#Model 5
#AI<-exp(5.059*market$coupon_stock+.012*market[,5])
#Model 6 -as estimated over 12 weeks
#AI<-exp(5.779*market$coupon_stock+.009*market[,5])
#Model 6 - as estimated over 20 weeks
#AI<-exp(4.965*market$coupon_stock+.011*market[,5])
# Model 2 
#AI<-exp(0*market$coupon_stock+0*market[,5])
#Model 1
#AI<-exp(0*market$coupon_stock+0*market[,5])

for (i in 1:52) {
#Model 6
gammaj[i]<-1-.966*(1-exp(-1.367*(i)))
#Model 3
#gammaj[i]<-1-.851*(1-exp(-1.144*(i)))
#Model 4
#gammaj[i]<-1-1.*(1-exp(-300.*(i)))
#Model 5
# gammaj[i]<-1-.912*(1-exp(-300.*(i)))
#Model 6 - as estimated over 12 weeks
# gammaj[i]<-1-1.*(1-exp(-2.353*(i)))
#Model 6 - as estimated over 20 weeks
#gammaj[i]<-1-1.*(1-exp(-1.747*(i)))
# Model 2
#gammaj[i]<-1-.750*(1-exp(-300.*(i)))
#Model 1
#gammaj[i]<-1-1.*(1-exp(-300.*(i)))
}

for (k in 1:100) {
  for (j in 1:1300) {
    t<-0
    t_sum<-0
    t_sum_plus<-1
    t_sum_real<-0
    lambda[1]<-rgamma(1,shape=r,scale=1/alpha)
    u0<-runif(52,min=0,max=1)
    u1<-runif(52,min=0,max=1)
    
    for (i in 2:52) {
      if (u0[i]<=gammaj[i-1]) {
        lambda[i]<-rgamma(1,shape=r,scale=1/alpha)
      }
      else {
        lambda[i]<-lambda[i-1]
      }
    }
    
    for (i in 1:52) {
      if (t_sum>=52) {
        t1[i,j]<-0
      }
      else {
        prob[,i]<-1-exp(-lambda[i,1]*c(rep(0,t_sum),cumsum(AI[t_sum_plus:52])))
        if (u1[i]<=prob[52,i]) {
          t<-sum(prob[t_sum_plus:52,i]>u1[i])
         
          if (t==52-t_sum) {
#           t1[i,j]<-u1[i]/prob[t_sum_plus,i]
#          print(t1[i, j])
          t1[i, j]<- log(1-u1[i])/log(1-prob[t_sum_plus,i])
#          print(t1[i, j])
#          t1[i, j]<- -log(1-u1[i])/lambda[i, 1]
#          print(t1[i, j])
            t_sum_real<-sum(t1[1:i,j])
          }
          else {
#           t1[i,j]<-(u1[i]-prob[52-t,i]+(52-t-t_sum_real)*prob[52-t+1,i])/prob[52-t+1,i]
#            print(t1[i, j]) 
           t1[i, j]<- (52-t-t_sum)*log(1-u1[i])/log(1-prob[52-t,i])
#          print(t1[i, j])
#         t1[i, j]<- -log(1-u1[i])/lambda[i, 1]
#         print(t1[i, j])
            t_sum_real<-sum(t1[1:i,j])
          }
        }
        else {
          t1[i,j]<-0
        }
      }
      t_sum<-floor(sum(t1[1:i,j]))
      t_sum_plus<-t_sum+1
    }
  }
  t1<-t1[,order(t1[1,],decreasing = TRUE)]
  t11<-t1[,1:(1300-sum(t1[1,]==0)+1)]
  
  for (j in 1:1499) {
    t<-0
    t_sum<-0
    t_sum_plus<-1
    t_sum_pplus<-53
    t_sum_real<-0
    lambda[1]<-rgamma(1,shape=r,scale=1/alpha)
    u0<-runif(52,min=0,max=1)
    u1<-runif(52,min=0,max=1)
    
    for (i in 2:52) {
      if (u0[i]<=gammaj[i-1]) {
        lambda[i]<-rgamma(1,shape=r,scale=1/alpha)
      }
      else {
        lambda[i]<-lambda[i-1]
      }
    }
    
    for (i in 1:52) {
      if (t_sum>=52) {
        t2[i,j]<-0
      }
      else {
        prob[,i]<-1-exp(-lambda[i,1]*c(rep(0,t_sum),cumsum(AI[t_sum_pplus:104])))
        if (u1[i]<=prob[52,i]) {
          t<-sum(prob[t_sum_plus:52,i]>u1[i])
          
          if (t==52-t_sum) {
            t2[i,j]<-u1[i]/prob[t_sum_plus,i]
#         t2[i, j]<- log(1-u1[i])/log(1-prob[t_sum_plus,i])
            t_sum_real<-sum(t2[1:i,j])
          }
          else {
#            t2[i,j]<-(u1[i]-prob[52-t,i]+(52-t-t_sum_real)*prob[52-t+1,i])/prob[52-t+1,i]
         t2[i, j]<- (52-t-t_sum)*log(1-u1[i])/log(1-prob[52-t,i])
            t_sum_real<-sum(t2[1:i,j])
          }
        }
        else {
          t2[i,j]<-0
        }
      }
      t_sum<-floor(sum(t2[1:i,j]))
      t_sum_plus<-t_sum+1
      t_sum_pplus<-t_sum+53
    }
    
  }
  t2<-t2[,order(t2[1,],decreasing = TRUE)]
  t21<-t2[,1:(1499-sum(t2[1,]==0)+1)]
  
  mod_t1<-data.frame(t11)
  mod_t1[mod_t1==0]<-NA
  cum_t1<-cumsum(mod_t1)
  cum_t1[is.na(cum_t1)]=0
  cum_t1<-ceiling(cum_t1)
  mod_t2<-data.frame(t21)
  mod_t2[mod_t2==0]<-NA
  cum_t2<-cumsum(mod_t2)
  cum_t2[is.na(cum_t2)]=0
  cum_t2<-ceiling(cum_t2)
  
  #summary
  total_transaction<-sum(cum_t1<53&0<cum_t1)+sum(0<cum_t2&cum_t2<53)
  nb_buyers<-
    sum(cum_t1[1,]<53&0<cum_t1[1,])+sum(cum_t2[1,]<53&0<cum_t2[1,])
  one_time_buyer<-
    sum(cum_t1[1,]<53&0<cum_t1[1,]&(cum_t1[2,]<1 | cum_t1[2,]>52))+
    sum(cum_t2[1,]<53&0<cum_t2[1,]&(cum_t2[2,]<1| cum_t2[2,]>52))
  two_time_buyer<-
    sum(cum_t1[2,]<53&0<cum_t1[2,]&(cum_t1[3,]<1 | cum_t1[3,]>52))+
    sum(cum_t2[2,]<53&0<cum_t2[2,]&(cum_t2[3,]<1 | cum_t2[3,]>52))
  three_time_buyer<-
    sum(cum_t1[3,]<53&0<cum_t1[3,]&(cum_t1[4,]<1 | cum_t1[4,]>52))+
    sum(cum_t2[3,]<53&0<cum_t2[3,]&(cum_t2[4,]<1 | cum_t2[4,]>52))
  four_time_buyer<-
    sum(cum_t1[4,]<53&0<cum_t1[4,]&(cum_t1[5,]<1 | cum_t1[5,]>52))+
    sum(cum_t2[4,]<53&0<cum_t2[4,]&(cum_t2[5,]<1 | cum_t2[5,]>52))
  result_simulation[k,1]<-total_transaction
  result_simulation[k,2]<-1-nb_buyers/2799
  result_simulation[k,3]<-one_time_buyer/2799
  result_simulation[k,4]<-two_time_buyer/2799
  result_simulation[k,5]<-three_time_buyer/2799
  result_simulation[k,6]<-four_time_buyer/2799
  
  #for (i in 1:52) {
  #  transaction_per_week[i,k]<-sum(cum_t1==i)+sum(cum_t2==i)
  #  result_transaction_times[i,k]<-
  #    sum(cum_t1[i,]<53&0<cum_t1[i,]&(cum_t1[i+1,]<1 | cum_t1[i+1,]>52))+
  #    sum(cum_t2[i,]<53&0<cum_t2[i,]&(cum_t2[i+1,]<1 | cum_t2[i+1,]>52))
  #  result_cum_week[i,6*(k-1)+1]<-sum(cum_t1<=i&cum_t1>0)+sum(cum_t2<=i&cum_t2>0)
  #  result_cum_week[i,6*(k-1)+2]<-2799-(sum(cum_t1[1,]<=i&0<cum_t1[1,])+sum(cum_t2[1,]<=i&0<cum_t2[1,]))
  #  result_cum_week[i,6*(k-1)+3]<-sum(cum_t1[1,]<=i&0<cum_t1[1,]&(cum_t1[2,]<1 | cum_t1[2,]>52))+
  #    sum(cum_t2[1,]<=i&0<cum_t2[1,]&(cum_t2[2,]<1 | cum_t2[2,]>52))
  #  result_cum_week[i,6*(k-1)+4]<-sum(cum_t1[2,]<=i&0<cum_t1[2,]&(cum_t1[3,]<1 | cum_t1[3,]>52))+
  #    sum(cum_t2[2,]<=i&0<cum_t2[2,]&(cum_t2[3,]<1 | cum_t2[3,]>52))
  #  result_cum_week[i,6*(k-1)+5]<-sum(cum_t1[3,]<=i&0<cum_t1[3,]&(cum_t1[4,]<1 | cum_t1[4,]>52))+
  #    sum(cum_t2[3,]<=i&0<cum_t2[3,]&(cum_t2[4,]<1 | cum_t2[4,]>52))
  #  result_cum_week[i,6*(k-1)+6]<-sum(cum_t1[4,]<=i&0<cum_t1[4,]&(cum_t1[5,]<1 | cum_t1[5,]>52))+
  #    sum(cum_t2[4,]<=i&0<cum_t2[4,]&(cum_t2[5,]<1 | cum_t2[5,]>52))
  #}
  for (i in 1:52) {
    transaction_per_week[i,k]<-sum(cum_t1==i)+sum(cum_t2==i)
    result_transaction_times[i,k]<-
      sum(cum_t1[i,]<53&0<cum_t1[i,]&(cum_t1[i+1,]<1 | cum_t1[i+1,]>52))+
      sum(cum_t2[i,]<53&0<cum_t2[i,]&(cum_t2[i+1,]<1 | cum_t2[i+1,]>52))
    result_cum_week[i,6*(k-1)+1]<-sum(cum_t1<=i&cum_t1>0)+sum(cum_t2<=i&cum_t2>0)
    result_cum_week[i,6*(k-1)+2]<-2799-(sum(cum_t1[1,]<=i&0<cum_t1[1,])+sum(cum_t2[1,]<=i&0<cum_t2[1,]))
    result_cum_week[i,6*(k-1)+3]<-sum(cum_t1[1,]<=i&0<cum_t1[1,]&(cum_t1[2,]>i|cum_t1[2,]==0))+
      sum(cum_t2[1,]<=i&0<cum_t2[1,]&(cum_t2[2,]>i|cum_t2[2,]==0))
    result_cum_week[i,6*(k-1)+4]<-sum(cum_t1[2,]<=i&0<cum_t1[2,]&(cum_t1[3,]>i|cum_t1[3,]==0))+
      sum(cum_t2[2,]<=i&0<cum_t2[2,]&(cum_t2[3,]>i|cum_t2[3,]==0))
    result_cum_week[i,6*(k-1)+5]<-sum(cum_t1[3,]<=i&0<cum_t1[3,]&(cum_t1[4,]>i|cum_t1[4,]==0))+
      sum(cum_t2[3,]<=i&0<cum_t2[3,]&(cum_t2[4,]>i|cum_t2[4,]==0))
    result_cum_week[i,6*(k-1)+6]<-sum(cum_t1[4,]<=i&0<cum_t1[4,]&(cum_t1[5,]>i|cum_t1[5,]==0))+
      sum(cum_t2[4,]<=i&0<cum_t2[4,]&(cum_t2[5,]>i|cum_t2[5,]==0))
  }
  
}

plot(1:52,rowMeans(transaction_per_week),type="l",col=30,ylab=expression(transaction_per_week),xlab=expression(week))
plot(1:52,cumsum(rowMeans(transaction_per_week)),type="l",col=30,ylab=expression(cum_transaction),xlab=expression(week))
output<-matrix(c(colMeans(result_simulation),sd(result_simulation[,1]),
                 sd(result_simulation[,2]),sd(result_simulation[,3]),
                 sd(result_simulation[,4]),sd(result_simulation[,5]),
                 sd(result_simulation[,6])),nrow=2,byrow=TRUE)
dimnames(output)<-list(c("mean","S.D."),c("total_transactions","non_buyers","one-time buyer","two-time buyer","three-time buyer","four-time buyer"))
setwd("C:/Users/bemmaor/Downloads/kiwibubbles")
print(output)
print(rowMeans(transaction_per_week))
print(rowMeans(result_transaction_times))
#result_cum_week<-data.frame(result_cum_week)
#result<-cbind(rowMeans(result_cum_week[grep('^total',names(result_cum_week))]),rowMeans(result_cum_week[grep('^non',names(result_cum_week))]),rowMeans(result_cum_week[grep('^one',names(result_cum_week))]),rowMeans(result_cum_week[grep('^two',names(result_cum_week))]),rowMeans(result_cum_week[grep('^three',names(result_cum_week))]),rowMeans(result_cum_week[grep('^four',names(result_cum_week))]))
#print(rowMeans(result_cum_week))

# The following output file shows: (i) the cumulative # of transactions, (ii) the cumulative # of nonbuyers (out of 2,799 households), (iii) the cumulative # of one-time buyers, (iv) the cumulative # of two-time buyers, 
#((v) the cumulative # of three-time buyers, and (vi) the cumulative # of four-time buyers over 52 weeks for 100 samples of households. 
write.table(result_cum_week,"C:/Users/bemmaor/Downloads/kiwibubbles/result_cum_week_Model_6_quart.txt")
# The following output file shows: (i) the total # of transactions, (ii) the % of nonbuyers, (iii) the % of one-time buyers, (iv) the % of two-time buyers, (v) the % of three-time buyers, and (vi) the % of four-time buyers over 52 weeks for 100 samples of households. 
write.table(result_simulation,"C:/Users/bemmaor/Downloads/kiwibubbles/result_simulation_Model_6_quart.txt")
# The following output file shows the # of one-time buyers, two-time buyers,..., 52-time buyers over a period of 52 weeks for 100 samples of households. One can easily obtain the # of nonbuyers for each sample. NA is equivalent to 0.   
write.table(result_transaction_times,"C:/Users/bemmaor/Downloads/kiwibubbles/result_transaction_times_Model_6_quart.txt")
# The following output file shows the mean and SD of (i) the total # of transactions, (ii) the % of nonbuyers, (iii) the % of one-time buyers, (iv) the % of two-time buyers, (v) the % of three-time buyers, and (vi) the % of four-time buyers over 52 weeks across 100 samples of households.The size of each sample equals 2,799.   
write.table(output,"C:/Users/bemmaor/Downloads/kiwibubbles/summary_Model_6_quart.txt")
# The following output file shows the weekly # of transactions over 52 weeks for 100 samples of households.  
write.table(transaction_per_week,"C:/Users/bemmaor/Downloads/kiwibubbles/transaction_Model_6_quart.txt")

#write.table(result_cum_week,"C:/Users/bemmaor/Downloads/kiwibubbles/result_cum_week_Model_3_quart.txt")
#write.table(result_simulation,"C:/Users/bemmaor/Downloads/kiwibubbles/result_simulation_Model_3_quart.txt")
#write.table(result_transaction_times,"C:/Users/bemmaor/Downloads/kiwibubbles/result_transaction_times_Model_3_quart.txt")
#write.table(output,"C:/Users/bemmaor/Downloads/kiwibubbles/summary_Model_3_quart.txt")
#write.table(transaction_per_week,"C:/Users/bemmaor/Downloads/kiwibubbles/transaction_Model_3_quart.txt")

#write.table(result_cum_week,"C:/Users/bemmaor/Downloads/kiwibubbles/result_cum_week_Model_4_quart.txt")
#write.table(result_simulation,"C:/Users/bemmaor/Downloads/kiwibubbles/result_simulation_Model_4_quart.txt")
#write.table(result_transaction_times,"C:/Users/bemmaor/Downloads/kiwibubbles/result_transaction_times_Model_4_quart.txt")
#write.table(output,"C:/Users/bemmaor/Downloads/kiwibubbles/summary_Model_4_quart.txt")
#write.table(transaction_per_week,"C:/Users/bemmaor/Downloads/kiwibubbles/transaction_Model_4_quart.txt")

#write.table(result_cum_week,"C:/Users/bemmaor/Downloads/kiwibubbles/result_cum_week_Model_5_quart.txt")
#write.table(result_simulation,"C:/Users/bemmaor/Downloads/kiwibubbles/result_simulation_Model_5_quart.txt")
#write.table(result_transaction_times,"C:/Users/bemmaor/Downloads/kiwibubbles/result_transaction_times_Model_5_quart.txt")
#write.table(output,"C:/Users/bemmaor/Downloads/kiwibubbles/summary_Model_5_quart.txt")
#write.table(transaction_per_week,"C:/Users/bemmaor/Downloads/kiwibubbles/transaction_Model_5_quart.txt")

#write.table(result_cum_week,"C:/Users/bemmaor/Downloads/kiwibubbles/result_cum_week_Model_6_12_weeks_quart.txt")
#write.table(result_simulation,"C:/Users/bemmaor/Downloads/kiwibubbles/result_simulation_Model_6_12_weeks_quart.txt")
#write.table(result_transaction_times,"C:/Users/bemmaor/Downloads/kiwibubbles/result_transaction_times_Model_6_12_weeks_quart.txt")
#write.table(output,"C:/Users/bemmaor/Downloads/kiwibubbles/summary_Model_6_12_weeks_quart.txt")
#write.table(transaction_per_week,"C:/Users/bemmaor/Downloads/kiwibubbles/transaction_Model_6_12_weeks_quart.txt")

#write.table(result_cum_week,"C:/Users/bemmaor/Downloads/kiwibubbles/result_cum_week_Model_6_20_weeks_quart.txt")
#write.table(result_simulation,"C:/Users/bemmaor/Downloads/kiwibubbles/result_simulation_Model_6_20_weeks_quart.txt")
#write.table(result_transaction_times,"C:/Users/bemmaor/Downloads/kiwibubbles/result_transaction_times_Model_6_20_weeks_quart.txt")
#write.table(output,"C:/Users/bemmaor/Downloads/kiwibubbles/summary_Model_6_20_weeks_quart.txt")
#write.table(transaction_per_week,"C:/Users/bemmaor/Downloads/kiwibubbles/transaction_Model_6_20_weeks_quart.txt")

#write.table(result_cum_week,"C:/Users/bemmaor/Downloads/kiwibubbles/result_cum_week_Model_2_quart.txt")
#write.table(result_simulation,"C:/Users/bemmaor/Downloads/kiwibubbles/result_simulation_Model_2_quart.txt")
#write.table(result_transaction_times,"C:/Users/bemmaor/Downloads/kiwibubbles/result_transaction_times_Model_2_quart.txt")
#write.table(output,"C:/Users/bemmaor/Downloads/kiwibubbles/summary_Model_2_quart.txt")
#write.table(transaction_per_week,"C:/Users/bemmaor/Downloads/kiwibubbles/transaction_Model_2_quart.txt")


#write.table(result_cum_week,"C:/Users/bemmaor/Downloads/kiwibubbles/result_cum_week_Model_1_quart.txt")
#write.table(result_simulation,"C:/Users/bemmaor/Downloads/kiwibubbles/result_simulation_Model_1_quart.txt")
#write.table(result_transaction_times,"C:/Users/bemmaor/Downloads/kiwibubbles/result_transaction_times_Model_1_quart.txt")
#write.table(output,"C:/Users/bemmaor/Downloads/kiwibubbles/summary_Model_1_quart.txt")
#write.table(transaction_per_week,"C:/Users/bemmaor/Downloads/kiwibubbles/transaction_Model_1_quart.txt")


