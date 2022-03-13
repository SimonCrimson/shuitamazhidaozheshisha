univariate_Cox<-function(db,gene) {
  expr<-get(paste(db,"_exp",sep=""))
  pati<-get(paste(db,"_patients",sep=""))
  
  patients_intersect<-intersect(colnames(expr),rownames(pati))
  
  expr<-expr[,patients_intersect]
  pati<-pati[patients_intersect,]
  
  if (gene %in% rownames(expr)==F) {
    print("该数据库无基因:")
    print(gene)
    return (NA)
  }
  
  pati$TEST_GENE<-as.vector(t(expr[gene,]))
  
  library(survival)
  

  mycox<-coxph(Surv(OS,Survival)~TEST_GENE,data=pati)
  mycox
  
  return (mycox)
}




multivariate_Cox<-function(db,gene,covariate) {
  expr<-get(paste(db,"_exp",sep=""))
  pati<-get(paste(db,"_patients",sep=""))
  
  patients_intersect<-intersect(colnames(expr),rownames(pati))
  
  expr<-expr[,patients_intersect]
  pati<-pati[patients_intersect,]
  
  if (gene %in% rownames(expr)==F) {
    print("该数据库无基因:")
    print(gene)
    return (NA)
  }
  
  pati$TEST_GENE<-as.vector(t(expr[gene,]))
  
  library(survival)
  
  fustr<-"Surv(OS,Survival)~TEST_GENE"
  
  for (v in covariate) {
    fustr<-paste(fustr,"+",v,sep="")
  }
  
  fu<-as.formula(fustr)
  
  
  mycox<-coxph(fu,data=pati)
  mycox
  
  return (mycox)
}


Surv_diff_median<-function(db,gene) {
  expr<-get(paste(db,"_exp",sep=""))
  pati<-get(paste(db,"_patients",sep=""))
  
  patients_intersect<-intersect(colnames(expr),rownames(pati))
  
  expr<-expr[,patients_intersect]
  pati<-pati[patients_intersect,]
  
  if (gene %in% rownames(expr)==F) {
    print("该数据库无基因:")
    print(gene)
    return (NA)
  }
  
  pati$TEST_GENE<-as.vector(t(expr[gene,]))
  
  su<-summary(pati$TEST_GENE)

  TEST_GENE_MEDIAN<-as.numeric(su['Median'])
  
  pati$TEST_GENE_GROUP<-"High"
  pati$TEST_GENE_GROUP[which(pati$TEST_GENE<=TEST_GENE_MEDIAN)]<-"Low"
  
  library(survival)
  
  mydiff<-survdiff(Surv(OS,Survival)~TEST_GENE_GROUP,data=pati)
  mydiff
  
  library(survminer)
  ggsurvplot(survfit(Surv(OS,Survival)~TEST_GENE_GROUP,data=pati),pval=T)
  
  return (mydiff)
}


Surv_diff_cutoffs<-function(db,gene,cutoffs) {
  expr<-get(paste(db,"_exp",sep=""))
  pati<-get(paste(db,"_patients",sep=""))
  
  patients_intersect<-intersect(colnames(expr),rownames(pati))
  
  expr<-expr[,patients_intersect]
  pati<-pati[patients_intersect,]
  
  if (gene %in% rownames(expr)==F) {
    print("该数据库无基因:")
    print(gene)
    return (NA)
  }
  
  pati$TEST_GENE<-as.vector(t(expr[gene,]))
  
  su<-summary(pati$TEST_GENE)
  
  i=1
  
  pati$TEST_GENE_GROUP<-"Group1"
  
  for (cu in cutoffs) {

    if (i>1) {
      pati$TEST_GENE_GROUP[which(pati$TEST_GENE>cu)]<-paste("Group",i,sep="")
    }
    i=i+1
  } 
  
  library(survival)
  
  mydiff<-survdiff(Surv(OS,Survival)~TEST_GENE_GROUP,data=pati)
  mydiff
  
  library(survminer)
  ggsurvplot(survfit(Surv(OS,Survival)~TEST_GENE_GROUP,data=pati),pval=T)
  
  return (mydiff)
}

Surv_diff_best_cutoffs<-function(db,gene) {
  expr<-get(paste(db,"_exp",sep=""))
  pati<-get(paste(db,"_patients",sep=""))
  
  patients_intersect<-intersect(colnames(expr),rownames(pati))
  
  expr<-expr[,patients_intersect]
  pati<-pati[patients_intersect,]
  
  if (gene %in% rownames(expr)==F) {
    print("该数据库无基因:")
    print(gene)
    return (NA)
  }
  
  pati$TEST_GENE<-as.vector(t(expr[gene,]))
  
  sc<-surv_cutpoint(data=pati,time="OS",event="Survival",variables = c("TEST_GENE"))
  cu<-sc$cutpoint$cutpoint[1]
  
  print("BEST CUTOFF:")
  print(cu)
  
  pati$TEST_GENE_GROUP<-"High"
  pati$TEST_GENE_GROUP[which(pati$TEST_GENE<=cu)]<-"Low"
  
  library(survival)
  
  mydiff<-survdiff(Surv(OS,Survival)~TEST_GENE_GROUP,data=pati)
  mydiff
  
  library(survminer)
  ggsurvplot(survfit(Surv(OS,Survival)~TEST_GENE_GROUP,data=pati),pval=T)
  
  return (mydiff)
}
