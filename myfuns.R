# functions for sentiment analysis NASDAQ

############fun for fitting panel data expectile model with random effects for diff tau


fit_eeffects<-function(formul,#formula as lm formula
                       dat,#dataframe
                       cluster=~Ticker + Date,#how to cluster the panel data
                       taus=seq(0.05,0.95,by=0.05),#expectile levels
                       iters=100#iterations
                       ){
  coefs0<-ses0<-NULL
  mod0<-list(0)
  r2<-NULL
  datt<-dat
  library(sandwich)
  environment(formul) <- environment()#for the lm function: env of formula

  for(i in 1:length(taus)){
    tau<-taus[i]
    w<-rep(0.5,nrow(datt))
    conv<-F
    iter<-0
    #combined sentiment
    datt$comb.tone<-tau*datt$LM.TONE+(1-tau)*datt$HIV4.TONE
    while(!conv){
      mod<-lm(formula=formul,data=datt,weights=w)
      wn<-ifelse(as.numeric(resid(mod))>0,tau,1-tau)
      if(all(wn==w)|iter>iters){
        conv<-T
      }
      iter<-iter+1
      w<-wn

    }

    mod<-lm(formula=formul,data=datt,weights=w)
    mod0[[i]]<-mod
    coefs0<-cbind(coefs0,coef(mod))
    ses0<-cbind(ses0, sqrt(diag(vcovPL(mod,cluster=cluster))))
    r2<-c(r2,summary(mod)$r.squared)
  }
  return(list(coefs=coefs0,ses=ses0, mods=mod0,r2=r2))
}


##########function for plotting the coefficients of the ERRE models with asymptotic $95\%$-CI for varying $\tau$.-->

plot_eefects<-function(mod,#fitted model
                       al=0.05,#alpha (CI error)
                       which_ones,#which coefs to plot
                       lbls=NA #labels fro the coefs on the plot
                       ){
  tau_axis<-seq(0.1,0.9)
  coefs<-mod$coefs[which_ones,]
  ses<-mod$ses[which_ones,]

  if(is.na(lbls[1])){
    lbls<-rownames(coefs)
  }
  par(mfrow=c(ceiling(length(which_ones)/2),2))
  for(j in 1:nrow(coefs)){
    plot(coefs[j,],type="l",main=lbls[j],xlab=bquote(tau), ylab=bquote(beta[tau]),col=4,lwd=2
         ,ylim=c(min(c(coefs[j,]- qnorm(1-al/2)*ses[j,],0)),max(c(coefs[j,]+ qnorm(1-al/2)*ses[j,]),0)),axes=F)
    axis(2)
    axis(1,at=seq(2,ncol(coefs),2),labels=seq(0.1,0.9,0.1))
    lines(coefs[j,] - qnorm(1-al/2)*ses[j,],lty=2,col=4)
    lines(coefs[j,] + qnorm(1-al/2)*ses[j,],lty=2,col=4)
    abline(h=0,lwd=1,col=2)
  }
}




##################make the coefficient table function-->

make_knitr_table<-function(mod, #fitted model
                           sent="lm",#sentiment measure used
                           which_tau=c(1,2,5,10,15,18,19), #taus to put in the table
                           caption="Coefficients and $R^2$ of the ERRE models with LM sentiment and different $\\tau$s",
                           taus=seq(0.05,0.95,0.05),
                           row.names=c("Constant","$\\log(size_{i,t})$","$\\alpha_{i,t}$","$turn_{i,t}$",
                                       "$nq_{i,t}$","$btm_{i,t-1}$","$news_{i,t}$",paste0("$",sent,"\\_tone_{i,t}$"),"$R^2$")
                           #,
                           #id="id=\"table1\""
){
  library(knitr)
  library(kableExtra)
  #https://stackoverflow.com/questions/53341155/coloring-rows-with-kableextra-based-on-cell-values
  #https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf

  #coefficient table
  options(scipen=999)
  tab<-round(mod$coefs[,which_tau],5)
  #pvals and stars
  pvals<-2*(1-pt(abs(mod$coefs/mod$ses)[,which_tau], mod$mods[[1]]$df.residual))
  stars<-matrix(rep("",nrow(mod$coefs)),dim(pvals)[1],dim(pvals)[2])
  stars[pvals<0.001]<-"***"; stars[pvals>0.001&pvals<0.01]<-"**";  stars[pvals>0.01&pvals<0.05]<-"*"
  #add stars
  tab[tab>0]<-paste0(" ",tab[tab>0])
  tab<-matrix(paste0(tab,stars,sep=""),dim(pvals)[1],dim(pvals)[2])
  tab<-rbind(tab,round(mod$r2[c(1,2,5,10,15,18,19)],5))
  rownames(tab)<-row.names

  coln<-paste0("$\\beta_{",taus[which_tau],"}$")

  knitr::kable(tab, align = "lllllll",
               col.names = coln,
               row.names = TRUE,
               #table.attr = "id=\"table1\"",
               digits = 5,
               caption = caption
  )%>%kable_styling() %>%
    row_spec(9, bold = T, color = "black", background = "gray")%>%
    row_spec(6:8, bold = T, color = "black", background = "yellow")
}

make_knitr_table_bm<-function(mod, #fitted model
                           sent="lm",#sentiment measure used
                           which_tau=c(1,2,5,10,15,18,19), #taus to put in the table
                           caption="Coefficients and $R^2$ of the ERRE models with LM sentiment and different $\\tau$s",
                           taus=seq(0.05,0.95,0.05),
                           row.names=c("Constant","$\\log(size_{i,t})$","$\\alpha_{i,t}$","$turn_{i,t}$",
                                       "$nq_{i,t}$","$btm_{i,t-1}$","$news_{i,t}$",paste0("$",sent,"\\_tone_{i,t}$"),"$R^2$")
                           #,
                           #id="id=\"table1\""
){
  library(knitr)
  library(kableExtra)
  #https://stackoverflow.com/questions/53341155/coloring-rows-with-kableextra-based-on-cell-values
  #https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf

  #coefficient table
  options(scipen=999)
  tab0<-round(mod$coefs[,which_tau],5)
  #pvals and stars
  pvals<-2*(1-pt(abs(mod$coefs/mod$ses)[,which_tau], mod$mods[[1]]$df.residual))
  stars<-matrix(rep("",nrow(mod$coefs)),dim(pvals)[1],dim(pvals)[2])
  stars[pvals<0.001]<-"***"; stars[pvals>0.001&pvals<0.01]<-"**";  stars[pvals>0.01&pvals<0.05]<-"*"
  #add stars
  tab<-tab0
  tab[tab>0]<-paste0(" ",tab[tab>0])
  tab<-matrix(paste0(tab,stars,sep=""),dim(pvals)[1],dim(pvals)[2])
  tab<-rbind(tab,round(mod$r2[c(1,2,5,10,15,18,19)],5))
  rownames(tab)<-row.names

  coln<-paste0("$\\beta_{",taus[which_tau],"}$")

  knitr::kable(tab, align = "lllllll", "latex",
               col.names = coln,
               row.names = TRUE,
               #table.attr = "id=\"table1\"",
               digits = 5,
               caption = caption,
               escape = FALSE
  )%>%kable_styling() %>%
    row_spec(9, bold = T, color = "black",
             background = "#DCDCDC"
               # cut(c(0,as.numeric(tab[9,])), c(0,0.00005,0.0001,0.001,0.01,0.1,10),
               #        #c("#F0F0F0", "#E8E8E8","#E0E0E0","#DCDCDC","#D3D3D3", "#C8C8C8")))%>%
               #     c("#666666", "#999999", "#BBBBBB","#666666", "#999999", "#BBBBBB"))
             )%>%
    #row_spec(6:8, bold = T, color = "black", background = "yellow")%>%
    row_spec(6, bold = T,color = ifelse(c(0,tab0[6,]) >= 0,  "blue","red"))%>%
    row_spec(7, bold = T,color = ifelse(c(0,tab0[7,]) >= 0,  "blue","red"))%>%
               #spec_color(x=c(0,tab0[7,]), alpha=0.1, begin=0.1,  end = 0.8,
               #                      option="C"))

    row_spec(8, bold = T,color = ifelse(c(0,tab0[8,]) >= 0,  "blue","red"))%>%
    kable_styling(latex_options = "scale_down")
}


