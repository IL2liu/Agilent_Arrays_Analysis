###################################################
##extracts correction factors from sliding window##
###################################################

local.norm<-function(f1,data,win=16)
{
 out=vector()
  #f1=data
  mr=max(data[,1])
  r=seq(f1[1]-win,f1[1]+win,1)
  r[which(r<1)]=abs(r[which(r<1)])+max(r)+1
  r[which(r>mr)]=min(r)-abs(r[which(r>mr)]-mr)
  r=r[order(r)]
  r=r[-which(r==f1[1])]
  mc=max(data[,2])
  c=seq(f1[2]-win,f1[2]+win,1)
  c[which(c<1)]=abs(c[which(c<1)])+max(c)+1
  c[which(c>mc)]=min(c)-abs(c[which(c>mc)]-mc)
  c=c[order(c)]
  c=c[-which(c==f1[2])]
  
  set=which((data[,1] %in% r) & (data[,2] %in% c))
  if (length(set)==0)
  {
   out[1]=NA
  }
  else
  {
   out[1]=median(data[set,3],na.rm=TRUE)
  }
  return(out)
}

local.norm.spikes<-function(f1,data,win=16)
{
 out=vector()
 #out1=vector()
  #f1=data
  mr=max(data[,1])
  #print(row.names(data))
  r=seq(f1[1]-win,f1[1]+win,1)
  r[which(r<1)]=abs(r[which(r<1)])+max(r)+1
  r[which(r>mr)]=min(r)-abs(r[which(r>mr)]-mr)
  r=r[order(r)]
  r=r[-which(r==f1[1])]
  mc=max(data[,2])
  c=seq(f1[2]-win,f1[2]+win,1)
  c[which(c<1)]=abs(c[which(c<1)])+max(c)+1
  c[which(c>mc)]=min(c)-abs(c[which(c>mc)]-mc)
  c=c[order(c)]
  c=c[-which(c==f1[2])]
  
  set=which((data[,1] %in% r) & (data[,2] %in% c) & (grepl('EC', row.names(data))==TRUE))
  #cat(length(set),"\n")
  if (length(set)==0)
  {
   out[1]=1
  }
  else
  {
   out[1]=median(data[set,3],na.rm=TRUE)
  }
  out[1]=paste(out[1],length(set),sep=';') 
  return(out)
}


####################################################
##old non vectorised script#########################
####################################################

#local.norm.old<-function(data, win=16,type=c("RG","MA"))
#{
# out=data.frame()
# step=1000
# for(i in 1:nrow(data))
# {
#  f1=data[i,]
#  mr=max(data$genes$Row)
#  r=seq(f1$genes$Row-win,f1$genes$Row+win,1)
#  r[which(r<1)]=abs(r[which(r<1)])+max(r)+1
#  r[which(r>mr)]=min(r)-abs(r[which(r>mr)]-mr)
#  r=r[order(r)]
#  mc=max(data$genes$Column)
#  c=seq(f1$genes$Column-win,f1$genes$Column+win,1)
#  c[which(c<1)]=abs(c[which(c<1)])+max(c)+1
#  c[which(c>mc)]=min(c)-abs(c[which(c>mc)]-mc)
#  c=c[order(c)]
#  
#  set=which((data$genes$Row %in% r) & (data$genes$Column %in% c))
#  set=set[-which(set==i)] 
#  if (length(set)==0)
#  {
#   out[i,1]=NA
#   out[i,2]=0
#  }
#  else if (type=="RG")
#  { 
#   out[i,1]=median(data$other[[4]][set],na.rm=TRUE)
#   out[i,2]=length(which(is.na(data$other[[4]][set])==TRUE))
#  }
#  else if (type=="MA")
#  {
#   out[i,1]=median(10^(data$M[set]),na.rm=TRUE)
#   out[i,2]=length(which(is.na(data$M[set])==TRUE))
#  }
# if(i>step)
# {
# cat(paste(i,"\n"))
# step=i+1000
# }
# }
#  
# if (type=="RG")
# {
#  out[,1]=data$other[[4]]/out[,1]
# }
# else if (type=="MA")
# {
#  out[,1]=(10^(data$M))/out[,1]
# } 
# 
# return(out)
# cat("\n")
#}
######################################################################
######################################################################
######################################################################

######################################################################
############for batch analysis, not adapted to normalisation yet######
######################################################################
arrayQCbatch<-function(directory=getwd(),tag=NULL, SD_cutoff_low=80, SD_cutoff_high=95, spikes=FALSE, nc=TRUE, kem=TRUE, tiling=FALSE, print.gpr=FALSE, local=FALSE)
{
 
 cat("Processing gpr files in ",directory,"\n\n")
 directory.hold=getwd()
 setwd(directory)
 toDo=dir(pattern='\\.gpr$')
 if (length(grep('processed\\.gpr', toDo)) !=0) 
 {
  toDo=toDo[-grep('processed\\.gpr', toDo)]
 }
 
 for (i in 1:length(toDo))
 {
  cat("Processing: ",toDo[i],sep='') 
  arrayQC(file=toDo[i], SD_cutoff_low=SD_cutoff_low, SD_cutoff_high=SD_cutoff_high, spikes=spikes, nc=nc, kem=kem, tiling=tiling, verbose=FALSE, print.file=TRUE, print.log=TRUE, print.gpr=print.gpr, local=local)
  cat("Finished.\n")
 }

 setwd(directory.hold)
}

########################################################################
##########for filtering and normalisation###############################
########################################################################

arrayQC<-function(file, SD_cutoff_low=80, SD_cutoff_high=95, tag=NULL, spikes=FALSE, nc=TRUE, kem=TRUE, tiling=FALSE, verbose=TRUE, print.file=TRUE, print.log=TRUE, print.gpr=FALSE, return.table=FALSE, background="subtract", global="none",local="none", window=15)
{

 err=try(library(limma),TRUE)

 if (class(err)=="try-error")
 {
  cat(
  "\n",
  "Please install the package Limma on your computer:","\n",
  sep='')
  return(0)
 }

 err=try(library(marray),TRUE)

 if (class(err)=="try-error")
 {
  cat(
  "\n",
  "Please install the package marray on your computer:","\n",
  sep='')
  return(0)
 }

 wd=getwd()
 cat("\n",sep='')

#####Loads data as an RG object and as a table. Loads gpr header. Records ratio formulation#####################

 info=try(read.delim(file, nrow=50, stringsAsFactor=F,header=F),TRUE)

 if (class(info)=="try-error")
 {
  cat(
  "########","\n",
  "Can't find file: ",file,"\n", 
  "Please check spelling and make sure that the file you want to analyse in the directory bellow:","\n",
  wd,"\n",
  "then, give it another go.","\n",
  "########","\n",
  sep='')
  return(0)
 }

 ratio.line=grep('RatioFormulations', info[,1])

 block=grep("Block",info[,1])-1
 info=info[1:block,]

 if (grepl('\\(532\\/635\\)',info[ratio.line,1]))
 {
  ratio.form="(532/635)"
  exp="GREEN"
  ctr="RED"
 }

 else if (grepl('\\(635\\/532\\)',info[ratio.line,1]))
 {
  ratio.form="(635/532)"
  exp="RED"
  ctr="GREEN"
 }

 ratio.form=paste("Median of Ratios",ratio.form)

 gpr=read.maimages(file, source="genepix", other.columns=c("Flags","% > B635+2SD2","% > B532+2SD2",ratio.form,"X","Y","Normalize","F532 % Sat.","F635 % Sat."), verbose=FALSE)

 out.file=sub('\\.gpr','',file)
 if(!is.null(tag))
 {
  out.file=paste(out.file,tag,sep='.')
 }

#####Records initial Flags status#############################

 flags.bad=length(which(gpr$other$Flags==-100))
 flags.notfound=length(which(gpr$other$Flags==-50))
 flags.before=flags.bad+flags.notfound

#####Change Flags#####################################

 gpr$other$Flags[gpr$other$Flags==-50]=-100
 
 kem.probes="Included"
 nc.probes="Included"
 spikes.probes="Included"

 if (kem==FALSE)
 {
  gpr$other$Flags[grep('PK_', gpr$genes$ID)]=-100  # Only present in 4x44k design, from Holstege's lab
  kem.probes="Flaged bad"
 }

 if (nc==FALSE)
 {
  gpr$other$Flags[grep('I$', gpr$genes$Name)]=-100        # remove intronic probes - not active at the moment
  gpr$other$Flags[grep('_AS', gpr$genes$Name)]=-100       # remove antisense probes - not active at the moment
  nc.probes="Flaged bad"
 }

 if (spikes==FALSE)
 {
  gpr$other$Flags[grep('EC', gpr$genes$Name)]=-100        #Spikes
  gpr$other$Flags[grep('EC', gpr$genes$ID)]=-100	  #Spikes
  spikes.probes="Flaged bad"
 }

 if (tiling==TRUE)
 {
  kem.probes=NA
  nc.probes=NA
  spikes.probes=NA
 }
 
 if (spikes==TRUE)
 {
  gpr$genes$Name[grep('EC', gpr$genes$ID)]=gpr$genes$ID[grep('EC', gpr$genes$ID)]	        #Spikes have no Name on the 8x15k array -> get them back.
  gpr$R[grep('EC', gpr$genes$ID)]=2000	        #Spikes have no Name on the 8x15k array -> get them back.
  gpr$G[grep('EC', gpr$genes$ID)]=2000	        #Spikes have no Name on the 8x15k array -> get them back.
  gpr$other[[4]][grep('EC', gpr$genes$ID)]=1	        #Spikes have no Name on the 8x15k array -> get them back.
 } 

 gpr$other$Flags[grep('NULL', gpr$genes$Name)]=-100   	   	# Bad spots from first generation of 4x44k arrays 
 gpr$other$Flags[grep('\\(\\+\\)',gpr$genes$ID)]=-100 		# control flags from Agilent arrays
 gpr$other$Flags[grep('\\(\\-\\)',gpr$genes$ID)]=-100 		# control flags from Agilent arrays
 gpr$other$Flags[grep('DarkCorner',gpr$genes$Name)]=-100 	# control flags from Agilent arrays
 gpr$other$Flags[grep('NegativeControl',gpr$genes$Name)]=-100 	# control flags from Agilent arrays
 gpr$other$Flags[grep('BrightCorner',gpr$genes$Name)]=-100 	# control flags from Agilent arrays
 gpr$other$Flags[grep('GE',gpr$genes$Name)]=-100 		# control flags from Agilent arrays
 gpr$other$Flags[grep('DCP',gpr$genes$Name)]=-100 		# control flags from Agilent arrays
 gpr$other$Flags[grep('ETG',gpr$genes$Name)]=-100 		# control flags from Agilent arrays
 gpr$other$Flags[grep('EQC',gpr$genes$Name)]=-100 		# control flags from Agilent arrays
 gpr$other$Flags[grep('RC',gpr$genes$Name)]=-100 		# control flags from Agilent arrays
 gpr$other$Flags[which(is.na(gpr$genes$Name))]=-100 		# NA...
 gpr$other$Flags[which(gpr$genes$Name=='')]=-100 		# empty names
 gpr$other$Flags[which(gpr$genes$Name=='RC')]=-100 		# empty names

 
# The probes below do not behave in self-self experiments - given -50 flags

 bad=c("SPAC18G6.09c","SNORNA.03","SNORNA.17","SNORNA.20","SNORNA.27","SNORNA.34","SNORNA.43","SNORNA.44","SNORNA.48","SNORNA.52","SPNCRNA.213","SPNCRNA.214","SPNCRNA.475_AS","SPNCRNA.450_AS","SPCC63.03","SPCC63.05","SPAC22F3.05c","SPAC8C9.11","SPAC1486.05I","SPAC6F6.03cI")

 gpr$other$Flags[which((gpr$genes$Name %in% bad)&(gpr$other$Flags==0))]=-50
 gpr$other$Flags[which((grepl('SPMIT', gpr$genes$Name))&(gpr$other$Flags==0))]=-50 
 gpr$other$Flags[which((grepl('RRNA', gpr$genes$Name))&(gpr$other$Flags==0))]=-50 

#flags spots with < SD_cutoff_low % pixels above 2SD of background in one or both channels
 
 SD_filtered=length(which(((gpr$other[[2]] < SD_cutoff_low) | (gpr$other[[3]] < SD_cutoff_low)) & ((gpr$other[[2]] < SD_cutoff_high)&(gpr$other[[3]] < SD_cutoff_high)) & (gpr$other$Flags==0)))
 
 gpr$other$Flags[which(((gpr$other[[2]] < SD_cutoff_low) | (gpr$other[[3]] < SD_cutoff_low)) & ((gpr$other[[2]] < SD_cutoff_high)&(gpr$other[[3]] < SD_cutoff_high)))]=-100

#falgs spots with saturated pixels################################

saturated=length(which((gpr$other$"F532 % Sat." != 0)|(gpr$other$"F635 % Sat." != 0)))
gpr$other$Flags[which((gpr$other$"F532 % Sat." != 0)|(gpr$other$"F635 % Sat." != 0))]=-100


######Create filtered RG object changing flaged probes to NA######

 gpr.filtered=gpr
 gpr.filtered$R[which(gpr.filtered$other$Flags != 0)]=NA
 gpr.filtered$G[which(gpr.filtered$other$Flags != 0)]=NA
 gpr.filtered$Rb[which(gpr.filtered$other$Flags != 0)]=NA
 gpr.filtered$Gb[which(gpr.filtered$other$Flags != 0)]=NA
 gpr.filtered$other[[2]][which(gpr.filtered$other$Flags != 0)]=NA
 gpr.filtered$other[[3]][which(gpr.filtered$other$Flags != 0)]=NA
 gpr.filtered$other[[4]][which(gpr.filtered$other$Flags != 0)]=NA



#####################################################################################################################
######Apply Normalisation############################################################################################
#####################################################################################################################



###initialise for qc plots###########################################
##
MA.global=normalizeWithinArrays(gpr.filtered, method="none",bc.method=background)
MA.global.only=MA.global
gpr.filtered.only=gpr.filtered
norm.ratios=2^MA.global.only$M
gpr.filtered.only$other$Normalize=norm.ratios
###

###
###global normalisation, done on signals######################################################################
###

if (global !="none")
{
 cat(paste("\ndoing ",'"',global,'"'," global normalisation with ",'"',background,'"'," as a background correction method...\n",sep=''),sep='')
 MA.global=normalizeWithinArrays(gpr.filtered, method=global,bc.method=background)
 norm.ratios=2^MA.global$M
 gpr.filtered$other$Normalize=norm.ratios

 ##keeps data with global normalisatio only for qc plots#####################################
 ##
 MA.global.only=MA.global
 gpr.filtered.only$other$Normalize=norm.ratios
 ##

 ##Adds local normalisation on top of global normalisation################################### 
 if (local=="inHouse")
 {
  cat("\ndoing in house local normalisation on normalised data...\n",sep='')
  data.local=cbind(gpr.filtered.only$genes$Row,gpr.filtered.only$genes$Column,(2^(MA.global$M)))
  rn=paste(gpr.filtered.only$genes$Name,sample(x=seq(1,10000,1)),sep='.')
  row.names(data.local)=rn
  #data.local=cbind(gpr.filtered.only$genes$Row,gpr.filtered.only$genes$Column,(2^(MA.global$M)))
  if (spikes==TRUE)
  {
  norm.fac=apply(data.local,1,local.norm.spikes,data=data.local, win=window)
  }
  else
  {
  norm.fac=apply(data.local,1,local.norm,data=data.local, win=window)
  }
####
return(norm.fac)
####
  gpr.filtered$other$Normalize=gpr.filtered$other$Normalize/norm.fac
  MA.global$M=log2(gpr.filtered$other$Normalize)
 }
 
 ##doesn't allow ma2D if previous global normalisation has been applied######################
 if (local=="ma2D")
 {
  stop("\nthis combination is not implemented...\n",sep='')
 }

}

###
###local normalisation without global normalisation, done on ratios imported from GenePix##########################################
###
else if ((global=="none") & (local=="inHouse"))
{
 cat("\ndoing in house local normalisation on raw data...\n",sep='')
 data.local=cbind(gpr.filtered$genes$Row,gpr.filtered$genes$Column,gpr.filtered$other[[4]])
 norm.fac=apply(data.local,1,local.norm,data=data.local, win=window)
 gpr.filtered$other$Normalize=gpr.filtered$other[[4]]/norm.fac
 MA.global$M=log2(gpr.filtered$other$Normalize) 
}

else if ((global=="none") & (local=="ma2D"))
{
 cat("\ndoing ma2D local normalisation on raw data...\n",sep='')
 layout=c(gpr$printer$ngrid.r,gpr$printer$ngrid.c,gpr$printer$nspot.r,gpr$printer$nspot.c)
 cor.fac=ma2D(x=gpr.filtered$genes$Row, y=gpr.filtered$genes$Column, z=log2(gpr.filtered$other[[4]]), g=layout)
 gpr.filtered$other$Normalize=2^(log2(gpr.filtered$other[[4]])-cor.fac)
 MA.global$M=log2(gpr.filtered$other[[4]])-cor.fac
}

else if ((global=="none") & (local=="none"))
{
 cat("\nno normalisation step added...\n",sep='')
 norm.ratios=rep(NA,nrow(gpr.filtered))
 gpr.filtered$other$Normalize=norm.ratios
 MA.global$M=norm.ratios
}
  

print("done with normalisation")

###################################################################################################
###################################################################################################
###################################################################################################



######exports filtered data##########################################################################

 filtered=data.frame(gpr$genes[4],gpr$genes[5], (gpr$R-gpr$Rb),(gpr$G-gpr$Rb),gpr$other[[4]],gpr.filtered$other$Normalize,gpr$other$Flags)
 colnames(filtered)=c("ID","Name","F635-B635","F532-B532",ratio.form,"Normalised","Flags")

 filtered=filtered[which(filtered$Flags==0),] 

 write.table(filtered, file=paste(out.file,"filtered",sep='.'),quote=F,row.names=FALSE,col.names=TRUE,sep="\t")
 
######qc plots#######################################################################################


 la=getLayout(gpr$gene)
 la2=getLayout(gpr.filtered$gene)
 g=which(gpr$other$Flags==0)
 norm.all=gpr$other[[4]]/median(gpr$other[[4]][g],na.rm=TRUE)
 gpr.median=normalizeWithinArrays(gpr.filtered, method="median")
 norm.good=2^(gpr.median$M)
 R=gpr$R-gpr$Rb
 R[which(R <= 0)]=0.25
 G=gpr$G-gpr$Gb
 G[which(G <= 0)]=0.25


########
 if (verbose==TRUE)
 {
  if((global != "none") | (local == "inHouse") | (local == "ma2D"))
  {
   par(mfcol=c(3,4))
  }
  else
  {
   par(mfrow=c(2,3))
  }
#
  plot(log2(R),log2(G),xlim=c(-2,20),ylim=c(-2,20),pch=20,col="grey",xlab="log2 R-Rb",ylab="log2 G-Gb", main="log2 signals plot not centered")
  lines(log2(R[g]),log2(G[g]),type="p", pch=20, col="orange")
  legend(x=0,y=18,legend=c("raw","filtered"),fill=c("grey","orange"), bty="n")
  abline(log2(2),1)
  abline(-log2(2),1)
  abline(0,1)
#
  imageplot(log2(gpr.filtered$other[[4]]),la2,zlim=c(-log2(2),log2(2)), main="Spatial distribution of filtered log2 ratios")
#
  imageplot(gpr$other$Flags,la,zlim=c(-1,1), main="Spatial distribution of flags after filtering")
#
  plotMA(gpr.median,main="MA plot of filtered centered data",ylim=c(-4,4))
  abline(h=c(-log2(2),log2(2)),col="red", lwd=1.5)
#
  plot(gpr$other$X,log2(norm.all),ylim=c(-4,4), pch=20,cex=0.7,col="grey", xlab="X position on the array", ylab="Centered median of ratios")
  lines(gpr$other$X,log2(norm.good),pch=20,cex=0.7,col="blue",type="p")
  abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
#
  plot(gpr$other$Y,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="Y position on the array", ylab="Centered median of ratios")
  lines(gpr$other$Y,log2(norm.good),pch=20,cex=0.7,col="blue",type="p")
  abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
#  
  
  if((global != "none") | (local == "inHouse") | (local == "ma2D"))
  {
   plotMA(MA.global.only,main="MA plot after global Nomalisation", ylim=c(-4,4))
   abline(h=c(log2(2),-log2(2)), col="red",lwd=1.5)
   plot(gpr$other$X,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="X position on the array", ylab="Normalised log2 ratios")
   lines(gpr.filtered.only$other$X,log2(gpr.filtered.only$other$Normalize),pch=20,cex=0.7,col="blue", type="p")
   abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
   plot(gpr$other$Y,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="Y position on the array", ylab="Normalised log2 ratios")
   lines(gpr.filtered.only$other$Y,log2(gpr.filtered.only$other$Normalize),pch=20,cex=0.7,col="blue", type="p")
   abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
   
   plotMA(MA.global,main="MA plot after local Normalisation", ylim=c(-4,4))
   abline(h=c(log2(2),-log2(2)), col="red",lwd=1.5)
   plot(gpr$other$X,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="X position on the array", ylab="Normalised log2 ratios")
   lines(gpr.filtered$other$X,log2(gpr.filtered$other$Normalize),pch=20,cex=0.7,col="blue", type="p")
   abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
   plot(gpr$other$Y,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="Y position on the array", ylab="Normalised log2 ratios")
   lines(gpr.filtered$other$Y,log2(gpr.filtered$other$Normalize),pch=20,cex=0.7,col="blue", type="p")
   abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
  }


 }
######Compute gene names stats on good genes ########################################################

 if (tiling==FALSE)
 {
  #name=unique(gpr$genes$Name[which(gpr$other$Flags != -100)])
  
  if(spikes==TRUE)
  {
  name=unique(gpr$genes$Name[which((gpr$other$Flags == 0) & (grepl('EC', gpr$genes$Name)==FALSE))])
  }
  else
  {
  name=unique(gpr$genes$Name[which(gpr$other$Flags == 0)])
  }
  if (verbose==TRUE)
  {
   cat("\n",length(name)," genes to process, be patient it will take a little while","\n", sep='')
  }
  
  ratio=list()
  ratio.good=list()
  all=vector()
  good=vector()

  for(i in 1:length(name))
  {
   if ((global !="none")|(local !="none"))
   {
   ratio[[i]]=gpr.filtered$other$Normalize[which(gpr.filtered$genes$Name==name[i])]
   ratio.good[[i]]=gpr.filtered$other$Normalize[which((gpr.filtered$genes$Name==name[i]) & (gpr.filtered$other$Flags==0))]
   }
   else
   {
   ratio[[i]]=gpr.filtered$other[[4]][which(gpr.filtered$genes$Name==name[i])]
   ratio.good[[i]]=gpr.filtered$other[[4]][which((gpr.filtered$genes$Name==name[i]) & (gpr.filtered$other$Flags==0))]
   }
   all[i]=length(ratio[[i]])
   good[i]=length(ratio.good[[i]])
   
   if(length(ratio.good[[i]])==0)
   {
    ratio.good[[i]]=0
   }
######
#if(length(ratio[[i]]) > 2)
#{
#print(name[i])
#}
######
   names(ratio)[i]=name[i]
   names(ratio.good)[i]=name[i]
  }


  s=sapply(ratio,sd,na.rm=TRUE)
  m=sapply(ratio,mean,na.rm=TRUE)
  h=sapply(ratio,max,na.rm=TRUE)
  l=sapply(ratio,min,na.rm=TRUE)
  cv=s/m
  md=h/l
  s.good=sapply(ratio.good,sd,na.rm=TRUE)
  m.good=sapply(ratio.good,mean,na.rm=TRUE)
  h.good=sapply(ratio.good,max,na.rm=TRUE)
  l.good=sapply(ratio.good,min,na.rm=TRUE)
  cv.good=s.good/m.good
  md.good=h.good/l.good
  m.good[which(m.good==0)]=NA
  md.good[which(md.good=="NaN")]=NA
  md[which(all==1)]=NA
  md.good[which(good==1)]=NA


  out=data.frame(round(all,2),round(good,2),round(m,2),round(md,2),round(cv,2),round(m.good,2),round(md.good,2),round(cv.good,2),stringsAsFactors=FALSE)
  row.names(out)=name
  colnames(out)=c("all_probes","good_probes","mean_all","max_diff_all","cv_all","mean_good","max_diff_good","cv_good")


########Compute stats############################

  cv30=length(which(out[,"cv_good"] > 0.3))
  md2=length(which(out[,"max_diff_good"] > 2))
  max.probe=max(out[,"all_probes"])
  min.probe=min(out[,"all_probes"])
 }

 flags.filter=length(which(gpr$other$Flags < 0))-flags.before

########Display stats############################
 if (verbose==TRUE)
 {
  cat(
  "\n",
  "File analysed: ",file,"\n",
  info[grep('Gal',info[,1]),1],"\n",
  info[grep('Pixel',info[,1]),1],"\n",
  info[grep('Barcode',info[,1]),1],"\n",
  info[grep('ArrayName',info[,1]),1],"\n",
  "\n",
  "The experiment is in ",exp,"\n",
  "The control is in ",ctr,"\n",
  "\n",
  "Probes from the Holstege lab: ",kem.probes,"\n",
  "Probes against Introns and antisense: ",nc.probes,"\n",
  "Probes against Control spikes: ", spikes.probes,"\n",
  "\n",
  "Number of spots: ",nrow(gpr$genes),"\n",
  "Number of spots with no names: ",length(which(gpr$genes$Name=='')),"\n",
  "Number of bad spots: ",flags.bad,"\n",
  "Number of spots not found: ",flags.notfound,"\n",
  "Number of spots flaged before applying filter: ",flags.before,"\n",
  "Number of spots flaged during filtering: ",flags.filter,"\n",
  "Number of spots with saturated pixels in either channel (discarded): ",saturated," (",round((saturated/nrow(gpr$genes)*100),2),"% of total)","\n",
  "Number of spots discarded only bacause they had less than ",SD_cutoff_low,"% of pixels > 2SD of background: ",SD_filtered,"\n",
  "Final number of flaged spots: ", sum(flags.before,flags.filter)," (",(round((sum(flags.before,flags.filter)/nrow(gpr$genes)*100),1)),"%)","\n",
  "\n",
  sep=''
  )

  if (tiling==FALSE)
  {
   cat(
   "Number of genes with 0 good probe(s): ",length(which(out[,"good_probes"]==0)),"\n",
   sep=''
   )

   for (i in min.probe:max.probe)
   {
    cat(
    "Number of genes with ",i," good probe(s): ",length(which(out[,"good_probes"]==i)),"\n",
    sep='' 
    )
   }

   cat(
   "\n",
   "Number of genes with a cv of probes ratios over 30%: ",cv30,"\n",
   "Number of genes with a maximum ratio fold change between probes over 2: ", md2,"\n","\n",
   sep=''
   )
  }
 }
##########Create output files if requested###########################

 if (print.file==TRUE)
 {
  png(file=paste(out.file,"QC.png",sep='.'), height=210, width=297, unit="mm", res=300)

########

  if((global != "none") | (local == "inHouse") | (local == "ma2D"))
  {
   par(mfcol=c(3,4))
  }
  else
  {
   par(mfrow=c(2,3))
  }
#
  plot(log2(R),log2(G),xlim=c(-2,20),ylim=c(-2,20),pch=20,col="grey",xlab="log2 R-Rb",ylab="log2 G-Gb", main="log2 signals plot not centered")
  lines(log2(R[g]),log2(G[g]),type="p", pch=20, col="orange")
  legend(x=0,y=18,legend=c("raw","filtered"),fill=c("grey","orange"), bty="n")
  abline(log2(2),1)
  abline(-log2(2),1)
  abline(0,1)
#
  imageplot(log2(gpr.filtered$other[[4]]),la2,zlim=c(-log2(2),log2(2)), main="Spatial distribution of filtered log2 ratios")
#
  imageplot(gpr$other$Flags,la,zlim=c(-50,50), main="Spatial distribution of flags after filtering")
#
  plotMA(gpr.median,main="MA plot of filtered centered data",ylim=c(-4,4))
  abline(h=c(-log2(2),log2(2)),col="red", lwd=1.5)
#
  plot(gpr$other$X,log2(norm.all),ylim=c(-4,4), pch=20,cex=0.7,col="grey", xlab="X position on the array", ylab="Centered median of ratios")
  lines(gpr$other$X,log2(norm.good),pch=20,cex=0.7,col="blue",type="p")
  abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
#
  plot(gpr$other$Y,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="Y position on the array", ylab="Centered median of ratios")
  lines(gpr$other$Y,log2(norm.good),pch=20,cex=0.7,col="blue",type="p")
  abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
#  
  
  if((global != "none") | (local == "inHouse") | (local == "ma2D"))
  {
   plotMA(MA.global.only,main="MA plot after global Nomalisation", ylim=c(-4,4))
   abline(h=c(log2(2),-log2(2)), col="red",lwd=1.5)
   plot(gpr$other$X,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="X position on the array", ylab="Normalised log2 ratios")
   lines(gpr.filtered.only$other$X,log2(gpr.filtered.only$other$Normalize),pch=20,cex=0.7,col="blue", type="p")
   abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
   plot(gpr$other$Y,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="Y position on the array", ylab="Normalised log2 ratios")
   lines(gpr.filtered.only$other$Y,log2(gpr.filtered.only$other$Normalize),pch=20,cex=0.7,col="blue", type="p")
   abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
   
   plotMA(MA.global,main="MA plot after local Normalisation", ylim=c(-4,4))
   abline(h=c(log2(2),-log2(2)), col="red",lwd=1.5)
   plot(gpr$other$X,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="X position on the array", ylab="Normalised log2 ratios")
   lines(gpr.filtered$other$X,log2(gpr.filtered$other$Normalize),pch=20,cex=0.7,col="blue", type="p")
   abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
   plot(gpr$other$Y,log2(norm.all),ylim=c(-4,4),pch=20,cex=0.7,col="grey", xlab="Y position on the array", ylab="Normalised log2 ratios")
   lines(gpr.filtered$other$Y,log2(gpr.filtered$other$Normalize),pch=20,cex=0.7,col="blue", type="p")
   abline(h=c(-log2(2),log2(2)),col="red", lwd=2)
  }
########
  dev.off()
 }

 if (print.log==TRUE)
 {
  cat(
  "\n",
  "File analysed: ",file,"\n",
  info[grep('Gal',info[,1]),1],"\n",
  info[grep('Pixel',info[,1]),1],"\n",
  info[grep('Barcode',info[,1]),1],"\n",
  info[grep('ArrayName',info[,1]),1],"\n",
  "\n",
  "The experiment is in ",exp,"\n",
  "The control is in ",ctr,"\n",
  "\n",
  "Number of spots: ",nrow(gpr$genes),"\n",
  "Number of spots with no names: ",length(which(gpr$genes$Name=='')),"\n",
  "Number of bad spots: ",flags.bad,"\n",
  "Number of spots not found: ",flags.notfound,"\n",
  "Number of spots flaged before applying filter: ",flags.before,"\n",
  "Number of spots flaged during filtering: ",flags.filter,"\n",
  "Number of spots with saturated pixels in either channel (discarded): ",saturated," (",round((saturated/nrow(gpr$genes)*100),2),"% of total)","\n",
  "Number of spots discarded only bacause they had less than ",SD_cutoff_low,"% of pixels > 2SD of background: ",SD_filtered,"\n",
  "Final number of flaged spots: ", sum(flags.before,flags.filter)," (",(round((sum(flags.before,flags.filter)/nrow(gpr$genes)*100),1)),"%)","\n",
  "\n",
  file=paste(out.file,".log",sep=''),
  sep='')

  if (tiling==FALSE)
  {
   cat(
   "Number of genes with 0 good probe(s): ",length(which(out[,"good_probes"]==0)),"\n",
   file=paste(out.file,".log",sep=''),
   append=TRUE,
   sep=''
   )

   for (i in min.probe:max.probe) 
   {
    cat(
    "Number of genes with ",i," good probe(s): ",length(which(out[,"good_probes"]==i)),"\n",
    file=paste(out.file,".log",sep=''),
    append=TRUE,
    sep=''
    )
   }

   cat(
   "\n",
   "Number of genes with a cv of probes ratios over 30%: ",cv30,"\n",
   "Number of genes with a maximum ratio fold change between probes over 2: ", md2,"\n","\n",
   file=paste(out.file,".log",sep=''),
   append=TRUE,
   sep=''
   )
  }
 }

 if (print.gpr==TRUE)
 { 
  out.file=paste(out.file,'_processed.gpr',collapse='',sep='')
  skip <-  grep("Row", readLines(file, n=100)) - 1
  info=readLines(file,n=skip+1)
  write.table (info,file=out.file, col.names=F,row.names=F, sep='\t',quote=F)
  gpr.table=read.delim(file=file,skip=skip,stringsAsFactors=FALSE,header=T,check.names=F)
  gpr.table[,"Flags"]=gpr$other$Flags
  gpr.table[,"Normalize"]=gpr.filtered$other$Normalize
  write.table (gpr.table,file=out.file, col.names=F,row.names=F, sep='\t',quote=F,append=T)
 }

 gc()

 if (return.table==TRUE)
 {
  return(out)
 }

 else
 {
  return("finished!")
  #return(gpr.filtered)
 }

}




































