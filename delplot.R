#https://rdrr.io/github/PredictiveEcology/Require/src/R/setLibPaths.R (add MitoSAlt bin to paths)
setLibPaths <- function(libPaths, standAlone = TRUE) {
  oldLibPaths <- .libPaths()
  #libPaths <- checkPath(normPath(libPaths), create = TRUE)#, mustWork = TRUE)

  shim_fun <- .libPaths
  shim_env <- new.env(parent = environment(shim_fun))
  if (isTRUE(standAlone)) {
    shim_env$.Library <- tail(.libPaths(), 1)
  } else {
    shim_env$.Library <- .libPaths()
  }
  shim_env$.Library.site <- character()

  environment(shim_fun) <- shim_env
  shim_fun(unique(libPaths))
}

setLibPaths("bin/",standAlone=FALSE) #add MitoSAlt bin to paths

#check Biostring library
check.biostring<-"Biostrings" %in% installed.packages()

if(check.biostring=="FALSE"){
    cat("Biostrings not installed, installing now!\n")
    
    rversion<-as.numeric(gsub(".*version\\s(.{3}).*","\\1",R.version.string))
    if(rversion<=3.5){
        source("http://bioconductor.org/biocLite.R")
        biocLite("Biostrings")
    }
    else if(rversion>3.5){
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("Biostrings",lib="bin")
    }
} 

#LOAD LIBRARIES
library(plotrix)
library(RColorBrewer)
library(Biostrings)

#Let's read the parameters and see what you put there!
args <- commandArgs(trailingOnly = TRUE)

mlength<-as.numeric(args[1])
ohs<-as.numeric(args[2])
ohe<-as.numeric(args[3])
ols<-as.numeric(args[4])
ole<-as.numeric(args[5])
del.limit<-as.numeric(args[6])
clsfile<-args[7]
bpfile<-args[8]
filename<-args[9]
hp.limit<-as.numeric(args[10])
genome.file<-args[11]
flank<-as.numeric(args[12])

#PARAMETERS TO IDENTIFY FLANKING SEQUENCES
mt.fa = paste(readDNAStringSet(genome.file))
mlength<-nchar(mt.fa)

endcheck.max<-mlength-20
endcheck.min<-20
random.shift<-550
del.limit<-10000
mat<-nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)

#FUNCTION
set.radius<-function(dat=NULL,radius=NULL,radius.diff=NULL,decreasing=NULL){
    dat<-dat[order(dat$delsize,decreasing=decreasing),]
    dat$radius<-1000
    for(i in 1:nrow(dat)){
        if(dat$radius[i]!=1000){next}
        
        start<-dat$start[i]
        end<-dat$end[i]
        deg1<-dat$deg1[i]
        deg2<-dat$deg2[i]
        arc<-ifelse(deg1<deg2,'clk','aclk')
        radius<-radius-radius.diff
        dat$radius[i]<-radius
        
        for(j in 1:nrow(dat)){
            if(dat$radius[j]!=1000){next}
            cstart<-dat$start[j]
            cend<-dat$end[j]
            cdeg1<-dat$deg1[j]
            cdeg2<-dat$deg2[j]
            carc<-ifelse(cdeg1<cdeg2,'clk','aclk')
            
            if(start<cend & end>cstart & arc=='clk' & carc=='clk'){}
            else if((start<cstart | end>cend) & arc=='clk' & carc=='aclk'){}
            else if((cstart<start | cend>end) & arc=='aclk' & carc=='clk'){}
            else if(start<cend & end>cstart & arc=='aclk' & carc=='aclk'){}
            else{
                dat$radius[j]<-radius
                tmp<-dat[dat$radius==radius,]
                accept<-'yes'
                
                for(k in 1:nrow(tmp)){
                    dstart<-tmp$start[k]
                    dend<-tmp$end[k]
                    ddeg1<-tmp$deg1[k]
                    ddeg2<-tmp$deg2[k]
                    darc<-ifelse(ddeg1<ddeg2,'clk','aclk')
                    
                    if(cstart==dstart & cend == dend){next}
                    else if(cstart<dend & cend>dstart & carc=='clk' & darc=='clk'){accept='no'}
                    else if((cstart<dstart | cend>dend) & carc=='clk' & darc=='aclk'){accept='no'}
                    else if((dstart<cstart | dend>cend) & carc=='aclk' & darc=='clk'){accept='no'}
                    else if(cstart<dend & cend>dstart & carc=='aclk' & darc=='aclk'){accept='no'}
                }
                if(accept=='no'){dat$radius[j]<-1000}
            }
        }
    }
    dat
}

#TO CHECK REPEAT SEQUENCES AROUND EACH PAIR OF BREAKPOINTS
align.bp<-function(start=NULL,end=NULL,mt.fa=NULL,mat=NULL,nb=NULL){
    
    bp<-apply(data.frame(start=start-nb,end=start+nb),1,function(x){substr(mt.fa,x[1],x[2])})
    bp_1<-apply(data.frame(start=start-nb,end=start),1,function(x){substr(mt.fa,x[1],x[2])})
    bp_2<-apply(data.frame(start=start+1,end=start+nb),1,function(x){substr(mt.fa,x[1],x[2])})
    bp_res<-apply(data.frame(bp_1,bp_2),1,function(x){paste(x[1],x[2],sep="*")})
    
    bp1<-apply(data.frame(start=end-nb,end=end+nb),1,function(x){substr(mt.fa,x[1],x[2])})
    bp1_1<-apply(data.frame(start=end-nb,end=end),1,function(x){substr(mt.fa,x[1],x[2])})
    bp1_2<-apply(data.frame(start=end+1,end=end+nb),1,function(x){substr(mt.fa,x[1],x[2])})
    bp1_res<-apply(data.frame(bp1_1,bp1_2),1,function(x){paste(x[1],x[2],sep="*")})
    
    bp<-gsub("N","A",bp)
    bp1<-gsub("N","A",bp1)  
    bp.df<-data.frame(a=bp,b=bp1,seq1=bp_res,seq2=bp1_res,seq=NA,score=NA,sbp1=NA,ebp1=NA,sbp2=NA,ebp2=NA,stringsAsFactors=F)
    
    for(j in 1:nrow(bp.df)){
        #cat(j,"\n")
        a<-bp.df$a[j]
        b<-bp.df$b[j]
        
        if (nchar(a)<31|nchar(b)<31){
            next
        }
        
        size <- 15
        found<-'no'
        for(i in 1:13){
            if(found=='yes'){break}
            size<-size-1
            posa<-1
            posb<-1+size-1
            
            while(posb <= 31){
                if(found=='yes'){break}
                #cat(size,posa,posb,"\n")
                if(posa<=17 & posb>=15){
                    tmp<-substr(a,posa,posb)
                    tmpm<-matchPattern(tmp,b)
                    
                    start1<-start(tmpm)
                    end1<-end(tmpm)
                    
                    if (length(start1)>=1){
                        if(length(start1)>1){
                            for(k in 1:length(start1)){
                                if(start1[k]<=17 & end1[k]>=15){found<-'yes'}
                            }
                            start1<-paste(start1,collapse=",")
                            end1<-paste(end1,collapse=",")
                        }
                        
                        else if(start1<=17 & end1>=15){
                            #cat(tmp,posa,posb,start1,end1,"\n")
                            found<-'yes'
                        }
                    }
                }
                posa=posa+1
                posb=posb+1
            }
        }
        if(found=='no'){
            tmp<-NA
            posa<-NA
            posb<-NA
            start1<-NA
            end1<-NA 
        }
        
        bp.df$seq[j]<-tmp
        bp.df$sbp1[j]<-posa-1
        bp.df$ebp1[j]<-posb-1
        bp.df$sbp2[j]<-start1
        bp.df$ebp2[j]<-end1
    }
    bp.df$score<-nchar(bp.df$seq)
    bp.df$a<-NULL
    bp.df$b<-NULL
    bp.df
}


plot.col<-brewer.pal(12,"Paired")

#SET COLOURS AND OUTPUT FILES
reds<-colorRampPalette(c('lavenderblush','lightpink','lightcoral','indianred','firebrick','#2C0E10FF'),space="rgb")
blues<-colorRampPalette(c('lightcyan','lightblue1','lightskyblue','royalblue','navyblue','#0f0e2cf0'),space="rgb")
plotfile<-paste("plot/",filename,".pdf",sep="")
textfile<-paste("indel/",filename,".tsv",sep="")

#SET COLOUR KEY
dat.col<-data.frame(value=1:10000/100,value.log=log(1:10000/100))
dat.col$del.col<-blues(1000)[cut(dat.col$value.log,breaks = 1000)]
dat.col$dup.col<-reds(1000)[cut(dat.col$value.log,breaks = 1000)]
dat.del.col<-data.frame(value=dat.col$value,col=dat.col$del.col)
dat.dup.col<-data.frame(value=dat.col$value,col=dat.col$dup.col)

#UPLOAD BREAKPOINT CLUSTERS
if(file.info(clsfile)$size>0){
    res <- read.delim(clsfile, header=FALSE,stringsAsFactors=F)
    colnames(res)<-c("cluster","read","del.start","del.end","lfstart","lfend","nread","tread","perc")
    res<-res[!is.na(res$cluster),]
    res$sample <- gsub(".cluster","",clsfile)
    res$del.start.median <- apply(res,1,function(x){median(as.numeric(unlist(strsplit(as.character(x[3]),","))))})
    res$del.end.median <- apply(res,1,function(x){median(as.numeric(unlist(strsplit(as.character(x[4]),","))))})
    res$del.start.min <- apply(res,1,function(x){min(as.numeric(unlist(strsplit(as.character(x[3]),","))))})
    res$del.start.max <- apply(res,1,function(x){max(as.numeric(unlist(strsplit(as.character(x[3]),","))))})
    res$del.end.min <- apply(res,1,function(x){min(as.numeric(unlist(strsplit(as.character(x[4]),","))))})    
    res$del.end.max <- apply(res,1,function(x){max(as.numeric(unlist(strsplit(as.character(x[4]),","))))})    
    res$del.start.range <- paste(res$del.start.min,"-",res$del.start.max)
    res$del.end.range <- paste(res$del.end.min,"-",res$del.end.max)
    list.reads<-strsplit(res$read,",",fixed=T)
    length.list.reads<-sapply(list.reads,length)
    res.read<-data.frame(sample=rep(res$sample,length.list.reads),cluster=rep(res$cluster,length.list.reads),read=unlist(list.reads))
    
    #UPLOAD RAW BREAKPOINTS
    bp <- read.delim(bpfile, header=FALSE,stringsAsFactors=F)[,c(2,4,5,10)]
    colnames(bp)<-c("read","del.start","del.end","dloop")
    
    #FILTER CLUSTERS
    res<-res[!is.na(res$cluster),]
    bp<-bp[!is.na(bp$read),]
    
    res.read.bp<-merge(res.read,bp,by="read")
    res.read.bp1<-unique(res.read.bp[,c(2,3,6)])
    
    res<-merge(res,res.read.bp1,by=c("sample","cluster"))
    res$delsize<-res$del.end.median-res$del.start.median
    res[res$dloop=="yes",]$delsize<-(mlength-res[res$dloop=="yes",]$del.end.median)+res[res$dloop=="yes",]$del.start.median

    
    #FILTER BY PREDICTED HETEROPLASMY
    if(nrow(res[res$perc>=hp.limit,])>0){
        res<-res[res$perc>=hp.limit,]
        res$final.event<-'del'
        
        #go through OriH and breakpoints to choose the deletions which can be flipped into duplications
        for(i in 1:nrow(res)){
            dloop<-res$dloop[i]
            Rs<-ohs
            Re<-ohe
            Ds<-res$del.start.median[i]
            De<-res$del.end.median[i]
            Dsr<-res$del.start.range[i]
            Der<-res$del.end.range[i]
            if(dloop == 'yes'){
                Ds<-res$del.start.median[i]
                De<-res$del.end.median[i]
                Dsr<-res$del.start.range[i]
                Der<-res$del.end.range[i]
                
            }

            if(Re>=Rs){
                #essential region NOT covering pos 0
                if(((Ds >= Rs) & (Ds <= Re)) | ((De >= Rs) & (De <= Re))){ #either start or end of the deletion overlaps with region
                    res$final.event[i]<-'dup'
                }
                else if((De > Ds) & (Ds <= Rs) & (De >= Re)){ #deletion NOT covering pos 0 AND completely overlaps region
                    res$final.event[i]<-'dup'
                }
                else if((De < Ds) & ((De >= Re) | (Ds <= Rs))){ #deletion IS covering pos 0 AND completely overlaps region
                    res$final.event[i]<-'dup'
                }
            }
            else{
                #essential region IS covering pos 0
                if((Ds >= Rs) | (Ds <= Re) | (De >= Rs) | (De <= Re)){ #either start or end of the deletion overlaps with region
                    res$final.event[i]<-'dup'
                }
                else if((De < Ds)){ #deletion IS covering pos 0 so MUST overlap the region
                    res$final.event[i]<-'dup'
                }
            }
        }
        
        #go through OriL and breakpoints to choose the deletions which can be flipped into duplications
        for(i in 1:nrow(res)){
            dloop<-res$dloop[i]
            Rs<-ols
            Re<-ole
            Ds<-res$del.start.median[i]
            De<-res$del.end.median[i]
            if(dloop == 'yes'){
                Ds<-res$del.end.median[i]
                De<-res$del.start.median[i]
            }
            
            if(Re>=Rs){
                #essential region NOT covering pos 0
                if(((Ds >= Rs) & (Ds <= Re)) | ((De >= Rs) & (De <= Re))){ #either start or end of the deletion overlaps with region
                    res$final.event[i]<-'dup'
                }
                else if((De > Ds) & (Ds <= Rs) & (De >= Re)){ #deletion NOT covering pos 0 AND completely overlaps region
                    res$final.event[i]<-'dup'
                }
                else if((De < Ds) & ((De >= Re) | (Ds <= Rs))){ #deletion IS covering pos 0 AND completely overlaps region
                    res$final.event[i]<-'dup'
                }
            }
            else{
                #essential region IS covering pos 0
                if((Ds >= Rs) | (Ds <= Re) | (De >= Rs) | (De <= Re)){ #either start or end of the deletion overlaps with region
                    res$final.event[i]<-'dup'
                }
                else if((De < Ds)){ #deletion IS covering pos 0 so MUST overlap the region
                    res$final.event[i]<-'dup'
                }
            }
        }


        
        #CLASSIFY POTENTIAL DUPLICATIONS
        dat<-data.frame(chr="MT",start=res$del.start.median,end=res$del.end.median,value=res$perc,dloop=res$dloop,delsize=res$delsize,final.event=res$final.event)
    
        #ADD DEGREES
        dat$deg1<-90+(358*dat$start/mlength)
        dat$deg2<-90+(358*dat$end/mlength)
        if(nrow(dat[dat$dloop=='no' & dat$final.event=='dup',])>0){dat[dat$dloop=='no' & dat$final.event=='dup',]$deg1<-360+dat[dat$dloop=='no' & dat$final.event=='dup',]$deg1}
        
        #SPLIT PREDICTED DELETIONS/DUPLICATIONS AND ADD COLORS
        dat1<-dat[dat$final.event=='del',]
        dat2<-dat[dat$final.event=='dup',]
        if(nrow(dat2)>0){dat2$delsize<-mlength-dat2$delsize}
        
        if(nrow(dat1)>0){
            dat1$value<-round(dat1$value,2)
            dat1<-merge(dat1,dat.del.col,by='value')
            dat1$col<-as.character(dat1$col)
        }
        if(nrow(dat2)>0){
            dat2$value<-round(dat2$value,2)
            dat2<-merge(dat2,dat.dup.col,by='value')
            dat2$col<-as.character(dat2$col)
        }
        
        dat<-rbind(dat1,dat2)
        
        #GET PLOTTING PARAMETERS
        if(nrow(dat)<=100){
            radius.diff<-4
            lwd.arc<-1.5
        } else if (nrow(dat)>100 & nrow(dat)<=250){
            radius.diff<-2
            lwd.arc<-0.5
        } else if (nrow(dat)>250 & nrow(dat)<=400){
            radius.diff<-1
            lwd.arc<-0.2
        } else if (nrow(dat)>400 & nrow(dat)<=800){
            radius.diff<-0.5
            lwd.arc<-0.1
        } else if (nrow(dat)>800){
            radius.diff<-0.2
            lwd.arc<-0.05
        }
        
        #ASSIGN RADII TO DELETIONS/DUPLICATIONS TRYING TO FIT NON OVERLAPPING ONES IN THE SAME LEVEL AND KEEPING THE LARGEST ONES NEAR THE PERIPHERY
        dat<-set.radius(dat,400,radius.diff,T)
        
        #BUILD MT AXIS
        mt.axis<-data.frame(name=0:floor(mlength/1000),position=seq(0,mlength,1000))
        mt.axis$deg.axis<-90+(358*mt.axis$position/mlength)
        mt.axis$name<-as.character(mt.axis$name)
        
        #PLOT DELETIONS/DUPLICATIONS
        pdf(plotfile)
        par(mar=c(1,1,1,1),xpd=TRUE)
        plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
        draw.circle(400,400,410,lwd=3,border="gray50")
        for(i in 1:nrow(mt.axis)){
            draw.radial.line(410, 415, center=c(400,400),deg=mt.axis$deg.axis[i],col="gray40")
            radialtext(mt.axis$name[i], center=c(400,400),deg=mt.axis$deg.axis[i],start=416,nice=T,col="gray30")
        }
        with(dat[order(dat$value),],draw.arc(400,400,radius,deg1=deg1,deg2=deg2, lwd=lwd.arc,col=col))
        
        dev.off()
        
        
        #Change deletion direction base on dloop value
        for(i in 1:nrow(res)){
            dloop<-res$dloop[i]
            if(dloop=="yes"){
                del.start.median <- res$del.start.median[i]
                del.end.median <- res$del.end.median[i]
                del.start.range <- res$del.start.range[i]
                del.end.range <- res$del.end.range[i]
                del.start <- res$del.start[i]
                del.end <- res$del.end[i]
                lfstart <- res$lfstart[i]
                lfend <- res$lfend[i]
                res$del.start.median[i] <- del.end.median
                res$del.end.median[i] <- del.start.median
                res$del.start[i] <- del.end
                res$del.end[i] <- del.start
                res$lfstart[i] <- lfend
                res$lfend[i] <- lfstart
                res$del.start.range[i] <- del.end.range
                res$del.end.range[i] <- del.start.range

                
            }
        }        
        res$dloop <- NULL
        res$del.start <- NULL
        res$del.end <- NULL
        res$lfstart <- NULL
        res$lfend <- NULL
        res$read <- NULL
        res$del.start.min <- NULL
        res$del.start.max <- NULL
        res$del.end.min <- NULL
        res$del.end.max <- NULL
        
        res$perc <- format(round(res$perc, 4), nsmall = 4)
        
        res <- within(res, {
            del.start.median <- ifelse(final.event == "del", del.start.median+1, del.start.median+1)
            final.event.size <- ifelse(final.event == "del", delsize, (mlength-delsize))
            final.end <- ifelse(final.event == "del", del.end.median, del.start.median-1)
            final.start <- ifelse(final.event == "del", del.start.median, del.end.median+1)
        })
        res <- within(res, {
            del.start.median <- ifelse(del.start.median == mlength+1, 1, del.start.median)
            final.start <- ifelse(final.start == mlength+1, 1, final.start)
        })
        
        #get flanking sequences and length of repeats near flanks
        res.adl<-align.bp(res$final.start-1,res$final.end,mt.fa,mat,flank)
        res.final<-as.data.frame(cbind(res,res.adl$seq1, res.adl$seq2, res.adl$seq))
        colnames(res.final) <- c("sample","cluster.id","alt.reads","ref.reads","heteroplasmy","del.start.median","del.end.median","del.start.range","del.end.range","del.size","final.event","final.start","final.end","final.size", "seq1", "seq2", "seq")
        
        res.final$del.start.median <- NULL
        res.final$del.end.median <- NULL
        
        #SAVE RESULTS
        write.table(res.final,file=textfile,sep="\t",quote=F,row.names=F)
    }
}


#____END____#
#HUMAN
#mlength<-16569
#ohs<-16024
#ohe<-191
#ols<-5734
#ole<-5759
#del.limit<-10000
#clsfile<-"D.cluster"
#bpfile<-"D.breakpoint"
#filename<-"test_a"
#hp.limit<-0.01

















