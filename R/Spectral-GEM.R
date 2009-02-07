installGEMcore<-function(machine="linux") 
{
  url="http://www.stat.cmu.edu/~jwu/GEMfolder/Spectral-GEM/"
  if (machine=="linux") {
     download.file(paste(url,"linuxGEMfiles/Spectral-GEM_Rv1.2.f90",sep=""),
                   "Spectral-GEM_Rv1.2.f90",mode="wb")
     download.file(paste(url,"linuxGEMfiles/GEM_sub.f",sep=""),"GEM_sub.f")
  } else {
     if (machine=="windows") {
     download.file(paste(url,"winGEMfiles/Spectral-GEM_Rv1.2.exe",sep=""),
         "Spectral-GEM_Rv1.2.exe",mode="wb")
     } else {
       cat("Error: must specify machine type\n")
       return()
    }
  }
  return()
}

SpectralGEM<-function(InputFile="matching_input.txt",machine="linux",CM="CM") 
{
  if (machine=="linux") {
     main="mainl"
     if (!file.exists("./Spectral-GEM")) {
         cat("You need to download fortran executables from\n")
         cat("http://www.stat.cmu.edu/~jwu/GEMfolder/Spectral-GEM/.\n")
         cat("You need to have ifort compiler otherwise try our windows version.\n")
         cat("Proceed? Enter \"Y\" (yes) or \"N\" (no) and hit enter twice.\n")
         get=scan(what="character")
         if (get=="Y") {
            installGEMcore(machine="linux")
          } else {
            return()
          }
      }
  } else {
     if (machine=="windows") { 
         main="mainw"
         if (!file.exists("./Spectral-GEM_Rv1.2.exe")) {
              cat("You need to download fortran executables from\n")
              cat("http://www.stat.cmu.edu/~jwu/GEMfolder/Spectral-GEM/.\n")
              cat("Proceed to download?\n")
              cat("Enter \"Y\" (yes) or \"N\" (no) and hit enter twice.\n")
              get=scan(what="character")
              if (get=="Y") {
                installGEMcore(machine="windows")
              } else {
                return()
              }
            }
       } else {
        cat("machine must be \"windows\" or \"linux\"\n")
        return()
     }
  } 
  # stage 1
  # The program runs main Fortran program.
  # input: matching_input.txt
  # enter C or M        :: C   

  U=NULL
  cl=NULL
  d=NULL
  if (CM=="C" || CM=="CM") {
    write(c(InputFile,"C"),file="tmp00001.txt",ncolumns=1)
    .C(main, PACKAGE = "SpectralGEM")
     ext=getVersion(InputFile)
    cat("Current version: ext=",ext,";\n")
      data=read.table(paste("clusters_",ext,".txt",sep=""),
            header=FALSE)
      cl=data[,c(1,5)]
      names(cl)<-c("sampleId","clusterId")
      U=data[,c(1,4,6:dim(data)[2])]
      names(U)<-c("sampleId","group",paste("U",(3:dim(U)[2])-3,sep=""))
      lam<-read.table(paste("significant_eigenvalues_",ext,".txt",sep=""),
             header=FALSE)
      lam=1-lam[,2]
      names(lam)<-paste("lambda",(1:length(lam))-1,sep="")
    
      # The program extracts the file extension
      # and generates an exlusion list by checking the case-control distances
      # User chooses using simulation or visual inspection
      # to select cut off 
      ext=getVersion(InputFile);
      x=getRecord(ext)
      excludeFile=getExcludeFile(InputFile);
      if (x$n.case>0 & x$n.cntrl>0) {
         updateExcludeFileDstr(excludeFile,ext)      
      }   
      updateInputFile(oldInputFile=InputFile,
          newInputFile=InputFile,stage="M",excludeFile)
      ext=getVersion(InputFile);
   }
   # stage 2

   # the program runs the main Fortran program
   # input: matching_input.txt; 
   # enter C or M        :: M
  ma=NULL
  if (CM=="M" || CM=="CM") {
      write(c(InputFile,"M"),file="tmp00001.txt",ncolumns=1)
     .C(main, PACKAGE = "SpectralGEM") 
      ext=getVersion(InputFile);
     cat("Current version: ext=",ext,";\n")
           
    data=read.table(paste("significant_eigenvector_",ext,".txt",sep=""),
            header=FALSE)
     U=data[,c(1,4,5:dim(data)[2])]
     names(U)<-c("sampleId","group",paste("U",(3:dim(U)[2])-3,sep=""))
     lam<-read.table(paste("significant_eigenvalues_",ext,".txt",sep=""),
             header=FALSE)
     lam=1-lam[,2]
     names(lam)<-paste("lambda",(1:length(lam))-1,sep="")
      x=getRecord(ext)      
      if (x$n.case==0 || x$n.cntrl==0) {
          cat("Missing case or control, cannot perform matching cases 
to controls.\n")
          return(list(cluster=cl,U=U,lambda=lam,d=d,ma=ma))
      }     
      ma=full_matching(ext);

     d=read.table(paste("distance_",ext,".matrix",sep=""),header=T)
  }
  return(list(cluster=cl,U=U,lambda=lam,d=d,ma=ma))
}

updateExcludeFileDstr<-function(excludeFile,ext) 
{
 data=read.table(paste("minimum_distance_",ext,".txt",sep=""),header=FALSE)
 n_dx=3
 n_dst=5
 x1=data[data[,n_dx] ==2, n_dst]
 x2=data[data[,n_dx] ==1, n_dst]
 maxX=max(max(x1),max(x2))
 par(mfrow=c(2,1))
 hist(x1,xlim=c(0,maxX),breaks=20,main="Distance of Cases to 
Controls",xlab="Minimum Distance")
 hist(x2,xlim=c(0,maxX),breaks=20,main="Distance of Controls to 
Cases",xlab="Minimum Distance")
 cat("\nEnter the cutoff to use for excluding individuals (hit enter twice) :\n")
 max_dist= scan(what='numeric')
 exclude=data[data[,n_dst] > as.numeric(max_dist),]
 #(pause, select cut off, say 0.05)
 if (dim(exclude)[1]>0) {
   write.table(c(exclude,ext),append=T,
   file=excludeFile,row.names=FALSE,col.names=FALSE,quote=FALSE)           
    cat(excludeFile,"was updated.\n")
  } else {
    cat(excludeFile,"was not updated.\n")
  }  
  return
}

getRecord<-function(ext)
{
 logfile=paste("GEM_log_",ext,".txt",sep="")
 x=scan(logfile,what="character",multi.line=T,quiet=T)
 xl=length(x)
 a=c("COUNTS", "AFTER", "EXCLUSION", "OF", "INDIVIDUALS")
 al=length(a) 
 for (i in 1:(xl-al)) {
   m=0
   for (j in 1:al) {
    if (x[i+j]==a[j]) m=m+1;
   }
   if (m==al) break 
 }
 i0=i+j 
 a=c("number","of", "cases")
 al=length(a) 
 for (i in i0:(xl-al)) {
   m=0
   for (j in 1:al) {
    if (x[i+j]==a[j]) m=m+1;
   }
   if (m==al) break 
 }
 n.case=as.numeric(x[i+al+2])
 a=c("number","of", "controls")
 al=length(a) 
 for (i in i0:(xl-al)) {
   m=0
   for (j in 1:al) {
    if (x[i+j]==a[j]) m=m+1;
   }
   if (m==al) break 
 }
 n.cntrl=as.numeric(x[i+al+2])
 a=c("NUMBER", "OF", "SIGNIFICANT", "EIGENVALUES")
 al=length(a)
 for (i in i0:(xl-al)) {
   m=0
   for (j in 1:al) {
    if (x[i+j]==a[j]) m=m+1;
   }
   if (m==al) break 
 }
 n.dim=as.numeric(x[i+al+2])+1

 return(list(n.dim=n.dim,n.case=n.case,n.cntrl=n.cntrl))
}



loadInputFile<-function(InputFile) 
{
 task1<-read.table(InputFile,sep="?")  
 inputlist=as.vector(task1[,1])
 return(inputlist)
}

getVersion<-function(InputFile) 
{
 task1<-read.table(InputFile,sep="?")
 inputlist=as.vector(task1[,1]) 
 ext=inputlist[1]
 exth=substr(ext,1,5)
 return(exth)
}

getExcludeFile<-function(InputFile)
{
 task1<-read.table(InputFile,sep="?")  
 inputlist=as.vector(task1[,1])  
 excludeFile=inputlist[4]
 excludeFile=unlist(strsplit(excludeFile," "))[1]
 excludeFile=unlist(strsplit(excludeFile,"\t"))[1]
 cat("excludeFile=",excludeFile,";\n")
 return(excludeFile)
}

InFile<-function(identifier="smal",stage="1",directory="./",
                    MMfile="MMprime.txt",excludefile="exclude.txt",
                    idlength=8,mincluster=10,logtype=0,outfile="matching_input.txt")
{
  if (nchar(identifier)!=4) {
      cat("Error: identifier must have 4 letters.\n")
      return()
  }
  stage=as.character(stage)
  if (nchar(stage)!=1) {
      cat("Error: stage must be 1 letter.\n")
   }
  x=c(paste(identifier,stage,sep=""),directory,MMfile,excludefile,idlength,mincluster,logtype)
  write(x,file=outfile,ncolumns=1)
}

MMfile<-function(H=H,sampleInfo=id.info,outfile="MMprime.txt")
{
  nsample=dim(H)[1]
  x1=paste(nsample,":: number of individuals")
  ntag=1000  #estimated number of tag SNPs
  x2=paste(ntag,":: number of tag SNPs")
  write(x1,file=outfile)
  write(x2,file=outfile,append=T)
  MM=cbind(sampleInfo,H)
  write(t(MM),file=outfile,ncolumns=dim(MM)[2],append=T)
}
         

updateInputFile<-function(oldInputFile,newInputFile,stage,
  excludeFile)
{
  inputlist=loadInputFile(oldInputFile)
  ext=inputlist[1]
  ext=paste(substr(ext,1,4),stage,sep="");
  inputlist[1]=ext
  x=inputlist[2]
  n=nchar(x)
  for (i in seq(n,1,-1)) {
     if (substr(x,i,i)!=" " && substr(x,i,i)!="\t" ) break
  }
  inputlist[2]=substr(x,1,i)
  x=inputlist[3]
  n=nchar(x)
  for (i in seq(n,1,-1)) {
     if (substr(x,i,i)!=" " && substr(x,i,i)!="\t" ) break
  }
  inputlist[3]=substr(x,1,i)
  inputlist[4]=excludeFile
  cat("\"",inputlist[1],"\"\n",file=newInputFile,sep="",append=F)
  for (i in 2:4) {      
    cat("\"",inputlist[i],"\"\n",file=newInputFile,sep="",append=T)
  }
  for (i in 5:length(inputlist)) {      
    cat(inputlist[i],"\n",file=newInputFile,append=T)
  }
    return()
}

pc_graphs_GEMpClusters<-function(ext) {
  lam<-read.table(paste("significant_eigenvalues_",ext,".txt",sep=""),
             header=FALSE)  
  lam=1-lam[,2]
  lambda=lam[-1]
  data=read.table(paste("clusters_",ext,".txt",sep=""),header=FALSE)
  data=data[,-6]
  n.pc=dim(data)[2]-5
  grp=data[,4]
  color=c("black","red","green","blue","orange","pink","purple","dark green","light red")

  type=data[,5]

  symbol=c("a","b","c","d","e","f","g","h","i","j",
           "k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
  SYMBOL=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N",
           "O","P","Q","R","S","T","U","V","W","X","Y","Z")
  symbol=c(symbol,SYMBOL)

  n_pair=floor((n.pc-1)/2)

  filename1=paste("clusters_",ext,".pdf",sep="")
  pdf(file=filename1)
  if (n_pair>0) {
    for (i.pair in 1:n_pair)
	{
	i.pc=(i.pair-1)*2+1
        j.pc=i.pc+1

	ax.1=(i.pc+5)
	ax.2=(j.pc+5)

 	x.lab=sprintf(" PC %2i",i.pc)
 	y.lab=sprintf(" PC %2i",j.pc)
 	plot(data[,ax.1]*sqrt(lambda[i.pc]),
             data[,ax.2]*sqrt(lambda[j.pc]),
           col=color[grp],pch=symbol[type],cex=.50,xlab=x.lab,ylab=y.lab)
	}
}
if((n_pair*2) < n.pc)
	{
		i.pc=n.pc-1
		j.pc=n.pc
		ax.1=(i.pc+5)
		ax.2=(j.pc+5)
		x.lab=sprintf(" PC %2i",i.pc)
		y.lab=sprintf(" PC %2i",j.pc)
 		plot(data[,ax.1]*sqrt(lambda[i.pc]),
                     data[,ax.2]*sqrt(lambda[j.pc]),
             col=color[grp],pch=symbol[type],cex=.50,xlab=x.lab,ylab=y.lab)
      }
dev.off() 

filename2=paste("clusters_",ext,".ps",sep="")
postscript(file=filename2)
if (n_pair>0) {
   for (i.pair in 1:n_pair)
	{
	i.pc=(i.pair-1)*2+1
      j.pc=i.pc+1

	ax.1=(i.pc+5)
	ax.2=(j.pc+5)

 	x.lab=sprintf(" PC %2i",i.pc)
 	y.lab=sprintf(" PC %2i",j.pc)
 	plot(data[,ax.1]*sqrt(lambda[i.pc]),
             data[,ax.2]*sqrt(lambda[j.pc]),
        col=color[grp],pch=symbol[type],cex=.50,xlab=x.lab,ylab=y.lab)
	}
}
if((n_pair*2) < n.pc)
	{
		i.pc=n.pc-1
		j.pc=n.pc
		ax.1=(i.pc+5)
		ax.2=(j.pc+5)
		x.lab=sprintf(" PC %2i",i.pc)
		y.lab=sprintf(" PC %2i",j.pc)
 		plot(data[,ax.1]*sqrt(lambda[i.pc]),
                     data[,ax.2]*sqrt(lambda[j.pc]),
             col=color[grp],pch=symbol[type],cex=.50,xlab=x.lab,ylab=y.lab)
      }
dev.off() 
cat("Figures are in ",filename1,"and",filename2,".\n")

  return()
}


pc_graphs_GEMp<-function(ext) {
  lam<-read.table(paste("significant_eigenvalues_",ext,".txt",sep=""),
             header=FALSE)  
  lam=1-lam[,2]
  lambda=lam[-1]
  data=read.table(paste("significant_eigenvector_",ext,".txt",sep=""),
    header=FALSE)
  data=data[,-5]
  color=c("black","red","green","blue","orange","pink","purple",
     "dark green","light red")

  n.pc=dim(data)[2]-4
grp=data[,4]
dx=data[,3]

n_pair=floor((n.pc-1)/2)

filename1=paste("scatter_",ext,".pdf",sep="")
pdf(file=filename1)
if (n_pair>0) {
  for (i.pair in 1:n_pair)
	{
	i.pc=(i.pair-1)*2+1
      j.pc=i.pc+1

	ax.1=(i.pc+4)
	ax.2=(j.pc+4)

 	x.lab=sprintf(" PC %2i",i.pc)
 	y.lab=sprintf(" PC %2i",j.pc)
	plot(data[,ax.1]*sqrt(lambda[i.pc]),
             data[,ax.2]*sqrt(lambda[j.pc]),
         col=color[grp],pch=dx,cex=.50,xlab=x.lab,ylab=y.lab)
  }
}

if((n_pair*2) < n.pc)
	{
		i.pc=n.pc-1
		j.pc=n.pc
		ax.1=(i.pc+4)
		ax.2=(j.pc+4)
		x.lab=sprintf(" PC %2i",i.pc)
		y.lab=sprintf(" PC %2i",j.pc)
		plot(data[,ax.1]*sqrt(lambda[i.pc]),
                     data[,ax.2]*sqrt(lambda[j.pc]),
                col=color[grp],pch=dx,cex=.50,xlab=x.lab,ylab=y.lab)
      }

dev.off() 

filename2=paste("scatter_",ext,".ps",sep="")
postscript(file=filename2)
if (n_pair>0) {
  for (i.pair in 1:n_pair)
	{
	i.pc=(i.pair-1)*2+1
      j.pc=i.pc+1

	ax.1=(i.pc+4)
	ax.2=(j.pc+4)

 	x.lab=sprintf(" EV %2i",i.pc)
 	y.lab=sprintf(" EV %2i",j.pc)
	plot(data[,ax.1]*sqrt(lambda[i.pc]),
             data[,ax.2]*sqrt(lambda[j.pc]),
         col=color[grp],pch=dx,cex=.50,xlab=x.lab,ylab=y.lab)
  }
}

if((n_pair*2) < n.pc)
	{
		i.pc=n.pc-1
		j.pc=n.pc
		ax.1=(i.pc+4)
		ax.2=(j.pc+4)
		x.lab=sprintf(" PC %2i",i.pc)
		y.lab=sprintf(" PC %2i",j.pc)
		plot(data[,ax.1]*sqrt(lambda[i.pc]),
                     data[,ax.2]*sqrt(lambda[j.pc]),
              col=color[grp],pch=dx,cex=.50,xlab=x.lab,ylab=y.lab)
      }

dev.off() 
cat("Figures are in ",filename1,"and",filename2,".\n")
  
return()
}



full_matching<-function(ext)
{ 
  idc=read.table(paste("distance_",ext,".matrix",sep=""),nrows=1,colClass="character",header=F)
  data=read.table(paste("distance_",ext,".matrix",sep=""),skip=1,row.names=1)
  names(data)<-idc
  idr=row.names(data)
  matches=fullmatch(as.matrix(data))
  list_dist=tapply(names(matches), matches, FUN = function(x, dmat)
	 {
 	# min(
	dmat[match(x, dimnames(dmat)[[1]]), match(x, dimnames(dmat)[[2]])]    
	#, na.rm=TRUE)
	 }, dmat = as.matrix(data))
 dist_list=unlist(list_dist)
 index=which(is.na(dist_list))
 short_list=dist_list[-index]

 # par(mfrow=c(1,1))
 #hist(short_list,breaks=20,main="Histogram of Distances in the Full 
 #Matched set",xlab="Distance")
 idm<-names(matches)
 nr=length(idr)
 nc=length(idc)
 cc=c(rep(2,nr),rep(1,nc))
 id=unlist(c(idr,idc))
 cc1=cc[order(id)]
 matches1=cbind(matches[order(idm)],cc1)
 m0=matches1[!is.na(matches1[,1]),]
 id1=as.data.frame(row.names(m0))
  m01=cbind(id1,m0)
 names(m01)<-c("sampleId","stratum","case/control")
  return(m01)
}

