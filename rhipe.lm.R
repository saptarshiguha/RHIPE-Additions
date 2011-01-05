## inputfile="/voip/modified.jitter.traffic.rate.database/part-r-00072/data"

count.levels <- function(inputfolder,type='sequence',factors=NULL,converter=NULL){
  if(is.null(factors)) return(list())
  map <- as.expression(bquote({
    factor.names <- .(FACNAME)
    if(is.null(.(CONVERTER))) {
      for(fname in factor.names){
        rhstatus(sprintf("Finding Unique Levels for factor:%s",fname))
        p <- table(unlist(lapply(map.values,function(df) df[,fname])))
        names.p <- names(p)
        counts.p <- as.vector(p)
        for(i in 1:length(names.p)){
          rhcollect(c(fname),names.p[i])
        }}}},list(FACNAME=factors,CONVERTER=converter)))
  reduce <- expression(
      pre={ coll = c() },
      reduce={ coll <- unique(append(coll, unlist(reduce.values)))},
      post = {rhcollect(reduce.key, coll)})
  .tmp <- tempfile(pattern="rhipe.lm.")
  on.exit({rhdel(sprintf("/tmp/%s",.tmp))})
  z <- rhmr(map=map,reduce=reduce, combiner=TRUE
            ,ifolder=inputfolder
            ,ofolder=sprintf("/tmp/%s",.tmp)
            ,inout=c(type,"sequence"))
  z.result <- rhex(z)
  results <- rhread(sprintf("/tmp/%s",.tmp),type='sequence')
}

rhlm <- function(fml,data,type='sequence',factors=NULL,transform=NULL,compfac=NULL,mapred=NULL,drop.na=TRUE,debug=FALSE,...){
  ##specify empty list for no factors
  m <- proc.time()[3]
  if(is.null(factors)) .faclevel <- factors else{
    .faclevel <- if(is.list(factors)) factors else{
      if(is.null(factors)) stop("Specify the names of the factors")
      cat(sprintf("Counting Levels of %s factor%s for input: %s\n",length(factors)
                  ,if(length(factors)>1) "s" else "", data))
      count.levels(data,type=type,factors=factors,converter=compfac)
    }
  }
  xtras <- list(...)
  .tmp <- sprintf("rhipe.lm%s.Rdata",paste(sample(letters,15),collapse=''))
  .tmp2 <- sprintf("/tmp/rhipe.lm%s",paste(sample(letters,15),collapse=''))
  en <- new.env();en$.tmp <- .tmp
  rhsave(.faclevel,file=sprintf("/tmp/%s",.tmp),envir=en)
  betahat <- c(NA)
  attr(betahat,"fac.levels") <- .faclevel
  on.exit(
           tryCatch(rhdel(c(.tmp2,sprintf("/tmp/%s",.tmp))),error=function(e) print(e),finally=betahat)
  )
  setup <- as.expression(bquote({load(.(tmp))},list(tmp=.tmp)))
  map <- as.expression(bquote({
    a <- do.call("rbind",map.values)
    rhcounter("rhlm","NROWS",nrow(a))
    tryCatch({
        for(facl in .faclevel){
          fn <- facl[[1]]
          a[,fn] <- factor(a[,fn], levels=facl[[2]])
        }
      if(!is.null(.(mapmod)))
        a <- .(mapmod)(a)
      if(!is.null(a) && nrow(a)>0){
        mm <- model.matrix(.(FORMULA), data=a,contrasts.arg=if('model.matrix' %in% names(.(xtras))) .(xtras)$model.matrix else NULL)
        mf <- model.frame(.(FORMULA), data=a);
        mf1 <- mf[,1]
        xpx <-  crossprod(mm) ## (t(mm) %*% mm) , see "Least Squares
                              ## Calculations in R", by Douglas Bates, R News,
                              ## 2004 - http://cran.r-project.org/doc/Rnews/Rnews_2004-1.pdf
        xpy <-  t(crossprod(mm, mf1) ## t(t(mm) %*% mf1)
        ypy <- sum(mf1 * mf1)    # sum of y^2
        ys <- sum(mf1)       # sum of y
        rhcollect(0L,xpx)
        rhcollect(1L,xpy)
        rhcollect(2L,as.numeric(c(ys,ypy)))
      }
    },error=function(e){
      rhcounter("R_ERRORS",as.character(e),1)
      rhcollect("error", head(a))
    })
  },list(mapmod=transform,FORMULA=fml,xtras=xtras)))
  reduce <- expression(
      pre={sums <- 0;} ,
      reduce = {
        if(reduce.key[[1]]==2L) sums <- sums+apply(do.call("rbind",reduce.values),2,sum)
        else for(i in reduce.values) sums <- sums+i
      },
      post = {
        rhcollect(reduce.key, sums) }
      )
  if(is.null(mapred)) mapred=list()
  if(!"rhipe_map_buff_size" %in% names(mapred)){
    mapred[["rhipe_map_buff_size"]] <- 500
  }
  if(!"mapred.reduce.tasks" %in% names(mapred)){
    mapred[["mapred.reduce.tasks"]] <- 2
  }
  cat(sprintf("Starting Regression for %s\n",data))
  z <- rhmr(map=map,reduce=reduce,combiner=TRUE,inout=c(type,"sequence")
            ,ifolder=data
            ,ofolder=.tmp2
          ,shared=sprintf("/tmp/%s",.tmp)
          ,mapred=mapred
          ,setup=list(map=setup),jobname=paste(deparse(fml),collapse=""))
  z.result=rhex(z)
  z.read <- rhread(.tmp2,type='sequence')
  which.is.xpx <- 
  xpx <- z.read[unlist(lapply(z.read,function(r) if (r[[1]]==0L) TRUE else FALSE))][[1]][[2]] ## 0L == xpx
  xpy <- z.read[unlist(lapply(z.read,function(r) if (r[[1]]==1L) TRUE else FALSE))][[1]][[2]] ## 1L == xpy
  others <- z.read[unlist(lapply(z.read,function(r) if (r[[1]]==2L) TRUE else FALSE))][[1]][[2]] ## 1L == xpy
  if(debug==2){
    y <- list(xpx=xpx,xpy=xpy,others=others)
    attr(y,"rhlm") <- z.result$counters
    return(y)
  }
  nro <- z.result$counters$rhlm['NROWS']
  so.xpx <- solve(xpx)
  betahat <- so.xpx %*% t(xpy)
  RSS <- others[2] - t(betahat) %*% xpx %*% betahat
  df <-nro - ncol(xpx)
  sigma.hat <- sqrt(RSS/df)
  stderr <- sqrt(diag(as.numeric(sigma.hat^2) * so.xpx))
  t.value <- betahat/stderr
  
  t.pr <- pt(abs(t.value), df=df, lower.tail=FALSE)*2
  r.square <- 1 - RSS/(others[2] - nro*(others[1]/nro)^2)
  betahat <- data.frame(Estimate=betahat, "Std. Error"=stderr,"t value"=t.value, "Pr(>|t|)"=t.pr)
  attr(betahat,"stats") <-c(sigmahat=sigma.hat, r.sq=r.square,df=as.numeric(df),n=as.numeric(nro))
  attr(betahat,'fac.levels') <- .faclevel
  attr(betahat,"counters") <- z.result$counters
  attr(betahat,"call") <- match.call()
  attr(betahat,"elapsed") <- proc.time()[3]-m
  if(debug==1) attr(betahat,"proj") <- list(xpx=xpx,xpy=xpy,others=others)
  return(betahat)
}

## rs <- rhlm(jitter~traffic.rate*rm.site, data=inputfile,factors=list(list("rm.site",c("Paris","Zurich","Brussels")))
##            ,type='map'
##            ,transform=function(a){
##              a$traffic.rate <- a$traffic.rate/1e6
##              a
##            })
## a <- rhread(inputfile,max=20)
## map.values <- lapply(a,"[[",2)
## a <- do.call("rbind",map.values)
## fn <- "rm.site";a[,fn] <- factor(a[,fn], attr(rs,"fac.levels")[[1]][[2]])
## b=a[1:3,]
## mm <- model.matrix(jitter~traffic.rate*rm.site, data=b)
## mf <- model.frame(jitter~traffic.rate*rm.site, data=b)
## xpx <-  (t(mm) %*% mm) 
## xpy <-  t(t(mm) %*% mf[,1])

make.words <- function(N,dest,cols=5,p=5,factor=1,local=FALSE){
  ## p is how long the word will be, longer more unique words
  ## factor, if equal to 1, then exactly N rows, otherwise N*factor rows
  ## cols how many columns per row
  map <- as.expression(bquote({
    P <- .(P)
    COLS <- .(COLS)
    F <- .(F)
    lapply(map.values,function(r){
      for(i in 1:F){
        f <- sapply(1:COLS, function(n) paste(sample(letters,P ),collapse=""))
        rhcollect(NULL,f)
      }
    })
  },list(COLS=cols,P=p,F=factor)))
  mapred <- list()
  if (local) mapred$mapred.job.tracker <- 'local'
  mapred[['mapred.field.separator']]=" "
  mapred[['mapred.textoutputformat.usekey']]=FALSE
  mapred$mapred.reduce.tasks=0
  z <- rhmr(map=map, N=N,ofolder=dest,inout=c("lapp","text"),
       mapred=mapred)
  rhex(z)
}


N <- 10000000
filename <- sprintf("/tmp/%s",paste(sample(letters,6),collapse=""))
make.words(N,filename,factor=10)

map <- expression({
  f <- table(unlist(strsplit(map.values," ")))
  n <- names(f)
  p <- as.numeric(f)
  sapply(seq_along(n),function(r) rhcollect(n[r],p[r]))
})
reduce <- expression(
    pre={ total <- 0},
    reduce = { total <- total+sum(unlist(reduce.values)) }
    post = { rhcollect(reduce.key,total) }
    )
z <- rhmr(map=map,reduce=reduce, inout=c("text","sequence")
          ,ifolder=filename
          ,ofolder=sprintf("%s-out",filename))
job.result <- rhstatus(rhex(z,async=TRUE),mon.sec=2)
if(job.result$status != "SUCCEEDED")
  error("THE JOB DID NOT WORK")
results <- rhread(sprintf("%s-out",filename))
results <- data.frame(words=unlist(lapply(results,"[[",1)), count = =unlist(lapply(results,"[[",2)))
message(sprintf("Did the job return the correct result? %s", sum(results[,2])== N*5))



         
         
