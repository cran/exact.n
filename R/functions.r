
#' Image plot of power
#'
#' Creates image plot of exact power for range of values of n0 and n1.
#' Optionally, the image plot will show where power exceeds a provided power
#' target.
#'
#' data will be one of the 85 saved objects corresponding to a selected value
#' of alpha and delta.
#'
#' @aliases power.image
#' @param data a matrix of powers with columns named beta p0.vals, n0 and n1.
#' @param p a scalar value for p0.vals that is used to select the subset of data.
#' @param binary If TRUE only the binary indicator of power exceeding beta is displayed.
#' @param beta target value of power
#' @return No return value and if the 'binary' parameter is TRUE then an image will be displayed
#' @author Chris J. Lloyd
#' @importFrom graphics image
#' @examples
#' oldpar <- graphics::par()
#' # Load toy version of power library for alpha=0.025, delta=0.20.
#' rdata_file = system.file('files', 'LIB.a025.d20.Rdata', package = 'exact.n')
#' load(rdata_file)
#' graphics::par(mfrow=c(1,2))
#' power.image(LIB.a025.d20,p=.5,binary=FALSE)
#' power.image(LIB.a025.d20,p=.2,beta=.7,binary=TRUE)
#' suppressWarnings(graphics::par(oldpar))
#' @export power.image
power.image=function(data,p=0.2,binary=TRUE,beta=.75){
  N0=max(data[,4])
  N1=max(data[,3])
  p0=data[,2]
  image.data=data[p0==p,]
  image.data=image.data[order(image.data[,4]),]
  image.data=image.data[order(image.data[,3]),]
  image.data=matrix(image.data[,1],ncol=(N0-19))
  if(!binary){graphics::image(20:N0,20:N1,t(image.data),
                   xlab=expression("n"[0]) ,ylab=expression("n"[1]))}
  if(binary){graphics::image(20:N0,20:N1,t(image.data>beta),
                   xlab=expression("n"[0]) ,ylab=expression("n"[1]))}
  #abline(0,1,lty=3,col="gray")
  NULL
}


#' Find smallest value of n1 that achieves target power
#'
#' Function calculates minimum value of n1 that achieves power beta. If there
#' is no solution less than 500, it models observed powers as a function of n1
#' and then extrapolates. It returns an infinite value if the power is
#' unattainable.
#'
#' This function is called by n1.get.vector and will likely never be run by the
#' user.
#'
#' The data matrix will be a subset of one of the 85 main databases. Supplying
#' alpha and delta loads the appropriate database and selecting a value of p0
#' further subsets this data base. The resulting matrix will have 4 columns and
#' 481^2 rows corresponding to all values of n0 and n1 from 20 to 500 inclusive
#' and is suitable for input into n1.get.
#'
#' If type=1, the smallest value N1 of n1 so that power > beta is returned. If
#' type=2, the smallest value N1 of n1 so that power > beta for all n1>=N1 is
#' returned.
#'
#' @aliases n1.get
#' @param data a matrix with column names beta, p0.vals, n0, n1 and 481^2 rows
#' (See details)
#' @param n0 scalar value of n0 between 20 and 500 inclusive
#' @param beta scalar target for power
#' @param delta value of clinically relevant difference
#' @param alpha value of nominal size of test
#' @param type type of minimum solution (See details)
#' @param plt if true produce diagnostic plots
#' @return a single scalar value of n1. An integer indicates a solution was
#' found in the database. A non-integer indicates an extrapolated solution. An
#' Inf value indicates no extrapolated solution was found.
#' @author Chris J. Lloyd
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to
#' size and power constraints. Preprint. \doi{10.13140/RG.2.2.11828.94085}
#' @importFrom stats qnorm lm pnorm
#' @importFrom graphics lines abline plot
#' @examples
#'
#' # Load toy version of power library for alpha=0.025, delta=0.20.
#' # Alternatively, load the full library using fetch.data(alpha=0.025,delta=0.20)
#' rdata_file = system.file('files', 'LIB.a025.d20.Rdata', package = 'exact.n')
#' load(rdata_file)
#' data=LIB.a025.d20[LIB.a025.d20[,2]==0.5,] # select subset with p0=.5
#' # For given value of n0, what minimum value of n1 ensures power at least 0.7?
#' n1.get(data, n0=70, beta=.7, delta=.2, alpha=0.025,
#'        type = 1, plt = TRUE) # Explicit solution 63 found in data base
#' n1.get(data, n0=50, beta=.7, delta=.2, alpha=0.025,
#'        type = 1, plt = TRUE) # Approximate solution 131 extrapolated from toy data base
#' # You can check the accuracy of this extrapolated result:
#' # POWER(n0=50,n1=131,alpha=0.025,delta=0.2,p0=0.5,type="elr")
#' # The minimum power at p0=0.5 is 0.699, slightly less than 0.7.
#' # With the full library n1.get returns the correct answer of n1=136.
#' n1.get(data, n0=30, beta=.7, delta=.2, alpha=0.025,
#'        type = 1, plt = TRUE) # Extrapolated solution is infinite
#'
#'
#' @export n1.get
n1.get=function(data,n0,beta,delta,alpha,type=1,plt=FALSE){
  # Function calculates minimum value of n1 that achieves power beta.
  # If there is no solution less than 500, it models observed
  # powers as a function of n1 and then extrapolates.
  # It returns an infinite value if the power is unattainable.
  #
  # ARGUMENTS:
  # data: a matrix typically selected for a subset of p0 (with 481^2 rows)
  # n0, delta, alpha, beta: scalar values as defined in paper
  # type: (1) find smallest value N1 of n1 so that power > beta
  #       (2) find smallest value N1 of n1 so that power > beta for all n1>=N1
  # plt: if true you will see a plot of power against n1 with model
  #
  N0=max(data[,4])
  pwr=data[data[,3]==n0,1]
  names(pwr)=NULL
  test=(pwr>=beta)
  if(sum(test)>0){
    if(type==1){value=min((20:N0)[as.logical(cummax(test))])}
    if(type==2){value=min((20:N0)[as.logical(rev(cummin(rev(test))))])}
  }
  if(sum(test)==0){
    y=delta^2/(stats::qnorm(pwr)+stats::qnorm(1-alpha))^2
    n.inv=1/(20:N0)
    wgt=(20:N0)
    wgt=wgt*(wgt>N0/2)
    params=summary(stats::lm(y~n.inv,weights=wgt))$coef[,1]
    # cat(params)
    # In pathalogical examples param[2] can be negative. In this case
    # fit the unweighted regression and the value should be negative
    if(params[2]<0){
      wgt=0*wgt+1
      params=summary(stats::lm(y~n.inv,weights=wgt))$coef[,1]
    }
    y=delta/(stats::qnorm(beta)+stats::qnorm(1-alpha))
    value=params[2]/(y^2-params[1])
  }
  if(plt){
    graphics::plot(20:N0, pwr,pch=".",ylab="power",xlab="n1")
    if(max(pwr)<beta){
      graphics::lines(20:N0,stats::pnorm(delta/sqrt(params[1]+params[2]/(20:N0))-stats::qnorm(1-alpha)),
            col="blue",lty=3)
    }
    graphics::abline(h=beta,v=value,lty=3)
  }
  names(value)=NULL
  if(value<0){value=Inf}
  value
}


#' Find minimal n1 achieving target power for all values of n0
#'
#' For values of n0 covered in the matrix data (typically from 20 to 500) the
#' function calls n1.get and finds the minimal value of n1 achieving the target
#' power.
#'
#' This function is called by n1.get.solution and will likely never be run by
#' the user.
#'
#' The data matrix will be one of the 85 main databases.
#'
#' If type=1, the smallest value N1 of n1 so that power > beta is returned. If
#' type=2, the smallest value N1 of n1 so that power > beta for all n1>=N1 is
#' returned.
#'
#' @param data power database for selected value of alpha and delta (See
#' details)
#' @param beta scalar target for power
#' @param p0 scalar value of baseline probability
#' @param delta value of clinically relevant difference
#' @param alpha value of nominal size of test
#' @param type type of minimum solution (See details)
#' @return a list with element x (giving range of values of n0) and y (giving
#' minimal solutions for n1 of power > beta). integer indicates
#' @author Chris J. Lloyd
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to
#' size and power constraints. Preprint. \doi{10.13140/RG.2.2.11828.94085}
#' @keywords ~htest ~design
#' @keywords internal
n1.get.vector=function(data,beta,p0,delta,alpha,type=1){
  sub.data=data[data[,2]==p0,]
  N0=max(sub.data[,4])
  n1.values=NULL
  n0.values=N0:20
  for(n0 in n0.values){
    last.three=0
    if(n0<(N0-3)){
      no.solution=is.infinite(n1.values)
      last.three=sum(no.solution[length(no.solution)+1-1:3])
      #    cat(c(n0,last.three,length(n1.values)),"\n")
    }
    if(last.three==3){val=Inf}
    if(last.three<3){
      val=n1.get(data=sub.data,n0=n0,delta=delta,
                 alpha=alpha,beta=beta,type=type,plt=F)
      # cat(val)
    }
    n1.values=c(n1.values,val)
  }
  n1.values=n1.values[order(n0.values)]
  n0.values=sort(n0.values)
  list(x=n0.values,y=n1.values)
}


#' Find minimal n1 achieving target power for range of values of p0
#'
#' This function calls n1.get.vector, runs if for all values of p0 within the
#' supplied range and then takes the worst (i.e. largest) solution for n1
#'
#'
#' @param data power database for selected value of alpha and delta
#' @param beta scalar target for power
#' @param p0 single value or range of values for baseline probability
#' @param delta value of clinically relevant difference
#' @param alpha value of nominal size of test
#' @param type type of minimum solution (See details)
#' @return vector of solutions for n1 with name vector equal to range of n0
#' values
#' @author Chris J. Lloyd
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to
#' size and power constraints. Preprint. \doi{10.13140/RG.2.2.11828.94085}
#' @examples
#'
#' # Load toy version of power library for alpha=0.025, delta=0.20.
#' rdata_file = system.file('files', 'LIB.a025.d20.Rdata', package = 'exact.n')
#' load(rdata_file)
#' # n0 solutions when p0=0.5
#' n1.get.solution(LIB.a025.d20,beta=.7,p0=0.5,delta=0.2,alpha=0.025,type=1)
#' # n0 solutions for p0 between 0.4 ad 0.5
#' n1.get.solution(LIB.a025.d20,beta=.7,p0=c(0.4,0.5),delta=0.2,alpha=0.025,type=1)
#'
#'
#' @export n1.get.solution
n1.get.solution=function(data,beta,p0,delta,alpha,type=1){
  N0=max(data[,4])
  p0.vals=(1:99)/100
  if(length(p0)==1){P0=p0}
  if(length(p0)>1){P0=p0.vals[p0.vals<=max(p0)&p0.vals>=min(p0)]}
  vals=NULL
  for(p in P0){
    vals=rbind(vals,
               n1.get.vector(data=data,beta=beta,p0=p,
                             alpha=alpha,delta=delta,type=type)$y)
  }
  vals=apply(vals,2,max)
  names(vals)=20:N0
  vals
}


#' Provide sample size solutions for target size and power.
#'
#' Function gives smallest values for n1 as function of n0 that achieve target
#' size and power.
#'
#' If out > 0 all solutions (including n1=Inf) are returned. If out=0, infinite
#' values are suppressed. If out < 0, only output satisfying the balance
#' criterion are output.
#'
#' @param alpha value of nominal size of test
#' @param delta value of clinically relevant difference
#' @param beta scalar target for power
#' @param p0 single value or range of values for baseline probability
#' @param type type of maximisation (see n1.get documentation)
#' @param plt If TRUE, plot n1 solutions versus
#' @param out More solutions output if out > 0 than out < 0 (see details)
#' @param b.lim maximum imbalance of sample sizes
#' @param prin If TRUE, error messages will be printed.
#' @return list with elements n0 and n1
#' @note The appropriate data file needs to have been downloaded corresponding
#' to the desired value of alpha and delta. This can be done with the
#' fetch.data() function.
#' @author Chris J. Lloyd
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to
#' size and power constraints. Preprint. \doi{10.13140/RG.2.2.11828.94085}
#' @importFrom graphics par title plot abline
#' @examples
#'
#' # We are interested in designs with power at least 0.75 when exact size
#' # 0.025 and delta=0.20. Therefore, you would need to have downloaded
#' # LIB.a025.d20 using fetch(0.015,0.20). The example below instead uses
#' # the toy data that comes with the package. The baseline probability is
#' # assumed to be between 0.3 and 0.5.
#' rdata_file = system.file('files', 'LIB.a025.d20.Rdata', package = 'exact.n')
#' load(rdata_file)
#' #' main.function(.025,0.20,p0=c(0.3,0.5),beta=0.75,plt=TRUE)
#' # The value of the function is the minimum value of n1 for a range
#' # of values of n0. The sample size ratio is limited to 5 by default.
#'
#' @export main.function

main.function=function(alpha,delta,beta=.75,p0=.5,type=2,plt=FALSE,
                       out=(-1),b.lim=5, prin = TRUE){
  # reset "par" using on.exit()
  init_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(init_par))
  # Check allowable values of delta
  delta.vals=(5:20)/100
  if(!is.element(delta,delta.vals)){
    if (prin) cat("Allowable values of delta are 0.05, 0.06, ... , 0.19, 0.20","\n")
    stop()
  }
  # Check allowable values of alpha
  alpha.vals=c(0.01,0.02,0.025,0.05,0.10)
  if(!is.element(alpha,alpha.vals)){
    if (prin) cat("Allowable values of alpha are 0.01, 0.02, 0.025, 0.05, 0.10","\n")
    stop()
  }
  #
  # Local function
  LIB.name=function(alpha,delta){
    digit1=floor(alpha*10)
    digit2=floor(alpha*100)-digit1*10
    digit3=floor(alpha*1000)-digit1*100-digit2*10
    char.alpha=paste(digit1,digit2,digit3,sep="")
    digit1=floor(delta*10)
    digit2=floor(delta*100)-digit1*10
    char.delta=paste(digit1,digit2,sep="")
    object.name=paste("LIB.a",char.alpha,".d",char.delta,sep="")
    return(object.name)
  }
  #
  # options(warn=-1)
  full.data=get(LIB.name(alpha,delta))
  if(plt){
    graphics::par(mfrow=c(1,2))
    #
    p0.vals=(1:99)/100
    if(length(p0)==1){P0=p0}
    if(length(p0)>1){P0=p0.vals[p0.vals<=max(p0)&p0.vals>=min(p0)]}
    P0.mean=round(100*mean(P0),0)/100
    power.image(data=full.data,p=P0.mean,beta=beta)
    graphics::title(bquote(italic(beta)*phantom() > .(beta) ~ ',' ~ italic(p)[0] %in% .(sprintf('[%s]', toString(p0)))))  }

  result=n1.get.solution(data=full.data,beta=beta,p0=p0,
                         delta=delta,alpha=alpha,type=type)
  n0=as.numeric(names(result))
  names(result)=NULL
  if(out==0){
    n0=n0[is.finite(result)]
    result=result[is.finite(result)]
  }
  if(out<0){
    balance=exp(abs(log(result/n0)))
    n0=n0[balance<=5]
    result=result[balance<=b.lim]
  }
  # options(warn=0)
  if(plt){
    graphics::plot(n0,result,pch=20,
                   xlab=expression("n"[0]),ylab=expression("n"[1]),)
    graphics::abline(0,1)
    graphics::abline(0,2,lty=2,col="gray")
    graphics::abline(0,3,lty=2,col="gray")
    graphics::abline(0,4,lty=2,col="gray")
    graphics::abline(0,1/2,lty=2,col="gray")
    graphics::abline(0,1/3,lty=2,col="gray")
    graphics::abline(0,1/4,lty=2,col="gray")
  }
  list(n0=n0,n1=result)
}


#' Exact Power of Test for Selected Sample Sizes.
#'
#' For known values of the sizes n0, this function computes the
#' exact probability of rejecting the null as a function of baseline
#' probability.
#'
#' p0 values must be between max(0,-delta) and min(1,1-delta)
#'
#' @param n0 control sample size
#' @param n1 treatment sample size
#' @param alpha value of nominal size of test
#' @param delta value of clinically relevant difference
#' @param psi null value of risk difference p1-p0
#' @param type either "lr" for approximate or "elr" for quasi-exact test
#' @param p0 baseline probability. If missing, a grid of values is created for
#' plotting. A scalar value can also be supplied.
#' @param sided (1 or 2 sided test)
#' @param obj Optional object with all possible p-values. Must be a list with
#' elements y0, y1, P (typically output of lr.rd or ESTEP.rd). If not supplied
#' then object is generated from n0, n1 and psi.
#' @return list with element x (containing values of baseline probability)
#' and element y (containing corresponding exact powers)
#' @author Chris J. Lloyd
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to
#' size and power constraints. Preprint. \doi{10.13140/RG.2.2.11828.94085}
#' @importFrom stats dbinom
#' @examples
#' oldpar <- graphics::par()
#' # Typical usage
#' #       POWER(n0=65,n1=82,psi=0,type="lr") # Exact size of approximate lr test
#' #       POWER(n0=65,n1=82,psi=0,type="elr",delta=.2) # Exact power of quasi exact test
#' # To make examples run faster, the package includes objects that contain
#' # all possible values of various tests when n0=65, n1=82.
#' load(system.file('files', "lr.stats.Rdata", package = 'exact.n'))
#' # All possible values of LR statistic p-values for testing p1-p0>0
#' load(system.file('files', "elr.stats.Rdata", package = 'exact.n'))
#' # All possible values of ELR statistic p-values for testing p1-p0>0
#' load(system.file('files', "elr10.stats.Rdata", package = 'exact.n'))
#' # Object contains all exact p-values for testing if p1-p0>0.1
#' # All possible values of ELR statistic p-values for testing p1-p0>0.1
#' #
#' graphics::par(mfrow=c(1,2))
#' # When delta=0 this gives type 1 error. The first plot is for the approximate
#' # lr based p-value, the second is for the quasi-exact e-p-value (alpha=0.05)
#' plot(POWER(n0=65,n1=82,alpha=0.05,psi=0,delta=0),type="l",
#'      xlab=expression("p"[0]),ylab="exact size")
#' abline(h=.05,lty=2)
#' plot(POWER(obj=elr.stats,alpha=0.05,delta=0),type="l",
#'      xlab=expression("p"[0]),ylab="exact size")
#' abline(h=.05,lty=2)
#' #
#' # For these sample sizes, power curve is calculated below for
#' # values of delta=0.1, 0.12, 0.14, 0.16, 0.18, 0.20. Power
#' # is poor for detecting a difference of 0.1 (see red).
#' plot(POWER(obj=lr.stats,alpha=0.05,delta=0.1),type="l",
#'      xlab=expression("p"[0]),ylab="exact power")
#' TITLE=expression('Exact power of LR test of p'[1]*'-p'[0]*'>0.')
#' title(main=TITLE,cex.main=0.8)
#' lines(POWER(obj=lr.stats,alpha=0.05,delta=0.12))
#' lines(POWER(obj=lr.stats,alpha=0.05,delta=0.14))
#' lines(POWER(obj=lr.stats,alpha=0.05,delta=0.16))
#' lines(POWER(obj=lr.stats,alpha=0.05,delta=0.18))
#' lines(POWER(obj=lr.stats,alpha=0.05,delta=0.20))
#' lines(POWER(obj=lr.stats,alpha=0.05,delta=0.10),col="red")
#' #
#' # The results below are for testing p1-p0>0.1.
#' plot(c(0,.9),c(0,1),type="n",
#'      xlab=expression("p"[0]),ylab="Pr(reject null)")
#' lines(POWER(obj=elr10.stats,alpha=0.05,psi=0.1,delta=0.1)) # Note delta=psi
#' abline(h=0.05,lty=2)
#' lines(POWER(obj=elr10.stats,alpha=0.05,delta=0.25),col="blue")
#' TITLE=expression('Exact size and power of test of p'[1]*'-p'[0]*'>0.1.')
#' title(main=TITLE,cex.main=0.8)
#' legend(.4,.4,lty=c(1,1),col=c("black","blue"),box.col="white",
#'     legend=c("size: p1-p0=0.1","power: p1-p0=0.25"),cex=.7)
#' # When using the package the above plots would be generated by
#' # lines(POWER(n0=65,n1=82,alpha=0.05,psi=0.1,delta=0.10,type="elr"))
#' # lines(POWER(n0=65,n1=82,alpha=0.05,psi=0.1,delta=0.25,type="elr"))
#' suppressWarnings(graphics::par(oldpar))
#'
#' @export POWER
POWER=function(n0,n1,alpha=0.05,delta=0,psi=0,type="lr",sided=1,p0=NULL,obj=NULL){
  # n0, n1 - sample sizes of control and treatment group
  # delta = fixed value (the null value is typically negative but non-null value could be positive).
  # type - either "lr" or "elr"
  # p0 - possible vector so p1=p0+delta
  # obj - must have elements y0, y1, P
  #
  # value - element x=nuisance parameter, y=significance
  #       - errors describe (sig-alpha)/alpha in % terms.
  #
  # Local function
  grid=function(B=100,theta=.5){
    b<-round(B/3)+1
    x<-(0:b)/b
    y<-.5+.5*x^theta
    grid.vals<-round(10000*c(1-y,x,y))/10000
    sort(unique(grid.vals))
  }
  #
  if(missing(obj)){
    if(type=="lr"){obj=lr.rd(n0,n1,psi=psi,sided=sided)}
    if(type=="elr"){obj=ESTEP.rd(n0,n1,psi=psi,sided=sided)}
  }
  y0<-obj$y0
  y1<-obj$y1
  n0<-max(y0)
  n1<-max(y1)
  #
  if(missing(p0)){
    p.min<-max(0,-delta)
    p.max<-min(1,1-delta)
    p0<-p.min+(p.max-p.min)*grid(B=100,theta=.1)
  }
  K<-length(p0)
  p<-NULL
  tail<-(obj$P<=alpha)
  if(sum(tail)==0){sig<-0*p0}
  if(sum(tail)>0){
    for(k in 1:K){
      p<-cbind(p,stats::dbinom(y1[tail],n1,p0[k]+delta)*stats::dbinom(y0[tail],n0,p0[k]))
    }
    sig<-apply(p,2,sum)
  }
  error<-error.exceed<-100*(sig-alpha)/alpha
  error.exceed[error<0]<-0
  #
  list(x=p0,y=sig)
}


#' Likelihood ratio statistics for one-sided tests of risk difference.
#'
#' Calculates all possible values of the (signed) LR statistic for testing
#' p1-p0 greater than a provided null value psi.
#'
#' If Y0,Y1 is not supplied then all possible values are output. If specific
#' values for Y0, Y1 are supplied then only these outcomes of the SRLR are
#' calculated.
#'
#' @param n0 control sample size
#' @param n1 treatment sample size
#' @param Y0 number of successes for control (see details)
#' @param Y1 number of successes for treatment (see details)
#' @param psi null value of risk difference p1-p0
#' @param sided (1 or 2 sided test)
#' #' @param sided 1-sided or 2-sided test
#' @param mod A very small adjustment to account for 0*log0 in certain
#' likelihood calculations. Should not need adjustment.
#' @param dec.places number of decimal places of t-values and p-value.
#' @return A list with elements \item{y0,y1}{data sets (scalar or vector)}
#' \item{T, P}{quasi-exact P-value(s) and equivalent T-value(s)}
#' \item{pmle}{profile ML estimates of baseline probability}
#' \item{index}{consistent code to select a single outcome} \item{psi}{scalar
#' null value of p1-p0}
#' @author Chris J. Lloyd
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to
#' size and power constraints. Preprint. \doi{10.13140/RG.2.2.11828.94085}
#' @importFrom stats pnorm
#' @keywords internal
lr.rd=function(n0, n1,Y0=NULL,Y1=NULL,psi=0,sided=1,mod=.0000001,dec.places=10){
  # FUNCTION COMPUTES LR STATISTICS FOR TESTING THE RISK DIFFERENCE IN A 2X2 TABLE FOR ALL POSSIBLE DATA SETS.
  # THE CASE PSI=0 CORRESPONDS TO A TEST OF INDEPENDENCE.
  #
  # ARGUMENTS
  #      psi - null value of p1-p0
  #
  logit=function(x){log(x/(1-x))}
  logL<-function(y0,n0,y1,n1,p0,p1,tol=0.00000000001){
    P0<-(p0+tol)/(1+2*tol)
    P1<-(p1+tol)/(1+2*tol)
    y0*logit(P0)+y1*logit(P1)+n0*log(1-P0)+n1*log(1-P1)
  }
  profile.mle.rd=function(y0,n0,y1,n1,psi,dec.places=7){
    # Calculates restrcited MLE of p0 for fixed value of probability
    # different psi=p1-p0. Thsi involves solving a cubic.
    # Can be used for vector of psi-values for fixed data, or for vector
    # of data sets for fixed psi-value. Inputs y0 and y1 should have same length.
    #
    if(abs(max(psi))>1){warning("Value of psi is outside [-1,1]")}
    #
    zero<-0*y0*y1*psi
    # Used for dimensioning coefficients
    L3<-n1+n0+0*y0*y1+zero
    L2<-n0*(1-2*psi)+n1*(1-psi)+y0+y1+zero
    L2<-(-L2)
    L1<-y0*(1-2*psi)+y1-n1*psi-n0*psi*(1-psi)+zero
    L0<-psi*(1-psi)*y0+zero
    #
    q=L2^3/(3*L3)^3-L1*L2/(6*L3^2)+L0/(2*L3)
    p=L2^2/(3*L3)^2-L1/(3*L3)
    p=sqrt(p)*(2*(q>=0)-1)
    check=q/p^3
    check[check>1]=1
    check[check<0]=0
    a=(pi+acos(check))/3
    rmle=2*p*cos(a)-L2/(3*L3)
    rmle[rmle>1-psi]=psi
    rmle[rmle>1]=1
    rmle[rmle<0]=0
    rmle[rmle<(-psi)]=(-psi)
    rmle
  }
  #
  value <- NULL
  # All possible data sets
  y0 <- rep(0:n0, n1 + 1)
  y1 <- rep(0:n1, rep(n0 + 1, n1 + 1))
  index<-1:((n0+1)*(n1+1))
   # Optional subset of all possible data sets
  if(!missing(Y0)){
    index<-index[y0==Y0]
    y1<-y1[y0==Y0]
    y0<-y0[y0==Y0]
  }
  if(!missing(Y1)){
    index<-index[y1==Y1]
    y0<-y0[y1==Y1]
    y1<-y1[y1==Y1]
  }
  #
  value$y0 <- y0
  value$y1 <- y1
  value$T<-0*y0
  #
  y1[y1==0]<-mod
  y1[y1==n1]<-n1-mod
  y0[y0==0]<-mod
  y0[y0==n0]<-n0-mod
  #
  p1<-y1/(n1)
  p0<-y0/(n0)
  p1[p1<0]<-0
  p0[p0>1]<-1
  #
  p0.p<-profile.mle.rd(y0,n0,y1,n1,psi=psi,dec.places=dec.places+2)
  #    p0.p<-(p0.p+mod)-2*mod*(p0.p-a)/(b-a)
  p1.p<-p0.p+psi
  #
  lmax<-logL(y0,n0,y1,n1,p0,p1)
  lpart<-logL(y0,n0,y1,n1,p0.p,p1.p)
  value$T<-sqrt(2*abs(lmax-lpart))*sign(p1-p0-psi)^(sided)
  value$T<-round(10^dec.places*value$T)/10^dec.places
  value$P<-sided*stats::pnorm(-value$T)
  value$P<-round(10^dec.places*value$P)/10^dec.places
  value$index<-index
  value$pmle<-p0.p
  value$psi<-psi
  value
}


#' Calculate E-P-values based on LR statistic.
#'
#' Calculates all possible values of the E-P-value for testing p1-p0 greater
#' than a provided null value psi.
#'
#' This function can take a long time for larger sample sizes. The computation
#' time is of order (n0*n1)^2.
#'
#' @param n0 control sample size
#' @param n1 treatment sample size
#' @param psi null value of risk difference p1-p0
#' @param J index of single data set if desired
#' @param sided (1 or 2 sided test)
#' @param dec.places number of decimal places of output t-values and p-values.
#' @param prin outputs expected time and progress of calculation
#' @return A list with elements \item{y0,y1}{data sets (scalar or vector)}
#' \item{oldP}{approximate p-value(s) before E-step} \item{T, P}{quasi-exact
#' P-value(s) and equivalent T-value(s)} \item{pmle}{profile ML estimates of
#' baseline probability} \item{index}{consistent code to select a single
#' outcome} \item{psi}{scalar null value of p1-p0}
#' @author Chris J. Lloyd
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to
#' size and power constraints. Preprint. \doi{10.13140/RG.2.2.11828.94085}
#' @importFrom stats dbinom qnorm
#' @keywords internal
ESTEP.rd=function(n0,n1,psi=0,J=NULL,sided=1,dec.places=10,prin=FALSE){
  # FUNCTION PRODUCES PLUG-IN OR ESTIMATED EXACT P-VALUES FROM APPROXIMATE P-VALUES
  #
  # ARGUMENTS:
  #        J - index of datasets to consider (default is all)
  #
  profile.mle.rd=function(y0,n0,y1,n1,psi,dec.places=7){
    if(abs(max(psi))>1){warning("Value of psi is outside [-1,1]")}
    #
    zero<-0*y0*y1*psi
    # Used for dimensioning coefficients
    L3<-n1+n0+0*y0*y1+zero
    L2<-n0*(1-2*psi)+n1*(1-psi)+y0+y1+zero
    L2<-(-L2)
    L1<-y0*(1-2*psi)+y1-n1*psi-n0*psi*(1-psi)+zero
    L0<-psi*(1-psi)*y0+zero
    #
    q=L2^3/(3*L3)^3-L1*L2/(6*L3^2)+L0/(2*L3)
    p=L2^2/(3*L3)^2-L1/(3*L3)
    p=sqrt(p)*(2*(q>=0)-1)
    check=q/p^3
    check[check>1]=1
    check[check<0]=0
    a=(pi+acos(check))/3
    rmle=2*p*cos(a)-L2/(3*L3)
    rmle[rmle>1-psi]=psi
    rmle[rmle>1]=1
    rmle[rmle<0]=0
    rmle[rmle<(-psi)]=(-psi)
    rmle
  }
  #
  stat=lr.rd(n0,n1,psi=psi,sided=sided)
  index<-1:((n0+1)*(n1+1))
  #
  phat<-profile.mle.rd(stat$y0,n0,stat$y1,n1,psi)
  obj<-function(p0,psi,n0,n1,y0,y1){
    sum(stats::dbinom(y0,n0,p0)*stats::dbinom(y1,n1,p0+psi))}
  if(missing(J)){J<-index}
  #
  flag<-is.null(stat$T)
  # Calculate estimated Pvalues
  est.P<-NULL
  how.often=round(sqrt(max(index)))*5
  if(prin){e.time=2+round(exp(-3.03183+(n0+n1)*0.0282))
  print(paste("Expected computation time is",e.time,"seconds."))
  }
  for(j in J){
    if(prin){
      if(j==how.often*round(j/how.often)){print(round(100*j/max(index)))}
      }
    if(flag){est.P<-c(est.P,obj(phat[j],psi,n0,n1,stat$y0[stat$P<=stat$P[j]],stat$y1[stat$P<=stat$P[j]]))}
    if(!flag){est.P<-c(est.P,obj(phat[j],psi,n0,n1,stat$y0[stat$T>=stat$T[j]],stat$y1[stat$T>=stat$T[j]]))}
  }
  # Output as list
  tol<-10^(-dec.places)
  list(y0=stat$y0[J],y1=stat$y1[J],oldP=stat$P[J],
       P=round(10^dec.places*est.P)/10^dec.places,T=-stats::qnorm((est.P+tol)/(1+2*tol)),pmle=phat[J],index=index[J],psi=psi)
}


#' Approximate and exact tests for single data set
#'
#' For a single provided data sets (y0, n0, y1, n1) calculate approximate and quasi-exact test statistics of the null value p1-p0=psi.
#'
#'
#' @param y0 number of successes for control
#' @param n0 control sample size
#' @param n1 treatment sample size
#' @param y1 number of successes for treatment
#' @param psi null value of risk difference p1-p0
#' @param sided 1-sided or 2-sided test
#' @param dec.places decimal places of output T and P values
#' @return List with elements \item{y}{the data vector y0,y1} \item{n}{the
#' sample size vector n0,n1} \item{approx}{approximate T-value and p-value}
#' \item{quasi.exact}{quasi.exact T-value and p-value} \item{pmle}{profile
#' maximum likelihood estimate of baseline probability} \item{psi}{null value
#' of p1-p0}
#' @author Chris J. Lloyd
#' @references C.J. Lloyd (2022) Exact samples sizes for clinical trials subject to
#' size and power constraints. Preprint. \doi{10.13140/RG.2.2.11828.94085}
#' @examples
#'
#' y0=25
#' y1=41
#' n0=65
#' n1=82
#' # Non-inferiority test of p1-p0>-0.1. Evidence is strong.
#' inference(y0,n0,y1,n1,psi=-0.1)
#'
#'
#' @export inference
inference=function(y0,n0,y1,n1,psi=0,sided=1,dec.places=4){
  lr=lr.rd(n0,n1,y0,y1,psi,sided=sided)
  elr=ESTEP.rd(n0,n1,psi,lr$index)
  approx=c(lr$P,lr$T)
  names(approx)=c("p-value","t-value")
  quasi.exact=c(elr$P,elr$T)
  names(quasi.exact)=c("p-value","t-value")
  #
  list(y=c(y0,y1),n=c(n0,n1),approx=round(10^dec.places*approx)/10^dec.places,
       quasi.exact=round(10^dec.places*quasi.exact)/10^dec.places,pmle=lr$pmle,psi=psi)
}

