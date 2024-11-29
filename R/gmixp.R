#' @export
gmixp<-function(data, s, n,a,parameters,  conf=0.95){
  if(!is.numeric(data)|| is.character(data) || is.matrix(data))
    stop("data must be a numeric vector")
  if(!is.numeric(s)|| is.character(s) || is.matrix(s) )
    stop("s must be a numeric vector")
  if(length(s) != 1 )
    stop("the length of s must be equal to 1")
  if(!is.numeric(n)|| is.character(n) || is.matrix(n) )
    stop("n must be a numeric vector")
  if( length(data)>s) { stop("s must be greater than data length")}
  if( length(data)>n) { stop("n must be greater than data length")}
  if( s>n) { stop("n must be greater than s")}
  if(!is.numeric(parameters)|| is.character(parameters) || is.matrix(parameters))
    stop("parameters must be a numeric vector")
  if(length(parameters) == 0 )
    stop("the length of parameters must be greater than or equal to 1")

usolve=function(r,s){
  comb=function(l,m){
    return( factorial(l) / (factorial(m) * factorial(l-m)) )
  }
  j=0:(s-r-1)
  S=numeric(1)
  a=factorial(n)/(factorial(r-1)*factorial(n-s)*factorial(s-r-1))
  f=function(x){
    for(i in 0:(r-1)){
      S=S+sum((-1)^(i+j)*comb(r-1,i)*comb(s-r-1,j)*((n-s+j+1)*((n-r+i+1)+(x*(n-s+j+1))))^(-1))}
    return((1- conf)-a*S)}
  uniroot(f,c(0,1),extendInt = "yes")
}
r<-length(data)
  u=usolve(r,s)$root
  l1=parameters[1]
  b1=parameters[2]
  l2=parameters[3]
  b2=parameters[4]
    pr=runif(n,0,.99)
    x=sort(data)
    l=x[r]
    m=x[s]
    pmix=a*(Igamma(l1,x[r]/b1,lower = T))/Cgamma(l1)+(1-a)*(Igamma(l2,x[r]/b2,lower = T))/Cgamma(l2)
    p=1-((1-pmix)^(u+1))
    fun=function(t){
      a*Igamma(l1,(b2/b1)*Igamma.inv(l2,t*Cgamma(l2)))/Cgamma(l1)+(1-a)*t-p
    }
    uni=uniroot(fun,c(0,1),extendInt = "yes")
    slv=uni$root
    x_s=b2*Igamma.inv(l2,slv*Cgamma(l2),lower = T)
    z=x_s
  lower=l
  xvalue=m
  upper=z
  if(lower>= upper){stop("Lower bound can not be greater than upper bound")}
  xest=(lower+upper)/2

  interval<-c(lower,upper)
  point<- xest
  names(interval)<-c("lower","upper")
  names(point)<-c("predicted point")
  interval<-interval[c("lower", "upper")]
  point<-point[c("predicted point")]
  int<-list(interval=interval, point=point, lower=lower, upper=upper)
  return(structure(int, class = "gmixp"))

}
#' @export
print.gmixp <- function(x, ...) {
  if (!inherits(x, "gmixp"))
    stop("Use only with 'gmixp' objects")
  cat("Prediction  for the future observation based on mixture of two gamma distribution \n")
  cat("Point:\n")
  print(cbind.data.frame("predicted point(mean of interval)" = x$point), ...)
  cat("Interval:\n")
  print(cbind.data.frame("PI" = x$interval), ...)
}
