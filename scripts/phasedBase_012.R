# From McKinney et al., 2021: 10.1111/mec.15712
phasedBase_012<-function(bases){
  orgAlleles<-unique(bases)
  #create named vector to convert
  newAlleles<-c("0","1")
  names(newAlleles)<-orgAlleles
  #convert bases to numbers
  numBases<-str_replace_all(bases,newAlleles)
  #swap 0 and 1 if necessary so 0 is most common
  if(sum(numBases==0)>sum(numBases==1)){
    numBases<-chartr("01","10",numBases)
  }
  numBases<-as.numeric(numBases)
  return(numBases)
}