# Data-driven analytical algorithm for basophil activation data
# Written by: Sarita Patil, MD. Function published 6/8/2017

#' BATgate: An algorithm for automated flow cytometry analysis of basophil activation data
#'
#' @return This function allows for data-driven flow cytometry analysis of basophil activation data
#' 
#' @author Sarita Patil, \email{sarita.patil@@mgh.harvard.edu}
#' 
#' @param wkdir The location of the working directory, where all of the .fcs files for a BAT experiment are located
#' @param fluorophores List of fluorophores, ordered the same as in the fset
#' @param SSCno The number of clusters within the SSC distribution
#' @param medcontroltube The sample number containing medium only or the control condition in the BAT assay,
#'        used to calculate the MFI for definition of CD63hi basophils
#' 
#' @return A csv file containing basophil flow cytometry output including percentage of CD63hi basophils,
#'        A pdf file containing graphical gating of flow cytometry analysis for visual inspection
#' 
#' @keywords models, cluster
#' 
#' @import flowCore
#' @import flowQ
#' @import flowViz
#' @import flowStats
#' @import flowUtils
#' @import flowClust
#' @import Biobase
#' @import xtable
#' @import knitr
#' @examples
#' BATgate(wkdir, c("Time","FSC-A","FSC-H","SSC-A","CD63","CCR3"), 2, 1)
#' BATgate()
#' @export
BATgate <- function(workingdir, fluorophores, SSCno, medcontroltube) {
  setwd(workingdir)
  file.list <- (list.files())
  for(i in 1:length(file.list)) {
    full.file.path <- paste(workingdir, "/", file.list[i], sep="")
    ffiles <- list.files(full.file.path, pattern = "*.fcs")
    fset <- read.flowSet(ffiles, path = full.file.path, name.keyword="TUBE NAME",phenoData=list(name="TUBE NAME",Filename="$FIL"))
    
    #Rename the column fluorophores by the markers, except for FSC/SSC
    colnames(fset) <- fluorophores
    
    # transform the entire flowSet
    atrans<-logicleTransform()
    fsettransform<-transform(fset,`CD63`=atrans(`CD63`),`CCR3`=atrans(`CCR3`))
    
    # Remove boundary events
    bf <- boundaryFilter(filterId="myBoundaryFilter", x=c("FSC-H", "SSC-A"), side="upper")
    fsettransform <- Subset(fsettransform, bf)
    
    # SSC-A lo gate (depends on the number of groupings in SSC, depending on flow acquisition)
    if(SSCno=="2") { 
    kmSSCA<-split(fsettransform,kmeansFilter("SSC-A"=c("Low","High"),filterId="SSCAkmeans"))
    SSCAhi <- exprs(kmSSCA$High[[1]])[,"SSC-A"]
    SSCAlo <- exprs(kmSSCA$Low[[1]])[,"SSC-A"]
    SSCAline <- ((min(SSCAhi)) + (max(SSCAlo1))/2)
    SSCAgate <- rectangleGate(`SSC-A`=c(0, SSCAline), filterId="SSCAlo")
    fsettransform <-Subset(fsettransform, SSCAgate)
    }
    if(SSCno=="3"){
    kmSSCA<-split(fsettransform,kmeansFilter("SSC-A"=c("Low","Med","High"),filterId="SSCAkmeans"))
    SSCAhi <- exprs(kmSSCA$Med[[1]])[,"SSC-A"]
    SSCAlo <- exprs(kmSSCA$Low[[1]])[,"SSC-A"]
    SSCAline <- ((min(SSCAhi)) + (max(SSCAlo))/2)
    SSCAgate <- rectangleGate(`SSC-A`=c(0, SSCAline), filterId="SSCAlo")
    fsettransform <-Subset(fsettransform, SSCAgate)  
    }else{print("Number of SSC groups is invalid")}
   
    # Additional CCR3 high non-baso gate: Make 2 different gates 
    basonormscale<- norm2Filter("CD63", "CCR3", scale.factor=5.5)
    
    fsetnonbaso <- Subset(fsettransform, basonormscale)
    CCR3line <- max(exprs(fsetnonbaso[[1]])[,"CCR3"])
    CCR3gate <- rectangleGate(`CCR3`=c(CCR3line, 5), filterId="CCR3hi")
    fsetbaso <-Subset(fsettransform, CCR3gate)
    
    # Narrow down the baso gate using a Norm2filter distribution
    baso2normscale2<- norm2Filter("SSC-A", "CCR3", scale.factor=2)
    
    #followed by subsetting the filter:
    fsetbason<-Subset(fsetbaso, baso2normscale2)
    
    # CD63 GATING as activation marker  
    # Gating strategy: By MichealBulhman medium only strategy gate setting on 2-2.5%
    # Predefined: First tube of experiment is Medium only or Background control
    # Know that [[1]] of the fset is always the Background control
    
    #Find the cutoff above which there is 2.5% in the medium only condition
    medcontroltube <- as.numeric(medcontroltube)
    medn <- exprs(fsetbason[[medcontroltube]])[,"CD63"]
    
    #Define & subset the medium-only cutoff defined subset for CD63hi
    #Define the filter as ms
    ms <- qnorm(0.975,mean=mean(medn),sd=sd(medn))
    msgate <- rectangleGate(`CD63`=c(ms, Inf), filterId="CD63med")
    fsetbasomed<-Subset(fsetbason, msgate)
    
    # Define the characteristics to output in table form
    # Allevents = all the events collected in the flow file, after boundary gating
    # Totalnorm = all basophils identified
    # CD63himed = number of CD63hi basophils identified
    # Make matrices of data types (median = MFI, mean= mean fl, sd = standard deviation)
    # CV is calculated as standard deviation/mean
    # CD63hi_perc
    Allevents <- as.numeric(fsApply(fsettransform, nrow, use.exprs = TRUE))
    Totalnorm <- as.numeric(fsApply(fsetbason, nrow, use.exprs = TRUE))
    CD63himed <- as.numeric(fsApply(fsetbasomed, nrow, use.exprs = TRUE))
    d1 <- data.frame(fsApply(fsetbason1, each_col, median))
    m1 <- data.frame(fsApply(fsetbason1, each_col, mean))
    s1 <- data.frame(fsApply(fsetbason1, each_col, sd))
    c1 <- s1/m1
    
    #Create a data.frame for export and subsequent analysis
    tubename <- data.frame(pData(phenoData(fset)))
    tubename <- list(tubename$name)
    
    data <- data.frame("Condition" = tubename,
                       "BasoCount"= Totalnorm,
                       "CD63hiBasoCount"= CD63himed, 
                       "SSCAMFIb"=d1$SSC.A,
                       "SSCACVb"=c1$SSC.A,
                       "CCR3MFIb"=d1$CCR3,
                       "CCR3CVb"=c1$CCR3,
                       "CD63MFIb"=d1$CD63,
                       "CD63CVb"=c1$CD63,
                       "Allevents"=Allevents)
    
    data$CD63hi_perc <- ((data[,3]/data[,2])*100)
    
    ## There are 2 outputs: a csv file with the datatable, and a pdf of the critical graphs
    
    # Final step: save the csv file in the same folder
    cname <- paste(full.file.path, "autoBAT.csv", sep="")
    write.csv(data, cname)
    
    # PDF of the graph checks through this script
    ctitle <- paste(full.file.path, "Graphs.pdf", sep="")
    pdf(ctitle)
    print(xyplot(`SSC-A` ~ `CCR3`, fsettransform, main="Setting the SSC gate"))
    print(xyplot(`CCR3` ~ `CD63`, fsettransform, filter= basonormscale,main="SSClo, km2, basoneg"))
    print(xyplot(`SSC-A` ~ `CCR3`, data= fsetbaso, filter = baso2normscale2, main="Basophil"))
    print(densityplot(~ `CD63`,data=fsetbason, refline=ms, main="CD63hi histograms"))
    dev.off() 
  }
}
  
  
