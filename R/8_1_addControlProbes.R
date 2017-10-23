# inner function #
getControlProbes <- function(outputDir, alignmentType)
{
  stopifnot(alignmentType %in% c("cdf", "pgf"))
  if (alignmentType == "pgf")
  {
    stopifnot(file.exists(file.path(
      outputDir, "tmpProbesTableDir", "ProbesPGF.csv")))
    probesPGF <- read.csv(file.path(outputDir, 
                                      "tmpProbesTableDir/ProbesPGF.csv"), 
                            header = TRUE )
    NonMainProbes <- subset(probesPGF, pstype != "main")
    NonMainProbes2 <- NonMainProbes[ ,c("fid", "man_fsetid")]
    if(any(duplicated(NonMainProbes2)))
    {NonMainProbes2 = NonMainProbes2[!duplicated(NonMainProbes2), ]}
    return(NonMainProbes)
  }
  
  if (alignmentType == "cdf")
  {
    stopifnot(file.exists(file.path(
      outputDir, "tmpProbesTableDir", "ProbesCDF.csv")))
    probesCDF <- read.csv(file.path(outputDir, 
                                      "tmpProbesTableDir/ProbesCDF.csv"), 
                            header = TRUE )
    ControlProbes <- probesCDF[grep("AFFX", as.vector(probesCDF$man_fsetid)), ]
    ControlProbes2 <- ControlProbes[ ,c("fid", "man_fsetid")]
    if(any(duplicated(ControlProbes2)))
    {ControlProbes2 = ControlProbes2[!duplicated(ControlProbes2), ]}
    return(ControlProbes)
  }
}
# main function #

addControlProbes <- function(ProbeSetsClean, outputDir, alignmentType)
{
  stopifnot(setequal(colnames(ProbeSetsClean), c("probe_id", 
                                                 "transcript_id", 
                                                 "gene_id")))
  ConPro <- getControlProbes(outputDir, alignmentType)
  ConProSets <- data.frame(probe_id = as.vector(ConPro$fid),
                           transcript_id = as.vector(ConPro$man_fsetid), 
                           gene_id = as.vector(ConPro$man_fsetid))
  # res = rbind(ConProSets, ProbeSetsClean) 
  # will need contrls separatelly, so return them, not all together!
  return(ConPro)
}