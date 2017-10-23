##############################################################################
#################### get annotation files from directory #####################
##############################################################################

getAnnotationFilesFromDir <- function(annotation.dir)
{
    # check annotation directory exists #
    if(!dir.exists(annotation.dir))
    {stop(message = "Wrong path to the annotation directory! 
                Input existing directory.\n")}
    
    if(!(
        (    length(list.files(annotation.dir, pattern = ".pgf", 
                               ignore.case = TRUE,
                               full.names = TRUE))==1 &
             length(list.files(annotation.dir, pattern = ".cdf", 
                               ignore.case = TRUE,
                               full.names = TRUE))==0    ) |
        (    length(list.files(annotation.dir, 
                               pattern = ".pgf", 
                               ignore.case = TRUE,
                               full.names = TRUE))==0 &
             length(list.files(annotation.dir, pattern = ".cdf", 
                               ignore.case = TRUE,
                               full.names = TRUE))==1    ))
    )
    {
        stop(message = "Check annotation directory: must be one and only one file 
                 with extention '.CDF' or '.PGF'. \n")
    }
    
    if (    length(list.files(annotation.dir, 
                              pattern = ".pgf", 
                              ignore.case = TRUE,
                              full.names = TRUE))==1    )
    {
        (pgf <- list.files(annotation.dir, 
                           pattern = ".pgf", 
                           ignore.case = TRUE,
                           full.names = TRUE))
        (clf <- list.files(annotation.dir, 
                           pattern = ".clf",
                           ignore.case = TRUE,
                           full.names = TRUE))
        (prob <- list.files(annotation.dir, 
                            pattern = ".probeset.csv",
                            ignore.case = TRUE,
                            full.names = TRUE))
        mps <- list.files(annotation.dir, 
                          pattern = "mps$", 
                          ignore.case = TRUE,
                          full.names = TRUE)
        trans <- list.files(annotation.dir, 
                            pattern="transcript.file.csv|transcript.csv",
                            ignore.case = TRUE,
                            full.names=TRUE)
        
        if (! (length(pgf) == 1 &
                     length(clf) == 1 & 
                     length(prob) == 1 & 
                     length(mps) == 1 & 
                     length(trans) == 1) )
        {
            stop(message = "Check that all annotation files are present in the 
                     annotation directory, and each file is present only once! 
                     (PGF-based annotation) \n")
        }
        res = list(pgf, clf, prob, mps, trans)
        names(res) <- c("pgf", "clf", "prob", "mps", "trans")
        return(res)
    }
    
    if(    length(list.files(annotation.dir, 
                             pattern = ".cdf", 
                             ignore.case = TRUE,
                             full.names = TRUE))==1     )
    {
        (cdf <- list.files( annotation.dir, 
                            pattern = "cdf", 
                            ignore.case = TRUE,
                            full.names = TRUE))
        (cel <- list.files( annotation.dir, pattern = ".CEL", ignore.case = TRUE,
                                                full.names = TRUE))
        (tab <- list.files(annotation.dir, pattern = "probe_tab", 
                                             ignore.case = TRUE,
                                             full.names = TRUE))
        if(!( length(cdf) == 1 &
                    length(cel) == 1 &
                    length(tab) == 1)
        )
        {stop(message = 
                        "Check annotation files in the annotation directory: 
                    for CDF-based annotation there must be 1 .CDF file, 
                    1 .CEL file and tab-delimited probe file 'probe_tab'.\n")
        }
        res = list(cdf, cel, tab)
        names(res) <- c("cdf", "cel", "tab")
        return(res)
        }
    
}

## for tests later ##
# dir1 = 
# "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
# dir2 <-    "/Users/milchevs/Downloads/CD_DrosGenome1\ 2/Full/DrosGenome1/LibFiles/"
# getAnnotationFilesFromDir(dir1)
# getAnnotationFilesFromDir(dir2)

