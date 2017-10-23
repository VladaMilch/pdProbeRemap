## ---- message=FALSE------------------------------------------------------
library("pdProbeRemap")

## ---- eval = FALSE-------------------------------------------------------
#  #evaluation worked
#  dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
#  outDir = "/Users/milchevs/ownCloud/ProbeRemapSteps14/OutPutDir/"
#  dir.create(outDir)
#  
#  annotation2fasta(originalAnnotationDir = dir1,
#                   organism = "Mouse", species = "Mus musculus",
#                   author = "Vladislava Milchevskaya",
#                   email = "milchv@gmail.com",
#                   outputDir = outDir,
#                   reverse = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  # evaluation worked
#  path_to_STAR = "/path/to/STAR"
#  sequence_fasta = "/path/to/sequence.fasta"
#  genomeDir = "/path/to/reference/genome/"
#  organism = "Human"
#  outfile = "output_file."
#  
#  star.command.make(path_to_STAR, sequence_fasta, genomeDir, organism, outfile)

## ---- eval =  FALSE------------------------------------------------------
#  # evaluation worked
#  outputDirAnotationDro <- "/Users/milchevs/Downloads/NewPackages/DroGene1_1"
#  
#  require("refGenome")
#  beg <- ensemblGenome()
#  basedir(beg) <- "/Users/milchevs/databases/sequences/Ensembl/annotation/Dmel/"
#  ens_gtf <- "dmel-all-r6.09.gtf" # adjusted file, excessive white spaces deleted
#  read.gtf(beg, ens_gtf)
#  #tableAttributeTypes(beg)
#  #moveAttributes(beg,c("gene_name","transcript_id","transcript_name","exon_number"))
#  
#  ANN_all.list <- (as.list(beg@ev))
#  ANN_all <- ANN_all.list$gtf
#  head(ANN_all)
#  save(ANN_all, file = file.path(outputDirAnotationDro, paste0(ens_gtf, ".AnnotationDataFrame.RData")))

## ---- eval = FALSE-------------------------------------------------------
#  alignmentDir = "/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/MoGene2/"
#  load("/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/Annotations/AnnotationDataFrame.Mus_musculus.GRCm38.82.RData")
#  dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
#  outDir = "/Users/milchevs/ownCloud/ProbeRemapSteps14/OutPutDir/"
#  o1 = makeOriginalPackageObject(originalAnnotationDir = dir1,
#                                 organism = "Mouse",
#                                 species = "Mus musculus",
#                                 outputDir = outDir )
#  
#  # NEW_ParsedData <- alignments2parsedData <- function(alignment.dir,
#  #                                   Annotation,
#  #                                   outputDir,
#  #                                   level,
#  #                                   min_probe_number,
#  #                                   package_seed)
#  
#  
#  library("pdProbeRemap")
#  data(Alignments_class_example, package = "pdProbeRemap")
#  data(Annotation_example)
#  data(seed_example)
#  alignmentDir <-
#   unlist(strsplit(dir(system.file("extdata",package="pdProbeRemap"),
#                    pattern="example_drosophila.Aligned.out.sam",
#                    full.names=TRUE), split = "example_drosophila.Aligned.out.sam"))[1]
#  
#  
#  makeNewAnnotationPackage(alignment.dir = alignmentDir,
#                           Annotation = Annotation_example,
#                           outputDir = outDir,
#                           #outputDir = ".",
#                           level = "gene",
#                           min_probe_number = 1,
#                           pkgNameSUFFIX = ".example")
#  

## ---- eval = FALSE-------------------------------------------------------
#  library("pdProbeRemap")
#  
#  alignmentDir = "/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/MoGene2/"
#  load("/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/Annotations/AnnotationDataFrame.Mus_musculus.GRCm38.82.RData")
#  dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
#  outDir = "/Users/milchevs/ownCloud/ProbeRemapSteps14/OutPutDir/Probe_Sequences_Fasta/"
#  

## ---- eval = FALSE-------------------------------------------------------
#  o1 = pdProbeRemap::makeOriginalPackageObject(originalAnnotationDir = dir1,
#                                               organism = "Mouse",
#                                               species = "Mus musculus",
#                                               outputDir = outDir )

## ---- eval = FALSE-------------------------------------------------------
#  NEW_PD = makeNewAnnotationPackage(alignment.dir = alignmentDir,
#                                    Annotation = ANN_MUS_all,
#                                    outputDir = outDir,
#                                    level = "gene",
#                                    min_probe_number = 1,
#                                    #package_seed = o1,
#                                    quiet = FALSE,
#                                    pkgNameSUFFIX = ".example")
#  #

## ---- eval = FALSE-------------------------------------------------------
#  CelFilesPath = c("/path/to/celfile1.CEL", "/path/to/celfile2.CEL")
#  celfiles = read.celfiles(filenames=CelFilesPath,
#                                pkgname = "pd.mogene.2.0.st")

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("/path/to/new/pd/annotation/package/pd.NEWPACKAGE",
#                   repos=NULL, type="source")
#  library("pd.NEWPACKAGE")
#  
#  CelFilesPath = c("/path/to/celfile1.CEL", "/path/to/celfile2.CEL")
#  celfiles = read.celfiles(filenames=CelFilesPath,
#                                pkgname = "pd.NEWPACKAGE")

## ---- eval = FALSE-------------------------------------------------------
#  library("pdProbeRemap")
#  
#  outDir <- "/Users/milchevs/ownCloud/ProbeRemapStep3/R/testdata/DroGene1_1"
#  OriginalAnnotDir <-  "/Users/milchevs/Downloads/CD_DrosGenome1\ 2/Full/DrosGenome1/LibFiles/"
#  
#  FastaFilePath_dro <- annotation2fasta(originalAnnotationDir = OriginalAnnotDir,
#                                        organism = "Dmel",
#                                        species = "Drosophila melanogaster",
#                                        outputDir = outDir,
#                                        author = "Name Surname",
#                                        email = "e@mail",
#                                        reverse = FALSE )

## ---- eval = FALSE-------------------------------------------------------
#  ######### step 2 ###########
#  #        alignment         #
#  star_command_plainprobes =
#    star.command.make(sequence_fasta = "/var/local/milchevs/ALIGNED/DroGene1/FASTA/DrosGenome1_probe_sequences_to_align.fasta", # the path is different because for helios!
#                      genomeDir = "/var/local/milchevs/databases/sequences/Ensembl/star/Dmel/",
#                      species = "Dmel",
#                      outfile = "/var/local/milchevs/ALIGNED/DroGene1/DrosGenome1.",
#                      runThreadN = 6)
#  
#  ######### step 3 ###########
#  load("/Users/milchevs/Downloads/NewPackages/DroGene1_1/dmel-all-r6.09.gtf.AnnotationDataFrame.RData") # Annotation
#  alignmentDir = file.path("/Users/milchevs/ownCloud/ProbeRemapStep3/R/testdata/DroGene1_1/", "Alignments")
#  
#  makeNewAnnotationPackage(alignment.dir = alignmentDir,
#                           Annotation = ANN_all,
#                           outputDir = outDir,
#                           level = "gene",
#                           min_probe_number = 1,
#                           quiet=FALSE,
#                           pkgNameSUFFIX = ".NewAnnotation"
#  )
#  
#  install.packages("/Users/milchevs/ownCloud/ProbeRemapStep3/R/testdata/DroGene1_1/pd.drosgenome1.NewAnnotation/", repos = NULL, type = "source")

