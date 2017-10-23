# THESE TESTS WILL ONLY WORK ON MY LOCAL MACHINE #
#
# alignmentDir = file.path("/Users/milchevs/ownCloud/ProbeRemapStep3/R/testdata/DroGene1_1/", "Alignments")
# # alignmentDir = "/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/MoGene2/"
# load("/Users/milchevs/Downloads/NewPackages/DroGene1_1/dmel-all-r6.09.gtf.AnnotationDataFrame.RData") # Annotation 
# # load("/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/Annotations/AnnotationDataFrame.Mus_musculus.GRCm38.82.RData")
# # load the alignments 
# options(stringsAsFactors = FALSE)
# 
# 
# 
# #load("8apr2016.RData")
# 
# ###########  PGF ###############################################
# alignmentDir = "/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/MoGene2/"
# load("/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/Annotations/AnnotationDataFrame.Mus_musculus.GRCm38.82.RData")
# dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
# outDir = "./tmpDRO/"
# o1 = makeOriginalPackageObject(originalAnnotationDir = dir1, organism = "Mouse", species = "Mus musculus", outputDir = outDir )
# 
# undebug(alignments2parsedData)
# 
# NEW_PD = alignments2parsedData(alignment.dir = alignmentDir, 
#                                Annotation = ANN_MUS_all,
#                                outputDir = outDir,
#                                level = "gene",
#                                min_probe_number = 1,
#                                package_seed = o1)
# 
# test_that("New Parsed Data corresponds to the array format (PGF)",
#           {
#             if(o1[[3]] == "pgf")
#             {
#               expect_true(length(NEW_PD) == 10)
#               expect_true(setequal(names(NEW_PD), c("featureSet", "pmFeatures", "mmFeatures", "geometry",   "pmSequence",
#                                                     "mmSequence", "chrom_dict", "level_dict", "type_dict",  "core")))
#             }
#           })
# 
####
# test_that("Error when wrong inputs",
#           {
#             expect_error(alignments2parsedData(alignment.dir = tempdir(), # dir does not contain what is needed
#                                                Annotation = ANN_MUS_all,
#                                                outputDir = outDir,
#                                                level = "gene",
#                                                min_probe_number = 1,
#                                                package_seed = o1), regexp = "Wrong alignment directory")
#             # generate wrong annotation table #
#             Annot_wrong = subset(ANN_MUS_all, feature == "exon")
#             expect_error(alignments2parsedData(alignment.dir = alignmentDir,
#                                                Annotation = Annot_wrong, # Annotation does not contain feature like transcript or RNA
#                                                outputDir = outDir,
#                                                level = "gene",
#                                                min_probe_number = 1,
#                                                package_seed = o1), regexp = "transcript|RNA")
# 
#             expect_error(alignments2parsedData(alignment.dir = alignmentDir,
#                                                Annotation = ANN_MUS_all,
#                                                outputDir = "./outDir/nonexisting/path/", #non existing path to output directory
#                                                level = "gene",
#                                                min_probe_number = 1,
#                                                package_seed = o1), regexp = "outputDir")
# 
#             expect_error(alignments2parsedData(alignment.dir = alignmentDir,
#                                                Annotation = ANN_MUS_all,
#                                                outputDir = outDir,
#                                                level = "ggene", # typo in level
#                                                min_probe_number = 1,
#                                                package_seed = o1), regexp = "level")
# 
#             expect_error(alignments2parsedData(alignment.dir = alignmentDir,
#                                                Annotation = ANN_MUS_all,
#                                                outputDir = outDir,
#                                                level = "gene",
#                                                min_probe_number = 0,
#                                                package_seed = o1), regexp = "min_probe_number")
# 
#             expect_error(alignments2parsedData(alignment.dir = alignmentDir,
#                                                Annotation = ANN_MUS_all,
#                                                outputDir = outDir,
#                                                level = "gene",
#                                                min_probe_number = 1,
#                                                package_seed = "?"), # wrong package_seed input
#                          regexp = "package_seed")
# 
# 
#           })