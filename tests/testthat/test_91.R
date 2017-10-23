# THESE TESTS WILL ONLY WORK ON MY LOCAL MACHINE #


# options(stringsAsFactors = FALSE)
# ##########################################################
# ############# load or generate needed input objects ######
# # ##########################################################
# # alignmentDir = "/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/MoGene2/"
# # load("/Users/milchevs/Documents/Biology/PROJECTS/Paul_Bertone/ReANNot/Annotations/AnnotationDataFrame.Mus_musculus.GRCm38.82.RData")
# # dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
# # outDir = "./tmpDRO/"
# # o1 = makeOriginalPackageObject(originalAnnotationDir = dir1, organism = "Mouse", species = "Mus musculus", outputDir = outDir )
# # 
# # Annotation = ANN_MUS_all
# # outputDir = outDir
# # level = "gene"
# # min_probe_number = 1
# # package_seed = o1
# # 
# # ParsedData <- package_seed$parsedDataOLD
# # Als <- processAlignments(alignment.dir = alignmentDir)
# # PS_mus = alignments2probesetsRaw(Alignments = Als, Annotation = Annotation, level = level)
# # PS_mus_cleaned <- cleanProbeSets(ProbeSets = PS_mus,
# #                                  level = level,
# #                                  min_probe_number = min_probe_number)
# # 
# # PS_controls <- addControlProbes(ProbeSetsClean = PS_mus_cleaned, outputDir = outputDir, alignmentType = package_seed[[3]])
# # shift = max(as.numeric(PS_controls$fsetid))
# # PSD <- make_objects_for_reannotation(Probesets = PS_mus_cleaned, level = level, Annotation = Annotation, shift = shift)
# # 
# # rm(ANN_MUS_all)
# # save.image("./data/PGF_mouse_allData_to_test_91.Rda")
# 
# # ######## PGF ########
# rm(list = setdiff(ls(), lsf.str()))
# 
# load("./data/PGF_mouse_allData_to_test_91.Rda")
# library(devtools)
# library(roxygen2)
# library(testthat)
# library("pdProbeRemap")
# load_all("pdProbeRemap")
# 
# PDnew <- build_new_ParsedDataOld_object_PGF(level = level, ParsedDataOld = ParsedData, PS_controls = PS_controls, Annotation = Annotation, ProbeSetData = PSD)
# debug(build_new_ParsedDataOld_object_PGF)
# 
# test_that("test structure of the new Parsed Data object", 
#           {
#             expect_true(setequal(names(PDnew), names(ParsedData)))
#             expect_is(PDnew$featureSet, "data.frame")
#             expect_is(PDnew$pmFeatures, "data.frame")
#             expect_is(PDnew$mmFeatures, "data.frame")
#             expect_is(PDnew$geometry, "character")
#             expect_is(PDnew$pmSequence, "DataFrame")
#             expect_true(class(PDnew$mmSequence) == class(ParsedData$mmSequence))
#             expect_is(PDnew$chrom_dict, "data.frame")
#             expect_is(PDnew$level_dict, "data.frame")
#             expect_is(PDnew$type_dict, "data.frame")
#             expect_is(PDnew$core, "data.frame")
#             })
# 
# 
# test_that("test featureSet", 
#           {
#             expect_true(all(PS_controls$fsetid %in% PDnew$featureSet$fsetid)) # all control probe sets are there
#             expect_true(all(PSD$ProbeSetDataAnnotation$fsetid %in% PDnew$featureSet$fsetid)) # all new main-type probe sets are there
#           })
# 
# 
# test_that("test pmFeatures",
#           {
#             mm_pm = merge(ParsedData$pmFeatures[ ,c("fid", "atom", "x", "y")], PDnew$pmFeatures[ ,c("fid","atom","x", "y")])
#             expect_true( all(PDnew$pmFeatures$pmFeatures$fid %in% mm_pm$fid )) # no new probe fid appeared
#             expect_true( all(PS_controls$fid %in% PDnew$pmFeatures$fid)) # control probes are in the table
#             expect_true(  all(PSD$ProbeSetDataFrame$fid %in% PDnew$pmFeatures$fid)  ) # all probes from new probe sets
#             
#             expect_true(setequal(PDnew$pmFeatures$fsetid, union(PSD$ProbeSetDataFrame$fsetid, PS_controls$fsetid))) # all and only control + main new probe set ids are there
#           })
# 
# test_that("test geometry: did not change",
#           {expect_equal(PDnew$geometry, ParsedData$geometry)})
# 
# test_that("test mmFeatures: does not contain pm probes",
#           {expect_true(!any(PDnew$pmFeatures$fid %in% PDnew$mmFeatures$fid))})
# 
# test_that("test pmSequence",
#           {expect_true(setequal(PDnew$pmFeatures$fid, PDnew$pmSequence$fid))})
# 
# test_that("test mmSequence",
#           {expect_true(setequal(PDnew$mmFeatures$fid, PDnew$mmSequence$fid))})
# 
# test_that("test level_dict: did not change",
#           {expect_equal(PDnew$level_dict, ParsedData$level_dict)})
# 
# test_that("test type_dict: did not change",
#           {expect_equal(PDnew$type_dict, ParsedData$type_dict)})
# 
# test_that("test core:",
#           {
#             expect_true(all(PDnew$featureSet$fsetid %in% PDnew$core$fsetid))
#             expect_true(all(PDnew$core$fsetid %in% PDnew$featureSet$fsetid))
#           })
# 
# test_that("test chrom_dict:",
#           {
#             expect_true(all(PDnew$chrom_dict$chrom_id %in% PSD$ProbeSetDataAnnotation$chrom))
#             expect_true(all(PSD$ProbeSetDataAnnotation$chrom %in% PDnew$chrom_dict$chrom_id ))
#           })
# 
# 
# ######## CDF ########
# 
# # dir2 <-  "/Users/milchevs/Downloads/CD_DrosGenome1\ 2/Full/DrosGenome1/LibFiles/"
# # load("/Users/milchevs/Downloads/NewPackages/DroGene1_1/dmel-all-r6.09.gtf.AnnotationDataFrame.RData") # Annotation 
# # alignmentDir = file.path("/Users/milchevs/ownCloud/ProbeRemapStep3/R/testdata/DroGene1_1/", "Alignments")
# # options(stringsAsFactors = FALSE)
# # 
# # outDir = "./tmpDRO/"
# # o2 = makeOriginalPackageObject(originalAnnotationDir = dir2, organism = "Dmel", species = "Drosophila melanogaster", outputDir = outDir )
# # 
# # Annotation = ANN_all
# # outputDir = outDir
# # level = "gene"
# # min_probe_number = 1
# # package_seed = o2
# # 
# # ParsedData <- package_seed$parsedDataOLD
# # Als <- processAlignments(alignment.dir = alignmentDir)
# # PS_dro = alignments2probesetsRaw(Alignments = Als, Annotation = Annotation, level = level)
# # PS_dro_cleaned <- cleanProbeSets(ProbeSets = PS_dro,
# #                                  level = level,
# #                                  min_probe_number = min_probe_number)
# # 
# # PS_controls <- addControlProbes(ProbeSetsClean = PS_dro_cleaned, outputDir = outputDir, alignmentType = package_seed[[3]])
# # shift = max(as.numeric(PS_controls$fsetid))
# # PSD <- make_objects_for_reannotation(Probesets = PS_dro_cleaned, level = level, Annotation = Annotation, shift = shift)
# 
# # save.image("./data/CDF_dro_allData_to_test_91.Rda")
# rm(list = setdiff(ls(), lsf.str()))
# load("./data/CDF_dro_allData_to_test_91.Rda")
# 
# PDnew <- build_new_ParsedDataOld_object_CDF(level = level, ParsedDataOld = ParsedData, PS_controls = PS_controls, Annotation = Annotation, ProbeSetData = PSD)
# 
# test_that("featureSet: all control and new target probe sets are present",
#           {
#             expect_true(all(PS_controls$fsetid %in% PDnew$featureSet$fsetid)) # all control probe sets are there
#             expect_true(all(PSD$ProbeSetDataAnnotation$fsetid %in% PDnew$featureSet$fsetid)) # all new main-type probe sets are there
#           })
# 
# test_that("pmFeatures: all control and main probe sets and probes are present",
#           {
#             mm_pm = merge(ParsedData$pmFeatures[ ,c("fid", "x", "y")], PDnew$pmFeatures[ ,c("fid","x", "y")])
#             expect_true( all(PDnew$pmFeatures$pmFeatures$fid %in% mm_pm$fid )) # no new probe fid appeared
#             expect_true( all(PS_controls$fid %in% PDnew$pmFeatures$fid)) # control probes are in the table
#             expect_true(  all(PSD$ProbeSetDataFrame$fid %in% PDnew$pmFeatures$fid)  ) # all probes from new probe sets
#             
#             expect_true(setequal(PDnew$pmFeatures$fsetid, union(PSD$ProbeSetDataFrame$fsetid, PS_controls$fsetid))) # all and only control + main new probe set ids are there
#           })
# 
# test_that("pmSequences: reflect pmFeatures",
#           {expect_true(setequal(PDnew$pmSequence$fid, PDnew$pmFeatures$fid))})
# 
# test_that("test mmSequence",
#           {
#             expect_true(!any( PDnew$mmSequence$fid %in% PDnew$pmFeatures$fid))
#             expect_true(all(PDnew$mmSequence$fid %in% PDnew$mmFeatures$fid)) # obratno ne trebuetsa!
#           })
# 
# test_that("test mmFeatures",
#           {expect_true(!any( PDnew$mmFeatures$fid %in% PDnew$pmFeatures$fid))})
# 
# test_that("test geometry",
#           {
#             expect_equal(ParsedData$geometry, PDnew$geometry)
#           })

