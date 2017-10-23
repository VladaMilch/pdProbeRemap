require(testthat)

# ####################### these tests can be run only on my local machine ##########################
# 
# dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
# dir2 <-  "/Users/milchevs/Downloads/CD_DrosGenome1\ 2/Full/DrosGenome1/LibFiles/"
# outDir = "."
# 
# o1 = makeOriginalPackageObject(originalAnnotationDir = dir1, organism = "Mouse", species = "Mus musculus", outputDir = outDir)
# o2 = makeOriginalPackageObject(originalAnnotationDir = dir2, organism = "Dmel", species = "Drosophila melanogaster", outputDir = outDir )
# 
# expect_true( setequal(names(o1), c("object", "parsedDataOLD", "AnnotationType")))
# expect_is(  o1[[1]],  "AffyGenePDInfoPkgSeed" )
# expect_is(  o1[[2]],  "list" )
# expect_is(  o1[[3]],  "character" )
# 
# expect_true( setequal(names(o2), c("object", "parsedDataOLD", "AnnotationType")))
# expect_is(  o2[[1]],  "AffyExpressionPDInfoPkgSeed" )
# expect_is(  o2[[2]],  "list" )
# expect_is(  o2[[3]],  "character" )
# 
# expect_true(dir.exists(file.path(outDir, "tmpProbesTableDir")))
# expect_true(file.exists(file.path(outDir, "tmpProbesTableDir", "ProbesCDF.txt")))
# expect_true(file.exists(file.path(outDir, "tmpProbesTableDir", "ProbesPGF.csv")))
# 
# expect_error(makeOriginalPackageObject(originalAnnotationDir = dir1, organism = "Muse", species = "Mus musculus", outputDir = outDir), regexp = "Organism.*Human.*Mouse.*Dmel")
# expect_error(makeOriginalPackageObject(originalAnnotationDir = dir1, organism = "Mouse", species = "Mus musculus", outputDir = "slkdjfh"), regex = "Output")
# expect_message(makeOriginalPackageObject(originalAnnotationDir = ".", organism = "Mouse", species = "Mus musculus", outputDir = outDir), regex = "Error: Original Annotation Input Files are missing.")