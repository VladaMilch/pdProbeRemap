library(testthat)
# # ####################### these tests can be run only on my local machine ##########################
# dir1 = "/Users/milchevs/ownCloud/ProbeRemapSteps14/inst/extdata/mogene2_affyLibFiles_pgfBased/"
# outDir = "."
# FastaFilePath1 <- annotation2fasta(originalAnnotationDir = dir1, organism = "Mouse", species = "Mus musculus", outputDir = outDir, reverse = FALSE)
# expect_true(file.exists(FastaFilePath))
# expect_true(file.exists( file.path(outDir, "tmpProbesTableDir/ProbesPGF.txt")))
# expect_true(file.exists( file.path(outDir, "tmpProbesTableDir/seed3.RData")))
# 
# system(command = paste0("rm -r ", file.path(outDir, "tmpProbesTableDir/")))
# system(command = paste0("rm -r ", file.path(outDir, "Probe_Sequences_Fasta/")))
# rm(FastaFilePath1)
# 
# # 
# dir2 <-  "/Users/milchevs/Downloads/CD_DrosGenome1\ 2/Full/DrosGenome1/LibFiles/"
# FastaFilePath2 <- annotation2fasta(originalAnnotationDir = dir2, organism = "Dmel", species = "Drosophila melanogaster", outputDir = outDir, reverse = FALSE )
# expect_true(file.exists(FastaFilePath2))
# expect_true(file.exists( file.path(outDir, "tmpProbesTableDir/ProbesPGF.txt")))
# expect_true(file.exists( file.path(outDir, "tmpProbesTableDir/seed3.RData")))

# system(command = paste0("rm -r ", file.path(outDir, "tmpProbesTableDir/")))
# system(command = paste0("rm -r ", file.path(outDir, "Probe_Sequences_Fasta/")))
# rm(FastaFilePath2)
# rm(dir1, dir2, outDir)

############### tests that can run everywhere ########
test_that("Input parameters are correct",
          {
            outDir = "."
            
            expect_error(suppressWarnings(
                          annotation2fasta(originalAnnotationDir = "./no/dir/", # wrong
                                          organism = "Mouse", 
                                          species = "Mus musculus", 
                                          outputDir = outDir, 
                                          author = "V M", 
                                          email = "v@m",
                                          reverse = FALSE)), 
                         regexp = "Annotation directory")
            
            expect_error(suppressWarnings(
                          annotation2fasta(originalAnnotationDir = "./no/dir/", 
                                          organism = "Mouse", 
                                          species = "Mmus musculus", 
                                          outputDir = outDir,
                                          author = "V M", 
                                          email = "v@m",
                                          reverse = FALSE)), 
                         regexp = "species")
            expect_error(suppressWarnings(
              annotation2fasta(originalAnnotationDir = "dir1", 
                                          organism = "Mmouse", 
                                          species = "Mus musculus", 
                                          outputDir = outDir, 
                                          author = "V M", 
                                          email = "v@m",
                                          reverse = FALSE)), 
                         regexp = "Mouse")
            expect_error(suppressWarnings(
              annotation2fasta(originalAnnotationDir = "dir1", 
                                          organism = "Mouse", 
                                          species = "Mus musculus", 
                                          outputDir = "outDir", 
                                          author = "V M", 
                                          email = "v@m",
                                          reverse = FALSE)), 
                         regexp = "Input existing directory path")
            expect_error(annotation2fasta(originalAnnotationDir = "dir1", 
                                          organism = "Mouse", 
                                          species = "Mus musculus", 
                                          outputDir = "outDir", 
                                          author = "V M", 
                                          email = "v@m",
                                          reverse = 0), 
                         regexp = "TRUE or FALSE")
          })



