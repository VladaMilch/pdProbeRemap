test_make_objects_for_reannotation <- function()
{
  ProbeSets <- data.frame(
    probe_id = c(101152, 101153, 101154, 101155, 101156, 101157), 
    transcript_id = rep(NA, 6),
    gene_id = c("gene2", "gene3", "gene1", "gene4", "gene4", "gene4")
  )
  Annotation1 <- data.frame(
    id = c(1,2,3,4,5,6,7,8,9,10),
    seqid = c("2L", "2L", "3L", "3L", "3L", "3L", "3L", "3R", "3R", "3R"),
    source = rep("FlyBase", 10),
    feature = rep("mRNA", 10),
    start = rep(1, 10), 
    end = rep(100, 10), 
    score = rep (15, 10), 
    strand = c("+", "+", "+", rep("-", 7)),
    transcript_symbol = rep("TS", 10), 
    transcript_id = c(1:10), 
    gene_symbol = rep("GS", 10), 
    gene_id = c("gene1", "gene1", 
                 "gene2", 
                 "gene4", "gene4", "gene4", "gene4",
                 "gene3", "gene3", "gene3")
  )
  
  gene_level_obj = make_objects_for_reannotation(Probesets = ProbeSets, Annotation = Annotation1, level = "gene", shift = 0)
  
  expect_equal(length(names(gene_level_obj)),2)
  expect_equal(names(gene_level_obj)[1], c("ProbeSetDataAnnotation"))
  expect_equal(names(gene_level_obj)[2], c("ProbeSetDataFrame"))
  expect_is(gene_level_obj[[1]], "data.frame")
  expect_is(gene_level_obj[[2]], "data.frame")
  expect_equal(colnames(gene_level_obj[[1]]),c("fsetname", "fsetid", "strand", "chrom"))
  expect_equal(colnames(gene_level_obj[[2]]),c("fid", "fsetid"))
  expect_true(  setequal(unique(gene_level_obj[[1]]$strand), c("+", "-"))  )
  
  #expect_true( subset(gene_level_obj[[2]], fid == "101153")$fsetid == 8 )
  
  expect_error(make_objects_for_reannotation(Probesets = ProbeSets, Annotation = Annotation1, level = "ssdf"))
  expect_error(make_objects_for_reannotation(Probesets = ProbeSets, Annotation = Annotation1[ ,1:3], level = "gene"))
  expect_error(make_objects_for_reannotation(Probesets = ProbeSets[, 1:2], Annotation = Annotation1, level = "transcript"))

  message("All tests passed for make_objects_for_reannotation. \n")
}

test_make_objects_for_reannotation()


