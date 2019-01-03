OCCURRENCE_THRESHOLD = Inf
COVERAGE_THRESHOLD = 0.01
ENRICHMENT_THRESHOLD = 10000000000000000001.10

#getTfs <- function() {
#  tfs <- read.csv('../inventory/tf_use.csv', comment.char='#', stringsAsFactors=FALSE)
#  tfs <- tfs[tfs$run == 1,]
#  return(tfs)
#}
#
#getTfPairs <- function(tfs) {
#  tmp <- tfs[, c('family', 'tf', 'primer', 'motif', 'final')]
#  cross <- merge(tmp, tmp, by=NULL)
#  cross <- cross[cross$family.x != cross$family.y, ]
#  #print(nrow(cross))
#  #cross <- cross[!(((cross$family.x == 'homeodomain') & (cross$family.y == 'ETS')) | ((cross$family.y == 'homeodomain') & (cross$family.x == 'ETS'))), ]
#  #print(nrow(cross))
#  return(cross)
#}
#
#getConfigs <- function() {
#  configs <- read.csv(text="lf,rf
##2,0
#1,1
##0,2
#", comment.char='#')
#  return(configs)
#}
#
#getShapes <- function() {
#  #shapes <- c("MGW", "HelT", "ProT", "Roll")
#  shapes <- c("HelT", "ProT")
#  shapes <- c("MGW", "HelT", "ProT")
#  shapes <- c("MGW")
#  return(shapes)
#}
#
#
#getCycles <- function() {
#  cycles <- c(4)
#  return(cycles)
#}
#
#getTopResDir <- function(levels_type, en_th) {
#  path <- paste('/panfs/pan1/aptax/results', levels_type, en_th, sep='/')
#  return(path)
#}
#
#getEnrichmentThresholds <- function() {
#  en_ths <- c('1.01', '1.05')
#  en_ths <- c('1.20')
#  return(en_ths)
#}
#
#TOP_RES_DIR <- getTopResDir('publish', '1.10')
#
#getTopResPairDir <- function(family1, family2) {
#  pair <- paste(min(family1, family2), max(family1, family2), sep='_')
#  return(paste(TOP_RES_DIR, pair, sep='/'))
#}
#getTopResFamilyDir <- function(family) {
#  return(paste(TOP_RES_DIR, family, sep='/'))
#}

getSummaryFilePart <- function(tf1, primer1, motif1, cycle1, tf2, primer2, motif2, cycle2, shape) {
  if ((tf1 == tf2) && (primer1 == primer2) && (motif1 == motif2) && (cycle1 == cycle2)) {
    return(paste(tf1, primer1, shape, cycle1, motif1, sep='.'))
  }
  return(paste(tf1, primer1, motif1, cycle1, tf2, primer2, motif2, cycle2, shape, sep='.'))
}


assignStar <- function(x) {
#  if (x <= 0.0001) return('****')
  if (x <= 0.001) return('***')
  if (x <= 0.01) return('**')
  if (x <= 0.05) return('*')
  return('No')
}

assignSignificance <- function(x) {
  if (x <= 0.05) return('Yes')
  return('No')
}


getProbFile <- function(ctx, shape, motif, dist, lflank, rflank, shape_levels) {
  fname <- paste('prob', ctx, motif, shape, shape_levels, sep='.')
  fname <- paste(fname, paste0('d', dist), sep='.')
  if (ctx == 'fg') {
    fname <- paste(fname, paste0('l', lflank), paste0('r', rflank), sep='.')
  }
  fname <- paste(fname, 'txt', sep='.')
  return(fname)
}

getInfo <- function(datadir, cycle, ctx, shape, motif, dist, lflank, rflank, shape_levels) {
  #print(paste(datadir, cycle, ctx, shape, motif, lflank, rflank))
  prefix <- paste(cycle, 'non', motif, sep='.')
  prefix <- paste(prefix, paste0('d', dist), ctx, sep='.')
  if (ctx == 'fg') {
    prefix <- paste(prefix, paste0('l', lflank), paste0('r', rflank), sep='.')
  }
  prefix <- paste(prefix, shape, shape_levels, sep='.')
  fname <- paste(prefix, 'mer.cnt', sep='.')
  fpath <- paste(datadir, fname, sep='/')
  #print(paste("Reading", fpath))
  covname <- paste(prefix, 'mer.cov', sep='.')
  covpath <- paste(datadir, covname, sep='/')
  data_name <- paste0(ctx, cycle)
  count_col <- paste0(data_name, '.count')
  frac_col <- paste0(data_name, '.frac')

  #print(fpath)

  df <- read.table(fpath, header=FALSE, col.names = c(count_col, 'kmer'))
  df[, frac_col] <- df[, count_col] / sum(df[, count_col])

  cov_col <- paste0(data_name, '.cov')
  cov_count_col <- paste0(data_name, '.cov.count')
  cov_uniq_col <- paste0(data_name, '.cov.uniq')
  cov_uniq_count_col <- paste0(data_name, '.cov.uniq.count')
  cov_repeat_col <- paste0(data_name, '.repeat')

  #print(paste("Reading", covpath))
  col_names <- c('kmer', cov_uniq_count_col, cov_uniq_col, cov_count_col, cov_col, 'skip', cov_repeat_col, 'skip')
  col_classes <- c('character', rep('numeric', 4), 'NULL', 'numeric', 'NULL')
  df.cov <- read.csv(covpath, header=FALSE, col.names = col_names, colClasses = col_classes)

  df <- merge(df, df.cov, all.x=TRUE)
  return(df)
}

is.true <- function(x) { 
  !is.na(x) & x 
}


getEstimate <- function(ctx, shape, motif, dist, lflank, rflank, shape_levels, probability_dir) {
  #print(paste('probability_dir', probability_dir))
  file_estimate <- getProbFile(ctx, shape, motif, dist, lflank, rflank, shape_levels)
  path_estimate <- paste(probability_dir, file_estimate, sep='/')
  cn_estimated <- paste0(ctx, '0.estimated')
  df <- read.table(path_estimate, header=FALSE, col.names=c('kmer', 'skip', cn_estimated),
                   colClasses=c('character', 'NULL', 'numeric'))
  return(df)
}

getEnrichmentInfo <-function(topdir, probdir, tf, primer, family, cycle, ctx, shape, motif, dist, lflank, rflank, shape_levels) {
  #print(paste('probdir', probdir))
  datadir = paste(topdir, family, tf, primer, sep='/')
  fg <- getInfo(datadir, cycle, ctx, shape, motif, dist, lflank, rflank, shape_levels)

  datadir0 = paste(topdir, 'NoFamily', 'ZeroCycle', primer, sep='/')
  fg0 <- getInfo(datadir0, 0, ctx, shape, motif, dist, lflank, rflank, shape_levels)

  cn_kmer <- 'kmer'
  cn_fg0_estimated <- paste0(ctx, '0.estimated')

  fg0.est <- getEstimate(ctx, shape, motif, dist, lflank, rflank, shape_levels, probdir)

  cn_fg_frac <- paste0(ctx, cycle, '.frac')
  cn_fg_cov <- paste0(ctx, cycle, '.cov')
  cn_fg_en <- paste0(ctx, cycle, '.en')
  cn_fg_en_cov <- paste0(ctx, cycle, '.en.cov')
  cn_fg_en_act <- paste0(ctx, cycle, '.en.act')
  cn_fg_en_act_rank <- paste0(ctx, cycle, '.en.act.rank')

  #fg <- merge(fg, fg0, all.x=TRUE)
  fg <- merge(fg, fg0, all=TRUE)
  fg <- merge(fg, fg0.est, all.x=TRUE)
  fg[, cn_fg_en] <- fg[, cn_fg_frac]/fg[, cn_fg0_estimated]
  fg[, cn_fg_en_cov] <- fg[, cn_fg_cov]/fg$fg0.cov
  fg[, cn_fg_en_act] <- fg[, cn_fg_frac]/fg$fg0.count.frac
  fg[, cn_fg_en_act_rank] <- rank(-fg[, cn_fg_en_act], ties.method="first")

  return(fg)
}

#getPromiscuousShapes <- function(shape) {
#  promiscuous <- read.csv(text="shape,kmer
#MGW,HHHHHH
#MGW,HHHHHM
#MGW,HHHHMH
#MGW,MHHHHM
#MGW,HHMHHS
#", comment.char='#', stringsAsFactors=FALSE)
#  return(promiscuous[promiscuous$shape==shape, 'kmer'])
#}

getPromiscuousShapes <- function(shape) {
  promiscuous <- read.csv('exclude_shapemers.csv', comment.char='#', stringsAsFactors=FALSE)
  return(promiscuous[promiscuous$shape==shape, 'kmer'])
}


getCrossSummary <-function(top_data_dir, top_res_dir, probdir, tf1, primer1, family1, motif1, dist1, cycle1, tf2, primer2, family2, motif2, dist2, cycle2, shape, lflank, rflank, shape_levels) {

  #print(paste(c('top_data_dir', 'top_res_dir', 'tf1', 'primer1', 'family1', 'motif1', 'dist1', 'cycle1', 'tf2', 'primer2', 'family2', 'motif2', 'dist2', 'cycle2', 'shape', 'lflank', 'rflank', 'shape_levels'), c(top_data_dir, top_res_dir, tf1, primer1, family1, motif1, dist1, cycle1, tf2, primer2, family2, motif2, dist2, cycle2, shape, lflank, rflank, shape_levels)))

  dir.create(top_res_dir, showWarnings = FALSE, recursive = TRUE)

  file_part <- getSummaryFilePart(tf1, primer1, motif1, cycle1, tf2, primer2, motif2, cycle2, shape)

  file_shapes <- paste('shapes', file_part, lflank, rflank, 'csv', sep='.')
  file_detailed <- paste('detailed', file_part, lflank, rflank, 'csv', sep='.')
  file_enriched <- paste('enriched', file_part, lflank, rflank, 'csv', sep='.')
  file_fisher <- paste('fisher', file_part, lflank, rflank, 'csv', sep='.')

  path_shapes <- paste(top_res_dir, file_shapes, sep='/')
  path_detailed <- paste(top_res_dir, file_detailed, sep='/')
  path_enriched <- paste(top_res_dir, file_enriched, sep='/')
  path_fisher <- paste(top_res_dir, file_fisher, sep='/')

  fg <- getEnrichmentInfo(top_data_dir, probdir, tf1, primer1, family1, cycle1, 'fg', shape, motif1, dist1, lflank, rflank, shape_levels)
  bg <- getEnrichmentInfo(top_data_dir, probdir, tf2, primer2, family2, cycle2, 'bg', shape, motif2, dist2, lflank, rflank, shape_levels)

  #print(fg)
  #print(bg)

  dshapes <- ddply(rbind(data.frame(kmer=fg$kmer, context='fg'), data.frame(kmer=bg$kmer, context='bg')), .(kmer), summarise, context=makeBoth(context))

  #write.csv(dshapes, file=path_shapes, row.names=FALSE)

  total <- nrow(dshapes)   # required for Fisher's test

  cn_kmer <- 'kmer'
  cn_label <- paste0('label')

  cn_fg_frac <- paste0('fg', cycle1, '.frac')
  cn_fg_cov <- paste0('fg', cycle1, '.cov')
  cn_fg_en <- paste0('fg', cycle1, '.en')
  cn_fg_en_rank <- paste0('fg', cycle1, '.en.rank')

  cn_bg_frac <- paste0('bg', cycle2, '.frac')
  cn_bg_cov <- paste0('bg', cycle2, '.cov')
  cn_bg_en <- paste0('bg', cycle2, '.en')
  cn_bg_en_rank <- paste0('bg', cycle2, '.en.rank')


  fg <- fg[fg[, cn_fg_frac] >= OCCURRENCE_THRESHOLD | fg[, cn_fg_cov] >= COVERAGE_THRESHOLD, ]
  bg <- bg[bg[, cn_bg_frac] >= OCCURRENCE_THRESHOLD | bg[, cn_bg_cov] >= COVERAGE_THRESHOLD, ]

  df <- merge(fg, bg, all=TRUE)

  detailed_headers <- order(names(df))
  #write.csv(df[, detailed_headers], row.names=FALSE, file=path_detailed)

  df.enriched <- df[is.true(df[, cn_fg_en] >= ENRICHMENT_THRESHOLD | df[, cn_bg_en] >= ENRICHMENT_THRESHOLD), ]
  df.enriched[, cn_label] <- ifelse(is.true(df.enriched[, cn_fg_en] >= ENRICHMENT_THRESHOLD & df.enriched[, cn_bg_en] >= ENRICHMENT_THRESHOLD), 'both', 
                                 ifelse(is.true(df.enriched[, cn_fg_en] >= ENRICHMENT_THRESHOLD), 'fg', 'bg')) 

############## remove promiscuous ################
##print(df.enriched[, c('kmer', 'label4')])
#  df.enriched <- df.enriched[!(df.enriched$kmer %in% getPromiscuousShapes(shape)),]
##print(df.enriched[, c('kmer', 'label4')])

  df.enriched[, cn_fg_en_rank] <- rank(-df.enriched[, cn_fg_en], ties.method="first")
  df.enriched[, cn_bg_en_rank] <- rank(-df.enriched[, cn_bg_en], ties.method="first")
  df.enriched <- df.enriched[order(-df.enriched[, cn_fg_en], -df.enriched[, cn_bg_en]), ]

  important_headers <- c(cn_kmer,
                         cn_fg_en_rank,
                         cn_bg_en_rank,
                         cn_fg_en,
                         cn_bg_en,
                         cn_label,
                         cn_fg_frac,
                         'fg0.estimated',
                         cn_bg_frac,
                         'bg0.estimated',
                         cn_fg_cov,
                         'fg0.cov',
                         cn_bg_cov,
                         'bg0.cov')

  write.csv(df.enriched[, important_headers], row.names=FALSE, file=path_enriched)

  FSBN <- sum(df.enriched[, cn_label] == 'fg')
  FNBS <- sum(df.enriched[, cn_label] == 'bg')
  FSBS <- sum(df.enriched[, cn_label] == 'both')
  FNBN <- total - (FSBN + FSBS + FNBS)
  pvalue <- fisher.test(matrix(c(FSBS, FSBN, FNBS, FNBN), ncol=2), alternative='greater')$p.value
  dfisher <- data.frame(FSBN=FSBN, FSBS=FSBS, FNBN=FNBN, FNBS=FNBS, pvalue=pvalue)

  write.csv(dfisher, file=path_fisher, row.names=FALSE)

  #print(path_fisher)
  #return(dfisher)
}


#getSameFisherName <- function(tfs, config, shape) {
#  file <- paste('fisher', tfs$tf[1], tfs$primer[1], shape, tfs$final[1], tfs$motif[1], config$lf[1], config$rf[1], 'csv', sep='.')
#  path <- paste(getTopResFamilyDir(tfs$family[1]), tfs$tf[1], tfs$primer[1], file, sep='/')
#  return(path)
#}
#
#getCrossFisherName <- function(cross, config, shape) {
#  file <- paste('fisher', cross$tf.x[1], cross$primer.x[1], cross$motif.x[1], cross$final.x[1], cross$tf.y[1], cross$primer.y[1], cross$motif.y[1], cross$final.y[1],
#                    shape, config$lf[1], config$rf[1], 'csv', sep='.')
#  print(file)
#  path <- paste(getTopResPairDir(cross$family.x[1], cross$family.y[1]), file, sep='/')
#  return(path)
#}
#
#readSameFisher <- function(tfs, config, shape) {
#  file <- getSameFisherName(tfs, config, shape)
#  #print(file)
#  dfisher <- read.csv(file, stringsAsFactors=FALSE)
#  return(dfisher)
#}
#
#readCrossFisher <- function(tfs, config, shape) {
#  file <- getCrossFisherName(tfs, config, shape)
#  #print(file)
#  dfisher <- read.csv(file, stringsAsFactors=FALSE)
#  return(dfisher)
#}
#
#getCrossFisherTableName <- function(shape, config) {
#  file <- paste('dfisher_cross', shape, paste0('l', config$lf[1]), paste0('r', config$rf[1]), 'csv', sep='.')
#  path <- paste(TOP_RES_DIR, file, sep='/')
#  return(path)
#}
#
#getSameFisherTableName <- function(shape, config) {
#  file <- paste('dfisher_same', shape, paste0('l', config$lf[1]), paste0('r', config$rf[1]), 'csv', sep='.')
#  path <- paste(TOP_RES_DIR, file, sep='/')
#  return(path)
#}
#
#writeCrossFisherTable <- function(df, shape, config) {
#  file <- getCrossFisherTableName(shape, config)
#  write.csv(df, file=file, row.names=FALSE)
#}
#
#writeSameFisherTable <- function(df, shape, config) {
#  file <- getSameFisherTableName(shape, config)
#  write.csv(df, file=file, row.names=FALSE)
#}
#
#readCrossFisherTable <- function(shape, config) {
#  file <- getCrossFisherTableName(shape, config)
#  return(read.csv(file, stringsAsFactors=FALSE))
#}
#
#readSameFisherTable <- function(shape, config) {
#  file <- getSameFisherTableName(shape, config)
#  return(read.csv(file, stringsAsFactors=FALSE))
#}
#
#
#
#
#getSameEnrichedName <- function(tfs, config, shape) {
#  file <- paste('enriched', tfs$tf[1], tfs$primer[1], shape, tfs$final[1], tfs$motif[1], config$lf[1], config$rf[1], 'csv', sep='.')
#  path <- paste(getTopResFamilyDir(tfs$family[1]), tfs$tf[1], tfs$primer[1], file, sep='/')
#  return(path)
#}
#
#getCrossEnrichedName <- function(cross, config, shape) {
#  file <- paste('enriched', cross$tf.x[1], cross$primer.x[1], cross$motif.x[1], cross$final.x[1], cross$tf.y[1], cross$primer.y[1], cross$motif.y[1], cross$final.y[1],
#                    shape, config$lf[1], config$rf[1], 'csv', sep='.')
#  path <- paste(getTopResPairDir(cross$family.x[1], cross$family.y[1]), file, sep='/')
#  return(path)
#}
#
#readSameEnriched <- function(tfs, config, shape) {
#  file <- getSameEnrichedName(tfs, config, shape)
#  denriched <- read.csv(file, stringsAsFactors=FALSE)
#  names(denriched)[names(denriched) == paste0('label', tfs$fina[1])] <- 'label'
#  return(denriched)
#}
#
#readCrossEnriched <- function(tfs, config, shape) {
#  file <- getCrossEnrichedName(tfs, config, shape)
#  #print(file)
#  denriched <- read.csv(file, stringsAsFactors=FALSE)
#  return(denriched)
#}
#
#getCrossEnrichedTableName <- function(shape, config) {
#  file <- paste('denriched_cross', shape, paste0('l', config$lf[1]), paste0('r', config$rf[1]), 'csv', sep='.')
#  path <- paste(TOP_RES_DIR, file, sep='/')
#  return(path)
#}
#
#getSameEnrichedTableName <- function(shape, config) {
#  file <- paste('denriched_same', shape, paste0('l', config$lf[1]), paste0('r', config$rf[1]), 'csv', sep='.')
#  path <- paste(TOP_RES_DIR, file, sep='/')
#  return(path)
#}
#
#writeCrossEnrichedTable <- function(df, shape, config) {
#  file <- getCrossEnrichedTableName(shape, config)
#  write.csv(df, file=file, row.names=FALSE)
#}
#
#writeSameEnrichedTable <- function(df, shape, config) {
#  file <- getSameEnrichedTableName(shape, config)
#  print(file)
#  write.csv(df, file=file, row.names=FALSE)
#}
#
#readCrossEnrichedTable <- function(shape, config) {
#  file <- getCrossEnrichedTableName(shape, config)
#  return(read.csv(file, stringsAsFactors=FALSE))
#}
#
#readSameEnrichedTable <- function(shape, config) {
#  file <- getSameEnrichedTableName(shape, config)
#  return(read.csv(file, stringsAsFactors=FALSE))
#}
