## code to perform the segmented block bootstrap analysis for the LASSi regions of selection ##
## we mostly followed the analyses described in the official vignette here https://bioconductor.org/packages/release/bioc/vignettes/nullranges/inst/doc/bootRanges.html#Segmented_block_bootstrap
## the only modification we made was to remove very large selection regions to preserve the isochore structure ##



libs = c("AnnotationHub", "ensembldb", "EnsDb.Hsapiens.v86", "nullranges", "GenomeInfoDb", "DNAcopy", "nullrangesData", "plyranges", "tidyr", "ggridges", "purrr", "ggplot2", "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "gridExtra")
suppressPackageStartupMessages(invisible(sapply(libs, library, character.only = TRUE)))

## read in different vip_sets ##

## maximum length of LASSi selection region to keep 
length_keep = 5e5

## population we are analysing
pop = "CKB"

regions = list.files(pattern=".regions", path="/home/sam/work/COVID_selection/data/VIPs", full.names=T)

## read in VIP list
vip_list = lapply(regions, fread)
names(vip_list) = str_remove_all(basename(regions), "\\.regions")

## turn vips into GRanges object
vip_list = lapply(vip_list, function(x) {
  colnames(x) = c("chr", "start", "end", "gene")
  x[, chr := paste0("chr", chr)]
  x = makeGRangesFromDataFrame(x, seqnames.field = "chr")

  x@seqinfo@seqlengths = seqlengths(Hsapiens)[as.numeric(str_remove_all(x@seqinfo@seqnames, "chr"))]

  x = x %>%
      mutate(id = seq_along(.)) 
  })

### read in LASSI etc ###
lassi_file = paste0("/home/sam/work/COVID_selection/data/LASSI/",pop,"_Lassi_peaks.txt") 
lassi = fread(lassi_file) %>%
  dplyr::rename(chr = 1, start = 2, end = 3) %>%
  mutate(chr = paste0("chr", chr))

## turn lassi into GRanges and then remove any regions which are too long
colnames(lassi) = c("chr", "start", "end")
lassi = makeGRangesFromDataFrame(lassi, seqnames.field = "chr")
lassi_remove = lassi %>% plyranges::filter(width(lassi) >= length_keep)

## remove any segments which are too long ##
lassi = lassi %>% plyranges::filter(width(lassi) <= length_keep)
lassi@seqinfo@seqnames = sprintf("chr%d", 1:22)
lassi@seqinfo@seqlengths = seqlengths(Hsapiens)[1:22]

lassi = lassi %>%
  mutate(id = seq_along(.)) 

## get exclusions

ah = AnnotationHub()
# hg38.Kundaje.GRCh38_unified_Excludable
exclude_1 = ah[["AH107305"]]
# hg38.UCSC.centromere
exclude_2 = ah[["AH107354"]]
# hg38.UCSC.telomere
exclude_3 = ah[["AH107355"]]
# hg38.UCSC.short_arm
exclude_4 = ah[["AH107356"]]
# HLA removal 
HLA = GRanges("chr6", IRanges(20 * 1e6, width=20e6), x_id=1)
# combine them
suppressWarnings({
  exclude = trim(c(exclude_1, exclude_2, exclude_3, exclude_4, lassi_remove, HLA))
})

exclude = sort(GenomicRanges::reduce(exclude))
## remove MT chrX and chrY from the 
seqlevels(exclude, pruning.mode="coarse") = setdiff(seqlevels(exclude), c("MT", "chrX", "chrY"))


## make segmentations according to gene density
edb = EnsDb.Hsapiens.v86
filt = AnnotationFilterList(GeneIdFilter("ENSG", "startsWith"))
g = genes(edb, filter = filt)

g = keepStandardChromosomes(g, pruning.mode = "coarse")
# MT is too small for bootstrapping, so must be removed
seqlevels(g, pruning.mode="coarse") = setdiff(seqlevels(g), c("MT", "X", "Y"))

# normally we would assign a new style, but for recent host issues
# that produced vignette build problems, we use `paste0`
## seqlevelsStyle(g) = "UCSC" 
seqlevels(g) = paste0("chr", seqlevels(g))

genome(g) = "hg38"
g = sortSeqlevels(g)
g = sort(g)
table(seqnames(g))


set.seed(5)
exclude = exclude %>%
  plyranges::filter(width(exclude) >= 500)
L_s = 1e6
seg_cbs = segmentDensity(g, n = 3, L_s = L_s, exclude = exclude, type = "cbs")


set.seed(5) # for reproducibility
R = 1000
blockLength = 1e6 
boots = bootRanges(lassi, blockLength, R = R, seg = seg_cbs, exclude=exclude)


combined = lassi %>% 
  mutate(iter=0) %>%
  bind_ranges(boots) %>% 
  plyranges::select(iter)

stats = combined %>% 
  group_by(iter) %>%
  summarize(n = n()) %>%
  as_tibble()
head(stats)


interdist = function(dat) {
    x = dat[-1,]
    y = dat[-nrow(dat),]
    ifelse(x$seqnames == y$seqnames,
           x$start + floor((x$width - 1)/2) -
           y$start - floor((y$width - 1)/2), NA)
}


# just looking at first 3 iterations...
combined %>% 
  plyranges::filter(iter %in% 0:3) %>%
  mutate(iter = droplevels(iter)) %>%
  plyranges::select(iter) %>%
  as_tibble() %>% 
  nest(data = !iter) %>%
  mutate(interdist = map(data, interdist)) %>% 
  dplyr::select(iter, interdist) %>% 
  unnest(interdist) %>% 
  mutate(type = ifelse(iter == 0, "original", "boot"),
         interdist = pmax(interdist, 0)) %>%
  filter(!is.na(interdist)) %>%
  ggplot(aes(log10(interdist + 1), iter, fill=type)) +
  geom_density_ridges(alpha = 0.75) +
  geom_text(data = head(stats, 4),
            aes(x=1.5, y=iter, label=paste0("n=",n), fill=NULL),
            vjust=1.5)

vip_list = lapply(vip_list, function(x) {
  x %>% mutate(n_overlaps = count_overlaps(., lassi))
  })


## plot the bootstrap results 
boot_stats_plot_list = lapply(names(vip_list), function(x) {

  boot_stats = vip_list[[x]] %>% 
    join_overlap_inner(boots, maxgap=0) %>%
    group_by(id.x, iter) %>%
    summarize(n_overlaps = n()) %>%
    as_tibble() %>%
    complete(id.x, iter, fill=list(n_overlaps = 0)) %>%
    group_by(iter) %>%
    summarize(sumOverlaps = sum(n_overlaps))

  label = paste0("p=", sum(boot_stats$sumOverlaps > sum(vip_list[[x]]$n_overlaps)) / R)

   ggplot(boot_stats, aes(sumOverlaps)) +
      geom_histogram(binwidth=2)+
      geom_vline(xintercept = sum(vip_list[[x]]$n_overlaps), linetype = "dashed") +
      ggpp::geom_text_npc(aes(npcx = x, npcy = y, label=label), size=12, colour='red',
                  data = data.frame(x = 0.75, y = 0.75, label=label)) +
      ggtitle(x) +
      theme(plot.title = element_text(hjust = 0.5, size=24))
  }
)

n = length(boot_stats_plot_list)
nCol <- floor(sqrt(n))

pdf(paste0("/home/sam/work/COVID_selection/bootstrapping/",pop,"_bootstrap_res.pdf"), width=18.2, height=7.2)
do.call("grid.arrange", c(boot_stats_plot_list, ncol=nCol))
dev.off()


##  use bootstraps to calculate empirical p-vals and confidence invervals
dt_res = data.table(do.call(rbind.data.frame, lapply(names(vip_list), function(x) {

  boot_stats = vip_list[[x]] %>% 
    join_overlap_inner(boots, maxgap=10*1000) %>%
    group_by(id.x, iter) %>%
    summarize(n_overlaps = n()) %>%
    as_tibble() %>%
    complete(id.x, iter, fill=list(n_overlaps = 0)) %>%
    group_by(iter) %>%
    summarize(sumOverlaps = sum(n_overlaps))


  n_real_overlap = sum(vip_list[[x]]$n_overlaps)
  enrichment_intervals = n_real_overlap / sapply(c(0.975, 1-0.975), function(x) quantile(boot_stats$sumOverlaps, prob=x))
  median_enrichment = n_real_overlap / median(boot_stats$sumOverlaps)
  mean_enrichment = n_real_overlap / mean(boot_stats$sumOverlaps)

  list(
    class = x,
    n_real_overlap = n_real_overlap, 
    upper_c = enrichment_intervals[1], 
    lower_c = enrichment_intervals[2], 
    p = sum(boot_stats$sumOverlaps > sum(vip_list[[x]]$n_overlaps)) / R,
    median_enrichment = median_enrichment, 
    mean_enrichment = mean_enrichment, 
    n_vips = length(width(vip_list[[x]])))
})))

fwrite(dt_res, paste0("/home/sam/work/COVID_selection/bootstrapping/",pop,"_bootstrap_res.txt"))
