### code to perform bootstrapping analysis to infer enrichment of selection regions at GWAS hits for COVID-19

libs = c("AnnotationHub", "ensembldb", "EnsDb.Hsapiens.v86", "nullranges", "GenomeInfoDb", "DNAcopy", "nullrangesData", "plyranges", "tidyr", "ggridges", "purrr", "ggplot2", "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "gridExtra")
suppressPackageStartupMessages(invisible(sapply(libs, library, character.only = TRUE)))

## read in different vip_sets ##
length_keep = 5e8
pop = "CKB"
gwas_significant = 5e-8

## filter GWAS hits to keep only unique, GWAS significant, autosomes and ##
gwas_hits = fread("/home/sam/work/COVID_selection/data/GWAS_stuff/hits_extended.txt") %>% 
  janitor::clean_names() %>%
  filter(chromosome_b38 %in% 1:22) %>%
  filter(pval <= gwas_significant) %>%
  filter(!duplicated(rsid))


## convert the data to the correct genomic build. 

lifted_coords = fread("/home/sam/work/COVID_selection/data/GWAS_stuff/hits_b38_chr_pos_liftedb37.txt") %>% 
  mutate(b38_pos = str_split(V4, "-", simplify=T)[,2])

gwas_hits = gwas_hits %>% 
  mutate(position_b37 = lifted_coords$V2[match(position_b38, lifted_coords$b38_pos)]) %>%
  select(chromosome_b38, position_b37, rsid) %>%
  mutate(chr = chromosome_b38, start = position_b37, end = position_b37) %>%
  select(chr, start, end, rsid) %>%
  filter(!(chr == 6 & data.table::between(start, 20000000, 40000000))) %>%
  filter(!(chr == 6 & data.table::between(end, 20000000, 40000000))) %>%
  rbind(data.table(chr = c(15, 18, 20, 22), start=1, end=1, rsid = "pascal")) %>%
  mutate(chr = paste0("chr", chr))

gwas_hits$start[gwas_hits$rsid == "rs71175234"] = 46373058 
gwas_hits$end[gwas_hits$rsid == "rs71175234"] = 46373058 

gwas_hits = makeGRangesFromDataFrame(gwas_hits)
gwas_hits@seqinfo@seqlengths = seqlengths(Hsapiens)[1:22]
gwas_hits = gwas_hits %>%
    mutate(id = seq_along(.)) 

### read in LASSI etc ###
lassi_file = paste0("/home/sam/work/COVID_selection/data/LASSI/",pop,"_Lassi_peaks.txt") 

lassi = fread(lassi_file) %>%
  dplyr::rename(chr = 1, start = 2, end = 3) %>%
  mutate(chr = paste0("chr", chr))

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
# HLA
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

# g = g %>% filter(seqnames %in% paste0("chr", 1:22))
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
R = 10
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

gwas_hits = gwas_hits %>% mutate(n_overlaps = count_overlaps(., lassi))
  
boot_stats = gwas_hits %>% 
  join_overlap_inner(boots, maxgap=10*1000) %>%
  group_by(id.x, iter) %>%
  summarize(n_overlaps = n()) %>%
  as_tibble() %>%
  complete(id.x, iter, fill=list(n_overlaps = 0)) %>%
  group_by(iter) %>%
  summarize(sumOverlaps = sum(n_overlaps))

label = paste0("p=", sum(boot_stats$sumOverlaps > sum(gwas_hits$n_overlaps)) / R)

ggplot(boot_stats, aes(sumOverlaps)) +
    geom_histogram(binwidth=2)+
    geom_vline(xintercept = sum(gwas_hits$n_overlaps), linetype = "dashed") +
    ggpp::geom_text_npc(aes(npcx = x, npcy = y, label=label), size=12, colour='red',
                data = data.frame(x = 0.75, y = 0.75, label=label)) +
    theme(plot.title = element_text(hjust = 0.5, size=24))

ggsave(paste0("/home/sam/work/COVID_selection/bootstrapping/",pop,"_bootstrap_res.GWAS_hits.pdf"), width=18.2, height=7.2)

boot_stats = gwas_hits %>% 
  join_overlap_inner(boots, maxgap=10*1000) %>%
  group_by(id.x, iter) %>%
  summarize(n_overlaps = n()) %>%
  as_tibble() %>%
  complete(id.x, iter, fill=list(n_overlaps = 0)) %>%
  group_by(iter) %>%
  summarize(sumOverlaps = sum(n_overlaps))

n_real_overlap = sum(gwas_hits$n_overlaps)
enrichment_intervals = n_real_overlap / sapply(c(0.975, 1-0.975), function(x) quantile(boot_stats$sumOverlaps, prob=x))
median_enrichment = n_real_overlap / median(boot_stats$sumOverlaps)
mean_enrichment = n_real_overlap / mean(boot_stats$sumOverlaps)
P = sum(boot_stats$sumOverlaps > sum(gwas_hits$n_overlaps)) / R
