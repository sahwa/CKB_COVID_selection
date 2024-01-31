library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(scales)
library(purrr)
library(tidyr)

setwd("/well/ckb/users/aey472/projects/CKB_COVID_VIP/LassI")

### function to find the unique peaks in the LASSI results which are ##
### above a particular quantile ##
### and a particular minimum distance ###

get_peaks = function(dat, threshold, min_snp_distance) {
  print(paste0("min snp distance is: ", min_snp_distance, " and threshold is: ", threshold))
  print(paste0("dat has ", nrow(dat), " rows"))

  dat[, index := rowid(chr)]
  numeric_cols = colnames(dat)
  numeric_cols = numeric_cols[-which(numeric_cols == "chr")]
  dat = dat[chr != "chr"]
  dat[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols, ]
  dat$chr = as.numeric(str_remove_all(dat$chr, "chr"))
  dat = dat[order(chr, pos)]
  dat[, c("nSNPs") := NULL]
  dat_list = split(dat, dat$chr)
  dat_list = map(dat_list, function(x) {
    upper = quantile(x$a_L, threshold, na.rm = T)
    x = x[a_L >= upper]
    x[, diff_index := abs(index - shift(index, n = 1, type = "lead")), by = "chr"]
    x[, group := dplyr::lag(cumsum(diff_index > min_snp_distance))]
    x = na.omit(x)
  })
  dat_groups = rbindlist(dat_list)
  dat_groups$group[dat_groups$group %in% NA] = 0
  dat_groups[, group := paste0(chr, "_", group)]
  list(dat_groups, dat_groups[, list(start = min(start), end = max(end)), group])
}

get_overlap_50_pc = function(selected_regions, vip_regions, wiggle_room) {
  vip_regions$start = vip_regions$start - wiggle_room
  vip_regions$end = vip_regions$end + wiggle_room

  selected_regions[, chr := as.character(chr)]
  vip_regions[, chr := as.character(chr)]

  setkey(selected_regions, chr, start, end)
  setkey(vip_regions, chr, start, end)

  overlap_res = na.omit(foverlaps(vip_regions, selected_regions))
  overlap_res[, index := .I]
  overlap_res_list = split(overlap_res, overlap_res$index)

  overlap_res_list_50_perc = map(overlap_res_list, function(x) {
    sum_overlap = sum(with(x, seq(start, end) %in% seq(i.start, i.end)))

    length_SR = with(x, length(seq(start, end)))
    perc_sr_overlap = sum_overlap / length_SR

    length_vip = with(x, length(seq(i.start, i.end)))
    perc_vip_overlap = sum_overlap / length_vip

    perc_sr_overlap > 0.5 | perc_vip_overlap > 0.5
  })

  overlap_res_list_50_perc = unlist(overlap_res_list_50_perc)

  rbindlist(overlap_res_list[overlap_res_list_50_perc])
}

get_overlap = function(selected_regions, vip_regions, wiggle_room) {
  vip_regions$start = vip_regions$start - wiggle_room
  vip_regions$end = vip_regions$end + wiggle_room

  selected_regions[, chr := as.character(chr)]
  vip_regions[, chr := as.character(chr)]

  setkey(selected_regions, chr, start, end)
  setkey(vip_regions, chr, start, end)

  foverlaps(vip_regions, selected_regions)
}
