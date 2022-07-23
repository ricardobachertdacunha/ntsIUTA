
### Files ---------------------------------------------------------------------------------------------------

spt <- samplesTable(dt)
spt <- spt[c(3:5), ]

### Load target dim --------------------------------------------------------------------------------------------

dim <- fread("C:\\Users\\Ricardo\\Documents\\CodeProjects\\ntsIUTA\\inst\\extdata\\mix1_RTs.csv")

rtRange <- c(min(dim[, .(orb_rt, tof_rt)]), max(dim[, .(orb_rt, tof_rt)]))
rtRange <- c(300, 700)
dim <- dim[rt >= rtRange[1] & rt <= rtRange[2], ]

mzRange <- c(min(dim[, mz]), max(dim[, mz]))
mzRange <- c(100, 400)


### Loop for trim -------------------------------------------------------------------------------------------

for (f in spt$file) {
  msf <- mzR::openMSfile(f)
  hd <- mzR::header(msf)
  hd2 <- hd[hd$retentionTime >= rtRange[1] & hd$retentionTime <= rtRange[2], ]
  
  spec <- mzR::peaks(msf, scans = hd2$seqNum)
  for (s in seq_len(length(spec))) {
    if (hd2$msLevel[s] == 1) {
      spec[[s]] <- spec[[s]][spec[[s]][, 1] >= mzRange[1] & spec[[s]][, 1] <= mzRange[2], ]
    }
  }
  
  hd2$seqNum <- seq_len(nrow(hd2))
  mzR::close(msf)
  
  
  format_tag <- "mzml"
  if (grepl("XML", f)) {
    format_tag <- "mzxml"
  }
  
  mzR::writeMSData(object = spec, file = paste0(dirname(f), "/trim_", basename(f)), header = hd2, outformat = format_tag)
  
  rm(msf, hd, hd2, spec, s, format_tag)
}

rm(f)
rm(rtRange, mzRange, dim, spt)