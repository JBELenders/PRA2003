# PRA2003 - Week 5: asymmetry between wild type (WT) and mutant, overall and as a
#                   function of momentum, for all 6 strain pairs and all 10 data sets.
#
# WHAT IT DOES
#   1. For every bacterium of the 6 pairs it computes the momentum |p| and puts it
#      in one of 5 momentum bins, then counts WT and mutant per file and per bin.
#   2. Overall result per pair (Section 6.2 of the README):
#        difference per event  D = (N_WT - N_mutant) / number of events
#        asymmetry             A = (N_WT - N_mutant) / (N_WT + N_mutant)
#      calculated in EACH file, then mean +- standard deviation over the 10 files.
#   3. The same asymmetry A in each momentum bin, again mean +- SD over the 10 files.
#   4. Two simple checks of "does A depend on momentum?":
#        - constant fit: is the same A in all 5 bins a good description? (chi-squared)
#        - slope: a straight line through the 5 bins is fitted in each file, then
#          mean slope +- SD over the 10 files, and slope / SD = significance.
#   5. One plot with a panel per pair, and three csv files with all numbers.
#
# METHOD (same as Weeks 4-5 / the lecture: sub-sampling)
#   value = mean over the 10 files, uncertainty = standard deviation over the 10 files,
#   significance = value / uncertainty, "significant" if it is 3 or more.
#
# HOW TO RUN
#   Put output-Set1.txt ... output-Set10.txt in the folder "DataFiles" next to this
#   script, and sub_sample_results.csv (from Week 4) next to it as well.
#   Then run:  Rscript Week5_Momentum_Analysis.R
#   The first run reads the raw data (about 20 minutes) and saves the counts in
#   momentum_perfile_results.csv and momentum_bin_edges.csv. Later runs use those two
#   files and take a few seconds. Delete them if you want to recount from scratch.
#
# FILES WRITTEN
#   momentum_perfile_results.csv  counts per file, pair and bin (the saved counts)
#   momentum_bin_edges.csv        the 6 edges of the 5 momentum bins
#   overall_asymmetry.csv         D and A per pair, mean +- SD over the files
#   asymmetry_per_bin.csv         A per pair and momentum bin, mean +- SD over the files
#   momentum_fits.csv             constant fit and slope per pair
#   asymmetry_vs_momentum.png     the plot
#
# A NOTE ON THE COMMENTS
#   I (Julian) am still learning R, so the comments below explain in my own words what
#   I think each line does. Where I am not 100% sure about a technical detail I say so.


# ---------------------------------------------------------------- Settings ----
# These are the numbers that I might want to change later. They are all put at the
# top so that I do not have to search through the whole script to find them.
# (In R, "<-" means "put the value on the right into the name on the left".)
nFiles <- 10                          # how many data files there are (Set 1 to Set 10)
nBins <- 5                            # how many momentum bins we use
chunkSize <- 1000000                  # number of lines read from a file at a time
dataFolder <- "DataFiles"             # the folder where the 10 big txt files are
# paste0() glues text together without spaces. Because 1:nFiles is 1,2,...,10 this
# gives 10 file names in one go: "DataFiles/output-Set1.txt", ..., "DataFiles/output-Set10.txt"
filenames <- paste0(dataFolder, "/output-Set", 1:nFiles, ".txt")


# The 6 pairs. The two assigned to this analysis (2212 and 3122) come first.
# c(...) makes a list of values (a "vector") and the three vectors below belong
# together: the first name goes with the first ID, the second name with the second ID, etc.
pairName <- c("Pseudomonas aeruginosa", "Streptococcus pneumoniae", "E. coli",
              "Bacillus subtilis", "Mycobacterium tuberculosis", "Salmonella")
pairWT <- c(2212, 3122, 211, 321, 3312, 3334)   # ID of the wild type
pairMutant <- -pairWT                           # the mutant has the same ID with a minus
nPairs <- length(pairName)                      # length() counts how many items are in the list (= 6)
targetIDs <- c(pairWT, pairMutant)              # all 12 IDs we care about, in one list


# Names of the csv files where the counted numbers are saved, so that we do not
# have to read the huge raw files again every time (that takes about 20 minutes).
countsFile <- "momentum_perfile_results.csv"
edgesFile <- "momentum_bin_edges.csv"


# Week 4 output, needed only for the number of (non-empty) events in each file.
# The script looks in two places for the file: next to the script, or in a folder
# called PRA2003. file.exists() gives TRUE/FALSE for each, we keep only the ones that
# exist and then take the first one ([1]). If none exists we get NA (= "nothing").
eventsFile <- c("sub_sample_results.csv", "PRA2003/sub_sample_results.csv")
eventsFile <- eventsFile[file.exists(eventsFile)][1]
if (is.na(eventsFile)) {
  # stop() ends the script with an error message, so we do not continue without this file
  stop("Cannot find sub_sample_results.csv. Run the Week 4 sub-sampling script first.")
}




# ------------------------------------------------------- Reading one file ----
# This is a function: a block of code with a name that we can use again for each file.
# Returns the momentum |p| and the ID of every bacterium that belongs to one of
# the 6 pairs. The file is read in chunks so that it never has to fit in memory.
# A bacterium line has 4 numbers (px py pz ID) = 3 spaces. Header lines
# (event number and number of bacteria) have only 2 numbers, so they are skipped.
# (So the trick is: count the spaces in a line. 3 spaces = bacterium, otherwise it is a header.)
readMomentumAndID <- function(filepath) {
  con <- file(filepath, "r")                # open the file for reading ("r")
  on.exit(close(con))                       # always close the file again, even if something goes wrong


  # Empty lists where we collect the results of each chunk. At the end we glue them together.
  momentumList <- list()
  idList <- list()
  chunkNumber <- 0                          # counts which chunk we are at


  # repeat{} runs the code inside again and again until we say break
  repeat {
    # read the next 1,000,000 lines. Because the file stays open, the next time we
    # call this it continues where it stopped (I think it works like a bookmark).
    chunk <- readLines(con, n = chunkSize)
    if (length(chunk) == 0) {
      break                                 # end of the file: nothing left to read, so stop the loop
    }
    chunkNumber <- chunkNumber + 1


    # gregexpr finds every place where a space occurs in each line and lengths() counts them.
    # So nSpaces has, for every line, the number of spaces in that line.
    nSpaces <- lengths(gregexpr(" ", chunk, fixed = TRUE))
    # keep only the lines with exactly 3 spaces (the bacteria lines, not the headers)
    bacteriumLines <- chunk[nSpaces == 3]


    if (length(bacteriumLines) > 0) {
      # sub() with "^.* " removes everything up to and including the last space, so only the
      # ID is left. as.integer() turns that text into a real number.
      ids <- as.integer(sub("^.* ", "", bacteriumLines))   # everything after the last space
      # %in% checks for each ID whether it is in our list of 12. Result is TRUE/FALSE per line.
      keep <- ids %in% targetIDs                           # only the 12 IDs we need
      if (any(keep)) {                                     # any() = is at least one of them TRUE?
        # read.table turns the text lines into a table of numbers (4 columns: px, py, pz, ID)
        numbers <- read.table(text = bacteriumLines[keep], colClasses = "numeric")
        # momentum = length of the vector (px, py, pz) = sqrt(px^2 + py^2 + pz^2).
        # numbers[[1]] means "first column", [[2]] second column, [[3]] third column.
        momentumList[[chunkNumber]] <- sqrt(numbers[[1]]^2 + numbers[[2]]^2 + numbers[[3]]^2)
        idList[[chunkNumber]] <- ids[keep]
      }
    }
  }
  # unlist() glues the chunk results into one long list. The function gives back
  # both the momenta and the IDs, in the same order, so they still belong together.
  list(momentum = unlist(momentumList), id = unlist(idList))
}


# Bin edges: the 5 bins each hold about 20% of the bacteria of the given momenta.
# The first and last bin are open-ended so that every bacterium is counted.
makeBinEdges <- function(momentum) {
  # quantile() finds the momentum values that cut the data in equal-sized parts.
  # seq(0, 1, length.out = 6) gives 0, 0.2, 0.4, 0.6, 0.8, 1, so we get 6 edges = 5 bins.
  edges <- quantile(momentum, probs = seq(0, 1, length.out = nBins + 1), names = FALSE)
  edges[1] <- -Inf                          # first edge becomes "minus infinity" ...
  edges[nBins + 1] <- Inf                   # ... and last edge "plus infinity", so nothing falls outside
  edges
}




# ------------------------------------- Counts per file, pair and bin (cached) ----
# "Cached" means saved for later. Here the script checks if the saved csv files already
# exist and have the right size, so it can skip the slow reading of the raw data.
# nFiles * nPairs * nBins = 10 * 6 * 5 = 300 rows are expected in the counts file.
countsOK <- file.exists(countsFile) && nrow(read.csv(countsFile)) == nFiles * nPairs * nBins
edgesOK <- file.exists(edgesFile) && nrow(read.csv(edgesFile)) == nBins + 1


# There are three situations, handled by if / else if / else:
if (countsOK && edgesOK) {
  # Situation 1: everything is saved already, so just load it (fast)
  cat("Using the saved counts in", countsFile, "\n")
  binEdges <- read.csv(edgesFile)$Edge      # the $ picks the column called "Edge"


} else if (countsOK) {
  # Situation 2: counts are there, only the bin edges (used for the axis labels) are missing
  # so we only read Set 1 again to calculate the edges
  cat("Recomputing the bin edges from", filenames[1], "\n")
  binEdges <- makeBinEdges(readMomentumAndID(filenames[1])$momentum)
  write.csv(data.frame(Edge = binEdges), edgesFile, row.names = FALSE)


} else {
  # Situation 3: nothing is saved, so we have to read all the raw data (slow, ~20 minutes)
  # First check that all 10 files really exist, otherwise stop with a clear message.
  for (f in filenames) {
    if (!file.exists(f)) stop("Cannot find '", f, "' (working folder: ", getwd(), ")")
  }


  countRows <- data.frame()                 # an empty table; we add rows to it below
  for (k in 1:nFiles) {                     # k = the file number, goes 1, 2, ..., 10
    cat("Reading", filenames[k], "...\n")   # cat() prints a message so we can see progress
    fileData <- readMomentumAndID(filenames[k])   # uses the function from above


    # The bin edges are determined once, from Set 1, and then used for all files
    # (so that the bins are the same everywhere and can be compared)
    if (k == 1) {
      binEdges <- makeBinEdges(fileData$momentum)
    }
    # cut() sorts every momentum into a bin using the edges. It first gives labels,
    # as.integer() turns them into the numbers 1 to 5 (= bin number).
    bin <- as.integer(cut(fileData$momentum, breaks = binEdges, include.lowest = TRUE))


    # Count, for every pair and every bin, how many WT and how many mutant there are.
    for (i in 1:nPairs) {
      for (b in 1:nBins) {
        # sum() of TRUE/FALSE counts the TRUEs. So this counts bacteria that have the
        # right ID AND (&) are in bin b. rbind() adds the new row at the bottom of the table.
        countRows <- rbind(countRows, data.frame(
          Set = k,
          Strain = pairName[i],
          MomentumBin = b,
          CountWT = sum(fileData$id == pairWT[i] & bin == b),
          CountMutant = sum(fileData$id == pairMutant[i] & bin == b)
        ))
      }
    }
  }
  options(scipen = 999)                     # stops R from writing big numbers like 1e+06 in the csv
  write.csv(countRows, countsFile, row.names = FALSE)
  write.csv(data.frame(Edge = binEdges), edgesFile, row.names = FALSE)
  cat("Saved", countsFile, "and", edgesFile, "\n")
}


counts <- read.csv(countsFile)              # from here on we always work with the saved csv


# Axis labels of the bins, e.g. "<0.65", "0.65-1.32", ">5.87"
# (just text that describes each bin, used in the csv and the README table)
binLabel <- character(nBins)                # an empty text list of length 5 to fill in
for (b in 1:nBins) {
  if (b == 1) {
    binLabel[b] <- paste0("<", round(binEdges[2], 2))          # first bin: everything below edge 2
  } else if (b == nBins) {
    binLabel[b] <- paste0(">", round(binEdges[nBins], 2))      # last bin: everything above the last real edge
  } else {
    binLabel[b] <- paste0(round(binEdges[b], 2), "-", round(binEdges[b + 1], 2))   # middle bins: from-to
  }
}




# --------------------------------------------------- Number of events per file ----
# sub_sample_results.csv has no "number of events" column, but
# total_count / average_per_event = number of non-empty events.
# E. coli WT (ID 211) is used because its average (about 20) is so large that the
# rounding of the average to 5 decimals makes no noticeable difference.
weekFour <- read.csv(eventsFile)
# The file name contains the set number (e.g. "output-Set7.txt"). This line cuts out just
# the number with a regular expression (the "\\1" means "keep the part in brackets").
weekFour$set <- as.integer(sub("^.*Set([0-9]+)\\.txt$", "\\1", weekFour$sub_sample))
ecoli <- weekFour[weekFour$code == 211, ]   # keep only the rows of E. coli WT
ecoli <- ecoli[order(ecoli$set), ]          # sort them by set number 1 to 10
# safety check: we must have exactly the sets 1 to 10, otherwise something is wrong
if (!identical(ecoli$set, 1:nFiles)) {
  stop("sub_sample_results.csv does not contain exactly the sets 1 to ", nFiles, " for ID 211.")
}
# number of events per file = total count / average per event (rounded to a whole number)
nEvents <- round(ecoli$total_count / ecoli$average_per_event)
cat("Non-empty events per file:", nEvents, "\n\n")




# ------------------------------------------ Put the counts in 3D arrays [file, pair, bin] ----
# The csv is a long table, but for the calculations it is easier to have a "block" of numbers
# (an array) where WT[k, i, b] is the count for file k, pair i and bin b. Like a cube.
WT <- array(0, dim = c(nFiles, nPairs, nBins))    # starts as all zeros, size 10 x 6 x 5
MUT <- array(0, dim = c(nFiles, nPairs, nBins))
for (k in 1:nFiles) {
  for (i in 1:nPairs) {
    for (b in 1:nBins) {
      # find the one row in the table that matches this file, pair and bin
      row <- counts[counts$Set == k & counts$Strain == pairName[i] & counts$MomentumBin == b, ]
      if (nrow(row) != 1) stop("Expected one row for set ", k, ", ", pairName[i], ", bin ", b)
      WT[k, i, b] <- row$CountWT            # put the numbers in the right place of the cube
      MUT[k, i, b] <- row$CountMutant
    }
  }
}


# Sanity check: all bins together must give the Week 4 totals (a warning if not)
# (If my bins lose or double count bacteria the totals would not match Week 4.)
for (k in 1:nFiles) {
  for (i in 1:nPairs) {
    for (s in 1:2) {                        # s = 1 means WT, s = 2 means mutant
      id <- c(pairWT[i], pairMutant[i])[s]
      week4Total <- weekFour$total_count[weekFour$set == k & weekFour$code == id]
      myTotal <- sum(if (s == 1) WT[k, i, ] else MUT[k, i, ])   # add up the 5 bins
      if (length(week4Total) == 1 && week4Total != myTotal) {
        warning("Set ", k, ", ID ", id, ": momentum counts (", myTotal,
                ") differ from the Week 4 total (", week4Total, ")")
      }
    }
  }
}




# ------------------------------------- Overall difference and asymmetry per pair ----
# Done in each file first (so the shared fluctuation of WT and mutant cancels),
# then the mean and the standard deviation over the 10 files are taken.
# apply(..., c(1, 2), sum) adds up the 5 bins, so we keep one number per file and pair.
WTtotal <- apply(WT, c(1, 2), sum)        # [file, pair]: all bins added up
MUTtotal <- apply(MUT, c(1, 2), sum)
D <- (WTtotal - MUTtotal) / nEvents       # difference per event, one value per file
A <- (WTtotal - MUTtotal) / (WTtotal + MUTtotal) * 100   # asymmetry in % (times 100)


# Build a results table, one row per pair.
# mean() = average over the 10 files, sd() = standard deviation over the 10 files.
# Significance = value / uncertainty (how many "sigmas" away from zero it is).
overall <- data.frame()
for (i in 1:nPairs) {
  overall <- rbind(overall, data.frame(
    Strain = pairName[i],
    ID = pairWT[i],
    Difference = mean(D[, i]),              # D[, i] = all 10 files for pair i
    DifferenceUnc = sd(D[, i]),
    DifferenceSig = mean(D[, i]) / sd(D[, i]),
    AsymmetryPercent = mean(A[, i]),
    AsymmetryPercentUnc = sd(A[, i]),
    AsymmetrySig = mean(A[, i]) / sd(A[, i])
  ))
}
# "significant" if the asymmetry is at least 3 sigma (abs() so that a negative value also counts)
overall$Significant <- ifelse(abs(overall$AsymmetrySig) >= 3, "yes", "no")




# ------------------------------------------ Asymmetry per momentum bin ----
# Same asymmetry formula as above, but now for every bin separately, so we get a cube
# again: Abin[file, pair, bin]. R does the calculation for all the numbers at once.
Abin <- (WT - MUT) / (WT + MUT) * 100     # [file, pair, bin], in %; NaN if a bin is empty


perBin <- data.frame()
for (i in 1:nPairs) {
  for (b in 1:nBins) {
    # na.rm = TRUE tells mean/sd to ignore missing values (NaN), in case a bin was empty
    m <- mean(Abin[, i, b], na.rm = TRUE)   # mean over the 10 files
    u <- sd(Abin[, i, b], na.rm = TRUE)     # standard deviation over the 10 files = error bar
    perBin <- rbind(perBin, data.frame(
      Strain = pairName[i],
      MomentumBin = b,
      MomentumRange = binLabel[b],
      AsymmetryPercent = m,
      AsymmetryPercentUnc = u,
      Significance = m / u
    ))
  }
}




# ---------------------------------------- Does the asymmetry depend on momentum? ----
# Two tests per pair: (a) is one constant value good enough for all 5 bins?
# (b) is there a slope (does A go up or down with momentum)?
binNumber <- 1:nBins                       # x-values for the straight line fit: bin 1, 2, 3, 4, 5
fits <- data.frame()
for (i in 1:nPairs) {
  rows <- perBin[perBin$Strain == pairName[i], ]   # the 5 rows of this pair
  y <- rows$AsymmetryPercent                       # the 5 asymmetry values
  u <- rows$AsymmetryPercentUnc                    # their 5 error bars
  w <- 1 / u^2                                     # weight: precise bins count more


  # (a) Constant fit: the best single value for all 5 bins is the weighted mean.
  #     chi2 measures how far the 5 bins are from that value, in units of their
  #     uncertainties. Rule of thumb: chi2 close to ndf (= 4) means a constant fits well.
  constant <- sum(w * y) / sum(w)                  # weighted mean
  constantUnc <- 1 / sqrt(sum(w))                  # its uncertainty (formula from the lecture / README)
  chi2 <- sum(((y - constant) / u)^2)              # chi-squared
  ndf <- nBins - 1                                 # degrees of freedom = 5 points - 1 fitted number
  # p-value: how likely is a chi2 this large (or larger) if a constant really is correct?
  pValue <- pchisq(chi2, ndf, lower.tail = FALSE)


  # (b) Slope: fit a straight line A = a + s * bin through the 5 bins of EACH file,
  #     then take the mean and the SD of the 10 slopes (sub-sampling method again).
  slopes <- numeric(nFiles)                        # an empty list of 10 numbers to fill in
  for (k in 1:nFiles) {
    # lm() fits a straight line. "Abin[k, i, ] ~ binNumber" reads as "A depends on bin number".
    # coef(...)["binNumber"] picks out the slope from the result.
    slopes[k] <- coef(lm(Abin[k, i, ] ~ binNumber))["binNumber"]
  }


  fits <- rbind(fits, data.frame(
    Strain = pairName[i],
    Constant = constant,
    ConstantUnc = constantUnc,
    ConstantSig = constant / constantUnc,
    Chi2 = chi2,
    ndf = ndf,
    pValue = pValue,
    Slope = mean(slopes),
    SlopeUnc = sd(slopes),
    SlopeSig = mean(slopes) / sd(slopes)
  ))
}




# ------------------------------------------------------------ Print and save ----
# options(width = 250) makes the printed tables wide so they do not wrap around.
options(width = 250, scipen = 999)
cat("==== Overall result, mean +- SD over", nFiles, "files ====\n")
print(overall, row.names = FALSE, digits = 4)
cat("\n==== Asymmetry A (%) per momentum bin ====\n")
print(perBin, row.names = FALSE, digits = 4)
cat("\n==== Constant fit and slope (A in % per bin) ====\n")
print(fits, row.names = FALSE, digits = 4)


# save the three result tables as csv files (these are the ones for the README / Canvas)
write.csv(overall, "overall_asymmetry.csv", row.names = FALSE)
write.csv(perBin, "asymmetry_per_bin.csv", row.names = FALSE)
write.csv(fits, "momentum_fits.csv", row.names = FALSE)




# ------------------------------------------------------------------ Plot ----
# One panel per pair. The x-axis is the real momentum (linear scale). Each point sits in
# the middle of its momentum bin, and the grey horizontal bar shows the width of the bin.
# The bins are open-ended at both sides, so for the plot only:
#   - the first bin starts at 0 (a momentum cannot be negative)
#   - the last bin is drawn as wide as the bin before it (a choice for the picture only;
#     the real last bin extends further, and the counts are not affected)
# Each panel has its own y-axis because the uncertainties differ a lot (abundant strains
# are measured much more precisely).
# (The real edges contain -Inf and Inf, which cannot be drawn, so we make a copy for the plot.)
plotEdges <- binEdges
plotEdges[1] <- 0
plotEdges[nBins + 1] <- binEdges[nBins] + (binEdges[nBins] - binEdges[nBins - 1])
binStart <- plotEdges[1:nBins]             # left side of each bin
binEnd <- plotEdges[2:(nBins + 1)]         # right side of each bin
binCentre <- (binStart + binEnd) / 2       # middle of each bin = where the point is drawn


# png() opens an image file; everything plotted after it goes into that file until dev.off().
png("asymmetry_vs_momentum.png", width = 2000, height = 1200, res = 150)
# par() sets general plot options:
#   mfrow = c(2, 3): 2 rows and 3 columns of small plots
#   mar = c(bottom, left, top, right): white space around each plot, in lines of text.
#     FIX: the top margin was 4.5, now 5.5, so there is room for three lines of text
#     above each plot (title + the two lines with the fit results).
#   oma: extra space around the whole figure (for the big title at the top)
#   las = 1: numbers on the y-axis are written horizontally
par(mfrow = c(2, 3),
    mar = c(7, 5, 5.5, 1),
    oma = c(0, 0, 2.5, 0),
    las = 1)

for (i in 1:nPairs) {
  rows <- perBin[perBin$Strain == pairName[i], ]
  fit <- fits[fits$Strain == pairName[i], ]
  y <- rows$AsymmetryPercent
  low <- y - rows$AsymmetryPercentUnc      # bottom end of each error bar
  high <- y + rows$AsymmetryPercentUnc     # top end of each error bar


  # y-axis range: all error bars and the zero line, plus empty room on top for the legend
  yMin <- min(low, 0, fit$Constant - fit$ConstantUnc)
  yMax <- max(high, 0, fit$Constant + fit$ConstantUnc)
  yMax <- yMax + 0.4 * (yMax - yMin)       # add 40% extra space on top so the legend does not cover data


  # type = "n" draws an EMPTY plot (only axes). We add the real things afterwards.
  # FIX: main is now empty, because the title is written with mtext() below. Before, the
  # title was drawn on top of the two lines of text with the fit results.
  plot(binCentre, y, type = "n", xlim = c(0, plotEdges[nBins + 1]), ylim = c(yMin, yMax),
       xlab = "", ylab = "Asymmetry A (%)", main = "")
  # mtext() writes text in the margin. side = 1 is below the plot, 3 is above it.
  # "line" says how far away from the plot (higher number = further away).
  mtext(expression("Momentum |p| (" * 10^-20 * " kg m/s)"),
      side = 1, line = 4, cex = 0.8)


  # The two lines of text above the plot with the fit results.
  # sprintf() fills the numbers into the text: %.2f means "a number with 2 decimals".
  mtext(sprintf("constant: A = %.2f +/- %.2f %% (%.1f sigma), chi2/ndf = %.1f/%d",
                fit$Constant, fit$ConstantUnc, fit$ConstantSig, fit$Chi2, fit$ndf),
        side = 3, line = 1.2, cex = 0.65)
  mtext(sprintf("slope = %.3f +/- %.3f %% per bin (%.1f sigma)",
                fit$Slope, fit$SlopeUnc, fit$SlopeSig),
        side = 3, line = 0.3, cex = 0.65)
  # FIX: the name of the strain (the title) goes on line 2.6, above the two lines of text.
  # font = 3 means italic (nice for bacteria names), cex = 1.1 makes it a bit bigger.
  mtext(pairName[i], side = 3, line = 2.6, font = 3, cex = 1.1)


  # rect() draws the light red rectangle = the +- 1 sigma band around the constant fit.
  # rgb(1, 0, 0, 0.12) is red (1,0,0) with a lot of see-through (0.12 = almost transparent).
  rect(0, fit$Constant - fit$ConstantUnc, plotEdges[nBins + 1], fit$Constant + fit$ConstantUnc,
       col = rgb(1, 0, 0, 0.12), border = NA)                    # +-1 sigma band
  abline(h = 0, lty = 3, col = "grey40")                         # no asymmetry (dotted line at 0)
  abline(h = fit$Constant, col = "red", lwd = 2)                 # constant fit (thick red line)
  # arrows() with angle = 90 and code = 3 makes a line with a small flat cap on both ends,
  # which is how I get "error bars". First one: the horizontal grey bar = width of the bin.
  arrows(binStart, y, binEnd, y, angle = 90, code = 3, length = 0.03,
         col = "grey60")                                         # width of each bin
  # Second one: the vertical black error bars (from low to high).
  arrows(binCentre, low, binCentre, high, angle = 90, code = 3, length = 0.04)
  points(binCentre, y, pch = 19, col = "darkblue")               # the dots (pch = 19 is a filled circle)


  # Label each bin with its actual momentum range
binLabels <- character(nBins)

for (b in 1:nBins) {
  if (b == 1) {
    binLabels[b] <- paste0("<", round(binEdges[2], 2))
  } else if (b == nBins) {
    binLabels[b] <- paste0(">", round(binEdges[nBins], 2))
  } else {
    binLabels[b] <- paste0(
      round(binEdges[b], 2), "-", round(binEdges[b + 1], 2)
    )
  }
}

# Label momentum bins; lower only the second label
labelLines <- c(2, 3, 2, 2, 2)

mtext(binLabels,
      side = 1,
      at = binCentre,
      line = labelLines,
      cex = 0.55,
      col = "darkblue")

  # The legend (the little box that explains the symbols). The four entries in legend =
  # line up with the four pch / lty / lwd / col values below (first with first, etc.)
  legend("topleft", bty = "n", cex = 0.65,
         legend = c("Mean of 10 files +/- SD", "Momentum range of bin (bin 5 open-ended)",
                    "Constant fit (band: +/- 1 sigma)", "A = 0 (no asymmetry)"),
         pch = c(19, NA, NA, NA), lty = c(NA, 1, 1, 3), lwd = c(NA, 1, 2, 1),
         col = c("darkblue", "grey60", "red", "grey40"))
}
# The big title over the whole figure (outer = TRUE means in the outer margin set by oma)
mtext(paste0("Asymmetry A = (WT - mutant) / (WT + mutant) vs momentum, all ", nFiles, " sets"),
      outer = TRUE, side = 3, line = 0.5, font = 2)
dev.off()                                  # closes the png file, only now is the image really saved


cat("\nWrote overall_asymmetry.csv, asymmetry_per_bin.csv, momentum_fits.csv and",
    "asymmetry_vs_momentum.png\n")