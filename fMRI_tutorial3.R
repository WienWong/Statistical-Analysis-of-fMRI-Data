
# fMRI Tutorial 3 -- Threshold Statistical Images

getwd()
setwd("D:/Coursera_R/fRMI/MoAEpilot_preproc2")

# Install packages
install.packages('oro.nifti')
library('fmri')
library('AnalyzeFMRI')
library('oro.nifti')

# Compute the Statistical Map

# Note the dimensions of the data will be slightly different from those in the earlier tutorial, as we aren't using normalized
# data.
data = read.ANALYZE("./fM00223/srfM00223_", numbered = TRUE, picstart = 4, numbpic = 96)

hrf = fmri.stimulus(96, onsets = seq(7, 96, by=12), durations = 6, rt = 7)

M = read.table("./fM00223/motion.txt")

X2 = fmri.design(cbind(hrf, M))
# image(X2)

spm = fmri.lm(data, X2)

### Below were just software processing message ### 
# fmri.lm: entering function
# fmri.lm: calculating AR(1) model
# 0% . 10% . 20% . 30% . 40% . 50% . 60% . 70% . 80% . 
# fmri.lm: smoothing with (hmax): 3.52 
# fmri.lm: finished
# fmri.lm: re-calculating linear model with prewithened data
# 0% . 10% . 20% . 30% . 40% . 50% . 60% . 70% . 80% . 90% . 
# fmri.lm: calculating spatial correlation
# fmri.lm: determining df: 81.18843 
# fmri.lm: exiting function
# Warning message:
# In fmri.lm(data, X2) :
#  Local smoothness characterized by large bandwidth 4 check residuals for structure,Local smoothness characterized by large
#  bandwidth 4 check residuals for structure,Local smoothness characterized by large bandwidth 4 check residuals for structure
### Above were just software processing message ### 

T = spm$cbeta/sqrt(spm$var)

# The variable T now consists of the t-statistic computed at each voxel of the brain. Our goal is to find an appropriate
# threshold for this statistical map so that we can determine which voxels should be considered active.

# Threshold the Statistical Map 

# Let's begin by reading in the anatomical image rsM00223_002.img, contained in the directory sM00223. We will use this image
# to present our results. We can read in the image using the function readNIfTI from the oro.nifti package:
ana = readNIfTI('./sM00223/rsM00223_002.img')

# To study the anatomical image in closer detail:
image(ana) # which plots all the axial slices of the brain.

# Now we are ready to threshold the statistical image T. In this tutorial we will use three different types of threshold: 
# (i) no correction for multiple comparisons; (ii) an FDR-corrected threshold; and (iii) a Bonferroni-corrected threshold.

# If we ignore the multiple comparisons problem all together, and threshold using the value appropriate for a single test, 
# we could use the following threshold:
thr1 = qt(0.95, spm$df)

# which represents the 95% quantile for a t-distribution with df degrees of freedom. Here the degrees of freedom were computed
# by the function fmri.lm and incorporate correction for the auto-correlation present in the data. To overlay the results onto
# our anatomical image type:
overlay(ana, ifelse(T > thr1, T, NA))

# The voxels whose t-statistic lie above the threshold should appear in red or yellow. Note there are many false positives
# throughout the brain.

# To compute the appropriate Bonferroni threshold we can use the command Threshold.Bonferroni from the AnalyzeFMRI package. 
# To compute this threshold and overlay significant voxels on the anatomical image type the commands:
thr2 = Threshold.Bonferroni(0.05, 64^3, type = c("t"), df1 = spm$df)

overlay(ana, ifelse(T > thr2, T, NA))

# Note here we correct for the fact that we performed 64^3 different tests, one for each voxel. In reality we should remove
# all voxels outside of the brain before computing this threshold, but this gives you a genral idea of the procedure. 

# Finally, to compute the FDR-corrected threshold we can use the command Threshold.FDR from the AnalyzeFMRI package. 
# To compute this threshold and overlay significant voxels on the anatomical image type the commands:
thr3 = Threshold.FDR(T,q=0.05,cV.type=2,type='t',df1=spm$df)

overlay(ana, ifelse(T > thr3, T, NA))

