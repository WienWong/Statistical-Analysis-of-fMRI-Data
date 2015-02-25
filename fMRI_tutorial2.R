
# fMRI Tutorial 2 -- Simple GLM Analysis

getwd()
setwd("D:/Coursera_R/fRMI/MoAEpilot_preproc")

install.packages('fmri')
library('fmri')

# Reading the Data

# Let's begin by reading the functional scan contained in the subdirectory 'fM00223', i.e. the files labeled 'swrfM00223_****'
data = read.ANALYZE("./fM00223/swrfM00223_", numbered = TRUE, picstart = 4, numbpic = 96)

# This function reads in 96 files starting with 'swrfM00223_0004'. It returns a list of class "fmridata". The 4-dimensional 
# data (64x64x64x96) can be extracted from the list using the command:
raw.data = extract.data(data)

#In addition the list contains basic information about the size of the data, the voxel size, and complete header information. 

# Defining a Model

# Next, we want to set up the design matrix to use in the GLM analysis. The predicted BOLD response is given by the convolution
# of the hemodynamic response function and the task indicator function. In our data set we should use:
hrf = fmri.stimulus(96, onsets = seq(7, 96, by=12), durations = 6, rt = 7)

# This creates the predicted BOLD response for a stimulus consisting of 96 scans, with onset times starting at the 7th scan and
# continuing every 12th point. The duration of each stimulus is 6 scans, and the time between subsequent scans is 7 s. Finally,
# the stimulus is convolved with a canonical HRF shape.

# To view the predicted BOLD response type:
plot(hrf, type = 'l', main="Predicted BOLD Response Type")

# In order to create the design matrix needed in the GLM analysis we can type:
X = fmri.design(hrf)

# By default the design matrix will include polynomial drift terms up to quadratic order. Thus the resulting design matrix will
# have dimensions 96x4. Use the command image(X) to further explore the structure of the design matrix.

image(X, main="Design Matrix \n (dim. of 96x4)")

# If we had multiple stimuli we would need to create separate predicted BOLD responses for each stimuli based on their onsets 
# and durations (e.g., hrf1, hrf2, etc.), and create a design matrix using the commands cbind() and fmri.design() together (e.g.
# fmri.design(cbind(hrf1, hrf2, ...))).

# Similarly, suppose we want to include the motion parameters (or any other covariate of interest) as regressors in our design 
# matrix. This can be done as follows:
M = read.table("fM00223/motion.txt")

X2 = fmri.design(cbind(hrf, M))

image(X2, main="Design Matrix incl. Motion Parameters")

# The first command reads in the motion parameters and the second creates a design matrix where the first column is the 
# predicted BOLD response (hrf), columns 2-7 are the motion parameters (M), and 8-10 are the polynomial drift terms. Hence, 
# X2 is now a 96x10 matrix. To study the estimated motion parameters type:
par(mfrow = (c(1, 2)))
matplot(M[, 1:3], type = 'l', main = "Translation Parameters")
matplot(M[, 4:6], type = 'l', main = "Rotation Parameters")

# This plots the translation and rotation parameters in separate plots.

# Estimating the Model

# The GLM can now be estimated using the command:

spm = fmri.lm(data, X2)

## Below information came from my RStudio console
# fmri.lm: entering function
# fmri.lm: calculating AR(1) model
# 0% . 10% . 20% . 30% . 40% . 50% . 60% . 70% . 80% . 90% . 
# fmri.lm: smoothing with (hmax): 3.52 
# fmri.lm: finished
# fmri.lm: re-calculating linear model with prewithened data
# 0% . 10% . 20% . 30% . 40% . 50% . 60% . 70% . 80% . 90% . 
# fmri.lm: calculating spatial correlation
# fmri.lm: determining df: 81.18843 
# fmri.lm: exiting function
# Warning message:
#     In fmri.lm(data, X2) :
#     Local smoothness characterized by large bandwidth  4  check residuals for structure,Local smoothness characterized by 
#     large bandwidth  4  check residuals for structure,Local smoothness characterized by large bandwidth  4  check residuals 
#     for structure
## Above information came from my RStudio console

summary(spm)
# Object of class fmrispm created by
# fmri.lm(data, X2)
# Data Dimension               : 79 95 68 96 
# Range of estimated parameters: -823.8255 ... 956.123 
# File(s)                      : ./fM00223/swrfM00223_004.img..../fM00223/swrfM00223_099.img 
# 
# Design Dimension             : 96 10 
# Prewhitening performed with smoothed map
# of autocorrelation parameter in AR(1) model for time series!

class(spm)
# "fmridata" "fmrispm" 

######## U don't need to run below chunk of code here ######### 

# where data is the data object obtained using read.ANALYZE() , and X2 is the design matrix created above. The term spm will 
# contain all the relevant information about the analysis. For the case where we have multiple stimuli we can also use the 
# command fmri.lm() to define contrasts of interest. For example, consider the case where we have two stimuli of interest 
# (hrf1 and hrf2). We begin by creating the design matrix:

X = fmri.design(cbind(hrf1, hrf2))

# Next, we can potentially look at the following contrasts:
# A contrast for stimulus 1 only
spm1 = fmri.lm(data, X, contrast = c(1,0))

# A contrast for stimulus 2 only
spm2 = fmri.lm(data, X, contrast = c(0,1))

# A contrast for contrast comparing the two
spm3 = fmri.lm(data, X, contrast = c(1,-1))

######## U don't need to run above chunk of code here #########

# For the auditory data set, where we only have one stimulus, then the contrast by default corresponds to that variable (i.e. 
# hrf). 

# Once we have estimated the model, we can compute t-values for the contrast of interest in each voxel of the brain using the 
# following command:
T = spm$cbeta/sqrt(spm$var)

#how to view beta???
beta <- spm$beta
class(beta)
# "array"
dim(beta)
# 79 95 68 10

summary(T)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.8230 -0.4413  0.2744  0.3610  1.0290 18.8500 

# We can visualize the t-map corresponding to a particular slice (e.g., slice 30) using the command:
dev.off() # reset the plot Or
# par(mfrow = (c(1, 1)))

image(T[, , 30], main = "T-map w.r.t. Slice 30")

par(mfrow = (c(1, 2)))
image(T[, , 30], main = "T-map w.r.t.\n Slice 30")
image(T[, , 40], main = "T-map w.r.t.\n Slice 40")


par(mfrow = (c(1, 2)))
image(X, main="Design Matrix")
image(X2, main="Design Matrix incl.\n Motion Parameters")


