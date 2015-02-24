
# fMRI Tutorial 1

getwd()
setwd("D:/Coursera_R/fRMI/MoAEpilot")

# Install package 'AnalyzeFMRI'. 
install.packages('AnalyzeFMRI')

# Begin by reading and visualizing the structural scan contained in the subdirectory 'sM00223'¯. We can use the function 
# f.read.analyze.volume from the AnalyzeFMRI package.

# Load the package
library(AnalyzeFMRI)

# Read in the data
img.filename = "./sM00223/sM00223_002.img"
# img.filename = "D:/Coursera_R/fRMI/MoAEpilot/sM00223/sM00223_002.img"
img = f.read.analyze.volume(img.filename)

# Check the dimensions of the data
dim(img)
# 256 256  54   1

# The structural image is now contained in the variable entitled img. Its dimensions are 256x256x54. We can visualize the
# image using the R-function image. For example, if we want to study slice number 30:
image(img[, , 30, 1], axes=FALSE)

# Turn to the functional data consisting of 96 files (one for each time point), contained in the subdirectory 'fM00223'.
# Read in the data from a single time point
fimg.filename =  "./fM00223/fM00223_010.img"
# fimg.filename =  "D:/Coursera_R/fRMI/MoAEpilot/fM00223/fM00223_010.img"
fimg = f.read.analyze.volume(fimg.filename)

# Check the dimensions of the data
dim(fimg)
# 64 64 64  1

# The structural image is now contained in the variable entitled fimg. Its dimensions are 64x64x64. Note these dimensions are
# not the same as the structural scan which has a higher spatial resolution. We can again visualize the image using the 
# R-function image. For example, if we want to study slice number 40:
image(img[, , 40, 1], axes=FALSE)

# Typically, we want to work on data from all time points in the experiment. Therefore, let us now read in all the functional
# data for a certain slice.
# List all the files in the subdirectory that end with .img
files = dir("./fM00223", pattern = "*.img", full.names = TRUE)
#files = dir("D:/Coursera_R/fRMI/MoAEpilot/fM00223", pattern = "*.img", full.names = TRUE)

# Choose a slice.
slice = 40

# Define a data array that will consist of 96 different 64x64 images.
data = array(0, dim = c(64, 64, 96))

# Cycle through the files, read in the data and place it into the data array
for (t in 1: length(files)) {
    
    file = files[t];
    
    img = f.read.analyze.volume(file)[, , slice, 1]
    
    data[, , t] <- img
    
}

# Now the variable data consists of the 40th slice sampled at each of the 96 time points. To look at the slice at time 30:
image(data[, , 30])

# To study the time series from a specific voxel (say the one at x=20 and y=20) type:
plot(data[20, 20, ], type='l')






