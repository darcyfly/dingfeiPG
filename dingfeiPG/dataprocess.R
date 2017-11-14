# read in image and plot
library(jpeg)
img = readJPEG("./23.jpg")
img = img[,,1]
par(mar = c(0,0,0,0))
plot(0,0,type = "n")
rasterImage(img, -1, -1, 1, 1)

# Define the set of elements that are not censored
#get the dim of the img matrix
m = dim(img)
n = dim(img)
#generate the Omega set
set.seed(100)
index = sample(c(1:length(img)),round(length(img)*0.4))
Omega = setdiff(1:length(img),index)

#plot the censored plot using the Omege set
Pimg = img
Pimg[index] = 0
par(mar = c(0,0,0,0))
plot(0,0,type = "n")
rasterImage(Pimg, -1, -1, 1, 1)
Pcomplete = dingfeiPG(Pimg,lambda = 0.23,Omega = Omega,index = index,maxinteration = 100)
par(mar = c(0,0,0,0))
plot(0,0,type = "n")
rasterImage(Pcomplete,-1,-1,1,1)

mm=svd(Pimg)
mm$d
