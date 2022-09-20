#!/usr/bin/env Rscript
# from https://rpubs.com/JanpuHou/300168
# Non-negative matrix factorization for topic modeling
# Janpu Hou
# August 16, 2017
# this guy has quite alot of rpubs. not sure how great his stuff is though.
# the datas ets are scattered and often not available.

# Simple Algebra Example of NMF (uses artificial 5x4 matrix)

# Non-negative Matrix Factorization (NMF) is a state of the art feature extraction algorithm. NMF is useful when there are many attributes and the attributes are ambiguous or have weak predictability. By combining attributes, NMF can produce meaningful patterns, topics, or themes. The unlabeled document or text collections are becoming larger and larger which is common and obvious; mining such data sets are a challenging task. During model apply, an NMF model maps the original data into the new set of attributes (features) discovered by the model.

# http://meyer.math.ncsu.edu/Meyer/Talks/SIAMSEAS_NMF.pdf. http://csweb.cs.wfu.edu/~pauca/publications/SIAMDM03.pdf.

# Both NMF and pLSA are instances of multinomial PCA (Buntine, 2002). pLSA is NMF with KL-divergence (Gaussier and Goutte, 2005). NMF can help estimates the parameters of the pLSA model.

# data set http://cogsys.imm.dtu.dk/toolbox/nmf/.

library(NMF)
library(ggplot2)
library(Cairo)

# create an artificial 5xr45 matrix
x1 <- c(5,4,1,1)
x2 <- c(4,5,1,1)
x3 <- c(1,1,5,5)
x4 <- c(1,1,4,5)
x5 <- c(1,1,5,4)
R <- as.matrix(rbind(x1,x2,x3,x4,x5))

# So for NMF, we're looking at R ~= WH
# calling NMF is this simple, but you end up none the wiser.
res <- nmf(R, 4,"lee") # I suppose to try.
res <- nmf(R, 3,"lee") 
V.hat <- fitted(res) # give a stab at recovering the original matrix from nmf() output
# print(V.hat) 

##         [,1]      [,2]      [,3]      [,4]
## x1 4.9992998 4.0007002 1.0186974 0.9813026
## x2 4.0007001 4.9992998 0.9813026 1.0186974
## x3 0.9999985 1.0000015 4.9999999 5.0000000
## x4 0.9813035 1.0186965 4.4992999 4.5007001
## x5 1.0186982 0.9813018 4.5007002 4.4992998

# For R ~=WH, basis() gives you W and coef() gives you H, I think.
w <- basis(res) # An NMF func ... W  user feature matrix matrix
# dim(w) # n x r (n= 5  r = 4)
# [1] 5 3
# print(w)
##          [,1]        [,2]       [,3]
## x1 0.59480680 0.002452533 0.30131218
## x2 0.12644947 0.034013168 0.53288729
## x3 0.09387198 0.344962625 0.05246055
## x4 0.08366655 0.309876765 0.06100590
## x5 0.10120521 0.308694909 0.05233407

h <- coef(res) # NMF func ... actually part of basis(): H  movie feature matrix
# dim(h) #  r x p (r = 4 p = 4)
## [1] 3 4
# print(h) 
##           [,1]      [,2]       [,3]       [,4]
## [1,] 5.2473068 2.2735776  1.3363828  1.2250250
## [2,] 0.5236661 0.9446866 14.0350786 14.0508446
## [3,] 6.2290272 8.7817357  0.6285386  0.7241304

# Something less artificial
# Test on Real Webpages
# Using the example in previous example: http://rpubs.com/JanpuHou/299832.

# don't have this.
m <- read.csv(file="D:/R_Files/corpus/tdm.csv")
head(m)

##                X d1.txt d2.txt d3.txt
## 1     additional      1      0      0
## 2 administrative      1      0      0
## 3        affairs      1      1      0
## 4       affected      2      0      0
## 5      affecting      1      0      0
## 6      afternoon      1      0      0

rownames(m) <- m[,1]
m[,1] <- NULL

res <- nmf(m, 3,"KL") 

w <- basis(res) #  W  user feature matrix matrix
dim(w)

## [1] 622   3

df <- as.data.frame(w)
head(df,10)

##                          V1           V2       V3
## additional     2.220446e-16 2.220446e-16 16.15177
## administrative 2.220446e-16 2.220446e-16 16.15177
## affairs        2.220446e-16 1.364603e+01 16.15177
## affected       2.220446e-16 2.220446e-16 32.30354
## affecting      2.220446e-16 2.220446e-16 16.15177
## afternoon      2.220446e-16 2.220446e-16 16.15177
## also           2.220446e-16 4.093808e+01 16.15177
## although       2.220446e-16 2.220446e-16 16.15177
## amid           2.220446e-16 2.220446e-16 16.15177
## anantharaman   2.220446e-16 2.220446e-16 16.15177

df$total <- rowSums(df)
df$word<-rownames(df)
colnames(df) <- c("doc1","doc2","doc3","total","word")
df <-df[order(-df$total),] 
head(df,20)

##                    doc1         doc2         doc3     total       word
## taiwan     1.182388e+02 8.187616e+01 1.130624e+02 313.17734     taiwan
## august     2.220446e-16 2.046904e+02 2.220446e-16 204.69040     august
## said       2.220446e-16 4.093808e+01 1.615177e+02 202.45579       said
## power      2.220446e-16 2.220446e-16 1.938213e+02 193.82125      power
## chinese    1.970646e+01 8.187616e+01 2.220446e-16 101.58262    chinese
## foundation 9.853232e+01 2.220446e-16 2.220446e-16  98.53232 foundation
## heritage   9.853232e+01 2.220446e-16 2.220446e-16  98.53232   heritage
## taiwans    2.220446e-16 5.458411e+01 3.230354e+01  86.88765    taiwans
## relations  5.911939e+01 2.729205e+01 2.220446e-16  86.41145  relations
## government 3.941293e+01 1.364603e+01 3.230354e+01  85.36250 government
## president  3.941293e+01 1.364603e+01 3.230354e+01  85.36250  president
## air        2.220446e-16 8.187616e+01 2.220446e-16  81.87616        air
## blackout   2.220446e-16 2.220446e-16 8.075886e+01  80.75886   blackout
## director   7.882586e+01 2.220446e-16 2.220446e-16  78.82586   director
## min        5.911939e+01 2.220446e-16 1.615177e+01  75.27116        min
## read       5.911939e+01 2.220446e-16 1.615177e+01  75.27116       read
## security   5.911939e+01 2.220446e-16 1.615177e+01  75.27116   security
## aircraft   2.220446e-16 6.823013e+01 2.220446e-16  68.23013   aircraft
## defense    2.220446e-16 6.823013e+01 2.220446e-16  68.23013    defense
## caused     2.220446e-16 1.364603e+01 4.845531e+01  62.10134     caused

wordMatrix = as.data.frame(w)
wordMatrix$word<-rownames(wordMatrix)
colnames(wordMatrix) <- c("doc1","doc2","doc3","word")


# Topic 1
newdata <-wordMatrix[order(-wordMatrix$doc1),] 
head(newdata)

##                 doc1         doc2         doc3       word
## taiwan     118.23879 8.187616e+01 1.130624e+02     taiwan
## foundation  98.53232 2.220446e-16 2.220446e-16 foundation
## heritage    98.53232 2.220446e-16 2.220446e-16   heritage
## director    78.82586 2.220446e-16 2.220446e-16   director
## min         59.11939 2.220446e-16 1.615177e+01        min
## read        59.11939 2.220446e-16 1.615177e+01       read

d <- newdata
df <- as.data.frame(cbind(d[1:10,]$word,as.numeric(d[1:10,]$doc1)))
colnames(df)<- c("Word","Frequency")

# for ggplot to understand the order of words, you need to specify factor order

df$Word <- factor(df$Word, levels = df$Word[order(df$Frequency)])
ggplot(df, aes(x=Word, y=Frequency)) + 
  geom_bar(stat="identity", fill="lightgreen", color="grey50")+
  coord_flip()+
  ggtitle("Topic 1")

# Topic 2
newdata <-wordMatrix[order(-wordMatrix$doc2),] 
head(newdata)

##                  doc1      doc2         doc3     word
## august   2.220446e-16 204.69040 2.220446e-16   august
## taiwan   1.182388e+02  81.87616 1.130624e+02   taiwan
## air      2.220446e-16  81.87616 2.220446e-16      air
## chinese  1.970646e+01  81.87616 2.220446e-16  chinese
## aircraft 2.220446e-16  68.23013 2.220446e-16 aircraft
## defense  2.220446e-16  68.23013 2.220446e-16  defense

d <- newdata
df <- as.data.frame(cbind(d[1:15,]$word,as.numeric(d[1:15,]$doc2)))
colnames(df)<- c("Word","Frequency")
df$Word <- factor(df$Word, levels = df$Word[order(df$Frequency)])
ggplot(df, aes(x=Word, y=Frequency)) + 
  geom_bar(stat="identity", fill="lightgreen", color="grey50")+
  coord_flip()+
  ggtitle("Topic 2")

# Topic 3
newdata <-wordMatrix[order(-wordMatrix$doc3),] 
head(newdata)

##                  doc1         doc2      doc3     word
## power    2.220446e-16 2.220446e-16 193.82125    power
## said     2.220446e-16 4.093808e+01 161.51771     said
## taiwan   1.182388e+02 8.187616e+01 113.06240   taiwan
## blackout 2.220446e-16 2.220446e-16  80.75886 blackout
## caused   2.220446e-16 1.364603e+01  48.45531   caused
## corp     2.220446e-16 2.220446e-16  48.45531     corp

d <- newdata
df <- as.data.frame(cbind(d[1:15,]$word,as.numeric(d[1:15,]$doc3)))
colnames(df)<- c("Word","Frequency")
df$Word <- factor(df$Word, levels = df$Word[order(df$Frequency)])
ggplot(df, aes(x=Word, y=Frequency)) + 
  geom_bar(stat="identity", fill="lightgreen", color="grey50")+
  coord_flip()+
  ggtitle("Topic 3")

