## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  message = FALSE,
  warning = FALSE,
  fig.align = "right",     # aligns logo right
  out.width = "120px"      # adjust logo size
)
options(width = 80)

## ----echo=TRUE, fig.width=8, fig.height=3, out.width="100%", results='hide'----
library(BAYSFWATCH)

#Generate a Synthetic Dataset or use your own data
#A set of outliers have been injected in Ts p x n matrices
#outv allows to pick a value for the value to the outliers
#outc is  an integer (1 to 4) defining the type of outlier injection:
#   - 1: Adds `outv` to all values in the matrix.
#   - 2: Adds `outv` to a single element (1,1) in the matrix.
#   - 3: Adds `outv` to `outn` randomly chosen elements in the matrix.
#   - 4: Adds `outv` to `Rcn` columns and `Rrn` rows, randomly selected.
# outn is the number of outliers to inject (only used in options 3 and 4).
#Must be lower than `pxn`.
#see the documentation to set the other parameters
?GenData

set.seed(42)
Data<-GenData(outv = 5, outc = 3, outn = 6, Rcn = 3, Rrn = 2)

matrix1 <- Data$X3D[,,79]
matrix2 <- Data$X3D[,,80]
matrix3 <- Data$X3D[,,81]

matrix1_no_outliers <- Data$X3D0[,,79]
matrix2_no_outliers <- Data$X3D0[,,80]
matrix3_no_outliers <- Data$X3D0[,,81]

# Define the common color scale
zlim <- range(c(matrix1, matrix2, matrix3))


# Plot the three matrices
par(mfrow = c(1, 3),    # 1 row, 3 columns
    mar = c(2, 2, 3, 2), # margins: bottom, left, top, right
    oma = c(0, 0, 2, 0), # outer margins (optional)
    asp = 1,   # square aspect
    xaxs = "i", yaxs = "i")# force square aspect ratio

# Use blue or green scale (e.g., topographical blues)
blue_col <- colorRampPalette(c("white", "lightblue", "blue"))(100)

# Plot
image(matrix1, col = blue_col, main = "Matrix 79", zlim = zlim)
image(matrix2, col = blue_col, main = "Matrix 80", zlim = zlim) #outliers
image(matrix3, col = blue_col, main = "Matrix 81", zlim = zlim)

# Plot Outliers
zlim <- range(c(matrix1-matrix1_no_outliers,
                matrix2-matrix2_no_outliers,
                matrix3-matrix3_no_outliers))

image(matrix1-matrix1_no_outliers, col = blue_col, main = "Matrix 79", zlim = zlim)
image(matrix2-matrix2_no_outliers, col = blue_col, main = "Matrix 80 - Outliers", zlim = zlim) 
image(matrix3-matrix3_no_outliers, col = blue_col, main = "Matrix 81", zlim = zlim)

## ----results='hide', fig.width=8, fig.height=3, out.width="100%"--------------
#Example 1: Outlier Detection 
#see the documentation to set the other parameters
?OutlierDetection

results <- OutlierDetection(Data$X3D,  M = Data$M, Sl = Data$Sl, V =  Data$V)


## ----results='hide', fig.width=8, fig.height=3, out.width="100%"--------------
#Example 2: Outlier Detection with MLE and one single outlier entry (outc=2).
#Note that when we introduce one single outlier with moderate absolute value, 
#the uncertainty about the nature of the outlier is higher.

Data<-GenData(outv = 1, outc = 2, outn = 6, Rcn = 3, Rrn = 2)
results <- OutlierDetection(Data$X3D)

## ----results='hide', fig.width=8, fig.height=3, out.width="100%"--------------
#Benchmark Datasets: Three benchmark datasets are available:
#-"EUdata" array 3x11x250
#-"Tradedata" array 27x27x22
#-"Volatilitydata" array 50x50x145

data("Volatilitydata")
?Volatilitydata

