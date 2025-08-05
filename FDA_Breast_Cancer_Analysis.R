# Load required package
library(fda)

# Load the Breast Cancer dataset
data(Cancerrate)

# Plot the raw data
plot(Cancerrate, main = "Raw Plot of Breast Cancer Rate by Age Group")

# View the dataset (optional, for checking structure)
View(Cancerrate)

# Extract the Numerical matrix from the functional data object
cancer = as.matrix(Cancerrate$y)
dim(cancer)

# Define Age Range and create the B-Spline Basis
AgeRng = c(45, 85)  # Age groups from 45 to 85+
Age = seq(45, 85, 5) # Every 5 years interval
norder = 4
nbasis = length(Age) + norder - 2
cancerbasis = create.bspline.basis(AgeRng, nbasis, norder)

# Set up the harmonic acceleration operator
harmaccelLfd = int2Lfd(max(0, norder - 2))

# Choosing Smoothing Parameter lambda by Generalized Cross-Validation (GCV)
loglam = seq(0, 2, .25)  # test range of lambda
nlam = length(loglam)
dfsave = rep(NA, nlam)
names(dfsave) = loglam
gcvsave = dfsave

for (ilam in 1:nlam) {
  cat(paste('log10 lambda =', loglam[ilam], '\n'))
  lambda = 10^loglam[ilam]
  fdParobj = fdPar(cancerbasis, harmaccelLfd, lambda)
  smoothlist = smooth.basis(Age, cancer, fdParobj)
  dfsave[ilam] = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}

# Plot GCV values to find the best lambda
plot(loglam, gcvsave, type = 'o', lwd = 2, 
     xlab = expression(log[10](lambda)), 
     ylab = expression(GCV(lambda)), 
     main = "GCV Curve for Selecting Smoothing Parameter")

# Smooth the data using the best lambda (choose the lambda minimizing GCV)
lambda = 10^0.8  # Assuming the best based on reference pattern

cancerfdParobj = fdPar(cancerbasis, harmaccelLfd, lambda)
cancer.fit = smooth.basis(Age, cancer, cancerfdParobj)
cancer.fd = cancer.fit$fd
fdnames = list("Age Group", "Years", "Breast Cancer Rate")
cancer.fd$fdnames = fdnames

# Plot the smoothed functional data
plot(cancer.fd, lwd = 1.5, main = "Smoothed Functional Data for Breast Cancer Rates")

# --- Descriptive Statistics ---

# Mean function
meancancer = mean.fd(cancer.fd)
plot(meancancer, lwd = 2, 
     xlab = "Age Group", ylab = "Breast Cancer Rate", 
     main = "Mean Function of Breast Cancer Rate")

# Standard deviation function
stddevcancer = std.fd(cancer.fd)
plot(stddevcancer, lwd = 2, 
     xlab = "Age Group", ylab = "Breast Cancer Rate", 
     main = "Standard Deviation Function")

# --- Covariance Analysis ---

# Bivariate covariance surface
cancervar.bifd = var.fd(cancer.fd)
cancervar_mat = eval.bifd(Age, Age, cancervar.bifd)

# 3D Perspective Plot of Covariance
persp(Age, Age, cancervar_mat, theta = -75, phi = 20, r = 7, expand = 0.5,
      ticktype = "detailed",
      xlab = "Age Group", ylab = "Age Group", zlab = "Variance (Breast Cancer Rate)",
      main = "3D Variance-Covariance Surface",
      col = rainbow(4000))

# Contour Plot of Covariance
contour(Age, Age, cancervar_mat, 
        xlab = "Age Group", ylab = "Age Group",
        lwd = 2, labcex = 1,
        main = "Contour Plot of Bivariate Variance-Covariance Surface")

# --- Phase-plane Plot (optional but useful) ---

# Generate the Phase-Plane plot (Velocity vs. Acceleration)
phaseplanePlot(Age, meancancer, 
               labels = list(evalarg = Age, labels = c("45", "50", "55", "60", "65", "70", "75", "80", "85+")),
               xlab = "Breast Cancer Rate Velocity",
               ylab = "Breast Cancer Rate Acceleration",
               main = "Phase-Plane Plot of Breast Cancer Rate")
# Preset the number of principal components (harmonics)
nharm = 3

# Compute Functional PCA
cancer.pcalist = pca.fd(cancer.fd, nharm)

# Extract eigenvalues
cancereig = cancer.pcalist$values
neig = 6

# Scree plot: eigenvalues
x = matrix(1, neig - nharm, 2)
x[,2] = (nharm + 1):neig
y = as.matrix(cancereig[(nharm + 1):neig])
c = lsfit(x, y, int=FALSE)$coef

# Plot Scree Plot
op <- par(mfrow=c(1,1), cex=1.2)
plot(1:neig, cancereig[1:neig], type="b", 
     xlab="Eigenvalue Number", ylab="Eigenvalue",
     main="Scree Plot for FPCA of Breast Cancer Rates")
lines(1:neig, c[1] + c[2]*(1:neig), lty=2)

# Redo FPCA with 2 principal components
nharm = 3
cancer.pcalist = pca.fd(cancer.fd, nharm)

# Print eigenvalues and variance proportions
print(cancer.pcalist$values[1:3])
print(cancer.pcalist$varprop)

# Plot the two principal components
plot.pca.fd(cancer.pcalist, expand=20)

# Hit <Enter> to see the next plot
par(ask=FALSE)

# Apply VARIMAX rotation
cancer.rotpcalist = varmx.pca.fd(cancer.pcalist)

# Plot the rotated principal components
plot.pca.fd(cancer.rotpcalist, expand=20)
# Hit <Enter /Return > to see the next plot :

par( ask = FALSE )


# Print rotated eigenvalues and variance proportions
print(cancer.rotpcalist$values[1:3])
print(cancer.rotpcalist$varprop)
