# Dimension Reduction and MARS
 Improve the performance of MARS by using linear combinations of the covariates which achieve sufficient dimension reduction. 

Paper: [[2302.05790\] Dimension Reduction and MARS (arxiv.org)](https://arxiv.org/abs/2302.05790)

### Usage ###

In this paper, we show the performance of drMARS in estimating the SDR space, the estimating dimensionality of the SDR space, and the prediction performance. The corresponding function implementations are "drMARS", "drMARS.CV" and "drMARS.fit", respectively. they can be called in the following manner:

* B = drMARS(x, y, degree = NULL, Xscale=F, plus=F)$B

  The estimation accuracy of the SDR space is evaluated with the true dimensions of the SDR space.

* d=drMARS.CV(B, max.dim=5, nfold=10)$ndir
  We select the dimension of SDR space (d) by 10-fold cross-validation.

* predictions = drMARS.fit(x,y,xnew,degree = NULL,Xadd=T,Xnorm=F,Xscale=F,plus=F,iter=F,ndir="NoPreSel",max.dim = 5,max.iter=50)$predicted

  We provide a number of arguments to make predictions using drMARS, and users need to adjust the arguments to improve the prediction accuracy of drMARS based on the data. The usage of the parameters is described in the file "drMARS.R".

In addition, we also provide a three-dimensional graph of the predicted performance of drMARS. The following figure shows that the left side is the fitted graph and the right side is the real graph, the similarity between them is very high, which indicates that the drMARS fitting performance is good. See "example_scripts" file for specific usage and performance in the example. 

![three-dimensional graphs](plot_drMARS.jpg)

### Datasets ###

The 7 real data used in the paper are in the "DataSets" folder and also contain the data loading and pre-processing script "DataSets.R". All datasets are from the UCI repository [UCI Machine Learning Repository](https://archive.ics.uci.edu/ml/index.php) and the kaggle  [Kaggle: Your Machine Learning and Data Science Community](https://www.kaggle.com/), and they can available freely online. The detailed use of the data is described in the paper.

### Contribution guidelines ###

Any improvements or conversions to other code formats would be appreciated, send me an email if you would like to contribute / require assistance. 
Yu Liu: liuyuchina123@gmail.com.
