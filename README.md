# Power Study to Compare Mean Vectors from Two-Populations


In order to compare mean vectors of two independent populations, we make the following assumptions:

(1) $X_{1,1}, ..., X_{1,n_{1}}$ is a random sample of size $n_{1}$ from a p-variate population with mean vector $\mu_{1}$ and covariance matrix $\Sigma_{1}$


(2) $X_{2,1}, ..., X_{2,n_{2}}$ is a random sample of size $n_{2}$ from a p-variate population with mean vector $\mu_{2}$ and covariance matrix $\Sigma_{2}$

(3) $\{ X_{1,1}, ..., X_{1,n_{1}} \}$ is independent of the other  $\{ X_{2,1}, ..., X_{2,n_{2}} \}$


<p>&nbsp;</p>


We want to test for the equality of the mean vectors from two populations:


$H_{0}: \mu_{1} = \mu_{2}$ vs. $H_{a}: \mu_{1} \neq \mu_{2}$


<p>&nbsp;</p>


If we add two more assumptions:

(4) $X_{1,1}, ..., X_{1,n_{1}}$ $\sim$  MN($\mu_{1}$, $\Sigma_{1}$) and $X_{2,1}, ..., X_{2,n_{2}}$ $\sim$ MN($\mu_{2}$, $\Sigma_{2}$)

(5) $\Sigma_{1}$ = $\Sigma_{2}$


Then, the likelihood ratio test of $H_{0}$ is based on the square of the statistical distance, $T^{2}$, and is given by:

$T^{2} = (\bar{X_{1}} - \bar{X_{2}} - (\mu_{1} - \mu_{2}))'(\frac{1}{n_{1}} + \frac{1}{n_{2}}) S_{pooled}^{-1} (\bar{X_{1}} - \bar{X_{2}} - (\mu_{1} - \mu_{2}))$

where:

* $T^{2}$ is distributed as $\frac{(n_{1} + n_{2} - 2)p}{n_{1} + n_{2} - p - 1} F_{p, n_{1} + n_{2} - p - 1}(\alpha)$

* $S_{pooled} = \frac{\sum_{j=1}^{n_{1}} (X_{1j} - \bar{X_1})(X_{1j} - \bar{X_1})' + \sum_{j=1}^{n_{1}} (X_{2j} - \bar{X_2})(X_{2j} - \bar{X_2})'}{n_{1} + n_{2} - 2} = \frac{n_{1} -1}{n_{1} + n_{2} - 2}S_{1} + \frac{n_{2} -1}{n_{1} + n_{2} - 2}S_{2}$

* the critical distance $c^{2}$ is determined from the distribution of the two-sample $T^{2}$ statistic so that $P(T^{2} < c^{2} = \frac{(n_{1} + n_{2} - 2)p}{n_{1} + n_{2} - p - 1} F_{p, n_{1} + n_{2} - p - 1}(\alpha) ) = 1- \alpha$.


<p>&nbsp;</p>


Therefore, we reject $H_{0}$ if: $T^{2} > c^{2}$


We conclude that all $\mu_{1} - \mu_{2}$ within squared statistical distance $c^{2}$ of $\bar{x_{1}} - \bar{x_{2}}$ constitute the confidence region. 
 
 
<p>&nbsp;</p>


When the assumptions (4) and/or (5) is/are not met, for large samples, this structure is sufficient for making inferences about the p × 1 vector $H_{0}$: $\mu_{1} = \mu_{2}$. 

The approximate test that allow us to test the hypotheses $H_{0}: \mu_{1} = \mu_{2}$ vs $H_{a}: \mu_{1} \neq \mu_{2}$ is:

$T^{2} = (\bar{X_{1}} - \bar{X_{2}} - (\mu_{1} - \mu_{2}))'(\frac{1}{n_{1}}S_{1} + \frac{1}{n_{2}}S_{2})  (\bar{X_{1}} - \bar{X_{2}} - (\mu_{1} - \mu_{2}) \sim \chi^{2}_{p}$


We reject $H_{0}$ at level $\alpha$ if $T^{2} > \chi^{2}_{p, \alpha}$.

<p>&nbsp;</p>


We construct a power study that estimates the power and type-I error of the approximate test when assumptions 4 and 5 are not met.
