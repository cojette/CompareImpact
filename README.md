# CompareImpact
R Package for Measuring and Comparing Time-series Data with Bayesian Prediction Model


## CompareImpact

### An R package for Measuring time-series effect with bayesian prediction modeling

#### What does the package do?
This R package implements an approach to integrate multivariate testing and a Bayesian time series prediction model to measure casual inference based on [CausalImpact package](https://github.com/google/CausalImpact) .  The benefit offering processes have suffered from a lack of effective methods due to insufficient customer data, a wide variety of benefits, and most of all, requiring predictions of customer reactions.

#### How does it work?
Given control time series (e.g. user visits in each group) and a response time series (e.g., total revenue), the package constructs Bayesian structural time-series models for each control time series. And using given event time as intevention,  these models show predictions and differences between time series. i.e. what group is the activist for a promotion. 


#### How is the package structured?
Main Function: CompareImpact()
Sub-Functions:  CreateCompImpPlot()

## 1. Installing the package

To install `CompareImpact`, type the following commands into an R session:

```{r, eval=FALSE}
## install from github
devtools::install_github("google/Causalmpact")

devtools::install_github("cojette/CompareImpact")
```
`CausalImpact` package install is necessary from Github. 

Once installed, the package can be loaded in a given R session using:

```{r,eval=FALSE}
library(CompareImpact)
```

## 3. Running an analysis
1) Multiple user segments
To estimate a causal effect of each user segment, we begin by specifying which period in the data should be used for training the model (*pre-intervention period*) and which period for computing a counterfactual prediction (*post-intervention period*).

```{r,eval=FALSE}
data(promo_data)
pre.period <- c(1, 57)
post.period <- c(58, 90)
```

This says that time points 1 ... 57 will be used for training, and time points 58 ... 100 will be used for computing predictions. There is a big promotion at time points 57

To perform inference, we run the analysis using:
```{r,eval=FALSE}

ci1 <- CompareImpact(promo_data[,-1], c(1,57), c(58,90))
## check impact data
ci1
```

The return value is a list of `CausalImpact` objects. You can see the detail of each CausalImpact objects with CausalImpact package functions and you may use summary, plot functions in CausalImpact package individualy.

## 4. Plotting the results

The easiest way of visualizing the results is to use the `CreateCompImpPlot()` function that is part of the package:
```{r, fig.width=8, fig.height=6, eval=FALSE}

CreateCompImpPlot(ci1)

```

You can control details with ggplot2 package functions. You can refer plot() of CausalImpact package for the detailed interpretation of the plot. 

### Note
CompareImpact package will be updated. And I would appreciate to discuss any aspects of this work (cojette@gmail.com).


## Example Site
[CompareImpact Sample Shiny Site](https://cojette.shinyapps.io/CompareImpactDash)

## Further resources

#### CausalImpact
[CausalImpact package](https://github.com/google/CausalImpact)

#### Measuring the benefit effect for customers with Bayesian predictive modeling
[Strata+Hadoop World 2015 London presentation](http://strataconf.com/big-data-conference-uk-2015/public/schedule/detail/39592)
