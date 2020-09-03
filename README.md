---
title: "CopulaR"
author: "Toke Zinn"
date: "1/9/2020"
output:
  github_document:
    pandoc_args: --webtex
---



# Introduction

```CopulaR``` is a lightweight vectorised, but non-optimised, copula framework in R. Copulas are multivariate distribution functions which have uniform marginals. Our interest in copulas primarily comes from Sklar's theorem stating that any multivariate distribution $F$ can be decomposed into its marginals $(F_1, \dots F_n)$ and a copula $C$.

## Regular Use
```CopulaR``` automatically approximates the density given a copula function. Conversely, given a copula density it approximates the copula. We do this through initialisation, namely the ```S4``` class ```Copula``` takes either the density or the copula as an input, i.e. 


```r
# Initialise with density
Independence = Copula(density = function(u,v){1})

# Initialise with copula
Independence_Alternative = Copula(copula = function(u,v){u*v})
```


The initialisation is rather slow. However, once initialised the framework is rather fast at sampling and evaluating. In order to evaluate the density, copula, or sample we use the functions ```d```, ```p```, and ```r``` respectively. 



```r
#Evaluate density 
d(Independence)(0.5,0.5)
```

```
## [1] 1
```

```r
#Evaluate copula 
p(Independence)(0.5,0.5)
```

```
## [1] 0.25
```

```r
#Sample from the copula (here 1000 points) 
plot(r(Independence)(1000), xlab = "U", ylab = "V")
```

<img src="Vignette_files/figure-gfm/IndependenceFigures-1.png" style="display: block; margin: auto;" />


Sampling is *rather* fast, although it can be optimised. 


```r
microbenchmark::microbenchmark(r(Independence)(10000), times = 10)
```

```
## Unit: milliseconds
##                    expr      min       lq     mean   median       uq      max
##  r(Independence)(10000) 442.3189 456.5712 472.2826 468.5299 494.6909 513.8638
##  neval
##     10
```


Furthermore, a standard plotting method is available. We can specify whether we want to plot the copula $C(u,v)$, the density $c(u,v)$, the partial derivatives $\partial C/\partial u$,$\partial C/\partial v$ using the ```surface``` argument. Standard surface is $C(u,v)$.


```r
#plot the copula
plot(Independence)
```

<img src="Vignette_files/figure-gfm/PlottingTypes-1.png" style="display: block; margin: auto;" />

```r
#plot the density
plot(Independence, surface = dC)
```

<img src="Vignette_files/figure-gfm/PlottingTypes-2.png" style="display: block; margin: auto;" />

```r
#plot dC/dU
plot(Independence, surface = dCdU)
```

<img src="Vignette_files/figure-gfm/PlottingTypes-3.png" style="display: block; margin: auto;" />

```r
#plot dC/dV
plot(Independence, surface = dCdV)
```

<img src="Vignette_files/figure-gfm/PlottingTypes-4.png" style="display: block; margin: auto;" />


# Methods
We present several methods implemented for the copula object. Throughout this section let $\mathscr{C}$ denote the set of bi-variate copulas. Recall that a copula $C$ is Lipschitz, and thus by Rademacher's Theorem it is differentiable *almost everywhere*, i.e. the set of points where $C$ is not differentiable has Lebesgue measure $0$. 

## The asterisk product

Let $(U,W)\sim A$ and $(W,V) \sim B$ with $A,B\in \mathscr{C}$. We define 

$$C = A \ast B =  \int_0^1 \frac{dA}{dw}(u,w)\cdot\frac{dB}{dw}(w,v)dw.$$
Then $C$ is a copula, and the probabilistic interpretation is that $C$ is the copula of $(U,V)$ given that $(U,V)$ are conditionally independent given $W$. First we initialise two copulas; a Clayton copula with $\theta = 3$ and a Gaussian copula with $\rho = 0.9$.



```r
A = Copula(copula = function(u,v){
  theta  = 3

  A = u^(-theta) + v^(-theta) - 1

  (A*(A>0))^(-1/theta)

})


B = Copula(copula = function(u,v){
  apply(X = cbind(qnorm(u),qnorm(v)),MARGIN = 1,FUN =  pmvnorm, lower = -Inf, sigma = cbind(c(1,0.9),c(0.9,1)))
})
```

We then apply the $\ast$ operation using ```%ast%```


```r
C = A %ast% B

plot(C)
```

<img src="Vignette_files/figure-gfm/AsteriskOperation-1.png" style="display: block; margin: auto;" />

Sampling from the new copula is readily available 


```r
plot(r(C)(1000))
```

<img src="Vignette_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

## Convex Combinations

The space $\mathscr{C}$ is closed under convex combinations. That is, let $\alpha \in (0,1)$, then $\alpha A + (1-\alpha)B\in\mathscr{C}$ where addition for copulas is defined pointwise in $(u,v)\in [0,1]^2$. 



```r
D = Convex(A,B,0.5)

plot(D, surface = C)
```

<img src="Vignette_files/figure-gfm/ConvexCombination-1.png" style="display: block; margin: auto;" />

# Approximation Methods and Internals 

We briefly discuss how the internals of the package works. The following is rather loose treatment of the subject. 

## Approximation of Copula, Density, and Partial Derivatives 

Suppose we initialse the copula ```C``` using the ```copula``` argument. The function is then evaluated on a $n\times m$ mesh with $Z_{i,j} = C(u_i, v_j)$ where $u_i$ belongs to some partition of $[0,1]$ and $v_j$ to some, possibly different, partition of $[0,1]$. 

Since copulas are differentiable almost everywhere, we apply (bi)linear interpolation between the points on the grid, with the expectation that the approximation is *reasonable*. This may be subject to change to better approximations methods. Similarly, the density and partial derivatives are linearly interpolated. 

In term of calculating partial derivates we have implemented the function ```Partial_Derivative``` in accordance with the ```dplyr``` and ```tidyverse``` framework; the grid $Z$ is expanded into "long" format, i.e. a dataframe with each row being a tuple $(u_i, v_j, z_{i,j})$ for each $i \in \{1,2,\dots, n\}$ and $j\in \{1,2,\dots, m\}$. In order to take a partial derivative we, say $dZ/dU$ we group by $V$ and apply the standard numerical derivatives; central differencing in the interior and forwards or backwards differencing on the boundary. The code is provided below.


```r
Partial_Derivative = function(.Data, Z, X, ...,  Name){
  if(missing(Name)){Name = paste0("d",as_string(ensym(Z)),"d",as_string(ensym(X))) }else{ Name = as_string(ensym(Name))}
  
  # Identify the Boundary
  
  boundary = .Data %>% pull({{X}}) %>% {c(min(.),max(.))}
  
  # Differencing
  .Data %>%
    group_by_at(vars(...,-{{Z}},-{{X}})) %>%
    arrange({{X}}) %>%
    mutate(!!Name := case_when(
      {{X}} == boundary[1] ~ (lead({{Z}}) - {{Z}})/(lead({{X}})-{{X}}),
      {{X}} == boundary[2] ~ ({{Z}} - lag({{Z}}))/({{X}} - lag({{X}})),
      TRUE ~ (lead({{Z}})- lag({{Z}}))/(lead({{X}})-lag({{X}}))
    )) %>% ungroup()

}
```











