CopulaR (WiP)
================
# Introduction

`CopulaR` is a lightweight vectorised, but non-optimised, copula
framework in R. Copulas are multivariate distribution functions which
have uniform marginals.

## Regular Use

`CopulaR` automatically approximates the density given a copula
function. Conversely, given a copula density it approximates the copula.
We do this through initialisation, namely the `S4` class `Copula` takes
either the density or the copula as an input, i.e. 

``` r
# Initialise with density
Independence = Copula(density = function(u,v){1})

# Initialise with copula
Independence_Alternative = Copula(copula = function(u,v){u*v})
```

The initialisation is rather slow. However, once initialised the
framework is rather fast at sampling and evaluating. In order to
evaluate the density, copula, or sample we use the functions `d`, `p`,
and `r` respectively.

``` r
#Evaluate density 
d(Independence)(0.5,0.5)
```

    ## [1] 1

``` r
#Evaluate copula 
p(Independence)(0.5,0.5)
```

    ## [1] 0.25

``` r
#Sample from the copula (here 1000 points) 
plot(r(Independence)(1000), xlab = "U", ylab = "V")
```

<img src="README_files/figure-gfm/IndependenceFigures-1.png" style="display: block; margin: auto;" />

Sampling is *rather* fast, although it can be optimised.

``` r
microbenchmark::microbenchmark(r(Independence)(10000), times = 10)
```

    ## Unit: milliseconds
    ##                    expr      min       lq     mean  median       uq     max
    ##  r(Independence)(10000) 437.8137 461.9183 477.9726 476.672 489.8568 520.069
    ##  neval
    ##     10

Furthermore, a standard plotting method is available. We can specify
whether we want to plot the copula
![C(u,v)](https://latex.codecogs.com/png.latex?C%28u%2Cv%29 "C(u,v)"),
the density ![c(u,v)](https://latex.codecogs.com/png.latex?c%28u%2Cv%29
"c(u,v)"), the partial derivatives ![\\partial C/\\partial
u](https://latex.codecogs.com/png.latex?%5Cpartial%20C%2F%5Cpartial%20u
"\\partial C/\\partial u"),![\\partial C/\\partial
v](https://latex.codecogs.com/png.latex?%5Cpartial%20C%2F%5Cpartial%20v
"\\partial C/\\partial v") using the `surface` argument. Standard
surface is ![C(u,v)](https://latex.codecogs.com/png.latex?C%28u%2Cv%29
"C(u,v)").

``` r
#plot the copula
plot(Independence)
```

<img src="README_files/figure-gfm/PlottingTypes-1.png" style="display: block; margin: auto;" />

``` r
#plot the density
plot(Independence, surface = dC)
```

<img src="README_files/figure-gfm/PlottingTypes-2.png" style="display: block; margin: auto;" />

``` r
#plot dC/dU
plot(Independence, surface = dCdU)
```

<img src="README_files/figure-gfm/PlottingTypes-3.png" style="display: block; margin: auto;" />

``` r
#plot dC/dV
plot(Independence, surface = dCdV)
```

<img src="README_files/figure-gfm/PlottingTypes-4.png" style="display: block; margin: auto;" />

# Methods

We present several methods implemented for the copula object. Recall
that a copula ![C](https://latex.codecogs.com/png.latex?C "C") is
Lipschitz, and thus by Rademacher’s Theorem it is differentiable *almost
everywhere*, i.e. the set of points where
![C](https://latex.codecogs.com/png.latex?C "C") is not differentiable
has Lebesgue measure ![0](https://latex.codecogs.com/png.latex?0 "0").

## The asterisk product

Let ![(U,W)\\sim
A](https://latex.codecogs.com/png.latex?%28U%2CW%29%5Csim%20A
"(U,W)\\sim A") and ![(W,V) \\sim
B](https://latex.codecogs.com/png.latex?%28W%2CV%29%20%5Csim%20B
"(W,V) \\sim B") with ![A,B](https://latex.codecogs.com/png.latex?A%2CB
"A,B") being bi-variate copulas. We define

  
![C = A \\ast B = \\int\_0^1
\\frac{dA}{dw}(u,w)\\cdot\\frac{dB}{dw}(w,v)dw.](https://latex.codecogs.com/png.latex?C%20%3D%20A%20%5Cast%20B%20%3D%20%20%5Cint_0%5E1%20%5Cfrac%7BdA%7D%7Bdw%7D%28u%2Cw%29%5Ccdot%5Cfrac%7BdB%7D%7Bdw%7D%28w%2Cv%29dw.
"C = A \\ast B =  \\int_0^1 \\frac{dA}{dw}(u,w)\\cdot\\frac{dB}{dw}(w,v)dw.")  
Then ![C](https://latex.codecogs.com/png.latex?C "C") is a copula, and
the probabilistic interpretation is that
![C](https://latex.codecogs.com/png.latex?C "C") is the copula of
![(U,V)](https://latex.codecogs.com/png.latex?%28U%2CV%29 "(U,V)") given
that ![(U,V)](https://latex.codecogs.com/png.latex?%28U%2CV%29 "(U,V)")
are conditionally independent given
![W](https://latex.codecogs.com/png.latex?W "W"). First we initialise
two copulas; a Clayton copula with ![\\theta
= 3](https://latex.codecogs.com/png.latex?%5Ctheta%20%3D%203
"\\theta = 3") and a Gaussian copula with ![\\rho
= 0.9](https://latex.codecogs.com/png.latex?%5Crho%20%3D%200.9
"\\rho = 0.9").

``` r
A = Copula(copula = function(u,v){
  theta  = 3

  A = u^(-theta) + v^(-theta) - 1

  (A*(A>0))^(-1/theta)

})


B = Copula(copula = function(u,v){
  apply(X = cbind(qnorm(u),qnorm(v)),MARGIN = 1,FUN =  pmvnorm, lower = -Inf, sigma = cbind(c(1,0.9),c(0.9,1)))
})
```

We then apply the ![\\ast](https://latex.codecogs.com/png.latex?%5Cast
"\\ast") operation using `%ast%`

``` r
C = A %ast% B

plot(C)
```

<img src="README_files/figure-gfm/AsteriskOperation-1.png" style="display: block; margin: auto;" />

Sampling from the new copula is readily
available

``` r
plot(r(C)(1000))
```

<img src="README_files/figure-gfm/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

## Convex Combinations

The space of copulas is closed under convex combinations. That is, let
![\\alpha \\in
(0,1)](https://latex.codecogs.com/png.latex?%5Calpha%20%5Cin%20%280%2C1%29
"\\alpha \\in (0,1)"), then ![\\alpha A +
(1-\\alpha)B](https://latex.codecogs.com/png.latex?%5Calpha%20A%20%2B%20%281-%5Calpha%29B
"\\alpha A + (1-\\alpha)B") is a copula where addition for copulas is
defined pointwise in ![(u,v)\\in
\[0,1\]^n](https://latex.codecogs.com/png.latex?%28u%2Cv%29%5Cin%20%5B0%2C1%5D%5En
"(u,v)\\in [0,1]^n").

``` r
D = Convex(A,B,0.5)

plot(D, surface = C)
```

<img src="README_files/figure-gfm/ConvexCombination-1.png" style="display: block; margin: auto;" />

# Approximation Methods and Internals

We briefly discuss how the internals of the package works. The following
is rather loose treatment of the subject.

## Approximation of Copula, Density, and Partial Derivatives

Suppose we initialse the copula `C` using the `copula` argument. The
function is then evaluated on a ![n\\times
m](https://latex.codecogs.com/png.latex?n%5Ctimes%20m "n\\times m") mesh
with ![Z\_{i,j} = C(u\_i,
v\_j)](https://latex.codecogs.com/png.latex?Z_%7Bi%2Cj%7D%20%3D%20C%28u_i%2C%20v_j%29
"Z_{i,j} = C(u_i, v_j)") where
![u\_i](https://latex.codecogs.com/png.latex?u_i "u_i") belongs to some
partition of ![\[0,1\]](https://latex.codecogs.com/png.latex?%5B0%2C1%5D
"[0,1]") and ![v\_j](https://latex.codecogs.com/png.latex?v_j "v_j") to
some, possibly different, partition of
![\[0,1\]](https://latex.codecogs.com/png.latex?%5B0%2C1%5D "[0,1]").

Since copulas are differentiable almost everywhere, we apply (bi)linear
interpolation between the points on the grid, with the expectation that
the approximation is *reasonable*. This may be subject to change to
better approximations methods. Similarly, the density and partial
derivatives are linearly interpolated.

In term of calculating partial derivates we have implemented the
function `Partial_Derivative` in accordance with the `dplyr` and
`tidyverse` framework; the grid
![Z](https://latex.codecogs.com/png.latex?Z "Z") is expanded into “long”
format, i.e. a dataframe with each row being a tuple ![(u\_i, v\_j,
z\_{i,j})](https://latex.codecogs.com/png.latex?%28u_i%2C%20v_j%2C%20z_%7Bi%2Cj%7D%29
"(u_i, v_j, z_{i,j})") for each ![i \\in \\{1,2,\\dots,
n\\}](https://latex.codecogs.com/png.latex?i%20%5Cin%20%5C%7B1%2C2%2C%5Cdots%2C%20n%5C%7D
"i \\in \\{1,2,\\dots, n\\}") and ![j\\in \\{1,2,\\dots,
m\\}](https://latex.codecogs.com/png.latex?j%5Cin%20%5C%7B1%2C2%2C%5Cdots%2C%20m%5C%7D
"j\\in \\{1,2,\\dots, m\\}"). In order to take a partial derivative we,
say ![dZ/dU](https://latex.codecogs.com/png.latex?dZ%2FdU "dZ/dU") we
group by ![V](https://latex.codecogs.com/png.latex?V "V") and apply the
standard numerical derivatives; central differencing in the interior and
forwards or backwards differencing on the boundary. The code is provided
below.

``` r
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
