library(cubature)
library(purrr)
library(tidyverse)
library(rlang)
library(reshape2)
library(mvtnorm)

#### Functions ###
Partial_Derivative = function(.Data, Z, X, ...,  Name){
  if(missing(Name)){Name = paste0("d",as_string(ensym(Z)),"d",as_string(ensym(X))) }else{ Name = as_string(ensym(Name))}

  boundary = .Data %>% pull({{X}}) %>% {c(min(.),max(.))}

  .Data %>%
    group_by_at(vars(...,-{{Z}},-{{X}})) %>%
    arrange({{X}}) %>%
    mutate(!!Name := case_when(
      {{X}} == boundary[1] ~ (lead({{Z}}) - {{Z}})/(lead({{X}})-{{X}}),
      {{X}} == boundary[2] ~ ({{Z}} - lag({{Z}}))/({{X}} - lag({{X}})),
      TRUE ~ (lead({{Z}})- lag({{Z}}))/(lead({{X}})-lag({{X}}))
    )) %>% ungroup()

}
Make_Grid = function(n){
  u = seq(from = 0, to = 1, length.out = n)
  expand.grid(u,u) %>% rename(U = Var1, V = Var2)
}

#### Helpers ###
Debug_Init = function(.Object){
  #browser()
  # Test if slots are filled
  empty_slots = c(is.na(try(.Object@density(),silent = T)),
                  is.na(try(.Object@copula(),silent = T))
  )

  empty_partials = c(is.na(try(.Object@cond_U(),silent = T)),
                     is.na(try(.Object@cond_V(),silent = T))
  )


  if(!all(empty_slots)){
    which_missing = which(empty_slots)


    U = V = seq(from = 0, to = 1, length.out = 100)
    Grid = expand.grid(U,V)
    names(Grid) = c("U","V")


    if(which_missing == 1){

      Grid = Grid %>% mutate(C = map2(U,V, .Object@copula)) %>%
        mutate(C = unlist(C)) %>%
        group_by(V) %>% arrange(U) %>%
        mutate(dC1 = case_when(
          U == 0 ~ (lead(C) - C)/(lead(U)-U),
          U == 1 ~ (C - lag(C))/(U - lag(U)),
          TRUE ~ (lead(C)- lag(C))/(lead(U)-lag(U))
        )) %>% ungroup() %>% group_by(U) %>% arrange(V) %>%
        mutate(d2C2 = case_when(
          V == 0 ~ (lead(dC1) - dC1)/(lead(V)-V),
          V == 1 ~ (dC1 - lag(dC1))/(V - lag(V)),
          TRUE ~ (lead(dC1)- lag(dC1))/(lead(V)-lag(V))
        )) %>% select(U,V,d2C2)



      density = function(u,v){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = Grid, "U","V",u,v,digits = 4))

      }

      .Object@density <- density



    }else if (which_missing == 2){
      Density_Wrapper = function(x){
        .Object@density(x[1],x[2])
      }
      integral_wrapper = function(u,v){
        hcubature(Density_Wrapper, lowerLimit = c(0,0),upperLimit = c(u,v))$integral
      }

      Grid = Grid %>% mutate(C = unlist(map2(U,V,integral_wrapper)))

      copula = function(u,v){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = Grid, "U","V",u,v,digits = 4))

      }

      .Object@copula <- copula


    }


    if(any(empty_partials)){
      temp = Grid %>% mutate(C = .Object@copula(U,V)) %>% select(U,V,C)

      if(empty_partials[1]){

        temp1 = temp %>% Partial_Derivative(C,U,Name = dC) %>% select(U,V,dC)

        cond_U = function(u,v){

          suppressMessages(Linear_Interpolator_Wrapper(.Data = temp1, "U","V",u,v,digits = 4))

        }

        .Object@cond_U <- cond_U

      }

      if(empty_partials[2]){
        temp2 = temp %>% Partial_Derivative(C,V,Name = dC) %>% select(U,V,dC)

        cond_V = function(u,v){

          suppressMessages(Linear_Interpolator_Wrapper(.Data = temp2, "U","V",u,v,digits = 4))

        }

        .Object@cond_V <- cond_V

      }



    }

    if(is.na(try(.Object@sampler()))){

      temp3 = Grid %>% select(U,V) %>%
        mutate(C = .Object@cond_U(U,V)) %>%
        rename(W = V) %>%
        group_by(U) %>% nest() %>%
        mutate(V = map2(U, data, function(U,data){
          f = approxfun(data$W, data$C,method = "linear",yleft = 0,yright = 1)

          if(all(sign(f(0)-data$W) == sign(f(1)-data$W))){
            cat(U)
          }

          sapply(data$W, function(w){if(!(abs(w) < 1e-6 | abs(1-w) < 1e-6)) uniroot(function(x){f(x) - w},interval = c(0,1))$root else w})

        })) %>% unnest(cols = c(data,V)) %>% ungroup() %>% select(U,W,V)


      sampler = function(u,w){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = temp3, "U","W",u,w,digits = 4))

      }

      .Object@sampler <- sampler

    }

  }else{
    if(all(empty_slots)){
      stop("Error! No density or copula supplied")
    }
  }

  .Object
}

#### Classes ###
Copula <- setClass("Copula",
         representation(
           density = "function",
           copula = "function",
           cond_U = "function",
           cond_V = "function",
           sampler = "function",
           grid = "ANY"
         ),
         prototype(
           density = function() NA,
           copula = function() NA,
           cond_U = function() NA,
           cond_V = function() NA,
           sampler = function() NA,
           grid = 0
         ))

Debug_Copula <- setClass("Debug_Copula",
                         representation(
                           density = "function",
                           copula = "function",
                           cond_U = "function",
                           cond_V = "function",
                           sampler = "function",
                           grid = "ANY"
                         ),
                         prototype(
                           density = function() NA,
                           copula = function() NA,
                           cond_U = function() NA,
                           cond_V = function() NA,
                           sampler = function() NA,
                           grid = NA
                         ))


# Fill Blank Method
setGeneric("Fill_Blank", function(.Object){standardGeneric("Fill_Blank")})
setMethod("Fill_Blank", signature(.Object = "Copula"), function(.Object){
  #browser()
  # Test if slots are filled
  empty_slots = c(is.na(try(.Object@density(),silent = T)),
                  is.na(try(.Object@copula(),silent = T))
  )

  empty_partials = c(is.na(try(.Object@cond_U(),silent = T)),
                     is.na(try(.Object@cond_V(),silent = T))
  )

  U = V = seq(from = 0, to = 1, length.out = 100)
  Grid = expand.grid(U,V)
  names(Grid) = c("U","V")

  Grid_Obj = Grid



  if(!all(empty_slots)){
    if(any(empty_slots)) which_missing = which(empty_slots) else which_missing = 0


    if(which_missing == 1){

      Grid_d = Grid %>% mutate(C = map2(U,V, .Object@copula)) %>%
        mutate(C = unlist(C)) %>%
        group_by(V) %>% arrange(U) %>%
        mutate(dC1 = case_when(
          U == 0 ~ (lead(C) - C)/(lead(U)-U),
          U == 1 ~ (C - lag(C))/(U - lag(U)),
          TRUE ~ (lead(C)- lag(C))/(lead(U)-lag(U))
        )) %>% ungroup() %>% group_by(U) %>% arrange(V) %>%
        mutate(d2C2 = case_when(
          V == 0 ~ (lead(dC1) - dC1)/(lead(V)-V),
          V == 1 ~ (dC1 - lag(dC1))/(V - lag(V)),
          TRUE ~ (lead(dC1)- lag(dC1))/(lead(V)-lag(V))
        )) %>% select(U,V,C,d2C2) %>% rename(dC = d2C2)

      Grid_Obj = Grid_Obj %>% left_join(Grid_d, by = c("U","V"))

      Grid_d = Grid_d %>% select(-C)

      density = function(u,v){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = Grid_d, "U","V",u,v,digits = 4))

      }

      .Object@density <- density



    }else if (which_missing == 2){
      Density_Wrapper = function(x){
        .Object@density(x[1],x[2])
      }
      integral_wrapper = function(u,v){
        hcubature(Density_Wrapper, lowerLimit = c(0,0),upperLimit = c(u,v))$integral
      }

      Grid_c = Grid %>% mutate(C = unlist(map2(U,V,integral_wrapper))) %>% mutate(dC = unlist(map2(U,V,.Object@density)))

      Grid_Obj = Grid_Obj %>% left_join(Grid_c,by = c("U","V"))

      Grid_c = Grid_c %>% select(U,V,C)

      copula = function(u,v){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = Grid_c, "U","V",u,v,digits = 4))

      }

      .Object@copula <- copula


    }


    if(any(empty_partials)){
      temp = Grid %>% mutate(C = .Object@copula(U,V)) %>% select(U,V,C)

      if(empty_partials[1]){

        temp1 = temp %>% Partial_Derivative(C,U,Name = dCdU) %>% select(U,V,dCdU)

        Grid_Obj = Grid_Obj %>% left_join(temp1, by = c("U","V"))

        cond_U = function(u,v){

          suppressMessages(Linear_Interpolator_Wrapper(.Data = temp1, "U","V",u,v,digits = 4))

        }

        .Object@cond_U <- cond_U

      }

      if(empty_partials[2]){
        temp2 = temp %>% Partial_Derivative(C,V,Name = dCdV) %>% select(U,V,dCdV)

        Grid_Obj = Grid_Obj %>% left_join(temp2, by = c("U","V"))

        cond_V = function(u,v){

          suppressMessages(Linear_Interpolator_Wrapper(.Data = temp2, "U","V",u,v,digits = 4))

        }

        .Object@cond_V <- cond_V

      }



    }

    if(is.na(try(.Object@sampler()))){

      temp3 = Grid %>% select(U,V) %>%
        mutate(C = .Object@cond_U(U,V)) %>%
        rename(W = V) %>%
        group_by(U) %>% nest() %>%
        mutate(V = map2(U, data, function(U,data){
          f = approxfun(data$W, data$C,method = "linear",yleft = 0,yright = 1)

          if(all(sign(f(0)-data$W) == sign(f(1)-data$W))){
            cat(U)
          }

          sapply(data$W, function(w){if(!(abs(w) < 1e-6 | abs(1-w) < 1e-6)) uniroot(function(x){f(x) - w},interval = c(0,1))$root else w})

        })) %>% unnest(cols = c(data,V)) %>% ungroup() %>% select(U,W,V)


      sampler = function(u,w){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = temp3, "U","W",u,w,digits = 4))

      }

      .Object@sampler <- sampler

    }

  }else{
    if(!all(.Object@grid == 0)){
      Grid = .Object@grid

      Grid_Obj = Grid

      Grid_c = Grid

      copula = function(u,v){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = Grid_c, "U","V",u,v,digits = 4))

      }

      .Object@copula <- copula

      Grid_d = Grid %>%
        group_by(V) %>% arrange(U) %>%
        mutate(dC1 = case_when(
          U == 0 ~ (lead(C) - C)/(lead(U)-U),
          U == 1 ~ (C - lag(C))/(U - lag(U)),
          TRUE ~ (lead(C)- lag(C))/(lead(U)-lag(U))
        )) %>% ungroup() %>% group_by(U) %>% arrange(V) %>%
        mutate(d2C2 = case_when(
          V == 0 ~ (lead(dC1) - dC1)/(lead(V)-V),
          V == 1 ~ (dC1 - lag(dC1))/(V - lag(V)),
          TRUE ~ (lead(dC1)- lag(dC1))/(lead(V)-lag(V))
        )) %>% select(U,V,C,d2C2) %>% rename(dC = d2C2) %>% select(-C)

      Grid_Obj = Grid_Obj %>%  left_join(Grid_d, by = c("U","V"))

      Grid_d = Grid_d

      density = function(u,v){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = Grid_d, "U","V",u,v,digits = 4))

      }


      temp = Grid %>% mutate(C = .Object@copula(U,V)) %>% select(U,V,C)

      temp1 = temp %>% Partial_Derivative(C,U,Name = dCdU) %>% select(U,V,dCdU)

      Grid_Obj  = Grid_Obj %>% left_join(temp1, by = c("U","V"))


      cond_U = function(u,v){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = temp1, "U","V",u,v,digits = 4))

      }

      .Object@cond_U <- cond_U

      temp2 = temp %>% Partial_Derivative(C,V,Name = dCdV) %>% select(U,V,dCdV)

      Grid_Obj = Grid_Obj %>% left_join(temp2, by = c("U","V"))

      cond_V = function(u,v){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = temp2, "U","V",u,v,digits = 4))

      }

      .Object@cond_V <- cond_V


      temp3 = Grid_Obj %>% select(U,V,dCdU) %>%
        rename(C = dCdU) %>%
        rename(W = V) %>%
        group_by(U) %>% nest() %>%
        mutate(V = map2(U, data, function(U,data){
          f = approxfun(data$W, data$C,method = "linear",yleft = 0,yright = 1)

          if(all(sign(f(0)-data$W) == sign(f(1)-data$W))){
            cat(U)
          }

          sapply(data$W, function(w){if(!(abs(w) < 1e-6 | abs(1-w) < 1e-6)) uniroot(function(x){f(x) - w},interval = c(0,1))$root else w})

        })) %>% unnest(cols = c(data,V)) %>% ungroup() %>% select(U,W,V)


      sampler = function(u,w){

        suppressMessages(Linear_Interpolator_Wrapper(.Data = temp3, "U","W",u,w,digits = 4))

      }

      .Object@sampler <- sampler



    }else{
      stop("Cannot Initialise: No Copula, Density or Grid Approximation Supplied")
    }
  }

  .Object@grid <- Grid_Obj


  .Object
})


# Initialize Copula
setMethod("initialize", signature(.Object = "Copula"), function(.Object, ...){
  .Object <- callNextMethod(.Object, ...)
  .Object <- Fill_Blank(.Object)
  .Object
})

# Convex Combinations
setGeneric("Convex", function(C_1,C_2,alpha){standardGeneric("Convex")})
setMethod("Convex", signature(C_1 = "Copula", C_2 = "Copula", alpha = "numeric"),function(C_1,C_2,alpha = 0.5){
            #browser()
            local = function(a,b){alpha*a + (1-alpha)*b}
            A = C_1@grid %>% select(U,V,C)
            B = C_2@grid %>% select(U,V,C)
            Grid = left_join(A,B, by = c("U","V"))
            Grid = Grid %>% mutate(C = reduce(select(., starts_with("C")), local)) %>% select(U,V,C)


            Copula(grid = Grid)

          })

# Distribution Functions
setGeneric("d", function(.Object){standardGeneric("d")})
setGeneric("p", function(.Object){standardGeneric("p")})
setGeneric("r", function(.Object){standardGeneric("r")})
setMethod("d", signature(.Object = "Copula"), function(.Object){.Object@density})
setMethod("p", signature(.Object = "Copula"), function(.Object){.Object@copula})
setMethod("r", signature(.Object = "Copula"), function(.Object){local = function(n){
  U = runif(n)
  W = runif(n)

  cbind(U,.Object@sampler(U,W))
}})

#Plot Method
setMethod("plot", signature(x = "Copula", y = "missing"), function(x,surface = C){
  x@grid %>% select(U,V,{{surface}}) %>% gather(key = "Function",value = "Value",-U,-V) %>%
    ggplot(aes(x = U, y = V, fill = Value)) + geom_raster() + facet_wrap(~Function) + theme_bw() +
    scale_fill_viridis_c()
})

# Asterisk Product
setGeneric("%ast%", function(C_1,C_2){standardGeneric("%ast%")})
setMethod("%ast%", signature(C_1 = "Copula", C_2 = "Copula"), function(C_1,C_2){

  dAdV <- C_1@grid %>% select(U,V,dCdV) %>% rename(W = V)
  dBdU <- C_2@grid %>% select(U,V,dCdU) %>% rename(W = U)

  C <- full_join(dAdV, dBdU, by = "W")

  C = C %>% group_by(U,V) %>%
    arrange(W) %>%
    mutate(dC = dCdV*dCdU) %>%
    summarise(C = sum((lead(dC)+dC)*(lead(W)-W)/2,na.rm = T)) %>%
    ungroup() %>%
    mutate(C = (C-min(C))/max(C))

  return(Copula(grid = C))
})




















