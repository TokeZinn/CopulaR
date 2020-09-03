# Functions #
library(tictoc)
library(reshape2)
library(rlang)
library(progress)



n = 10
x = y = seq(from = 0, to = 1, length.out = n)
Z = matrix(0,n,n)

for(i in 1:n){
  for(j in 1:n){
    Z[i,j] = x[i]*y[j]
  }
}


# Bilinear Interpolator #

rownames(Z) = colnames(Z) = x
Z = Z %>% melt() %>% rename(x = Var1, y = Var2, z = value)

Z %>% ggplot(aes(x = x, y = y, fill = z)) + geom_raster()

k = 100
x_new = seq(from = 0, to = 1, length.out = k)
y_new = seq(from = 0, to = 1, length.out = k)


Make_Matrix = function(.Data, x, y, digits = 4){
  #browser()
  X = .Data %>% pull(x) %>% unique() %>% sort()
  Y = .Data %>% pull(y) %>% unique() %>% sort()

  M = dcast(.Data, paste(x,y,sep = "~")) %>%  .[,-1] %>%  as.matrix()

  rownames(M) = sprintf(paste0("%.",digits,"f"),X)
  colnames(M) = sprintf(paste0("%.",digits,"f"),Y)

  return(M)
}

Interpolater = function(x0,y0, x, y, M, digits = 4, pb = NA){
  #browser()

  if(!is.na(pb)){
    pb$tick()
  }

  if ( x0 == 1 & y0 != 1){

    I_ym = max(which(y <= y0))
    Iy = c(I_ym, I_ym + 1)

    y_grid = y[Iy]

    D = diff(y_grid)
    N = y0 - y_grid[1]

    drop((N/D)*diff(M[nrow(M),sprintf(paste0("%.",digits,"f"),y_grid)]) +M[nrow(M),sprintf(paste0("%.",digits,"f"),y_grid[1])])


  } else if  (x0 != 1 & y0 == 1){

    I_xm = max(which(x <= x0))
    Ix = c(I_xm, I_xm + 1)

    x_grid = x[Ix]

    D = diff(x_grid)
    N = x0 - x_grid[1]

    drop((N/D)*diff(M[sprintf(paste0("%.",digits,"f"),x_grid),ncol(M)]) +M[sprintf(paste0("%.",digits,"f"),x_grid[1]),ncol(M)])

  } else if (x0 != 1 & y0 != 1) {

    I_xm = max(which(x <= x0))
    Ix = c(I_xm,I_xm+1)

    I_ym = max(which(y <= y0))
    Iy = c(I_ym,I_ym+1)

    x_grid = x[Ix]
    y_grid = y[Iy]

    x_bar = c(x_grid[2]-x0, x0-x_grid[1])
    y_bar = c(y_grid[2]-y0, y0-y_grid[1])

    Q = M[sprintf(paste0("%.",digits,"f"),x_grid),sprintf(paste0("%.",digits,"f"),y_grid)]
    D = diff(x_grid)*diff(y_grid)

    drop((1/D)*t(x_bar)%*% Q %*% y_bar)

  } else{

    M[nrow(M),ncol(M)]
  }

}

linear_interpolator = function(.Data, U,V, u, v, digits = 4, progress = T){
  #browser()

  if(progress){
    pb = progress_bar$new(total = length(u),
                          format = "Compiling Function [:bar] :percent eta: :eta",
                          clear = T,
                          width = 60,)
  }else{
    pb = NA
  }

  M = Make_Matrix(.Data,x = U,y = V,digits = digits)

  U = .Data %>% pull(U) %>% unique() %>% sort()
  V = .Data %>% pull(V) %>% unique() %>% sort()

  Z = map2(u,v, Interpolater, x = U, y = V, M = M, digits = digits, pb = pb)

  res = tibble(u,v, C = unlist(Z))

  return(res)
}


grid = expand.grid(x_new,y_new)

test = linear_interpolator(Z, "x","y", grid$Var1, grid$Var2,progress = F)
test %>%
  ggplot(aes(x = u, y = v, fill = C)) + geom_raster()
