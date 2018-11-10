
##For a given matrix A, 
#find lower trigular matrix L and upper triangular U such that A = LU 
LUfactorization <- function(mat){
  
  num_rows = nrow(mat)
  num_cols = ncol(mat)
  
  U = mat
  L = diag(num_rows)
  
  num_iter = min(num_rows, num_cols) - 1
  sapply(1:num_iter, function(curr_col){
    M <- diag(num_rows)
    Minv <- diag(num_rows)
    curr_pivot_col <- U[((curr_col+1):num_rows),
                        curr_col]/U[curr_col, curr_col]
    M[((curr_col+1):num_rows),curr_col] <- curr_pivot_col * (-1)
    Minv[((curr_col+1):num_rows),curr_col] <- curr_pivot_col
    U <<- M %*% U
    L <<- L %*% Minv 
    
  })
  
  return(list("L" = L, "U" = U))
}

##For a given matrix A, 
#find lower trigular matrices L and M, and diagonal matrix D such that A = LDM^T 
LDMfactorization <- function(mat){
  
  num_rows <- nrow(mat)
  num_cols <- ncol(mat)
  
  LUfactors <- LUfactorization(mat)
  U <- LUfactors$U
  
  D <- matrix(data = 0, num_rows, num_rows)
  Mtran <- U
  
  sapply(1:num_rows, function(curr_row){
    D[curr_row, curr_row] <<- Mtran[curr_row, min(curr_row, num_cols)]
    Mtran[curr_row,] <<- Mtran[curr_row,]/D[curr_row, curr_row]
  })
  
  M <- t(Mtran)
  
  return(list("L" = LUfactors$L, "D" = D, "M" = M))
}

##For a given symmetric positive definite matrix A, 
#find lower trigular matrix L such that A = LL^T using LDM^T factorization
CholeskyFactorizationByLDM <- function(mat){
  
  num_rows <- nrow(mat)
  
  LDLTfactors <- LDMfactorization(mat)
  
  L <- LDLTfactors$L
  D <- LDLTfactors$D
  
  sqrtD <- matrix(data = 0, num_rows, num_rows)
  sapply(1:num_rows, function(curr_row){
    sqrtD[curr_row, curr_row] <<- sqrt(D[curr_row, curr_row])
  })
  
  L <- L %*% sqrtD
  
  return(L)
}


##For a given symmetric positive definite matrix A, 
#find lower trigular matrix L such that A = LL^T using outer product method
CholeskyByOuter <- function(mat){
  
  num_rows <- nrow(mat)
  L <- mat
  
  sapply(1:num_rows, function(curr_row){
    alpha <- L[curr_row, curr_row]
    beta <- sqrt(alpha)
    if(curr_row < num_rows){
      B <- L[(curr_row+1):num_rows, (curr_row+1):num_rows]
      v <- L[(curr_row+1):num_rows,curr_row]
      vbeta <- v/beta
      v_outer <- (v %*% t(v)) / alpha
      G <- B - v_outer
      L[curr_row, curr_row] <<- beta
      L[(curr_row+1):num_rows, curr_row] <<- vbeta
      L[curr_row, (curr_row+1):num_rows] <<- 0
      L[(curr_row+1):num_rows, (curr_row+1):num_rows] <<- G
    }else{
      L[curr_row, curr_row] <<- beta
    }
    
  })
  
  return(L)
}

##For a given matrix A, 
#find orthogonal matrix Q and upper triangular matrix R such that A = QR 
#using Householder reflections
QRfactorizationHouseholder <- function(mat) {
  
  num_cols <- ncol(mat)
  num_rows <- nrow(mat)
  
  R <- mat
  
  num_iter <- min(num_cols, num_rows)
  
  H <- lapply(1:num_iter, function(i){
    x <- R[i:num_rows, i]
    
    len_x <- length(x)
    
    x_ext <- matrix(c(rep(0, (i-1)),
                      x), ncol = 1)
    
    e <- matrix(c(rep(0, (i-1)),
                  1,
                  rep(0, (len_x - 1))), ncol = 1)
    u <- sqrt(sum(x^2))
    v <- x_ext + sign(x[1]) * u * e
    
    v_inner <- (t(v) %*% v)[1]
    
    Hi <- diag(num_rows) - 2 * (v %*% t(v)) / v_inner
    
    R <<- Hi %*% R
    
    return(Hi)
  })
  
  Q <- Reduce("%*%", H)
  res <- list('Q'=Q,'R'=R)
  
  return(res)
}


##For a given matrix A, 
#find orthogonal matrix Q and upper triangular matrix R such that A = QR 
#using Givens rotations
QRfactorizationGivens <- function(mat) {
  
  givens <- function(a, b){
    c <- 0
    s <- 0
    if (b == 0){
      c <- 1
      s <- 0
    }else{
      if(abs(b) > abs(a)){
        tau <- -1 * a/b
        s <- 1 / sqrt(1 + tau^2)
        c <- s * tau
      }else{
        tau <- -1 * b/a
        c <- 1 / sqrt(1 + tau^2)
        s <- c * tau
      }
    }
    
    return(list("c" = c, "s" = s))
  }
  
  num_cols <- ncol(mat)
  num_rows <- nrow(mat)
  
  R <- mat
  
  num_iter <- min(num_cols, num_rows) - 1
  
  G <- diag(num_rows)
  sapply(1:num_iter, function(j){
    sapply(num_rows:(j+1), function(i){
      givens_curr <- givens(R[(i-1), j], R[i, j])
      rotation_mat <- matrix(c(givens_curr$c, -givens_curr$s,
                               givens_curr$s, givens_curr$c), ncol = 2)
      R[(i-1):i, j:num_cols] <<- t(rotation_mat) %*% R[(i-1):i, j:num_cols]
      G_curr <- diag(num_rows)
      G_curr[(i-1):i, (i-1):i] <- rotation_mat
      G <<- G %*% G_curr
    })
  })

  Q <- G
  res <- list('Q'=Q,'R'=R)
  
  return(res)
}

#Inverse of an upper triangular matrix
InverseUpperTriangular <- function(U){
  
  if(any(diag(U) == 0)){
    stop("Upper triangular matrix not invertible, 0 found in diagonal")
  }
  num_rows <- nrow(U)
  UI <- cbind(U, diag(num_rows))
  
  sapply(num_rows:1, function(curr_row){
    curr_pivot <- UI[curr_row, curr_row]
    UI[curr_row, ] <<- UI[curr_row, ]/curr_pivot
    if(curr_row > 1){
      sapply((curr_row-1):1, function(i){
        curr_row_mult <- UI[i, curr_row]
        UI[i, ] <<- UI[i, ] - UI[curr_row, ] * curr_row_mult
      })
    }
    
  })
  Uinv <- UI[, (num_rows+1):(2*num_rows)]
  return(Uinv)
}

##For a given matrix A and vector b, 
#find least square solution x for Ax = b
#using QR decomposition
LeastSquareSolutionByQR <- function(A, b){
  
  y <- matrix(data = b, ncol=1)
  
  num_rows <- nrow(A)
  num_cols <- ncol(A)
  if(num_cols > num_rows){
    stop("Underdetermined system")
  }
  
  if(nrow(y) != num_rows){
    stop("Incompatible dimensions of A and b")
  }
  
  QRfactor <- QRfactorizationHouseholder(A)
  
  Q <- QRfactor$Q
  R1 <- QRfactor$R[1:num_cols,]
  
  Rinv <- InverseUpperTriangular(R1)
  Qty <- t(Q) %*% y
  
  c <- Qty[1:num_cols]
  
  x_LS <- Rinv %*% c
  
  return(x_LS)
}


