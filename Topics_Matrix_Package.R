#*************************
# Learn "Matrix" package
#*************************
#ref: Package ‘Matrix’
# Version 1.4-1
# Date 2022-03-21
# Title Sparse and Dense Matrix Classes and Methods

library(Matrix)


#========
# band()
#========

# Descption:
  # Returns a new matrix formed by extracting the lower triangle (tril) or
  # the upper triangle (triu)
  # or a general band relative to the diagonal (band), 
  # and setting other elements to zero

# Usage
#band(x, k1, k2, ...)
#tril(x, k = 0, ...)
#triu(x, k = 0, ...)


# Arguments
# x: a matrix-like object
# k,k1,k2: 
  # integers specifying the diagonal bands that will not be set to zero. 
  # relative to the main diagonal, which is k=0. 
  # A negative value of k indicates a diagonal below the main diagonal 
  # and a positive value indicates a diagonal
  # above the main diagonal.


# Value:
  # An object of an appropriate matrix class
  # result is of class sparseMatrix only if x is


# Methods
  #x = "CsparseMatrix" method for compressed, sparse, column-oriented matrices.
  #x = "TsparseMatrix" method for sparse matrices in triplet format.
  #x = "RsparseMatrix" method for compressed, sparse, row-oriented matrices.
  #x = "ddenseMatrix" method for dense numeric matrices, including packed numeric matrices.


set.seed(15-04-2022)
m <- matrix(0, 5, 5) # 25 elements
length(m) # 1
m[sample(length(m), size = 14, replace = T)] <-rep(1:9, length=14)
(m_Csp<- as(m, "CsparseMatrix"))
(m_Rsp<- as(m, "RsparseMatrix"))

tril(m_Csp)
triu(m_Csp)

tril(m_Csp, -1) # strict lower tri w/o  main diag
triu(m_Csp, 1) # stric upper tri w/o main diag

band(m_Csp, -1, 2)
# [1,] 1 8 2 . .
# [2,] 5 . 9 . .
# [3,] . . . . 4
# [4,] . . 3 1 5
# [5,] . . . 3 .



m5 <- matrix(rnorm(25), nrow = 5)
m5
tril(m5)  
tril(m5, -1)
triu(m5, 1)
band(m5, -1, 2)


m65 <- matrix(rnorm(30), nrow = 6)
# not square
triu(m65)
sym5 <- crossprod(m65)
sym5
band(sym5, -1, 1)
as(band(sym5, -1, 1), "sparseMatrix")
# prefer


#===========
# Block diag
#===========
# bdiag()

#Description
  #Build a block diagonal matrix given several building block matrices.

# Usage
  #bdiag(.)
  #.bdiag(lst)

# Arguments
  #... individual matrices or a list of matrices
  # non-empty list of matrice

# Value
  # A *sparse* matrix obtained by combining the arguments 
    # into a block diagonal matrix

# Note:
  # efficient for the case of relatively few block matrices which
    # are typically sparse themselves
  # inefficient for the case of many small dense block matrices
  # bdiag_m() function is an order of magnitude faster

bdiag(matrix(1:4, nrow = 2), diag(3))

## combine "Matrix" class and traditional matrices:
Diagonal(2)  # 
# 2 x 2 diagonal matrix of class "ddiMatrix"
      #[,1] [,2]
#[1,]    1    .
#[2,]    .    1


# traditional diag matrix
diag(3:2)

bdiag(Diagonal(2), matrix(1:3, 3,4), diag(3:2))


mlist <- list(1, 2:3, diag(x=5:3), 27, cbind(1,3:6), 100:101)
bdiag(mlist)


stopifnot(identical(bdiag(mlist), 
                    bdiag(lapply(mlist, as.matrix))))


??nMatrix # numerical matrix ?

ml <- c(as(matrix((1:24)%% 11 == 0, 6,4),"nMatrix"),
        rep(list(Diagonal(2, x=TRUE)), 3))

is(bdiag(ml), "lsparseMatrix")


## random (diagonal-)block-triangular matrices:

Rd_Blk_Tri<- function(nb, max.ni, lambda = 3) {
  .bdiag(replicate(
    nb, 
    expr = {
      n <- sample.int(max.ni, 1)
      tril(Matrix(rpois(n*n, lambda = lambda), n, n))
    }))
}

rd3 <- Rd_Blk_Tri(nb = 3, max.ni = 10)
image(rd3)

rd12 <- Rd_Blk_Tri(12, 20)
image(rd12)


Rd_blk_bnd <- function(nb, max.int, lambda, lbd, ubd) {
  .bdiag(replicate(nb, expr = {
    n <- sample.int(max.int, 1)
    band(Matrix(rpois(n * n, lambda = lambda), n, n), 
         lbd, ubd)
  }))
}
B12 <- Rd_blk_bnd(nb = 12, max.int = 20, 
                  lambda = 3, lbd = -2, ubd = 1)
image(B12)


# try matrix value from N(0, 1)
Rd_Blk_N <- function(nb, max.int) {
  .bdiag(
    replicate(nb, expr = {
      n <- sample.int(max.int, 1)
      Matrix(rnorm(n * n), n, n)
    }))
}

N1 <-Rd_Blk_N(nb = 10, max.int = 10)
image(N1)


#============
# bandSparse 
#============

# Construct Sparse Banded Matrix from (Sup-/Super-) Diagonals
# Construct a sparse banded matrix by specifying its non-zero sup- and super-diagonals

# Usage
# bandSparse(n, m = n, k, diagonals, symmetric = FALSE,
           # repr = "C", giveCsparse = (repr == "C"))

#Arguments
  # n,m the matrix dimension (n;m) = (nrow; ncol).
  # k integer vector, identical meaning as in band(*,k),
    # relative to the main diagonal, which is k=0.
  # diagonals: optional list of sub-/super- diagonals; 
    # if missing, the result will be a pattern matrix,

  # repr character string, one of "C", "T", or "R", 

# Value
  # a sparse matrix (of class CsparseMatrix) of dimension n m 
  # with diagonal “bands” as specified

diags <- list(1:30, 10*(1:20), 100*(1:20))
s1 <- bandSparse(13, k = -c(0:2, 6), diag = c(diags, diags[2]), symm=TRUE)
s1
s2 <- bandSparse(13, k = c(0:2, 6), diag = c(diags, diags[2]), symm=TRUE)
s2

identical(s1, t(s2))
stopifnot(identical(s1, t(s2)), is(s1,"dsCMatrix"))


## a pattern Matrix of *full* (sub-)diagonals:
bk <- c(0:4, 7,9)
(s3 <- bandSparse(30, k = bk, symm = TRUE))


n <- 1e4
bk <- c(0:5, 7,11)
bMat <- matrix(1:8, n, 8, byrow=TRUE)
bLis <- as.data.frame(bMat)
B <- bandSparse(n, k = bk, diag = bLis)
Bs <- bandSparse(n, k = bk, diag = bLis, symmetric=TRUE)

B [1:15, 1:30]
Bs[1:15, 1:30]

### can use a list *or* a matrix for specifying the diagonals:
identical(B, bandSparse(n, k = bk, diag = bMat))
# [1] TRUE

identical(Bs, bandSparse(n, k = bk, diag = bMat, symmetric = T))
# [1] TRUE



##=======
# colSums
##=======

(M <- bdiag(Diagonal(2), matrix(1:3, 3,4), diag(3:2)))
colSums(M)

d <- Diagonal(10, c(0,0,10,0,2,rep(0,5)))

MM <- kronecker(d, M)
dim(MM)
MM@x
length(MM@x)
MM <- drop0(MM)
length(MM@x)


cm <- colSums(MM)
sp_cm <- colSums(MM, sparseResult = T)

is(sp_cm, "sparseVector") # [1] TRUE
str(sp_cm)
identical(cm, as.numeric(sp_cm))
# [1] TRUE

rowSums(MM, sparseResult = T)
# sparse vector (nnz/length = 14/70) of class "dsparseVector"









