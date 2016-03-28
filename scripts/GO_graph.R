# testing 
test <- matrix(1:9, 3)
test
?eigen
eigen(test)

test <- matrix(1:9, 2)
test
eigen(test)

# other matrix decompositions for non-square matrices
# SVD
?svd
svd(test)

# next! DUV figure these out. 
# D vector min of n or p

library("qpgraph")

