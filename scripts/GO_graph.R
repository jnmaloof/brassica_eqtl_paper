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

# next! duv figure these out. 