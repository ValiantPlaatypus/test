from scipy.sparse import coo_matrix, bmat

A = coo_matrix([[1, 2], [3, 4]])
B = coo_matrix([[5], [6]])
C = coo_matrix([[7]])
print(' A: ')
print(A)
D = bmat([[A, A]])

print(' D: ')
print(D.toarray())


