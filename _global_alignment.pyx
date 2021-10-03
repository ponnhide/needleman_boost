import numpy as np 
cimport numpy as np
cimport cython
from cpython cimport array

@cython.boundscheck(False)
@cython.wraparound(False)
def align(char *x, char *y, int match, int mis_match, int gap_extension, int gap_start, int end_weight):
    #{"1":h, "2":v, "3":d, "4":hv, "5":hd, "6":dv, "7",hdv}
    cdef int lx
    cdef int ly
    cdef int m
    cdef int p
    cdef int q
    cdef int i
    cdef int j 
    cdef int k
    cdef int l 
    cdef int v
    cdef int s 
    cdef int t 
    cdef int n
    cdef int pre
    cdef int maxv
    cdef int[:,:] M
    cdef int[:,:] P
    cdef int[:,:] Q
    cdef array.array[char] arrx, arry, template = array.array('b')
    lx  = len(x) 
    ly  = len(y)
    arrx = array.clone(template, lx + ly, False)
    arry = array.clone(template, lx + ly, False)
    M = np.zeros((lx+1, ly+1), dtype=np.int32)
    P = np.zeros((lx+1, ly+1), dtype=np.int32)
    Q = np.zeros((lx+1, ly+1), dtype=np.int32)

    for i in range(1, lx+1):
        M[i,0] = -65535
        P[i,0] = -65535
        Q[i,0] = end_weight * (gap_start + i * gap_extension)
    
    for i in range(1, ly+1):
        M[0,i] = -65535
        P[0,i] = end_weight * (gap_start + i * gap_extension) 
        Q[0,i] = -65535

    for i in range(1, lx+1):
        for j in range(1, ly+1):
            if x[i-1] == y[j-1]:
                v = match 
            else:
                v = mis_match 

            M[i,j] = v + max(M[i-1,j-1], P[i-1,j-1], Q[i-1,j-1])
            P[i,j] = max(
                    gap_extension + gap_start + M[i,j-1],
                    gap_extension + P[i,j-1],
                    gap_extension + gap_start + Q[i,j-1]
            )

            Q[i,j] = max(
                    gap_extension + gap_start + M[i-1,j],
                    gap_extension + gap_start + P[i-1,j],
                    gap_extension + Q[i-1,j]
            )
    
    
    k = lx
    l = ly
    s = max(M[lx,ly], P[lx,ly], Q[lx,ly])
    
    if M[k,l] == s:
        t = 1
    elif P[k,l] == s:
        t = 2
    else:
        t = 3

    n = 0
    while k > 0 or l > 0:
        if s == M[k,l] and t == 1:
            if x[k-1] == y[l-1]:
                v = match 
            else:
                v = mis_match  
            if v + M[k-1,l-1] == s:
                s = M[k-1,l-1] 
                t = 1
            elif v + P[k-1,l-1] == s:
                s = P[k-1,l-1]
                t = 2
            elif v + Q[k-1,l-1] == s:
                s = Q[k-1,l-1]
                t = 3
            
            arrx[n] = x[k-1]
            arry[n] = y[l-1]
            k = k - 1
            l = l - 1

        elif s == P[k,l] and t == 2:
            if gap_extension + gap_start + M[k,l-1] == s:
                s = M[k,l-1]
                t = 1
            elif gap_extension + P[k,l-1] == s or k == 0:
                s = P[k,l-1]
                t = 2
            elif gap_extension + gap_start + Q[k,l-1] == s:
                s = Q[k,l-1]
                t = 3
            
            arrx[n] = "-"
            arry[n] = y[l-1]
            k = k
            l = l - 1

        elif s == Q[k,l] and t == 3:
            if gap_extension + gap_start + M[k-1,l] == s:
                s = M[k-1,l]
                t = 1
            elif gap_extension + gap_start + P[k-1,l] == s:
                s = P[k-1,l]
                t = 2
            elif gap_extension + Q[k-1,l] == s or l == 0:
                s = Q[k-1,l]
                t = 3 
            
            arrx[n] = x[k-1]
            arry[n] = "-"
            k = k - 1
            l = l
        n += 1
    arrx = arrx[:n][::-1]
    arry = arry[:n][::-1]
    return M, P, Q, arrx, arry
