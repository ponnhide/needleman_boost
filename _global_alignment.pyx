import numpy as np 
cimport numpy as np
cimport cython

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
    cdef int pre
    cdef int maxv
    cdef int[:,:] M
    cdef int[:,:] P
    cdef int[:,:] Q
    cdef int[:,:] S
    cdef int[:,:] T

    lx = len(x) 
    ly = len(y) 
    M = np.zeros((lx+1, ly+1), dtype=np.int32)
    P = np.zeros((lx+1, ly+1), dtype=np.int32)
    Q = np.zeros((lx+1, ly+1), dtype=np.int32)
    S = np.zeros((lx+1, ly+1), dtype=np.int32)
    T = np.zeros((lx+1, ly+1), dtype=np.int32)

    for i in range(1, lx+1):
        M[i,0] = -65535
        P[i,0] = -65535
        S[i,0] = -65535
        Q[i,0] = end_weight * (gap_start + i * gap_extension)
        #T[i,0] = 2
    
    for i in range(1, ly+1):
        M[0,i] = -65535
        P[0,i] = end_weight * (gap_start + i * gap_extension) 
        Q[0,i] = -65535
        S[0,i] = -65535
        #T[0,i] = 1

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
            S[i,j] = -65535
    
    S[lx,ly] = max(M[lx,ly], P[lx,ly], Q[lx,ly])
    for i in range(0, lx+1):
        k = lx - i
        for j in range(0, ly+1): 
            l = ly - j
            if S[k,l] == -65535:
                pass 
            else:
                if S[k,l] == M[k,l] and (T[k,l] == 1 or (i==0 and j ==0)):
                    if x[k-1] == y[l-1]:
                        v = match 
                    else:
                        v = mis_match  
                    if v + M[k-1,l-1] == S[k,l]:
                        S[k-1,l-1] = M[k-1,l-1] 
                        T[k-1,l-1] = 1
                    elif v + P[k-1,l-1] == S[k,l]:
                        S[k-1,l-1] = P[k-1,l-1]
                        T[k-1,l-1] = 2
                    elif v + Q[k-1,l-1] == S[k,l]:
                        S[k-1,l-1] = Q[k-1,l-1]
                        T[k-1,l-1] = 3

                elif S[k,l] == P[k,l] and (T[k,l] == 2 or (i==0 and j ==0)):
                    if gap_extension + gap_start + M[k,l-1] == S[k,l]:
                        S[k,l-1] = M[k,l-1]
                        T[k,l-1] = 1
                    elif gap_extension + P[k,l-1] == S[k,l] or k == 0:
                        S[k,l-1] = P[k,l-1]
                        T[k,l-1] = 2
                    elif gap_extension + gap_start + Q[k,l-1] == S[k,l]:
                        S[k,l-1] = Q[k,l-1]
                        T[k,l-1] = 3

                elif S[k,l] == Q[k,l] and (T[k,l] == 3 or (i==0 and j ==0)):
                    if gap_extension + gap_start + M[k-1,l] == S[k,l]:
                        S[k-1,l] = M[k-1,l]
                        T[k-1,l] = 1
                    elif gap_extension + gap_start + P[k-1,l] == S[k,l]:
                        S[k-1,l] = P[k-1,l]
                        T[k-1,l] = 2
                    elif gap_extension + Q[k-1,l] == S[k,l] or l == 0:
                        S[k-1,l] = Q[k-1,l]
                        T[k-1,l] = 3 

    return M, P, Q, S
