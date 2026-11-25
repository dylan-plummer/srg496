"""Find a 196-column Hadamard submatrix invariant under the order-seven action."""
from collections import Counter
from itertools import product
import numpy as np
from ortools.sat.python import cp_model
from family import HERE,make,power,mul8


def translate7(i):
    if i<16:
        return (i//8)*8+mul8(power(-1),i%8)
    p,x=divmod(i-16,8)
    if p<56:
        p=(p//28)*28 + 4*((p//4+1)%7)+p%4
    else:
        x=mul8(power(-1),x)
    return 16+8*p+x


def search(i,j):
    F=make(i,j)
    H=(F.T@F-14*np.eye(496,dtype=np.int64)).astype('int16')
    perm=[translate7(k) for k in range(496)]
    assert np.array_equal(H,H[np.ix_(perm,perm)])
    goals=[]
    visited=set()
    for k in range(496):
        if k in visited:
            continue
        orbit=[]
        z=k
        while z not in visited:
            visited.add(z)
            orbit.append(z)
            z=perm[z]
        if len(orbit)==7:
            goals.append(orbit)
    assert len(goals)==70
    idx=np.array(goals)
    model=cp_model.CpModel()
    var=[model.new_bool_var(str(a)) for a in range(70)]
    model.add(sum(var)==28)
    # Every selected vector must lie in the 28-eigenspace of the partial
    # frame operator. This tests the equality in the entire ambient space.
    gram=F.T@F
    for a,orbit in enumerate(goals):
        i0=orbit[0]
        summed=(F*gram[i0,:])[...,idx].sum(axis=2)
        for row,f in zip(summed,F[:,i0]):
            if np.any(row) or f:
                model.add(cp_model.LinearExpr.weighted_sum(var,row)==int(28*f)).only_enforce_if(var[a])
    solver=cp_model.CpSolver()
    solver.parameters.num_search_workers=2
    solver.parameters.max_time_in_seconds=45
    res=solver.solve(model)
    print(i,j,solver.status_name(res),flush=True)
    if res in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        X=idx[[a for a in range(70) if solver.value(var[a])]].ravel()
        M=H[np.ix_(X,X)]
        assert np.array_equal(M@M,196*np.eye(196))
        np.save(HERE/f'H196_{i}{j}.npy',M)
        print('HIT; regularization remains to be checked',flush=True)


if __name__=='__main__':
    for i,j in product(range(3),repeat=2):
        search(i,j)
