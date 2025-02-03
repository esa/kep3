import numpy as _np
import pykep as _pk

def primer_vector(DVi, DVj, Mji, Mjk):
    """This function computes the primer vector in a point k, relative
    to impulses given at i and j.
    
    Args:
        *DVi* (:class:`ndarray` - (3,)): the impulse at point i.
        
        *DVj* (:class:`ndarray` - (3,)): the impulse at point j.
        
        *Mji* (:class:`ndarray` - (6,6)): the state transition matrix from i to j (dxj = Mji dxi).
        
        *Mjk* (:class:`ndarray` - (6,6)): the state transition matrix from k to j (dxj = Mji dxi).
        
        Returns:
            :class:`tuple`: The primer vector, the Aik matrix, the Ajk matrix.
            
        Note:
            The impulse transfer matrix Anm is defined as those matrices that allow to compute the
            variation of the impulse at point n given the variation of the impulse at point m. In formal terms, 
            dDVn = Anm dDVm. All variations are such that the terminal state is kept fixed.
    """
    Aik = -(_np.linalg.inv(Mji[:3,3:]))@Mjk[:3,3:]
    Ajk = -(Mji[3:,3:]@Aik + Mjk[3:,3:])
    p = - Aik.T@DVi/_np.linalg.norm(DVi) - Ajk.T@DVj/_np.linalg.norm(DVj)
    return p, Aik, Ajk

def primer_vector_surrogate(DVk, Mki, Mkj):
    Aij = -(_np.linalg.inv(Mki[:3, 3:]) @ Mkj[:3, 3:])
    Akj = -(Mki[3:, 3:] @ Aij + Mkj[3:, 3:])
    B = Aij
    b =  - Akj.T@DVk / _np.linalg.norm(DVk)
    p_surrogate_norm, u_star = _pk.trajopt.minBu_bu_p(B, b)
    p_surrogate = p_surrogate_norm * u_star
    return p_surrogate, Aij, Akj
    
