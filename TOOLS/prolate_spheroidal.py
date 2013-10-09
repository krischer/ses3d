""" package for prolate spheroidal wave functions
"""

from pylab import *
from scipy import linalg

##########################################################################################

def PSMat(N,W):

    """ define the NxN matrix for the prolate spheroidal sequence

    S=PSMat(N,W)

    INPUT:
    N=length of the index series
    W=frequency bandwidth W in [-0.5, 0.5]

    OUTPUT:
    S=NxN matrix in the algebraic eigenvalue problem
    """

    s=zeros((N,N))

    for i in range(N):
        for j in range(max(0,i-1),1+min(i+1,N-1)):

            if j==i-1:
                s[i,j]=0.5*i*(N-i)
            elif j==i:
                s[i,j]=((0.5*(N-1)-i)**2)*cos(2*pi*W)
            elif j==i+1:
                s[i,j]=0.5*(i+1)*(N-i-1)

    S=matrix(s)

    return S

##########################################################################################

def PSSeq(N,W):

    """ compute prolate spheroidal sequences as solutions of the eigenvalue equation 118 in Slepian, 1977.
    The sequences are normalised to 1 in the frequency interval [-W,W].

    eigval,c = PSSeq(N,W)

    INPUT:
    N=length of the index series
    W=frequency bandwidth, W in [-0.5, 0.5]

    OUTPUT:
    eigval[k] is the eigenvalue of the k-th PSWF
    c[:,k] is the corresponding ith eigenvector
    """

    #- solve eigenvalue problem

    S=PSMat(N,W)

    l,c=linalg.eig(S)

    #- normalise and compute eigenvalues of the PSWFs

    const=zeros(N,dtype=float)
    eigval=zeros(N,dtype=float)

    for n in range(N):
        eigval_help=zeros(N,dtype=float)
        for m in range(N):
            const=const+c[n,:]*c[m,:]*sinc(2*W*(n-m))
            eigval_help=eigval_help+c[m,:]*sinc(2*W*(n-m))

        eigval=eigval+abs(eigval_help)**2

    const=sqrt(2*W*const)
    for k in range(N):
        c[:,k]=c[:,k]/const[k]

    eigval=4*W*W*eigval/(const**2)

    #- sort the eigenvalues in descending order

    sort_index=argsort(-eigval)

    return eigval[sort_index],c[:,sort_index]

##########################################################################################

def PSFreq(N,W,k_max,nu):

    """ compute frequency-domain prolate spheroidal wave function
    These prolate spheroidal wave functions are normalised to 1 in the interval [-W,W].

    eigval,psi=PSFreq(N,W,k,nu)

    INPUT:
    N=length of the index series
    W=frequency bandwidth, W in [-0.5, 0.5]
    k_max=maximum index of the prolate spheroidal wave function (between 0 and N-1)
    nu=frequency axis

    OUTPUT:
    eigval: eigenvalues 0 to k_max
    Psi[:,k]: frequency-domain prolate spheroidal wave functions with index k (k in [0, k_max])
    """

    #- compute eigenvalues and eigenvectors

    eigval,c=PSSeq(N,W)

    #- compute frequency-domain prolate spheroidal wave function

    Psi=zeros((len(nu),k_max+1),dtype=complex)

    for k in range(k_max+1):
        for n in range(N):
            Psi[:,k]=Psi[:,k]+c[n,k]*exp(-2*pi*1j*n*nu)

        Psi[:,k]=Psi[:,k]*exp(1j*pi*(N-1)*nu)

    return eigval[0:k_max+1],Psi

##########################################################################################

def PSIdx(N,W,k_max,n):

    """ compute index-domain prolate spheroidal wave function

    eigval,psi=PSIdx(N,W,k,n)

    INPUT:
    N=length of the index series
    W=frequency bandwidth, W in [-0.5, 0.5]
    k_max=maximum index of the prolate spheroidal wave function (between 0 and N-1)
    n=index axis

    OUTPUT:
    eigval=eigenvalues from index 0 to k_max
    psi[:,k]=k-th prolate spheroidal sequence normalised to 1 in [-W,W] in the frequency domain
    """

    #- compute eigenvectors

    eigval,c=PSSeq(N,W)

    #- compute index-domain prolate spheroidal wave function

    psi=zeros((len(n),k_max+1),dtype=complex)

    for k in range(k_max+1):
        for m in range(N):
            psi[:,k]=psi[:,k]+2*W*c[m,k]*sinc(W*(N-1-2*m+2*n))

    return eigval[0:k_max+1],psi

##########################################################################################

def PSPlot(N,W,k):

    """ plot prolate spheroidal wave functions

    INPUT:
    N=length of the index series
    W=frequency bandwidth, W in [-0.5, 0.5]
    k=index of the prolate spheroidal wave function (between 0 and N-1)
    """
    nu=arange(-0.5,0.5,0.002)
    n=arange(-3*N,3*N+1)

    eigval,Psi=PSFreq(N,W,N-1,nu)
    eigval,psi=PSIdx(N,W,N-1,n)

    #- plot results

    subplot(211)
    plot(nu,real(Psi[:,k]),color=[0.8,0.8,0.8],linewidth=2)
    plot(nu,imag(Psi[:,k]),color=[0.8,0.8,0.8],linewidth=2)
    plot(nu[abs(nu)<=W],real(Psi[abs(nu)<=W,k]),'k')
    plot(nu[abs(nu)<=W],imag(Psi[abs(nu)<=W,k]),'k--')
    xlabel('normalised frequency, nu')
    ylabel('\Psi(nu)')
    grid(True)

    subplot(212)
    plot(n,real(psi[:,k]),color=[0.7,0.7,0.7],linewidth=2)
    plot(n,real(psi[:,k]),'kx')
    xlabel('n')
    ylabel('psi_n')
    grid(True)

    show()

