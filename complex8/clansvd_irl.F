      subroutine clansvd_irl(which,jobu,jobv,m,n,dim,p,neig,maxiter,
     c     aprod,U,ldu,Sigma,bnd,V,ldv,tolin,work,lwork,cwork,lcwrk,
     c     iwork,liwork,soption,ioption,info,cparm,iparm)

c     DLANSVD_IRL: Compute the leading singular triplets of a large
c     and sparse matrix A by implicitly restarted Lanczos bidiagonalization 
c     with partial reorthogonalization.
c
c     Parameters:
c
c     WHICH: CHARACTER*1. Decides which singular triplets to compute. 
c            If WHICH.EQ.'L' then compute triplets corresponding to the K
c            largest singular values. 
c            If WHICH.EQ.'S' then compute triplets corresponding to the K
c            smallest singular values. 
c     JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
c           Otherwise the array U is not touched.
c     JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular 
c           vectors. Otherwise the array V is not touched.
c     M:    INTEGER. Number of rows of A.
c     N:    INTEGER. Number of columns of A.
c     DIM:  INTEGER. Dimension of the Krylov subspace.
c     P:    INTEGER. Number of shift per restart.
c     NEIG: INTEGER. Number of desired singular triplets.  
c           NEIG <= MIN(DIM-P,M,N)
c     MAXITER: INTEGER. Maximum number of restarts.
c     APROD: Subroutine defining the linear operator A. 
c            APROD should be of the form:
c
c           SUBROUTINE DAPROD(TRANSA,M,N,X,Y,DPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           COMPLEX X(*),Y(*),DPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y = A * X.
c           If TRANSA.EQ.'T' then the function should compute the matrix-vector
c           product Y = A^H * X, where A^H is the conjugate transpose (adjoint)
c           of A .
c           The arrays IPARM and DPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,KMAX+1): COMPLEX array. On return the first K columns of U
c               will contain approximations to the left singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH)  singular values of A.
c               On entry the first column of U contains the starting vector
c               for the Lanczos bidiagonalization. A random starting vector
c               is used if U is zero.
c     LDU: INTEGER. Leading dimension of the array U. LDV >= M.
c     SIGMA(K): REAL array. On return Sigma contains approximation
c               to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     BND(K)  : REAL array. Error estimates on the computed 
c               singular values. The computed SIGMA(I) is within BND(I)
c               of a singular value of A.
c     V(LDV,KMAX): COMPLEX array. On return the first K columns of V
c               will contain approximations to the right singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     TOLIN: REAL. Desired relative accuracy of computed singular 
c            values. The error of SIGMA(I) is approximately 
c            MAX( 16*EPS*SIGMA(1), TOLIN*SIGMA(I) )
c     WORK(LWORK): REAL array. Workspace of dimension LWORK.
c     LWORK: INTEGER. Dimension of WORK. 
c            If JOBU.EQ.'N' and JOBV.EQ.'N' then  LWORK should be at least
c            M + N + 10*KMAX + 2*KMAX**2 + 5 + MAX(M+N,4*KMAX+4).
c            If JOBU.EQ.'Y' or JOBV.EQ.'Y' then LWORK should be at least
c            M + N + 10*KMAX + 5*KMAX**2 + 4 + 
c            MAX(3*KMAX**2+4*KMAX+4, NB*MAX(M,N)), where NB>1 is a block 
c            size, which determines how large a fraction of the work in
c            setting up the singular vectors is done using fast BLAS-3 
c            operation. 
c     ZWORK: COMPLEX array of dimension ????.
c     IWORK: INTEGER array. Integer workspace of dimension LIWORK.
c     LIWORK: INTEGER. Dimension of IWORK. Should be at least 8*KMAX if
c             JOBU.EQ.'Y' or JOBV.EQ.'Y' and at least 2*KMAX+1 otherwise.
c     DOPTION: REAL array. Parameters for LANBPRO.
c        soption(1) = delta. Level of orthogonality to maintain among
c          Lanczos vectors.
c        soption(2) = eta. During reorthogonalization, all vectors with
c          with components larger than eta along the latest Lanczos vector
c          will be purged.
c        soption(3) = anorm. Estimate of || A ||.
c        soption(4) = min relgap. Smallest relgap allowed between any shift
c                     the smallest requested Ritz value.
c
c     IOPTION: INTEGER array. Parameters for LANBPRO.
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical Gram-Schmidt. IF CGS.EQ.0 then 
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c
c     INFO: INTEGER. 
c         INFO = 0  : The K largest or smallest (depending on the 
c                     value of WHICH) singular triplets were computed 
c                     succesfully.
c         INFO = J>0, J<K: An invariant subspace of dimension J was found.
c         INFO = -1 : K singular triplets did not converge within KMAX
c                     iterations.   
c     ZPARM: COMPLEX array. Array used for passing data to the APROD
c         function.   
c     IPARM: INTEGER array. Array used for passing data to the APROD
c         function.   
c
c
c     (C) Rasmus Munk Larsen, Stanford University, 2000,2004
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 which,jobu,jobv
      integer m,n,p,neig,maxiter,ldu,ldv,iter,liwork,lcwrk
      integer iwork(liwork),lwork,info,ioption(*)
      complex U(ldu,*),V(ldv,*), cwork(lcwrk),cparm(*)
      real Sigma(*),bnd(*),work(lwork)
      real tolin,soption(*)
      integer iparm(*)
      external aprod

c     %------------%
c     | Parameters |
c     %------------%
      real one, zero, FUDGE
      parameter(one = 1.0, zero = 0.0, FUDGE = 1.01)
            
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer k,i,ibnd,iwrk,ierr,ip,iq,nconv,lwrk,kold,dim
      integer ialpha,ibeta,ialpha1,ibeta1,ishift,nshft,lapinfo
      real eps,eps34,epsn2,epsn,sfmin,anorm,rnorm,tol
      real shift, relgap
      integer st,cnk,wst,wcnk, tid, nt

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external szero,izero,scopy,sbdsqr
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      real slamch,pscnrm2
      external pscnrm2,lsame
      external slamch
#ifdef _OPENMP
      integer omp_get_thread_num, omp_get_num_threads
      external omp_get_thread_num, omp_get_num_threads
#endif

c-------------------- Here begins executable code ---------------------
      
c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = slamch('e')
      eps34 = eps**(3.0/4.0)
      epsn = real(max(m,n))*eps/2.0
      epsn2 = sqrt(real(max(m,n)))*eps/2.0
      sfmin = slamch('s')
      

c     %--------------------------------%
c     | Guard against absurd arguments |
c     %--------------------------------%
      dim = min(dim,n+1,m+1)
      k = dim-p
      tol = min(one,max(16.0*eps,tolin))
      anorm = zero

c     %------------------------------%
c     | Set pointers into work array |
c     %------------------------------%
      ibnd = 1
      ialpha = ibnd + dim+1
      ibeta = ialpha + dim
      ialpha1 = ibeta + dim
      ibeta1 = ialpha1 + dim
      ishift = ibeta1 + dim
      ip = ishift + dim
      iq = ip + (dim+1)**2
      iwrk = iq + dim**2
      lwrk = lwork-iwrk+1
      call szero(8*dim + 3 + 2*dim**2,work,1)

c     %---------------------------------------------------------------%
c     | Set up random starting vector if none is provided by the user |
c     %---------------------------------------------------------------%
      rnorm = pscnrm2(m,U(1,1),1)
      if (rnorm.eq.zero) then
         call cgetu0('n',m,n,0,1,U,rnorm,U,ldu,aprod,
     c        cparm,iparm, ierr,ioption(1),anorm,cwork)     
      endif

      iter = 0
      info = 0
      nconv = 0
      kold = 0       
c     %------------------------------%
c     | Iterate until convergence... |
c     %------------------------------%
      do while (nconv.lt.neig .and. iter.lt.maxiter)


c     %---------------------------------------------------%
c     | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
c     %---------------------------------------------------%

         call clanbpro(m, n, kold, dim, aprod, U, ldu, V, ldv,
     c        work(ialpha),dim,rnorm,soption(1),ioption(1),
     c        work(iwrk), cwork, iwork, cparm, iparm, ierr)         
         kold = k
c     %---------------------------------------------%
c     | Compute and analyze SVD(B) and error bounds |
c     %---------------------------------------------%
         call scopy(dim, work(ialpha),1,work(ialpha1),1)
         call scopy(dim, work(ibeta),1,work(ibeta1),1)
         call szero(dim+1,work(ibnd),1)

         call sbdqr((dim.eq.min(m,n)),'N',dim,work(ialpha1),
     c        work(ibeta1),work(ibnd+dim-1),work(ibnd+dim),
     c        work(ip),dim+1)
         
         call sbdsqr('u',dim,0,1,0,work(ialpha1),work(ibeta1),work,1,
     c        work(ibnd),1,work,1,work(iwrk),lapinfo)
         nbsvd = nbsvd + 1

         if (dim.gt.5) then
            anorm = work(ialpha1)
         else
            anorm = max(anorm,work(ialpha1))
         endif
         do i=1,dim
            work(ibnd+i-1) = abs(rnorm*work(ibnd+i-1))
c            write (*,*) 'bnd(',i,') = ',work(ibnd+i-1)
         enddo

c     %---------------------------------------------%
c     | Refine error bounds using the "Gap theorem" |
c     %---------------------------------------------%
         if (lsame(which,'s')) then
            call srefinebounds(min(m,n),dim,work(ialpha1),work(ibnd),
     c           epsn*anorm,eps34) 
         else
            call srefinebounds(min(m,n), min(dim,neig),work(ialpha1),
     c           work(ibnd),epsn*anorm,eps34) 
         endif
c         do i=1,dim
c            write (*,*) 'bnd(',i,') = ',work(ibnd+i-1)
c         enddo
         
c     %----------------------------------------------------%
c     | Determine the number of converged singular values  |
c     %----------------------------------------------------%
c         do i=1,min(dim,neig)
c            write(*,*)  iter,i,work(ialpha1+i-1),work(ibnd+i-1)
c         enddo
         if (lsame(which,'s')) then
            i = dim-neig+1
            nconv = 0
            do while(i.le.dim)
c               write(*,*) 'iter = ',iter,'sigma = ',work(ialpha1+i-1),
c     c              'error = ',work(ibnd+i-1)
c               write(*,*)  iter,work(ialpha1+i-1),work(ibnd+i-1)
               if (work(ibnd+i-1).le.tol*work(ialpha1)) then
                  nconv = nconv + 1
                  sigma(nconv) = work(ialpha1+i-1)
                  bnd(nconv) = work(ibnd+i-1)
               endif
               i = i+1
            enddo
         else
            i = 1
            nconv = 0
            do while(i.le.min(dim,neig))
c               write(*,*)  iter,work(ialpha1+i-1),work(ibnd+i-1)
               if (work(ibnd+i-1).le.tol*work(ialpha1+i-1)) then
                  nconv = nconv + 1
                  sigma(nconv) = work(ialpha1+i-1)
                  bnd(nconv) = work(ibnd+i-1)
                  i = i+1
               else
                  i = k+1
               endif
            enddo
         endif
            
         
c     %-----------------------------------------------%
c     | Test if an invariant subspace have been found |
c     %-----------------------------------------------%
         if (ierr.lt.0) then
            if (dim.lt.k) then
               write(*,*) 'WARNING: Invariant subspace found.',
     c              ' Dimension = ',dim
               info = dim
            endif
            goto 50               
         endif
         if (nconv.lt.neig) then

c     %---------------------------------------------------------------------%
c     | Implicit restart:                                                   |
c     |   Apply  shifts mu_1, mu_2,...,mu_p to "back up" the          |
c     |   bidiagonalization (dim-k) steps to                                |
c     |                                                                     |
c     |      A V_{k}^{+} = U_{k+1}^{+} B_{k}^{+}                            |
c     | corresponding to starting vector                                    |
c     |      u_1^{+} = \prod_{i=1}^{p} (A A^T - mu_i^2) u_1             |
c     |                                                                     |
c     | We use exact shifts mu_i for which the relative gap between mu_i    |
c     | and the lower bound on the k'th Ritzvalue is larger than soption(4) |
c     %---------------------------------------------------------------------%
            call szero(dim-k,work(ishift),1)
            nshft = 0
            if (lsame(which,'s')) then
               do i=1,k
                  relgap = (work(ialpha1+i-1)-work(ibnd+i-1) -
     c                 work(ialpha1+dim-neig-1))
                  if (relgap.gt.soption(4)*work(ialpha1+dim-neig-1))
     c                 then
                     work(ishift + nshft) = work(ialpha1+i-1)
                  else
                     work(ishift + nshft) = work(ialpha1)
                  endif
c                  write (*,*) 'shift = ',work(ishift + nshft)
                  nshft = nshft + 1
               enddo      
            else
               do i=dim,k+1,-1
                  relgap = work(ialpha1+k-1) -
     c                 (work(ialpha1+i-1)+work(ibnd+i-1))
                  if (relgap.gt.soption(4)*work(ialpha1+k-1)) then
                     work(ishift + nshft) = work(ialpha1+i-1)
                  else
                     work(ishift + nshft) = zero
                  endif
                  nshft = nshft + 1
               enddo      
            endif

c     %--------------------------------------------------%
c     | Apply shifts and accumulate rotations such that  |
c     |   B_{dim}^{+} = P * B_{dim} * Q^T                |
c     %--------------------------------------------------%
            call szero((dim+1)*(dim+1),work(ip),1)
            call szero(dim*dim,work(iq),1)
            do i=1,dim+1
               work(ip+(i-1)*(dim+2)) = one
            enddo
            do i=1,dim
               work(iq+(i-1)*(dim+1)) = one
            enddo

            do i=dim,k+1,-1
               shift = work(ishift+dim-i)
               call sbsvdstep('y','y',dim+1,dim,i,shift,work(ialpha),
     c              work(ibeta),work(ip),dim+1,work(iq), dim)
            enddo

c     %---------------------------------------------------%
c     | Compute first k+1 left and first k right updated  |
c     | Lanczos vectors                                   |
c     |   U_{dim+1}^{+} = U_{dim+1} * P(:,1:k+1)          |
c     |   V_{dim}^{+} = V_{dim} * Q(:,1:k)                |
c     %---------------------------------------------------%
#ifdef _OPENMP
c$OMP PARALLEL private(tid,nt,cnk,st,wcnk,wst) 
            tid = omp_get_thread_num()
            nt = omp_get_num_threads()
#else
            tid = 0
            nt = 1
#endif
            wcnk = lcwrk/nt
            wst = tid*wcnk+1
            cnk = m/nt
            st = tid*cnk+1
            if (tid.eq.nt-1) then
               wcnk = lcwrk-wst+1
               cnk = m-st+1
            endif
            call csgemm_ovwr_left('n',cnk,k+1,dim+1,U(st,1),ldu,
     c           work(ip),dim+1,cwork(wst),wcnk)            
            cnk = n/nt
            st = tid*cnk+1
            if (tid.eq.nt-1) then
               cnk = n-st+1
            endif
            call csgemm_ovwr_left('n',cnk,k,dim,V(st,1),ldv,work(iq),
     c           dim,cwork(wst),wcnk)
#ifdef _OPENMP
c$OMP END PARALLEL
#endif
            rnorm = work(ibeta+k-1)
            nrestart = nrestart + 1
         endif         
         iter = iter + 1
      enddo

 50   if ((nconv.ge.neig .or. info.gt.0) .and. 
     c     (lsame(jobu,'y') .or. lsame(jobv,'y'))) then
c     %-----------------------------------------%
c     | Calculate singular vectors if requested %
c     %-----------------------------------------%
         call scopy(dim, work(ialpha),1,work(ialpha1),1)
         call scopy(dim, work(ibeta),1,work(ibeta1),1)
         lwrk = lwrk + dim**2 + (dim+1)**2
         call critzvec(which, jobu,jobv,m,n,nconv,dim,work(ialpha1),
     c        work(ibeta1),work(ialpha1),U,ldu,V,ldv,work(ip),
     c        lwrk,cwork,lcwrk,iwork)
      endif
      neig = nconv
      nlandim = dim
      end
