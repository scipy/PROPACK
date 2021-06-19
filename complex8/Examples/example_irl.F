c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c
      program example

      implicit none

      include 'matvec.h'

c     %------------%
c     | Parameters |
c     %------------%
      integer ldwrk,maxiter,liwrk,lcwrk
      parameter(ldwrk=mmax+nmax+14*kmax+8*kmax**2+32*mmax+9)
      parameter(liwrk = 8*kmax,lcwrk=mmax+nmax+32*mmax)
      parameter(maxiter=1000)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      character*1 jobu, jobv
      character*256 matfilename, sigmafilename,ufilename,vfilename
      integer dim,i,m,n,neig,p,ioption(10),iwork(liwrk)
      integer info,iparm(1), compuv, uvform
      complex U(mmax,kmax+1),V(nmax,kmax),cwork(lcwrk),cparm(1)
      real soption(10),Sigma(kmax),swork(ldwrk)
      real tol,bnd(kmax),eps
c     Save large arrays in a common block to avoid stack overflow on
c     certain machines.
      common/cexample/U,V,Sigma,bnd,swork,cwork,iwork

c     %--------------------%
c     | External Functions |
c     %--------------------%
      real slamch
      external slamch, matvecmult

c-------------------- Here begins executable code ---------------------

c      
c     Get input parameters from stdin.
c
      write (*,*) 'Format of matrix file ' //
     c     '(1=Harwell-Boeing, 2=dense ASCII, 3=dense binary,' //
     c     ' 4=coordinate ASCII, 5=coordinate binary, 6=diagonal' //
     c     'ASCII, 7=diagonal binary)?'
      read (5,*) matform 
      write(*,*) matform
      write (*,*) 'Name of matrix input file?'
      read (5,'(a)') matfilename
      write(*,*) matfilename
      write (*,*) 'Desired number of singular values?'
      read (5,*) neig
      write (*,*) neig
      write (*,*) 'Desired relative accuracy?'
      read (5,*) tol
      write (*,*) tol
      write (*,*) 'Maximal dimension of subspace?'
      read (5,*) dim
      write (*,*) dim
      write (*,*) 'Number of shifts per restart?'
      read (5,*) p
      write (*,*) p
      write (*,*) 'Name of singular values output file?'
      read (5,'(a)') sigmafilename
      write(*,*) sigmafilename
      write (*,*) 'Should singular vectors be computed (0=no, 1=yes)?'
      read (5,*) compuv
      write (*,*) compuv
      if (compuv.eq.1) then
         jobu='y'
         jobv='y'
         write (*,*) 'File format for singuler vectors' //
     c        ' (1=ASCII, 2=binary)?'
         read (5,*) uvform
         write (*,*) uvform
         write (*,*) 'Name of U output file?'
         read (5,'(a)') ufilename
         write(*,*) ufilename
         write (*,*) 'Name of V output file?'
         read (5,'(a)') vfilename
         write(*,*) vfilename
      else
         jobu='n'
         jobv='n'
      endif

c
c     Read matrix in Harwell-Boeing format.
c
      write (*,*) 'Reading matrix from file ',matfilename
      call readmatrix(matfilename,m,n)
      write (*,*) '...done'
      write (*,*) 'Dimensions of matrix: m,n = ',m,n
      if (m.gt.mmax .or. n.gt.nmax) then
         write (*,'(a, I6, a, i6)') 'Matrix kan be at most',mmax,' by ',
     c        nmax
         stop 'Matrix exceeded maximum size defined in matvec.h.'
      endif
         
c
c     Set the first column of U to zero. This will make dlansvd_irl
c     generate a random starting vector.      
c
      call czero(m,U,1)

c
c     Clear statistics counters and timers
c
      call clearstat

c
c     Set options for DLANSVD_IRL
c
      eps = slamch('e')
      ioption(1) = 0
      ioption(2) = 1      
      soption(1) = sqrt(eps)
      soption(2) = eps**(3.0/4.0)
      soption(3) = 0.0
      soption(4) = 0.002
c
c     Compute the neig largest singular values
c
      write (*,*) 'Computing the neig=',neig,' largest singular ' //
     c     'values and corresponding singular vectors'
      dim = min(dim,kmax)
      p = min(dim-neig,p)
      if (p.le.0) then
         write (*,*) 'ERROR: min(dim,kmax)-neig must be larger than 0.'
         stop 'Invalid dimensions or not enough workspace.'
      endif
      call clansvd_irl('l',jobu,jobv,m,n,dim,p,neig,maxiter,matvecmult,
     c     U,mmax,Sigma,bnd,V,nmax,tol,swork,ldwrk,cwork,lcwrk,iwork,
     c     liwrk,soption,ioption,info,cparm,iparm)

c
c     Print singular values and computed error bounds.
c
      open(10,file=sigmafilename,status='unknown')
      write(*,*)
      write(*,*) '  i          sigma(i)              bnd(i)'
      do i=1,neig
         write(*,'(1I5,1p,1e15.7,1p,1e12.3)') i,Sigma(i),bnd(i)
         write(10,'(1I5,1p,1e15.7,1p,1e12.3)') i,Sigma(i),bnd(i)
      enddo 
      close(10)
c
c     Write singular vectors.
c
      if (compuv.eq.1) then
         write (*,*) 'Writing left singular vectors U to ',ufilename
         call writedense(ufilename,uvform+1,m,neig,U,mmax)
         write (*,*) 'Writing right singular vectors V to ',vfilename
         call writedense(vfilename,uvform+1,n,neig,V,nmax)
      endif

c
c     Print runtine statistics 
c
      call printstat
      end

