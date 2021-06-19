c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c
      program example

      implicit none

      include 'matvec.h'

c     %------------%
c     | Parameters |
c     %------------%
      integer lwrk,liwrk
      parameter(lwrk=mmax+nmax+13*kmax+8*kmax**2+32*mmax+8)
      parameter(liwrk = 8*kmax)
      double precision dzero, done
      parameter(dzero = 0.0, done = 1.0)

      character*1 jobu, jobv
      character*256 matfilename, sigmafilename,ufilename,vfilename
      integer i,m,n,neig,ioption(10),iwork(liwrk)
      integer iparm(1),info, compuv, uvform
      double precision U(mmax,kmax+1),V(nmax,kmax),work(lwrk)
      double precision doption(10),Sigma(kmax),dparm(1)
      double precision tol,bnd(kmax), eps
c     Save large arrays in a common block to avoid stack overflow on
c     certain machines.
      common/cexample/U,V,Sigma,bnd,work,iwork

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch
      external dlamch, matvecmult

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
         write (*,'(a, I6, a, i6)') 'Matrix can be at most',mmax,' by ',
     c        nmax
         stop 'Matrix exceeded maximum size defined in matvec.h.'
      endif
         
c
c     Set the first column of U to zero. This will make dlansvd_irl
c     generate a random starting vector.      
c
      do i=1,m
         U(i,1) = 0d0
      enddo

c
c     Clear statistics counters and timers
c
      call clearstat

c
c     Set options for DLANSVD
c
      eps = dlamch('e')
      ioption(1) = 0
      ioption(2) = 1
      doption(1) = sqrt(eps)
      doption(2) = eps**(3.0/4.0)
      doption(3) = dzero
c
c     Compute the neig largest singular values
c
      write (*,*) 'Computing the neig=',neig,' largest singular ' //
     c     'values and corresponding singular vectors'
      call dlansvd(jobu,jobv,m,n,neig,kmax,matvecmult,U,mmax,Sigma,bnd,
     c     V,nmax,tol,work,lwrk,iwork,liwrk,doption,ioption,info,
     c     dparm,iparm)

c
c     Print singular values and computed error bounds.
c
      write(*,*)
      write(*,*) '  i          sigma(i)              bnd(i)'
      open(10,file=sigmafilename,status='unknown')
      do i=1,neig
         write(*,'(1I5,1p,1e25.15,1p,1e12.3)') i,Sigma(i),bnd(i)
         write(10,'(1I5,1p,1e25.15,1p,1e12.3)') i,Sigma(i),bnd(i)
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
c     Print runtime statistics 
c
      call printstat
      end

