c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c
      program example

      implicit none

      include 'matvec.h'

      character*128 sigmafilename1,ufilename1,vfilename1
      character*128 sigmafilename2,ufilename2,vfilename2
      integer i,k,k1,k2,m,n,m1,m2,n1,n2,uvform,idummy
      complex U1(mmax),V1(nmax)
      real Sigma1(kmax),bnd1(kmax)
      complex U2(mmax),V2(nmax)
      complex udiff(mmax), vdiff(nmax),dotprod
      real Sigma2(kmax),bnd2(kmax)
      real norm2,norm1,normd,eps,tol,delta
      real abserr,relerr,maxabserr,maxrelerr
      integer relidx, absidx, info
c     Save large arrays in a common block to avoid stack overflow on
c     certain machines.
      common/cexample/U1,V1,Sigma1,bnd1,U2,V2,Sigma2,bnd2

      complex cdotc
      real scnrm2,slamch
      external scnrm2,slamch,cdotc

c-------------------- Here begins executable code ---------------------
      
c
c     Get input parameters from stdin.
c
      
      write (*,*) 'Dimensions m,n of original matrix?'
      read (5,*) m,n
      write (*,*) m,n
      write (*,*) 'Tolerance on computed singular values?'
      read (5,*) tol
      write (*,*) tol      
      write (*,*) 'Level of orthogonality among Lanczos vectors?'
      read (5,*) delta
      write (*,*) delta    
      write (*,*) 'File format for singular vectors' //
     c     ' (1=ASCII, 2=binary)?'
      read (5,*) uvform
      write (*,*) uvform
      write (*,*) 'Name of file with reference singular values?'
      read (5,'(a)') sigmafilename1
      write(*,*) sigmafilename1
      write (*,*) 'Name of file with reference U?'
      read (5,'(a)') ufilename1
      write(*,*) ufilename1
      write (*,*) 'Name of file with reference V?'
      read (5,'(a)') vfilename1
      write(*,*) vfilename1
      write (*,*) 'Name of file with test singular values?'
      read (5,'(a)') sigmafilename2
      write(*,*) sigmafilename2
      write (*,*) 'Name of file with test U?'
      read (5,'(a)') ufilename2
      write(*,*) ufilename2
      write (*,*) 'Name of file with test V?'
      read (5,'(a)') vfilename2
      write(*,*) vfilename2

c
c     Check singular values and computed error bounds.
c
      eps = slamch('e')

      open(10,file=sigmafilename1,status='old')
      do k1=1,kmax
         read(10,*,end=10) idummy,Sigma1(k1),bnd1(k1)
      enddo 
 10   close(10)
      open(10,file=sigmafilename2,status='old')
      do k2=1,kmax
         read(10,*,end=20) idummy,Sigma2(k2),bnd2(k2)
      enddo 
 20   close(10)
      k1 = k1-1
      k2 = k2-1
      if (k1.ne.k2) then
         print *,'*** Comparison failed, different number of' //
     c        'singular values:'
         print *,sigmafilename1,' contains ',k1,' singular value.'
         print *,sigmafilename2,' contains ',k2,' singular value.'
         goto 1000
      endif
      print *,'Analyzing singular values'
      print *,'Sigma1(',1,') = ',Sigma1(1)
      print *,'Sigma2(',1,') = ',Sigma2(1)
      print *,'Sigma1(',k1,') = ',Sigma1(k1)
      print *,'Sigma2(',k1,') = ',Sigma2(k1)

      maxabserr = abs(Sigma1(1)-Sigma2(1))
      absidx=1
      maxrelerr = maxabserr / abs(Sigma1(1))
      relidx=1
      call testsigma(tol,eps,1,m,n,Sigma1,Sigma2)
      do k=2,k1
         call testsigma(tol,eps,k,m,n,Sigma1,Sigma2)
         abserr = abs(Sigma1(k)-Sigma2(k))
         relerr = abserr / (abs(Sigma1(k))+abs(Sigma2(k)))
         if (relerr.gt.maxrelerr) then
            maxrelerr = relerr
            relidx = k
         endif
         if (abserr.gt.maxabserr) then 
            maxabserr = abserr
            absidx = k
         endif
      enddo
      print *,'Maximum relative error occured in component ',relidx
      print *,'s_ref = ',Sigma1(relidx),' s = ',Sigma2(relidx)
      print *,'|s_ref-s|/|s_ref| = ',maxrelerr
      print *,'Maximum absolute error occured in component ',absidx
      print *,'s_ref = ',Sigma1(absidx),' s = ',Sigma2(absidx)
      print *,'|s_ref-s| = ',maxabserr


      print *,''
      print *,'Analyzing left singular vectors'
      open(10,file=ufilename1,status='old')
      open(11,file=ufilename2,status='old')
      read(10,*) m1,k1
      read(11,*) m2,k2
      if ((k1.ne.k2) .or. (m1.ne.m2)) then
         print *,'*** Comparison failed, dimensions of U in' //
     c        'reference and test files differ.'
         print *,ufilename1,' contains ',k1,' vectors of length',m1
         print *,ufilename2,' contains ',k2,' vectors of length',m2
         goto 1000
      endif
      absidx=1
      do k = 1,k1
         read(10,*) (U1(i),i=1,m1)
         read(11,*) (U2(i),i=1,m1)
         norm1 = scnrm2(m1,U1,1)
         norm2 = scnrm2(m2,U2,1)
         dotprod = cdotc(m1,U1,1,U2,1)
         call testu(delta,eps,k,m,n,k1,Sigma1,norm1,norm2,abs(dotprod))
      enddo
      close(10)
      close(11)



      print *,''
      print *,'Analyzing right singular vectors'
      open(10,file=ufilename1,status='old')
      open(11,file=ufilename2,status='old')
      read(10,*) n1,k1
      read(11,*) n2,k2
      if ((k1.ne.k2) .or. (n1.ne.n2)) then
         print *,'*** Comparison failed, dimensions of V in' //
     c        'reference and test files differ.'
         print *,ufilename1,' contains ',k1,' vectors of length',n1
         print *,ufilename2,' contains ',k2,' vectors of length',n2
         goto 1000
      endif
      absidx=1
      do k = 1,k1
         read(10,*) (V1(i),i=1,n1)
         read(11,*) (V2(i),i=1,n1)
         norm1 = scnrm2(n1,U1,1)
         norm2 = scnrm2(n2,U2,1)
         dotprod = cdotc(n1,V1,1,V2,1)
         call testu(delta,eps,k,m,n,k1,Sigma1,norm1,norm2,abs(dotprod))
      enddo
      close(10)
      close(11)
 1000 end


      subroutine testu(delta,eps,i,m,n,k,sigma,norm_ref,norm,dotprod)
      implicit none
      integer m,n,k,i
      real eps,sigma(*),sep,norm_ref,norm,dotprod
      real theta,delta,errbnd

      if (abs(norm-1.0).gt.10*max(m,n)*eps) then
         print *,'Warning: U(:,',i,') is not normalized.' //
     c        ' ||u|| = ',norm
      endif

      if (i.gt.k .or. k.eq.1 .or. (i.eq.k .and. k.ne.min(m,n)))
     c     return
      theta = acos(min(1.0,dotprod/sqrt(norm_ref*norm)))
      if (i.eq.1) then
         sep = sigma(1)-sigma(2)
      else if (i.lt.k ) then        
         sep = min(sigma(i-1)-sigma(i),sigma(i)-sigma(i+1))
      else
         sep = min(sigma(i-1)-sigma(i),sigma(i))
      endif
      errbnd = 10*max(max(m,n)*eps*sigma(1)/sep,sqrt(real(k))*delta)
      if (theta.gt.errbnd) then
         print *,'Warning: U(:,',i,') differs significantly from' //
     c        ' reference: theta = ',theta,' > errbnd = ',
     c        errbnd 
      endif
      end



      subroutine testsigma(tol,eps,i,m,n,Sigma1,Sigma2)
      implicit none 
      integer i,m,n
      real eps,Sigma1(*),Sigma2(*),errbnd,err,tol

      err = abs(Sigma1(i)-Sigma2(i))
      errbnd = 10*max(m,n)*eps*Sigma1(1)
c      print *,'err = ',err,' errbnd = ',errbnd
      if (err.gt.errbnd) then
         print *,'Warning: Sigma(',i,') differs significantly from' //
     c        ' reference:'
         print *,'Sigma_ref(',i,') = ',Sigma1(1)
         print *,'Sigma(',i,') = ',Sigma2(1)
      endif
      end
