c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

c 
c     DGETU0: Attempt to generate a pseudo-random vector in SPAN(Op(A)) 
c     orthogonal to span(U(:,1:j)), where Op(A) = A if transa='n' and
c     Op(A) = A^T if transa='t'.
c

      subroutine dgetu0(transa, m, n, j, ntry, u0, u0norm, U, ldu,
     c     aprod, dparm, iparm, ierr, icgs, anormest, work)

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 transa
      integer m, n, j, ntry, ldu, ierr,icgs
      integer iparm(*)
      double precision u0norm,anormest
      double precision u0(*),U(*),work(*),dparm(*)
      external aprod

c     %------------%
c     | Parameters |
c     %------------%
      double precision kappa
      parameter(kappa = 0.717)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer itry,idist,iseed(4),rsize,usize,index(3)
      double precision nrm

c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      double precision pdnrm2,rand
      external pdnrm2,lsame,rand

c-------------------- Here begins executable code ---------------------

      iseed(1) = 1
      iseed(2) = 3
      iseed(3) = 5
      iseed(4) = 7

      if (lsame(transa,'n')) then
c     %-------------------------%
c     | u0 is to be an m-vector |
c     %-------------------------%
         rsize = n
         usize = m
      else
c     %-------------------------%
c     | u0 is to be an n-vector |
c     %-------------------------%
         rsize = m
         usize = n
      endif

      idist = 2
      ierr = 0
      do itry=1,ntry
         call dlarnv(idist, iseed, rsize, work)
         nrm = pdnrm2(rsize,work,1)
         call aprod(transa,m,n,work,u0,dparm,iparm)
         nopx = nopx+1

         u0norm = pdnrm2(usize,u0,1)
         anormest = u0norm/nrm
         
         if (j.ge.1) then
            index(1) = 1
            index(2) = j
            index(3) = j+1            
            call dreorth(usize,j,U,ldu,u0,u0norm,index,kappa,
     c           work,icgs)
         endif

         if (u0norm.gt.0) goto 9999
      enddo
      ierr = -1
 9999 return
      end
      
