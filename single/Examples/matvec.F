c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c
      subroutine readmatrix(filename,m,n)
      implicit none
      character*256 filename
      integer m,n
      include 'matvec.h'

      if (matform.eq.1) then
         call readhb(filename,m,n)
      else if (matform.eq.2 .or. matform.eq.3) then
         call readdense(filename,matform,m,n,matrix,densemmax)
      else if (matform.eq.4 .or. matform.eq.5) then
         call readcoord(filename,matform,m,n)
      else if (matform.eq.6 .or. matform.eq.7) then
         call readdiag(filename,matform,m,n)
      endif
      
      end


      subroutine matvecmult(transa,m,n,x,y,dparm,iparm)
      implicit none
      
      character*1 transa
      integer m,n,iparm(*)
      real x(*),y(*),dparm(*)
      include 'matvec.h'

      if (matform.eq.1) then
         call hbmv(transa,m,n,x,y,dparm,iparm)
      else if (matform.eq.2 .or. matform.eq.3) then
         call densemv(transa,m,n,x,y,dparm,iparm)
      else if (matform.eq.4 .or. matform.eq.5) then
         call coordmv(transa,m,n,x,y,dparm,iparm)
      else if (matform.eq.6 .or. matform.eq.7) then      
         call diagmv(transa,m,n,x,y,dparm,iparm)
      endif
      end

ccccccccccccccccccccc Harwell-Boeing format cccccccccccccccccccccccccccc
c Note: this is also known as compressed column storage (CCS) format.
 
      subroutine readhb(filename,m,n)
      implicit none
c
c     Read matrix in Harwell-Boeing format
c
      include 'matvec.h'
      character*256 filename
       
      integer lunit,m,n
      parameter(lunit=10)
      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 , RHSTYP*3,
     1               PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     1               NELTVL, I,  NRHS  , NRHSIX
      logical lsame
      external lsame

      
      open(lunit,file=filename,status='old')

C    ------------------------
C     ... READ IN HEADER BLOCK
C     ------------------------
      read ( lunit, 1000 ) title , key   ,
     1                     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     2                     mxtype, nrow  , ncol  , nnzero, neltvl,
     3                     ptrfmt, indfmt, valfmt, rhsfmt
      if  ( rhscrd .gt. 0 )
     1    read ( lunit, 1001 ) rhstyp, nrhs, nrhsix
 1000 format ( a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20 )
 1001 format ( a3, 11x, 2i14 )

      write (*,1000) title , key   ,
     1                     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     2                     mxtype, nrow  , ncol  , nnzero, neltvl,
     3                     ptrfmt, indfmt, valfmt, rhsfmt
c     -------------------------
C     ... READ MATRIX STRUCTURE
C     -------------------------

      if (ncol.gt.nmax) stop 'ERROR in readHB: ncol > nmax'
      if (nrow.gt.mmax) stop 'ERROR in readHB: nrow > mmax'
      if (nnzero.gt.nnzmax) stop 'ERROR in readHB: nnzero > nnzmax'
      

      read ( lunit, ptrfmt ) ( colptr (i), i = 1, ncol+1 )
      read ( lunit, indfmt ) ( rowind (i), i = 1, nnzero )

      if  ( valcrd .gt. 0 )  then

c         ----------------------
c         ... read matrix values
c         ----------------------

          read ( lunit, valfmt ) ( values (i), i = 1, nnzero )

      endif

      m = nrow
      n = ncol
      close(lunit)
      end


      subroutine hbmv(transa,m,n,x,y,dparm,iparm)
      implicit none
      character*1 transa
      integer m,n,i,j,iparm(*)
      real x(*),y(*),dparm(*),xi,sum
#ifdef _OPENMP
      real ylocal(m)
      integer omp_get_max_threads
      external omp_get_max_threads
#endif
      logical lsame
      external lsame
      include 'matvec.h'
      
      if (lsame(transa,'n')) then
c
c     computes  y <- A*x
c
#ifdef _OPENMP
        if ( omp_get_max_threads().gt.1 ) then
c$OMP PARALLEL private(ylocal) private(i,j,xi) shared(y,x,/hb/)
CDIR$ LOOP COUNT(10000)
c$OMP DO
           do i=1,m
              y(i) = 0.0
           enddo
c$OMP ENDDO NOWAIT
CDIR$ LOOP COUNT(10000)
           do i=1,m
              ylocal(i) = 0.0
           enddo
CDIR$ LOOP COUNT(10000)
c$OMP DO 
           do i=1,n
              xi = x(i)
cDIR$ IVDEP
              do j=colptr(i),colptr(i+1)-1
                 ylocal(rowind(j)) = ylocal(rowind(j)) + values(j)*xi
              enddo
           enddo
c$OMP ENDDO
CDIR$ LOOP COUNT(10000)
c$OMP CRITICAL
           do i=1,m
              y(i) = y(i) + ylocal(i)
           enddo
c$OMP END CRITICAL
c$OMP END PARALLEL
        else
#endif
CDIR$ LOOP COUNT(10000)
           do i=1,m
              y(i) = 0.0
           enddo
CDIR$ LOOP COUNT(10000)
           do i=1,n
              xi = x(i)
CDIR$ IVDEP
              do j=colptr(i),colptr(i+1)-1
                 y(rowind(j)) = y(rowind(j)) + values(j)*xi
              enddo
           enddo
#ifdef _OPENMP
        endif
#endif
      else
c
c     computes  y <- A'*x
c
CDIR$ LOOP COUNT(10000)
c$OMP PARALLEL DO private(j,sum)
         do i=1,n
            sum = 0.0
CDIR$ IVDEP
            do j=colptr(i),colptr(i+1)-1
               sum = sum + values(j)*x(rowind(j))
            enddo
            y(i) = sum
         enddo
      endif
      end

ccccccccccccccccccccc Coordinate format cccccccccccccccccccccccccccc
c
c     For the purpose of running PROPACK, it would often be faster first
c     to convert coordinate format matrices into e.g. compressed column 
c     storage format (also known as the Harwell-Boeing format) and 
c     use that for computing matrix-vector multiplies, but we supply
c     the following routines for educational purposes.
c

      subroutine readcoord(filename, form, m, n)
      implicit none
      character*256 filename
      integer lunit,form,m,n,i
      parameter(lunit=10)
      include 'matvec.h'

      if (form.eq.4) then
         open(lunit,file=filename,status='old')
         read(lunit,*) m,n,nnzero
         do i=1,nnzero
            read(lunit,*) coordinate(i,1),coordinate(i,2),values(i)
         enddo
      else if (form.eq.5) then
         open(lunit,file=filename,status='old',form='unformatted')
         read(lunit) m,n,nnzero
         do i=1,nnzero
            read(lunit) coordinate(i,1),coordinate(i,2),values(i)
         enddo
      endif
      close(lunit)
      end

      subroutine coordmv(transa,m,n,x,y,dparm,iparm)
      implicit none
      character*1 transa
      integer m,n,i,j,k,iparm(*)
      real x(*),y(*),dparm(*),aij
#ifdef _OPENMP
      real ylocal(m+n)
#endif
      logical lsame
      external lsame
      include 'matvec.h'
      
#ifdef _OPENMP
      if (lsame(transa,'n')) then
c$OMP PARALLEL private(ylocal,i,j,aij)
c$OMP DO 
         do i=1,m
            y(i) = 0.0
         enddo
c$OMP ENDDO NOWAIT
         do i=1,m
            ylocal(i) = 0.0
         enddo
c$OMP DO 
         do k=1,nnzero
            i = coordinate(k,1)
            j = coordinate(k,2)
            aij = values(k)
            ylocal(i) = ylocal(i) + aij*x(j)
         enddo
c$OMP critical
         do i=1,m
            y(i) = y(i) + ylocal(i)
         enddo
c$OMP end critical
c$OMP END PARALLEL
      else
c$OMP PARALLEL private(ylocal,i,j,aij)
c$OMP DO 
         do i=1,m
            y(i) = 0.0
         enddo
c$OMP ENDDO NOWAIT
         do i=1,m
            ylocal(i) = 0.0
         enddo
c$OMP DO 
         do k=1,nnzero
            i = coordinate(k,1)
            j = coordinate(k,2)
            aij = values(k)
            ylocal(j) = ylocal(j) + aij*x(i)
         enddo
c$OMP critical
         do i=1,n
            y(i) = y(i) + ylocal(i)
         enddo
c$OMP end critical
c$OMP END PARALLEL
      endif
#else
      if (lsame(transa,'n')) then
         do i=1,m
            y(i) = 0.0
         enddo
         do k=1,nnzero
            i = coordinate(k,1)
            j = coordinate(k,2)
            aij = values(k)
            y(i) = y(i) + aij*x(j)
         enddo
      else
         do i=1,m
            y(i) = 0.0
         enddo
         do k=1,nnzero
            i = coordinate(k,1)
            j = coordinate(k,2)
            aij = values(k)
            y(j) = y(j) + aij*x(i)
         enddo
      endif
#endif
      end


ccccccccccccccccccccc dense format cccccccccccccccccccccccccccc
c
c     Note: We assume the matrix is always stored in major 
c     column order.
c
      subroutine readdense(filename, form, m, n, A, lda)
      implicit none
      character*256 filename
      integer lunit,form,m,n,i,j,lda
      real A(lda,*)
      parameter(lunit=10)
      include 'matvec.h'

      if (form.eq.2) then
         open(lunit,file=filename,status='old')
         read(lunit,*) m,n
         read(lunit,*) ((A(i,j),i=1,m),j=1,n) 
      else if (form.eq.3) then
         open(lunit,file=filename,status='old',form='unformatted')
         read(lunit) m,n
         read(lunit) ((A(i,j),i=1,m),j=1,n)         
      endif
      close(lunit)      
      end

      subroutine writedense(filename, form, m, n, A, lda)
      implicit none
      character*256 filename
      integer lunit,form,m,n,i,j,lda
      real A(lda,*)
      parameter(lunit=10)
      include 'matvec.h'

      if (form.eq.2) then
         open(lunit,file=filename,status='unknown')
         write(lunit,*) m,n
         do j=1,n
            write(lunit,'(e23.15)') (A(i,j),i=1,m)
         enddo
      else if (form.eq.3) then
         open(lunit,file=filename,status='unknown',form='unformatted')
         write(lunit) m,n
         write(lunit) ((A(i,j),i=1,m),j=1,n)         
      endif
      close(lunit)      
      end


      subroutine densemv(transa,m,n,x,y,dparm,iparm)
      implicit none
      character*1 transa
      integer m,n,iparm(*)
      real x(*),y(*),dparm(*)
      include 'matvec.h'
     
      call sgemv(transa,m,n,1.0,matrix,densemmax,x,1,0.0,y,1)

      end

ccccccccccccccccccccc diagonal format cccccccccccccccccccccccccccc

      subroutine readdiag(filename,form, m, n)
      implicit none
      character*256 filename
      integer lunit,m,n,i,form
      parameter(lunit=10)
      include 'matvec.h'

      if (form.eq.6) then
         open(lunit,file=filename,status='old')
         read(lunit,*) m,n
         read(lunit,*) (values(i),i=1,min(m,n))
      else if (form.eq.7) then
         open(lunit,file=filename,status='old',form='unformatted')
         read(lunit) m,n
         read(lunit) (values(i),i=1,min(m,n))
      endif
      
      end


      subroutine diagmv(transa,m,n,x,y,dparm,iparm)
      implicit none
      character*1 transa
      integer m,n,i,iparm(*)
      real x(*),y(*),dparm(*)
      logical lsame
      external lsame
      include 'matvec.h'
      
c$OMP PARALLEL DO 
      do i=1,min(m,n)
         y(i) = values(i)*x(i)
      enddo
      if (m.ne.n) then
         if (lsame(transa,'n')) then
c$OMP PARALLEL DO 
            do i=min(m,n)+1,m
               y(i) = 0.0
            enddo
         else
c$OMP PARALLEL DO 
            do i=min(m,n)+1,n
               y(i) = 0.0
            enddo
         endif
      endif
      end

