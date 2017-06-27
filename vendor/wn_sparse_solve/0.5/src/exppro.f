c-----------------------------------------------------------------------
      subroutine exppro (n, m, eps, tn, u, w, x, y, fpar,
     &                   ipar, indic, ierr)
      parameter (nf = 6, ni = 5)
      integer n, m, indic, ierr, ipar(ni)
      double precision eps, tn, u(n,m+1), w(n), x(n), y(n), fpar(nf)

c-----------------------------------------------------------------------
c     
c this subroutine computes an approximation to the vector
c
c        	w :=  exp( - A * tn ) * w
c
c where A is an arbitary matrix and w is a given input vector 
c uses a dynamic estimation of internal time advancement (dt) 
c----------------------------------------------------------------------- 
c THIS IS A REVERSE COMMUNICATION IMPLEMENTATION. 
c------------------------------------------------- 
c USAGE: (see also comments on indic below).
c------ 
c
c      indic = 0
c 1    continue
c      call exppro (n, m, eps, tn, u, w, x, y, indic)
c      if (indic .eq. 1) goto 2 <-- indic .eq. 1 means job is finished
c      call matvec(n, x, y)     <--- user's matrix-vec. product
c                                    with x = input vector, and
c                                     y = result = A * x.
c      goto 1
c 2    continue
c      .....
c
c-----------------------------------------------------------------------
c
c on entry:
c---------- 
c n	= dimension of matrix
c
c m	= dimension of Krylov subspace (= degree of polynomial 
c         approximation to the exponential used. )
c
c eps   = scalar indicating the relative error tolerated for the result. 
c         the code will try to compute an answer such that 
c         norm2(exactanswer-approximation) / norm2(w) .le. eps 
c
c tn	= scalar by which to multiply matrix. (may be .lt. 0)
c         the code will compute an approximation to exp(- tn * A) w
c         and overwrite the result onto w.
c
c u	= work array of size n*(m+1) (used to hold the Arnoldi basis )
c
c w	= real array of length n = input vector to  which exp(-A) is
c         to be applied. this is also an output argument 
c
c x, y  = two real work vectors of length at least  n each. 
c         see indic for usage.
c
c fpar  = formerly saved variables moved to array (BSM 2013-12-17):
c
c           fpar(1) = beta
c           fpar(2) = fnorm
c           fpar(3) = dtl
c           fpar(4) = tcur
c           fpar(5) = told
c           fpar(6) = errst
c   
c ipar  = formerly saved variables moved to array (BSM 2013-12-17):
c
c           ipar(1) = i1
c           ipar(2) = i
c           ipar(3) = job
c           ipar(4) = idebug
c           ipar(5) = m1
c
c indic = integer used as indicator for the reverse communication.
c         in the first call enter indic = 0. See below for more.
c
c on return:
c-----------
c 
c w     = contains the resulting vector exp(-A * tn ) * w when 
c         exppro has finished (see indic) 
c
c indic = indicator for the reverse communication protocole.
c       * INDIC .eq. 1  means that exppro has finished and w contains the
c         result. 
c       * INDIC .gt. 1 ,  means that exppro has not finished and that
c         it is requesting another matrix vector product before
c         continuing. The user must compute Ax where A is the matrix
c         and x is the vector provided by exppro, and return the 
c         result in y. Then exppro must be called again without
c         changing any other argument. typically this must be
c         implemented in a loop with exppro being called as long
c         indic is returned with a value .ne. 1.
c
c ierr  = error indicator. 
c         ierr = 1 means phipro was called with indic=1 (not allowed)
c         ierr = -1 means that the input is zero the solution has been 
c         unchanged.
c 
c NOTES:  m should not exceed 60 in this version  (see mmax below)
c-----------------------------------------------------------------------
c written by Y. Saad -- version feb, 1991.
c B. Meyer moved variables to in/out 2013-12-16 to make routines thread safe.
c----------------------------------------------------------------------- 
c For reference see following papers : 
c (1) E. Gallopoulos and Y. Saad: Efficient solution of parabolic 
c     equations by Krylov approximation methods. RIACS technical
c     report 90-14.
c (2) Y.Saad: Analysis of some Krylov subspace approximations to the
c     matrix exponential operator. RIACS Tech report. 90-14 
c-----------------------------------------------------------------------
c local variables 
c 
      integer mmax
      parameter (mmax=100) 
      double precision red, dabs, dble
      double precision hh(mmax+2,mmax+1), z(mmax+1)
      complex(8)   wkc(mmax+1) 
      integer ih
      logical verboz
c-----------------------------------------------------------------------
c indic = 3  means  passing through only with result of y= Ax to exphes
c indic = 2  means exphes has finished its job
c indic = 1  means exppro has finished its job (real end)/
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     Debugging
c-----------------------------------------------------------------------
      if (ipar(4) .eq. 1 ) then
        verboz = .true.
      else
        verboz = .false.
      endif
c-----------------------------------------------------------------------

      ierr = 0 
      ih = mmax 

      if (indic .eq. 3) goto 101
      if (indic .eq. 1) then 
         ierr = 1
         return
      endif
c----- 
      m  = min0(m,mmax) 
      fpar(4) = 0.0d0
      fpar(3) = tn-fpar(4)
      ipar(3) = -1
c-------------------- outer loop ----------------------------- 
 100  continue
      if(verboz) print *,'In EXPPRO, current time = ', tcur ,'---------'
c------------------------------------------------------------- 
c ---- call exponential propagator --------------------------- 
c------------------------------------------------------------- 
      fpar(5) = fpar(4)
 101  continue
c     if (told + dtl .gt. tn) dtl = tn-told
      call  exphes (n,m,eps,u,w,z,wkc,fpar,ipar,hh,ih,
     *              x,y,indic,ierr) 
c-----------------------------------------------------------------------
      if (ierr .ne. 0) return
      if (indic .ge. 3) return
      fpar(4) = fpar(5) + fpar(3) 
      if(verboz) print *, ' tcur now = ', fpar(4), ' dtl = ', fpar(3)
c
c     relative error 
c      if(verboz) print *, ' beta', fpar(1)
      fpar(6) = fpar(6) / fpar(1)
c---------
       if (fpar(4) .eq. tn) goto 102
cc 
c      if ((errst .le. eps) .and. ( (errst .gt. eps/100.0) .or.
c     *     (fpar(4) .eq. tn))) goto 102
cc
c     
c     use approximation :  [ new err ] = fact**m  * [cur. error]
c     
      red =  (0.5*eps / fpar(6))**(1.0d0 /dble(m) ) 
      fpar(3) = fpar(3)*red 
      if (dabs(fpar(5)+fpar(3)) .gt. dabs(tn) )  fpar(3) = tn-fpar(5)
      if(verboz) print *, ' red =',red,' , reducing dt to: ', fpar(3) 
c-------
      ipar(3) = 1 
      goto 101
c-------
 102  continue 
c
      call exppro_project(n,m,u,z,w)
c never go beyond tcur
      ipar(3) = 0
      fpar(3) = dmin1(fpar(3), tn - fpar(4))
      if (dabs(fpar(4)+fpar(3)) .gt. dabs(tn)) fpar(3) = tn-fpar(4) 
      if (dabs(fpar(4)) .lt. dabs(tn)) goto 100
      indic = 1
c
      return
      end
c----------end-of-expro-------------------------------------------------
c-----------------------------------------------------------------------
      subroutine exphes (n,m,eps,u,w,z,wkc,fpar,ipar,
     &                   hh,ih,x, y, indic,ierr) 
      parameter (nf = 6, ni = 5)
      integer n, m, ih, indic, ierr, ipar(ni)
      double precision hh(ih+2,m+1), u(n,m+1), w(n), z(m+1), x(n), y(n)
      complex(8) wkc(m+1) 
      double precision eps, fpar(nf)
c-----------------------------------------------------------------------
c this subroutine computes the Arnoldi basis and the corresponding 
c coeffcient vector in the approximation 
c 
c        	w  ::= beta  Vm  ym 
c               where ym = exp(- Hm *dt) * e1
c
c to the vector exp(-A dt) w where A is an arbitary matrix and 
c w is a given input vector. In case job = 0 the arnoldi basis 
c is recomputed. Otherwise the
c code assumes assumes that  u(*) contains an already computed 
c arnoldi basis and computes only the y-vector (which is stored in v(*))
c-----------------------------------------------------------------------
c on entry:
c---------- 
c n	= dimension of matrix
c
c m	= dimension of Krylov subspace (= degree of polynomial 
c         approximation to the exponential used. )
c
c dt	= scalar by which to multiply matrix. Can be viewed
c         as a time step. dt must be positive [to be fixed].
c
c eps   = scalar indicating the relative error tolerated for the result. 
c         the code will try to compute an answer such that 
c         norm2(exactanswer-approximation) / norm2(w) .le. eps 
c
c u	= work array of size n*(m+1) to contain the Arnoldi basis
c
c w	= real array of length n = input vector to  which exp(-A) is
c         to be applied. 
c
c job	= integer (ipar(3)). job indicator. If job .lt.  0 then the Arnoldi
c         basis is recomputed. If job .gt. 0 then it is assumed
c         that the user wants to use a previously computed Krylov
c         subspace but a different dt. Thus the Arnoldi basis and
c         the Hessenberg matrix Hm are not recomputed. 
c	  In that case the user should not modify the values of beta
c         and the matrices hh and u(n,*) when recalling phipro. 
c         job = -1 : recompute basis and get an initial estimate for 
c                    time step dt to be used.
c         job = 0  : recompute basis and do not alter dt.
c         job = 1  : do not recompute arnoldi basis. 
c
c z     = real work array of  size (m+1)
c wkc   = double complex work array of size (m+1) 
c
c fpar, ipar = parameter arrays (see exppro--BSM 2013-12-17)
c hh    = work array of size size at least (m+2)*(m+1)
c
c ih+2	= first dimension of hh as declared in the calling program.
c         ih must be .ge. m.
c 
c-----------------------------------------------------------------------
c on return:
c-----------
c w2	= resulting vector w2 = exp(-A *dt) * w
c beta  = real equal to the 2-norm of w. Needed if exppro will
c         be recalled with the same Krylov subspace and a different dt. 
c errst = rough estimates of the 2-norm of the error.
c hh	= work array of dimension at least (m+2) x (m+1)
c
c----------------------------------------------------------------------- 
c local variables 
c 
      integer ndmax
      parameter (ndmax=20) 
      double precision alp0, t, rm, ddot, dabs, dsqrt, dsign,dble
      complex(8) alp(ndmax+1), rd(ndmax+1)
      integer j, k, ldg, i0
      logical verboz
c------debug information---------------------------------------------
      if (ipar(4) .eq. 1 ) then
        verboz = .true.
      else
        verboz = .false.
      endif
c------use degree 14 chebyshev all the time -------------------------- 
c
c------input fraction expansion of rational function ------------------ 
c chebyshev (14,14) 
      ldg= 7
      alp0 =  0.183216998528140087E-11
      alp(1)=( 0.557503973136501826E+02,-0.204295038779771857E+03)
      rd(1)=(-0.562314417475317895E+01, 0.119406921611247440E+01)
      alp(2)=(-0.938666838877006739E+02, 0.912874896775456363E+02)
      rd(2)=(-0.508934679728216110E+01, 0.358882439228376881E+01)
      alp(3)=( 0.469965415550370835E+02,-0.116167609985818103E+02)
      rd(3)=(-0.399337136365302569E+01, 0.600483209099604664E+01)
      alp(4)=(-0.961424200626061065E+01,-0.264195613880262669E+01)
      rd(4)=(-0.226978543095856366E+01, 0.846173881758693369E+01)
      alp(5)=( 0.752722063978321642E+00, 0.670367365566377770E+00)
      rd(5)=( 0.208756929753827868E+00, 0.109912615662209418E+02)
      alp(6)=(-0.188781253158648576E-01,-0.343696176445802414E-01)
      rd(6)=( 0.370327340957595652E+01, 0.136563731924991884E+02)
      alp(7)=( 0.143086431411801849E-03, 0.287221133228814096E-03)
      rd(7)=( 0.889777151877331107E+01, 0.166309842834712071E+02)
c-----------------------------------------------------------------------

      if (indic .ge. 3) goto 60
c     
c     if ipar(3) .gt. 0 skip arnoldi process:
c     
      if (ipar(3) .gt. 0) goto 2
c------normalize vector w and put in first column of u -- 
      fpar(1) = dsqrt(ddot(n,w,1,w,1)) 
c-----------------------------------------------------------------------
      if(verboz) print *, ' In EXPHES, beta ', fpar(1)
      if (fpar(1) .eq. 0.0d0) then
         ierr = -1 
         indic = 1
         return
      endif
c
      t = 1.0d0/fpar(1)
      do 25 j=1, n
         u(j,1) = w(j)*t 
 25   continue
c------------------Arnoldi loop ------------------------------------- 
      ipar(1) = 1
 58   ipar(2) = ipar(1)
      ipar(1) = ipar(2) + 1
      do 59 k=1, n
         x(k) = u(k,ipar(2))
 59   continue
      indic = 3
      return
 60   continue
      do 61 k=1, n
         u(k,ipar(1)) = y(k) 
 61   continue
      i0 =1
c     
c switch  for Lanczos version 
c     i0 = max0(1, ipar(2)-1)
      call exppro_mgsr (n, i0, ipar(1), u, hh(1,ipar(2)))
      fpar(2) =
     &      fpar(2) + ddot(ipar(1), hh(1,ipar(2)),1, hh(1,ipar(2)),1)
      if (hh(ipar(1),ipar(2)) .eq. 0.0) m = ipar(2)
      if  (ipar(2) .lt. m) goto 58
c--------------done with arnoldi loop ---------------------------------
      rm = dble(m) 
      fpar(2) = dsqrt( fpar(2) / rm )
c-------get  : beta*e1 into z 
      ipar(5) = m+1 
      do 4 ii=1,ipar(5)
         hh(ii,ipar(5)) = 0.0
 4    continue
c
c     compute initial dt when  job .lt. 1 
c
      if (ipar(3) .ge. 0) goto 2 
c
c     t = eps / beta 
c     
      t = eps 
      do 41 k=1, m-1
         t = t*(1.0d0 - dble(m-k)/rm ) 
 41   continue
c
      t = 2.0d0*rm* (t**(1.0d0/rm) )  / fpar(2)
      if(verboz) print *, ' t, dt = ', t, fpar(3), fpar(2)
c remove fpar(2)*/
      t = dmin1(dabs(fpar(3)),t) 
      fpar(3) = dsign(t, fpar(3)) 
c     
 2    continue 
      z(1) = fpar(1)
      do 3 k=2, ipar(5)
         z(k) = 0.0d0 
 3    continue
c-------get  : exp(H) * beta*e1
      call exppro_hes(ldg,ipar(5),hh,ih,fpar(3),z,rd,alp,alp0,wkc)
c-------error estimate 
      fpar(6) = dabs(z(ipar(5)))
      if(verboz) print *, ' error estimate =', fpar(6)
c-----------------------------------------------------------------------
      indic = 2
      return
      end
c-----------------------------------------------------------------------
      subroutine exppro_mgsr (n, i0, i1, ss, r)
      integer n, i0, i1
      double precision ss(n,i1), r(i1)
c-----------------------------------------------------------------------
c modified gram - schmidt  with  partial  reortho. the vector ss(*,i1) is
c orthogonalized against the first i vectors  of ss  (which  are  already
c orthogonal).  the coefficients of the orthogonalization are returned in
c the array r
c------------------------------------------------------------------------
c local variables 
c 
      integer i, j, k, it
      double precision hinorm, tet, ddot, t, dsqrt
      data  tet/10.0d0/

      do 53 j=1, i1
         r(j) = 0.0d0
 53   continue
      i = i1-1
      it = 0
 54   hinorm = 0.0d0
      it = it +1
      if (i .eq. 0) goto 56
c     
      do 55 j=i0, i
         t = ddot(n, ss(1,j),1,ss(1,i1),1)
         hinorm = hinorm + t**2
         r(j) = r(j) + t
         call daxpy(n,-t,ss(1,j),1,ss(1,i1),1)
 55   continue
      t = ddot(n, ss(1,i1), 1, ss(1,i1), 1)
 56   continue
c     
c     test for reorthogonalization see daniel et. al.
c     two reorthogonalization allowed ---
c     
      if (t*tet .le. hinorm .and. it .lt. 2) goto 54
      t =dsqrt(t)
      r(i1)= t
      if (t .eq. 0.0d0) return
      t = 1.0d0/t
      do 57  k=1,n
         ss(k,i1) = ss(k,i1)*t
 57   continue
      return
      end
c----------end-of-exppro_mgsr--------------------------------------------------
c-----------------------------------------------------------------------
      subroutine exppro_project(n,m,u,v,w)
      integer n, m
      double precision u(n,m), v(m), w(n)
c
c     computes the vector w = u * v
c
c local variables 
c
      integer j, k

      do 1 k=1,n
         w(k) = 0.d0
 1    continue

      do 100 j=1,m
         do 99 k=1,n
            w(k) = w(k) + v(j) * u(k,j) 
 99      continue
 100  continue
      return
      end
c-----------------------------------------------------------------------     
      subroutine exppro_hes (ndg,m1,hh,ih,dt,y,root,coef,coef0,w2)
      integer ndg, m1, ih 
      double precision hh(ih+2,m1), y(m1)
      complex(8) coef(ndg), root(ndg), w2(m1)
      double precision dt, coef0
c--------------------------------------------------------------------  
c computes exp ( H dt) * y    (1)
c where H = Hessenberg matrix (hh)
c y	  = arbitrary vector.
c ----------------------------
c ndg	= number of poles as determined by getrat
c m1    = dimension of hessenberg matrix
c hh	= hessenberg matrix (real)
c ih+2	= first dimension of hh
c dt	= scaling factor used for hh (see (1)) 
c y	= real vector. on return exp(H dt ) y is computed
c         and overwritten on y.
c root  = poles of the rational approximation to exp as
c         computed by getrat
c coef, 
c coef0 = coefficients of partial fraction expansion 
c         
c  exp(t) ~ coef0 +  sum     Real [   coef(i) / (t - root(i)  ]
c                  i=1,ndg  
c
c valid for real t.
c coef0 is real, coef(*) is a double complex array.
c
c--------------------------------------------------------------------  
c local variables 
c
      integer m1max
      parameter (m1max=61) 
      complex(8) hloc(m1max+1,m1max), t, zpiv, dcmplx
      double precision yloc(m1max), dble
      integer i, j, ii
c     
c      if (m1 .gt. m1max) print *, ' *** ERROR : In HES, M+1 TOO LARGE'
c     
c     loop associated with the poles.
c     
      do 10 j=1,m1
         yloc(j) = y(j)
         y(j)    = y(j)*coef0
 10   continue
c     
      do 8 ii = 1, ndg
c     
c     copy Hessenberg matrix into temporary
c     
         do 2 j=1, m1
            do 1 i=1, j+1
               hloc(i,j) = dcmplx( dt*hh(i,j) )
 1          continue
            hloc(j,j) = hloc(j,j) - root(ii) 
            w2(j)     = dcmplx(yloc(j)) 
 2       continue 
c
c forward solve 
c 
         do 4 i=2,m1
            zpiv  = hloc(i,i-1) / hloc(i-1,i-1)
            do 3 j=i,m1
               hloc(i,j) = hloc(i,j) - zpiv*hloc(i-1,j)
 3          continue
            w2(i)     = w2(i) - zpiv*w2(i-1)
 4       continue 
c     
c     backward solve
c     
         do 6 i=m1,1,-1 
            t=w2(i)
            do 5 j=i+1,m1
               t = t-hloc(i,j)*w2(j)
 5          continue
            w2(i) = t/hloc(i,i)
 6       continue
c     
c     accumulate result in y.
c     
         do 7 i=1,m1
            y(i) = y(i) + dble ( coef(ii) * w2(i) ) 
 7       continue
 8    continue
      return
      end
c----------end-of-hes---------------------------------------------------


