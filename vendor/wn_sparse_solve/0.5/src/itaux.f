c-----------------------------------------------------------------------
      function distdot(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 distdot, x(*), y(*), ddot
      external ddot
      distdot = ddot(n,x,ix,y,iy)
      return
      end
c-----end-of-distdot
c-----------------------------------------------------------------------
c
      function afun (x,y,z)
      real*8 afun, x,y, z 
      afun = -1.0D0
      return 
      end
      
      function bfun (x,y,z)
      real*8 bfun, x,y, z 
      bfun = -1.0D0
      return 
      end
      
      function cfun (x,y,z)
      real*8 cfun, x,y, z 
      cfun = -1.0D0
      return 
      end
      
      function dfun (x,y,z)
      real*8 dfun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha
      dfun = gammax*exp(x*y)
      return 
      end
      
      function efun (x,y,z)
      real*8 efun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha
      efun = gammay*exp(-x*y) 
      return 
      end
      
      function ffun (x,y,z)
      real*8 ffun, x,y, z 
      ffun = 0.0D0
      return 
      end
      
      function gfun (x,y,z)
      real*8 gfun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha
      gfun = alpha 
      return 
      end
      
      function hfun (x,y,z)
      real*8 hfun, x,y, z, gammax, gammay, alpha
      common /func/ gammax, gammay, alpha
      hfun = alpha * sin(gammax*x+gammay*y-z)
      return 
      end
      

      function betfun(side, x, y, z)
      real*8 betfun, x, y, z
      character*2 side
      betfun = 1.0
      return
      end

      function gamfun(side, x, y, z)
      real*8 gamfun, x, y, z
      character*2 side
      if (side.eq.'x2') then
         gamfun = 5.0
      else if (side.eq.'y1') then
         gamfun = 2.0
      else if (side.eq.'y2') then
         gamfun = 7.0
      else
         gamfun = 0.0
      endif
      return
      end
c-----------------------------------------------------------------------
c     functions for the block PDE's 
c-----------------------------------------------------------------------
      subroutine afunbl (nfree,x,y,z,coeff)
      return
      end
c     
      subroutine bfunbl (nfree,x,y,z,coeff)
      return 
      end
      
      subroutine cfunbl (nfree,x,y,z,coeff)
c     
      return 
      end
      
      subroutine dfunbl (nfree,x,y,z,coeff)
      
      return
      end
c     
      subroutine efunbl (nfree,x,y,z,coeff)
      return 
      end
c     
      subroutine ffunbl (nfree,x,y,z,coeff)
      return 
      end
c     
      subroutine gfunbl (nfree,x,y,z,coeff)
      return 
      end
      
