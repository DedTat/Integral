program lab1

    INTEGER NOFUN
      REAL A,B,RELERR,ABSERR,RESULT,ERREST,FLAG,x
      DATA A/0./,B/1.57/,RELERR/1.E-06/,ABSERR/0.0/
      CHARACTER(15)      :: form
      integer, parameter:: N = 11, nl = 11
      integer :: k, i, j, Out = 0
      real :: XArray(N+1), FArray(N+1), xk(N-1), LgrArray(N), SplArray(N), Dif(N)
      real :: C(n), D(n), E(n), SEVAL

      i = 1
      do x = 0.5, 1.5, 0.1
      CALL QUANC8(FUN,A,B,ABSERR,RELERR,RESULT,ERREST,NOFUN,FLAG)
      XArray(i) = x
      FArray(i) = RESULT
      i = i + 1
      PRINT 1,RESULT,ERREST,NOFUN,FLAG
      enddo

    1 FORMAT(10X,'RESULT=',E14.7,3X,'ERREST=',E12.5/11X,'NOFUN=',I8,11X,'FLAG=',F10.3)

      call SPLINE(N,XArray,FArray,C,D,E)

      do k = 1, 10, 1
        xk(k) = 0.55 + 0.1*k
        SplArray(k) = SEVAL(N,xk(k),XArray,FArray,C,D,E)
        LgrArray(k) = Lgr(xk(k))
        Dif(k) = SplArray(k) - LgrArray(k)
      enddo

      print *, ' '
      do i = 1, 10, 1
      print *, "Spline f(x): ", SplArray(i), "Lagrange f(x): ", LgrArray(i)
      end do

      print *, ' '
      do i = 1, 10, 1
      print *, "Dif between Spline and Lagrange : ", Dif(i)
      end do

      contains

      real FUNCTION FUN(Z)
      REAL Z
      FUN=(Z/(SIN(Z)**2 + x*COS(Z)**2))
      RETURN
      END function FUN

      real function Lgr(x)
      real x, numerator, denominator
      integer j,i
      i=0
      Lgr=0
      do i=0, N, 1
        numerator=1
        denominator=1
            do j=0, N, 1
                if(i /= j) then
                numerator = numerator * (x - XArray(j))
                denominator = denominator * (XArray(i) - XArray(j))
                end if
            end do
        Lgr=Lgr+FArray(i)*(numerator/denominator) !
      end do
      end function Lgr


end program lab1
