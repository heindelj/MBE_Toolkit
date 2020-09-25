module dfpmin_mod

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dfpmin(n,p,gtol,tolx, stpmx,iter,fret, iprint, dfunc)
!.......................................................................
      implicit none
      integer,                        intent( in    ) :: n
      double precision, dimension(n), intent( inout ) :: p(n)
      double precision,               intent( in    ) :: gtol, tolx
      double precision,               intent( in    ) :: stpmx
      integer,                        intent(   out ) :: iter
      double precision,               intent(   out ) :: fret
      integer,                        intent( in    ) :: iprint
      interface
         subroutine dfunc(n, p, g, fret)
         implicit none
         integer,                        intent( in    ) :: n
         double precision, dimension(n), intent( in    ) :: p
         double precision, dimension(n), intent(   out ) :: g
         double precision,               intent(   out ) :: fret
         end subroutine dfunc
      end interface
      !...
      integer, parameter :: ITMAX=999999
      INTEGER i,j,its
      LOGICAL check
      double precision, dimension(:,:), allocatable :: hessin
      double precision, dimension(:), allocatable :: g, dg, hdg, pnew, xi
      double precision :: den,fac,fad,fae,fp,stpmax,sum0,sumdg,sumxi,temp,test
      double precision :: tmp, testx, testg,EPS

      allocate(hessin(n,n))
      allocate(g(n))
      allocate(dg(n))
      allocate(hdg(n))
      allocate(pnew(n))
      allocate(xi(n))
      EPS=TOLX/4.d0
      call dfunc(n, p, g, fp)
      sum0=0.d0 
      hessin(1:n, 1:n) = 0.d0
      do i=1,n
         hessin(i,i)=1.D0
         xi(i)=-g(i)
         sum0=sum0+p(i)**2
      enddo
      stpmax=STPMX*dmax1(dsqrt(sum0),dble(n))
      do 27 its=1,ITMAX
        iter=its
        call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,tolx, check, dfunc)
        fp=fret
        !Here is my attempt to constrain the oxygen positions JPH
        !pnew(n/3+1:n) = p(n/3+1:n)
        xi(1:n) = pnew(1:n) - p(1:n)
        p(1:n) = pnew(1:n)
        test=0.d0
        do 14 i=1,n/3
          temp=dabs(xi(i))/dmax1(dabs(p(i)),1.d0)
        !JPH
          if(temp.gt.test)test=temp
14      continue
        if(test.lt.TOLX) goto 1111
        dg(1:n)=g(1:n)
        !...................
        call dfunc(n, p, g, tmp)
        test=0.d0
        den=dmax1(fret,1.d0)
        do 16 i=1,n
          temp=dabs(g(i))*dmax1(dabs(p(i)),1.d0)/den
          if(temp.gt.test)test=temp
16      continue
        if(test.lt.gtol) goto 1111
        testg=test
        if (iprint>0.and.mod(its,iprint)==0)write(*,*)'iter=',iter, fp

  
        dg(1:n) = g(1:n) - dg(1:n)
        do i=1,n
          hdg(i)=0.d0
          do  j=1,n
            hdg(i)=hdg(i)+hessin(i,j)*dg(j)
          enddo
        enddo
        fac = sum(dg(1:n)*xi(1:n))
        fae = sum(dg(1:n)*hdg(1:n))
        sumdg=sum(dg(1:n)**2)
        sumxi=sum(xi(1:n)**2)
        if(fac.gt.dsqrt(EPS*sumdg*sumxi))then
          fac=1.d0/fac
          fad=1.d0/fae
          do 22 i=1,n
            dg(i)=fac*xi(i)-fad*hdg(i)
22        continue
          do 24 i=1,n
            do 23 j=i,n
              hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+  &
                                              fae*dg(i)*dg(j)
              hessin(j,i)=hessin(i,j)
23          continue
24        continue
        endif
        do 26 i=1,n
          xi(i)=0.d0
          do 25 j=1,n
            xi(i)=xi(i)-hessin(i,j)*g(j)
25        continue
26      continue
27    continue
      stop 'too many iterations in dfpmin'
1111  continue   ! return

      fret = fp

      deallocate(xi)
      deallocate(hessin)
      deallocate(g)
      deallocate(dg)
      deallocate(hdg)
      deallocate(pnew)
      END SUBROUTINE dfpmin
!.......................................................................
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,tolx, check, dfunc)
!.......................................................................
      implicit none
      LOGICAL check
      integer :: n
      double precision ::  f,fold,stpmax,g(n),p(n), x(n),xold(n), g0(n)
      double precision, parameter :: ALF=1.d-4
      EXTERNAL dfunc
      INTEGER i
      double precision :: a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum, &
                          temp, test,tmplam, tolx

      check=.false.
      sum=0.d0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=dsqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.d0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      if(slope.ge.0.d0) stop 'roundoff problem in lnsrch'
      test=0.d0
      do 14 i=1,n
        temp=dabs(p(i))/dmax1(dabs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.d0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        call dfunc(n, x, g0, f)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          goto 1111
        else if(f.le.fold+ALF*alam*slope)then
          goto 1111
        else
          if(alam.eq.1.d0)then
            tmplam=-slope/(2.d0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.d0*b)
            else
              disc=b*b-3.d0*a*slope
              if(disc.lt.0.d0)then
                tmplam=.5d0*alam
              else if(b.le.0.d0)then
                tmplam=(-b+dsqrt(disc))/(3.d0*a)
              else
                tmplam=-slope/(b+dsqrt(disc))
              endif
            endif
            if(tmplam.gt.0.5d0*alam)tmplam=0.5d0*alam
          endif
        endif
        alam2=alam
        f2=f
        alam=dmax1(tmplam,0.1d0*alam)
      goto 1
1111  return
      END SUBROUTINE lnsrch
end module dfpmin_mod

