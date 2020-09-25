!***************************************************************************!
!     George S. Fanourgakis, PNL Laboratory                                 !
!     george.fanourgakis@pnl.gov                                            !
!                                                                           !
!     Disclaimer: This code is provided "as is". No representation is made  !
!                 that is  suitable for any particular purpose              !
!                                                                           !
!     Reference: Fanourgakis, G. S. and Xantheas, S. S.                     !
!                 J. Phys. Chem. A. vol 110, page 4100-4106, year 2006      !
!                                                                           !
!***************************************************************************!
program test_ttm2f
!*** this program is an example of how to use the subroutine "ttm2f" for the
!*** energy and derivatives calculation
!use model_mod
!use ttm2f_mod
use potential_mod
use dfpmin_mod

implicit none
character(len=100) :: fconfig
integer :: iw, Nw, ia, ixyz
double precision, dimension(:), allocatable :: p, g
double precision, dimension(:,:), allocatable :: RR, dRR, dRR0
double precision :: Energy, Energy0
double precision :: delta, diff
double precision :: fret, gtol, tolx, stpmx
integer :: iter, neqs,  iprint
!***
!*** reads the filename of the configuration file
!***
debug = .false.
!print*,'Enter model (2: TTM2-F, 21:TTM21-F,   3:TTM3-F)'
!read(*,*)imodel
imodel=21
!print*,'Enter configuration filename'
IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
  WRITE(*,*)'ERROR, FILE SHOULD BE PROVIDED ON COMMAND LINE, STOPPING'
  STOP
ENDIF
CALL GET_COMMAND_ARGUMENT(1, fconfig)
!read(*,'(a100)')fconfig
open(10,file=fconfig,status='old')  ! 
call read_natoms(10, Nw)   ! read the number of water molecules (Nw)
allocate(RR(3, 3*Nw))      ! allocate the array for the coodinates  (RR)
allocate(dRR(3, 3*Nw))     ! allocate the array for the derivatives (dRR)
call readXYZ(10, Nw, RR)   ! read and sort properly the coordinates of atoms
close(10)
!***
! Now we are ready to call the potential
! Print out also the energy (in kcal/mol) and the derivatives(kcal/mol/A)
!***
!call ttm2f(Nw,RR, dRR, Energy)
call potential(Nw,RR, dRR, Energy)
write(*,'("Energy (kcal/mol) = ",f12.6)')Energy
write(*,*)
print*,'----------Derivatives (kcal/mol/A)  ---------------'
!do iw=1,3*Nw
!   write(*,'(i3,3x,3(f12.6,2x))')iw,dRR(1:3,iw)
!enddo
!the lines are stored with oxygens first and then the hydrogens
!but the hydrogens are swapped from the original positions in the file
do iw=1,Nw
   write(*,'(i3,3x,3(f12.6,2x))')iw,dRR(1:3,iw)
   write(*,'(i3,3x,3(f12.6,2x))')iw,dRR(1:3,Nw+2*iw)
   write(*,'(i3,3x,3(f12.6,2x))')iw,dRR(1:3,Nw+2*iw-1)
enddo
!***
!***
!***   Test numerical and analytical derivatives
!***
!print*,'----------Test analytical and numerical Derivatives (kcal/mol/A)  ---------------'
!allocate(dRR0(3, 3*nw))
!delta=1.d-8
!call potential(Nw,RR, dRR0, Energy0)
!do ia=1, 3*nw
!   do ixyz=1, 3
!      RR(ixyz, ia) = RR(ixyz, ia) + delta
!      call potential(Nw,RR, dRR, Energy)
!      diff = (Energy-Energy0)/delta
!      RR(ixyz, ia) = RR(ixyz, ia) - delta
!      write(*,'(i5,2x,i2, 5x,"anal= ", f12.6,2x,"num= ", f12.6)')ia, ixyz, dRR0(ixyz, ia), diff
!   enddo
!enddo
!deallocate(dRR0)
!***
!***
!***   Minimize energy
!***
!print*,'----------  GEOMETRY OPTIMIZATION -------------------'
!Neqs = 3*3*Nw
!allocate(p(neqs))
!allocate(g(neqs))
!p = reshape(RR, (/neqs/))
!call min_func(neqs, p, g, fret)
!print*,'fret=', fret
!gtol = 1.d-8
!tolx = 1.d-8
!stpmx = 0.01d0
!iter=1000
!iprint=10
!
!call dfpmin(neqs, p, gtol, tolx, stpmx, iter, fret, iprint, min_func)
!
!RR = reshape(p,(/3, 3*nw/))
!call potential(Nw,RR, dRR, Energy)
!print*,'Final minimized Energy = ', Energy
!open(15, file='minimum.xyz', status='unknown')
!write(15,*)3*nw
!write(15,*)Energy
!do iw=1, nw
!   write(15,'(a3,3(1x, f15.8))')'O', RR(1:3, iw)
!   write(15,'(a3,3(1x, f15.8))')'H', RR(1:3, iw+Nw)
!   write(15,'(a3,3(1x, f15.8))')'H', RR(1:3, iw+2*Nw)
!enddo
!close(15)

! deallocate all arrays
!***
deallocate(RR)
deallocate(dRR)

!Contains
!   subroutine min_func(n, p, g, fret)
!   implicit none
!   integer, intent(in) :: n
!   double precision, dimension(n), intent(in) :: p
!   double precision, dimension(n), intent(out) :: g
!   double precision, intent(out) :: fret
!   !...
!   integer :: Neqs, iopt, na, nw
!   double precision, dimension(:,:), allocatable :: R, dR
!
!   na = n/3
!   nw = na/3
!   allocate(R(3,na))
!   allocate(dR(3, na))
!   R =reshape(p,(/3, na/))
!
!   call potential(Nw,R, dR, Energy)
!
!   fret = Energy
!   g = reshape(dR,(/n/))
!
!   deallocate(r, dr)
!   end subroutine min_func

end program test_ttm2f
