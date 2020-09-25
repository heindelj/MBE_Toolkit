subroutine read_natoms(iunit, Nw)
!*** reads the number of atoms from a file
implicit none
integer, intent(in) :: iunit
integer, intent(out) :: Nw
integer :: Natoms

read(10,*)Natoms
read(10,*)
Nw = Natoms/3

end subroutine read_natoms
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine readXYZ(iunit, Nw, RR)
!*** reads the coordinates of atoms from a file and stores them properly in
!*** the array (RR)
implicit none
integer, intent(in) :: iunit
integer, intent(in) :: Nw
double precision, intent(out), dimension(3, 3*Nw) :: RR
double precision, parameter :: RcSQ_wat =2.2025d0
integer         :: i, j, i_top, i_end, atom_id, i_conn, idx, Natoms
double precision, dimension(3) :: tmpR, R_O
character(len=100) :: filename
character(len=2) :: atom_name
character(len=2) :: adj

Natoms = 3*Nw
i_top = 1
i_end = NAtoms
do i=1,NAtoms
   read(10,*)atom_name,tmpR(1:3)
    adj = adjustl(atom_name)
   if ((adj == "8") .or. (adj=="O") .or. (adj=="o") ) then
      RR(1:3, i_top) = tmpR(1:3)
      i_top = i_top + 1
   else if ((adj == "1") .or. (adj=="H") .or. (adj=="h"))  then
      RR(1:3, i_end) = tmpR(1:3)
      i_end = i_end-1
   endif
enddo

do i=1, Nw
   R_O(1:3) = RR(1:3, i)
   i_conn = 0
   do j=Nw+1, NAtoms
      if (dot_product(R_O(1:3)-RR(1:3,j), R_O(1:3)-RR(1:3,j)) < RcSQ_wat ) then
         i_conn=i_conn+1
         idx = Nw+ 2*(i-1) + i_conn
         tmpR(1:3)    = RR(1:3,j)
         if (i_conn<=2) then
            RR(1:3,j)    = RR(1:3, idx)
            RR(1:3, idx) = tmpR(1:3)
         endif
      endif
   enddo
   if (i_conn/=2) then
      print*, 'Error in the det,of the structure water=',i,' conn. H=',i_conn
      stop
   endif
enddo
return
end subroutine readXYZ
