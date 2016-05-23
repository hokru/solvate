
! process command line options
subroutine read_opts(infile)
use options
implicit none
integer maxarg,i
character(*) infile
character(120), allocatable :: arg(:)
character(120) aa,narg
real(8) s2r
integer s2i
logical fstr

! set defaults
buffer=5
tol=2.5
solvent_mass=18
p=1d0
isolv=1 ! water

do_cema=.false.
do_addion=.false.
do_cell=.false.


call getarg(1,infile)
if(infile=='-h ') call help
if(infile=='') stop 'no structure file!!'

maxarg=iargc()
if(maxarg.gt.0) then
allocate(arg(maxarg))

do i=1,maxarg
   call getarg(i,arg(i))
enddo

! loop over all arguments
do i=2,maxarg
  aa=arg(i)
  if(i/=maxarg) narg=trim(arg(i+1))


  if(fstr(aa,'-h ')) call help

  if(fstr(aa,'-water ')) isolv=1

  if(fstr(aa,'-cema ')) do_cema=.true.
  if(fstr(aa,'-cell ')) do_cell=.true.

  if(fstr(aa,'-b ').or.fstr(aa,'-buffer')) then
   buffer=s2r(narg)
  endif

  if(fstr(aa,'-mass')) then
   solvent_mass=s2r(narg)
  endif

  if(fstr(aa,'-tol')) then
   tol=s2r(narg)
  endif

  ! read a solvent file
  if(fstr(aa,'-solvent')) then
   isolv=-1
   solvfile=trim(narg)
  endif

  ! select from database
  if(fstr(aa,'-select')) then
   isolv=s2i((narg))
  endif


! future stuff
  if(fstr(aa,'-addion')) then
   do_addion=.true.
  endif

enddo
endif
end subroutine


subroutine help
implicit none
 print '(2x,a)', ''
 print '(2x,a)', ''
 print '(2x,a)',' *** H E L P ***'
 print '(2x,a)', ''
 print '(2x,a)', '-water              :   shortcut to "-select 1" '
 print '(2x,a)', ''
 print '(2x,a)', '-b/-buffer <float>  :   specify buffer region in A'
 print '(2x,a)', ''
 print '(2x,a)', '-cema               :   do center of mass transformation (default off)'
 print '(2x,a)', ''
 print '(2x,a)', '-cell               :   cell info only'
 print '(2x,a)', ''
 print '(2x,a)', '-mass <float>       :   specify solvent mass'
 print '(2x,a)', ''
 print '(2x,a)', '-tol                :   packmol packing tolerance in A (default 2.5)'
 print '(2x,a)', ''
 print '(2x,a)', '-solvent <file>     :   specify file with solvent molecule (formats: xyz,tmol)'
 print '(2x,a)', ''
 print '(2x,a)', '[-addion            :   not implemented yet]'
 print '(2x,a)', ''
 print '(2x,a)', '-select <int>       :   database solvent (pbe0)'
 print '(4x,a)', 'select from:'
 print '(4x,a)', '1 = water           / H2O'
 print '(4x,a)', '2 = di-methyl-ether / C2H6O'
 print '(4x,a)', '3 = ethane          / C2H6 '
 print '(4x,a)', '4 = formaldehyde    / CH2O'
 print '(4x,a)', '5 = formic-acid     / CH2O2'
 print '(4x,a)', '6 = hexane          / C6H14 '
 print '(4x,a)', '7 = DMSO            / C2H6OS'
 print '(2x,a)', ''
stop 

end subroutine
