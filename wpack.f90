! write a packmol input

subroutine wpack(cell,nwat,molname,string)
use options, only: solvfile
implicit none
real(8) cell(3)
character(*) molname,string
character(200) name2,packname
integer i,nwat,io

name2='solv-'//molname
packname='pack'//trim(string)//'.in'

!cell=cell*0.5d0


open(newunit=io,file=trim(packname))
  write(io,'(a)') '#automated packmol input'
  write(io,'(a)') ''
  write(io,'(a)') ' seed -1'
  write(io,'(a)') ''
  write(io,'(a)') ' tolerance 2.5'
  write(io,'(a)') ''
  write(io,'(a)') ' filetype xyz'
  write(io,'(a)') ' '
  write(io,'(a,a)') 'output ',name2
  write(io,'(a)') ''
  write(io,'(a,a)') 'structure ',molname
  write(io,'(a)') ' number 1'
  write(io,'(a)') ' center'
  write(io,'(a)') ' fixed 0. 0. 0. 0. 0. 0.'
  write(io,'(a)') 'end structure'
  write(io,'(a)') ''
  write(io,'(a)') 'structure '//trim(solvfile)
  write(io,'(a,I5)') 'number ',nwat
  write(io,'(a,6(F7.2,1x))') ' inside box ',-cell(1:3)*0.5d0,cell(1:3)*0.5d0
  write(io,'(a)') ' end structure'
close(io)


print*,' ->  writing  : ',trim(packname)

return
end subroutine

! read xyz file solvent
subroutine solvent_read(solvfile,smass,svol)
use atomdata, only: ams
implicit none
character(*) solvfile
integer io,i,j
integer nat
real(8), allocatable:: xyz(:,:)
integer, allocatable:: iat(:)
real(8) smass,cmax(3),cmin(3)
real(8) vol,svol,cell(3)

call tmolrd(nat,xyz,iat,solvfile,.false.,.true.)
allocate(xyz(3,nat),iat(nat))
call tmolrd(nat,xyz,iat,solvfile,.true.,.false.)

xyz=xyz*0.52917726d0

smass=0d0
do i=1,nat
  smass=smass+ams(iat(i))
enddo
print*, ' solvent mass',smass


cmax=-1d9
cmin=1d9
do i=1,nat
 do j=1,3
  if(xyz(j,i)>cmax(j)) then
    cmax(j)=xyz(j,i)
  endif
  if(xyz(j,i)<cmin(j)) then
   cmin(j)=xyz(j,i)
  endif
 enddo
enddo

! svol can become 0 if we have a flat structure
cell=cmax-cmin
vol=cell(1)*cell(2)*cell(3)
call substrate_volume(cmax,cmin,nat,iat,xyz,vol,svol)
!print*, ' est. solvent volume ',svol

end subroutine

! writes a xyz file from pre-defined solvents
subroutine solvent_db(isolv)
use options, only: solvfile
implicit none
integer isolv,io
! PBE0-D3(BJ)/def2-TZVP opt

write(*,'(a,$)') ' Providing DFT(PBE0) solvent file ='
select case(isolv)
 case(1)
   write(*,'(a)') ' WATER (wat.xyz) '
    solvfile='wat.xyz'
    open(newunit=io,file=solvfile)
    write(io,'(a)') ' 3'
    write(io,'(a)') ' water '
    write(io,'(a)')'O  -1.551007  -0.114520   0.000000' 
    write(io,'(a)')'H  -1.934259   0.762503   0.000000'
    write(io,'(a)')'H  -0.599677   0.040712   0.000000'
 case(2)
   write(*,'(a)') ' dimethylether (ether.xyz) '
    solvfile='ether.xyz'
    open(newunit=io,file=solvfile)
    write(io,'(a)')' 9 '
    write(io,'(a)')'  dimethylether '
    write(io,'(a)')'H     1.4734091    0.1727921    1.5447025 '
    write(io,'(a)')'C     1.0136192    0.0526327    0.5634890 '
    write(io,'(a)')'H     1.3699874    0.8558984   -0.0976041 '
    write(io,'(a)')'H     1.3315664   -0.9096565    0.1366547 '
    write(io,'(a)')'O    -0.3738635    0.1054901    0.7342934 '
    write(io,'(a)')'C    -1.0573440   -0.0402874   -0.4776021 '
    write(io,'(a)')'H    -0.7969680    0.7588678   -1.1867049 '
    write(io,'(a)')'H    -0.8350960   -1.0067386   -0.9528601 '
    write(io,'(a)')'H    -2.1253105    0.0110013   -0.2643684'
 case(3)
   write(*,'(a)') ' ETHANE (ethane.xyz) '
    solvfile='ethane.xyz'
    open(newunit=io,file=solvfile)
    write(io,'(a)')'  8'
    write(io,'(a)')'   ethane'
    write(io,'(a)')'H     1.1860047   -0.0042746    0.9843370 '
    write(io,'(a)')'C     0.7587257   -0.0226867   -0.0210586 '
    write(io,'(a)')'H     1.1679000    0.8304883   -0.5673351 '
    write(io,'(a)')'H     1.1167607   -0.9299977   -0.5133547 '
    write(io,'(a)')'C    -0.7587253    0.0226890    0.0210609 '
    write(io,'(a)')'H    -1.1678977   -0.8304647    0.5673717 '
    write(io,'(a)')'H    -1.1167625    0.9300233    0.5133116 '
    write(io,'(a)')'H    -1.1860056    0.0042232   -0.9843328'
 case(4)
   write(*,'(a)') ' FORMALDEHYDE (formal.xyz) '
    solvfile='formal.xyz'
    open(newunit=io,file=solvfile)
    write(io,'(a)')'   4'
    write(io,'(a)')'  formaldehyde'
    write(io,'(a)')'H    -0.9400284    0.0000000    0.5927836'
    write(io,'(a)')'C     0.0000000    0.0000000    0.0048954'
    write(io,'(a)')'H     0.9400284    0.0000000    0.5927836'
    write(io,'(a)')'O     0.0000000    0.0000000   -1.1904875'
 case(5)
   write(*,'(a)') ' FORMIC ACID (formac.xyz) '
    solvfile='formac.xyz'
    open(newunit=io,file=solvfile)
    write(io,'(a)')'5'
    write(io,'(a)')' formic acid'
    write(io,'(a)')'H     0.8829913   -1.1152319    0.0000000'
    write(io,'(a)')'C     0.3970836   -0.1293889    0.0000000'
    write(io,'(a)')'O     0.9695763    0.9185929    0.0000000'
    write(io,'(a)')'O    -0.9295061   -0.2804407    0.0000000'
    write(io,'(a)')'H    -1.3201451    0.6064685    0.0000000'
 case(6)
   write(*,'(a)') ' HEXANE (hexane.xyz) '
    solvfile='hexane.xyz'
    open(newunit=io,file=solvfile)
    write(io,'(a)')'20                                     '
    write(io,'(a)')' hexane '                           
    write(io,'(a)')'H     3.6929269   -0.3507066    1.6827636 '
    write(io,'(a)')'C     3.0898434   -0.2810925    0.7745685 '
    write(io,'(a)')'H     3.4429516    0.5843683    0.2066353 '
    write(io,'(a)')'H     3.2911382   -1.1719213    0.1728394 '
    write(io,'(a)')'C     1.6097018   -0.1593034    1.0914455 '
    write(io,'(a)')'H     1.2878360   -1.0217466    1.6866418 '
    write(io,'(a)')'H     1.4386160    0.7218340    1.7208196 '
    write(io,'(a)')'C     0.7416199   -0.0599623   -0.1518206 '
    write(io,'(a)')'H     1.0618406    0.8038036   -0.7482156 '
    write(io,'(a)')'H     0.9127937   -0.9408143   -0.7836253 '
    write(io,'(a)')'C    -0.7415692    0.0605347    0.1518283 '
    write(io,'(a)')'H    -1.0618743   -0.8034882    0.7478038 '
    write(io,'(a)')'H    -0.9127030    0.9411016    0.7840439 '
    write(io,'(a)')'C    -1.6095994    0.1605057   -1.0914233 '
    write(io,'(a)')'H    -1.2888344    1.0242506   -1.6853228 '
    write(io,'(a)')'H    -1.4372094   -0.7194856   -1.7220453 '
    write(io,'(a)')'C    -3.0899549    0.2797875   -0.7745921 '
    write(io,'(a)')'H    -3.4420162   -0.5871687   -0.2082936 '
    write(io,'(a)')'H    -3.2925385    1.1692560   -0.1712832 '
    write(io,'(a)')'H    -3.6929686    0.3502475   -1.6827681 '
 case(7)
   write(*,'(a)') ' DMSO (dmso.xyz) '
    solvfile='dmso.xyz'
    open(newunit=io,file=solvfile)
    write(io,'(a)')'  10 '
    write(io,'(a)')' dmso '
    write(io,'(a)')'C   -1.19694   0.61049  -0.19515'
    write(io,'(a)')'S    0.58922   0.81630  -0.24748'
    write(io,'(a)')'C    0.99104  -0.73980   0.56018'
    write(io,'(a)')'O    0.99374   0.67799  -1.66562'
    write(io,'(a)')'H   -1.62661   1.46643  -0.71528'
    write(io,'(a)')'H   -1.53998   0.59680   0.84082'
    write(io,'(a)')'H   -1.46561  -0.31199  -0.71282'
    write(io,'(a)')'H    0.61661  -0.73499   1.58529'
    write(io,'(a)')'H    2.07798  -0.81933   0.56427'
    write(io,'(a)')'H    0.56055  -1.56189  -0.01421'
 case default
  print*, 'select from:'
  print*, '1 = water           / H2O'
  print*, '2 = di-methyl-ether / C2H6O'
  print*, '3 = ethane          / C2H6 '
  print*, '4 = formaldehyde    / CH2O'
  print*, '5 = formic-acid     / CH2O2'
  print*, '6 = hexane          / C6H14 '
  print*, '7 = DMSO            / C2H6OS'
  print*, '8 ='
  print*, '9 ='
  print*, '10 ='
  stop 'error'
end select

close(io)
return
end subroutine



