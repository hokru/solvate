!**************************************************************************************
!    Copyright 2016 Holger Kruse                                                      *
!                                                                                     *
!    This file is part of solvate.                                                    *
!                                                                                     *
!    solvate is free software: you can redistribute it and/or modify                  *
!    it under the terms of the GNU Lesser General Public License as published by      *
!    the Free Software Foundation, either version 3 of the License, or                *
!    (at your option) any later version.                                              *
!                                                                                     *
!    solvate is distributed in the hope that it will be useful,                       *
!    but WITHOUT ANY WARRANTY; without even the implied warranty of                   *
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                    *
!    GNU Lesser General Public License for more details.                              *
!                                                                                     *
!    You should have received a copy of the GNU Lesser General Public License         *
!    along with solvate.  If not, see <http://www.gnu.org/licenses/>.                 *
!                                                                                     *
!                                                                                     *
!  Feel free to contact me: Holger Kruse (mail2holger@gmail.com)                      *
!                                                                                     *
!**************************************************************************************
!*************************************************
!* usage: Solvate <molecule> <buffer size>       *
!*                                               *
!* writes packmol input                          *
!* determines solute Volume by crude MC          *
!*                                               *
!*                                               *
!*  if cema instead of <buffer size> it          *
!*  will make a center of mass transformation    *
!*  output:cell.xyz                              *
!*                                               *
!*                                               *
!* author(s):                                    *
!* HK (mail2holger@gmail.com)                    *
!*************************************************




program getcellsize
use atomdata, only: ams
use options
implicit none
integer nat,i,j,k,l
integer, allocatable :: iat(:)

real(8) xmax(2),ymax(2),zmax(2),xmin(2),ymin(2),zmin(2),cell(3),vol,bohr,vander(100)
real(8) cmax(3),cmin(3)
real(8), allocatable :: xyz(:,:)

real(8) mass,nwat,molvol
real(8) smass, svol,com
character(120) infile,option


data vander/1.20d0,1.20d0,1.37d0,1.45d0,1.45d0,1.50d0,1.50d0, &
                 1.40d0,1.35d0,1.30d0,1.57d0,1.36d0,1.24d0,1.17d0, &
                 1.80d0,1.75d0,1.70d0,17*2.5d0,2.3d0,65*2.5d0/

bohr=0.52917726d0


call getarg(1,infile)
call read_opts(infile)

call tmolrd(nat,xyz,iat,infile,.false.,.true.)
allocate(xyz(3,nat),iat(nat))
call tmolrd(nat,xyz,iat,infile,.true.,.false.)

xyz=xyz*bohr

! mol mass
mass=0d0
do i=1,nat
  mass=mass+ams(iat(i))
enddo
print *,''
print *, 'molecular mass', mass
call composition(nat,iat)


if(isolv>0) call solvent_db(isolv)

print*,'using solvent file:', trim(solvfile)

if(do_cema) then
print*,''
print*,''
print*,' center of mass transformation'
 call getCOM(com,nat,xyz,iat,.true.)
 call getIntertia(nat,iat,xyz,.true.)
print*,''
print*,''
endif


xmax=0
ymax=0
zmax=0
xmin=0
ymin=0
zmin=0


do i=1,nat
 if(xyz(1,i)>xmax(1)) then
  xmax(1)=xyz(1,i)
  xmax(2)=i
 endif

 if(xyz(2,i)>ymax(1)) then
  ymax(1)=xyz(2,i)
  ymax(2)=i
 endif

 if(xyz(3,i)>zmax(1)) then
   zmax(1)=xyz(3,i)
   zmax(2)=i
 endif

 if(xyz(1,i)<xmin(1)) then
   xmin(1)=xyz(1,i)
   xmin(2)=i
 endif

 if(xyz(2,i)<ymin(1))then
   ymin(1)=xyz(2,i)
   ymin(2)=i
 endif

 if(xyz(3,i)<zmin(1)) then
  zmin(1)=xyz(3,i)
  zmin(2)=i
 endif
enddo

cmax(1)=xmax(1)
cmax(2)=ymax(1)
cmax(3)=zmax(1)

cmin(1)=xmin(1)
cmin(2)=ymin(1)
cmin(3)=zmin(1)

! debug
!print *,'max coord (x-y)',cmax
!print *,'min coord (x-y)',cmin
!print*,'Assuming water density 1g/ml)'

cell(1)=xmax(1)-xmin(1)
cell(2)=ymax(1)-ymin(1)
cell(3)=zmax(1)-zmin(1)

vol=cell(1)*cell(2)*cell(3)
print *,''
print*,'******  minimal box ************'
write(*,'(a,3F8.2)') ' cell / A   : ',cell(1:3)
print *,'vol  / A^3 : ',vol
print*,'****************************'
print *,''


print*,'******  VDW BOX ************'
! add vdw radii of corner atoms
cell(1)=xmax(1)-xmin(1)+vander(iat(int(xmax(2))))+vander(iat(int(xmin(2))))
cell(2)=ymax(1)-ymin(1)+vander(iat(int(ymax(2))))+vander(iat(int(ymin(2))))
cell(3)=zmax(1)-zmin(1)+vander(iat(int(zmax(2))))+vander(iat(int(zmin(2))))
vol=cell(1)*cell(2)*cell(3)
write(*,'(a,3F8.2)') ' cell+vdw buffer/ A   : ',cell(1:3)
print *,'vol(+vdW)  / A^3 : ',vol
print*,'****************************'

! for -cell we exit here
if(do_cell) return

print '(a)', ' estimating substrate volume..'
call substrate_volume(cmax,cmin,nat,iat,xyz,vol,molvol)
print*, ' est. V(solute) / A^3= ',molvol
print*,''
print*,''


print*,'******  minimal buffer BOX **********'
! add buffer region
cell(1)=xmax(1)-xmin(1)+buffer*2
cell(2)=ymax(1)-ymin(1)+buffer*2
cell(3)=zmax(1)-zmin(1)+buffer*2
vol=cell(1)*cell(2)*cell(3)
print*,' Buffer to each side: ', buffer
write(*,'(a,3F8.2)') ' cell / A   : ',cell(1:3)
print *,'vol(+buffer)  / A^3 : ',vol
print*, 'vol for water / A^3= ',vol-molvol

call solvent_read(solvfile,smass,svol)
nwat=(p*0.602*vol-molvol)/smass
print '(a,2x,I5)',' --> number of integer waters needed:', nint(nwat)
print '(a,2x,F7.3)', 'effective water density [g/mL]   : ', (nint(nwat)*smass)/(vol-molvol)/0.602
print '(a,2x,F7.3)', 'effective liquid density [g/mL  ]: ', (nint(nwat)*smass+mass)/vol/0.602
print*,'****************************'
call wpack(cell,int(nwat),trim(infile),'')

print*,''
print*,'******  equal size box **********'
cell=maxval(cell)
!cell=vol
vol=cell(1)*cell(2)*cell(3)
write(*,'(a,3F8.2)') ' cell / A   : ',cell(1:3)
print *,'vol(+buffer)  / A^3 : ',vol

nwat=(p*0.602*vol-molvol)/smass
print '(a,2x,F7.1)',' --> number of waters needed:', nwat

call wpack(cell,int(nwat),trim(infile),'cube')
print*,'****************************'

print*,'******  vdW buffer size **********'
cell(1)=vander(iat(int(xmax(2))))+vander(iat(int(xmin(2))))
cell(2)=vander(iat(int(ymax(2))))+vander(iat(int(ymin(2))))
cell(3)=vander(iat(int(zmax(2))))+vander(iat(int(zmin(2))))
write(*,'(a,3F8.2)') '  / A^3   : ',cell(1:3)
print*,'****************************'

print*, 'Suggested simulation box (buffer+vdw) '
cell(1)=xmax(1)-xmin(1)+vander(iat(int(xmax(2))))+vander(iat(int(xmin(2))))+buffer*2
cell(2)=ymax(1)-ymin(1)+vander(iat(int(ymax(2))))+vander(iat(int(ymin(2))))+buffer*2
cell(3)=zmax(1)-zmin(1)+vander(iat(int(zmax(2))))+vander(iat(int(zmin(2))))+buffer*2
write(*,'(a,3F8.2)') ' minimal / A^3   : ',cell(1:3)
cell=maxval(cell)
write(*,'(a,3F8.2)') ' cube / A^3   : ',cell(1:3)


print*,''



end program


subroutine substrate_volume(cmax,cmin,nat,iat,xyz,vol,molvol)
use atomdata, only: rvdw
implicit none
integer iat(nat),nat
integer(8) i,j
real(8) xyz(3,nat),molvol,cmax(3),cmin(3),randr,len,p(3),vol
real(8) old,hit,ii,nmax

! volume via simple monte carlo
! V(m)=V(box)*hits/tries
! nmax>1e9 will be REALLY SLOW!

!nmax=1e+8
nmax=nat*6*1000
if(nmax>1e7) nmax=1e7
hit=0
ii=0
molvol=0

! loop using floats
do 
 ii=ii+1
 old=molvol
 ! test point p
 p(1)=Randr(cmin(1),cmax(1))
 p(2)=Randr(cmin(2),cmax(2))
 p(3)=Randr(cmin(3),cmax(3))
 ! test if point p is within any atomic vdW radii of atom j
  do j=1,nat
  call veclen2(p,xyz(1,j),len)
  if(len<=rvdw(iat(j))) then
    hit=hit+1
    exit !exit at first hit to avoid double counting
  endif
 enddo
 if(ii>nmax) exit
enddo
  molvol=vol*hit/ii
!molvol=vol*dble(hit)/dble(nmax)
end subroutine


subroutine veclen(x,v)
implicit none
real(8) x(3),v
v=sqrt(dot_product(x,x))
end subroutine


subroutine veclen2(a,b,v)
implicit none
real(8) a(3),b(3),v,x(3)
x=a-b
v=sqrt(dot_product(x,x))
end subroutine




!************************************************************
!* reads a turbomole (bohr) or xmol (angst)rom file.        *
!* Tests if xmol starts with "number of atoms + blank" or   *
!* directly with the coordinates.                           *
!************************************************************
subroutine tmolrd(nat,xyz,iat,infile,echo,c_nat)
!use parm
implicit none

integer i,j,k,l,nat
real(8) xyz(3,*),energy,s2r
integer ifrez(nat),iat(nat)

character*2 cc,ff
character*80  atmp
character*(*) infile
real(kind=8) txyz(3,nat),xx(5)
real(kind=8) bohr
integer tiat(nat),nn,tifrez(nat),istat,iff
logical da,c_nat,echo
bohr=0.52917726d0
i=0
tifrez=0
iff=0

inquire(file=infile,exist=da)
select case (da)
case (.true.)
      if(echo) write(*,'('' reading...'',$)')

open(unit=33,file=infile)
! test for tmol or txyz file

 read(33,'(a)') atmp ! $coord
rewind(33)
if(index(atmp,'$coord').ne.0) then

 ! count number of atoms
 do while (da)
  read(33,'(a)',end=100) atmp ! $coord
   if(index(atmp,'$coord').ne.0) cycle
   if(index(atmp,'$').ne.0) exit
   i=i+1
  enddo
 nat=i
 100 continue
 if(c_nat) then
  close(33)
  return  ! just return number of atoms
 endif
 rewind(unit=33)

 ! read TMOL file
 read(33,*) atmp ! $coord
 do j=1,nat
    read(33,'(a)') atmp ! $coord
    backspace(33)
    if(index(atmp,' f ').ne.0) then
    read(33,*) txyz(1,j),txyz(2,j),txyz(3,j),cc,ff
     tifrez(j)=1
    iff=iff+1
    else ! default
     read(33,*) txyz(1,j),txyz(2,j),txyz(3,j),cc
   endif
   call elem(cc,tiat(j))
   txyz(1:3,j)=txyz(1:3,j)
  enddo
 if(echo) write(*,*) ' Turbomole file [bohr] :  ', trim(infile)

 close(33)

!**********
! PDB FILE
!**********
elseif(index(atmp,'REMARK').ne.0.or.index(atmp,'ATOM').ne.0) then
rewind(33)

 if(c_nat) then ! just count atoms
  do
   read (33,'(a)',end=777) atmp
   if(index(atmp,'ATOM').ne.0) i=i+1
  enddo
  777 continue
  nat=i
  close(33)
  return
  endif

   i=0
   do
     read (33,'(a)',end=778) atmp
     if(index(atmp,'ATOM').ne.0) then
      i=i+1
      txyz(1,i)=s2r(atmp(31:38))
      txyz(2,i)=s2r(atmp(39:46))
      txyz(3,i)=s2r(atmp(47:54))
     endif
   enddo
  txyz=0
 778 continue
   txyz=txyz*1d0/bohr
   close(33)
   if(echo) write(*,'(5x,'' PDB file [angst]: '',a)')  trim(infile)

else ! txyz file
       read(33,'(a)',end=101) atmp
! check for first two lines
       call readl(atmp,xx,nn)
        if(nn.gt.1) then   ! more than one argument found, assuming they are coords
           do
            nat=nat+1
            read(33,'(a)',end=123) atmp
           enddo
            if(c_nat) then
              close(33)
              return  ! just return number of atoms
            endif
          else
            nat=idint(xx(1))
            if(c_nat) then
             close(33)
             return  ! just return number of atoms
            endif
           read(33,'(a)',end=101) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
            read(33,'(a)') atmp
            call readl(atmp,xx,nn)
            call elem(atmp,tiat(i))
            txyz(1:3,i)=xx(1:3)*1d0/bohr
       enddo
 101  close(33)
      if(echo) write(*,'(5x,'' XYZ file [angst]: '',a)')  trim(infile)
      endif

if(maxval(tifrez,nat).eq.1) then
 if(echo) then
  write(*,'(a,x,I4,x,a)') '  found ',iff,' frozen cart. coordinates'
  if(iff.lt.50) then ! dont spam the output to much ...
   write(*,'(a,$)') '  atom nr: '
   do i=1,nat
     if(tifrez(i)==1) write(*,'(x,I2,$)') i
   enddo
   print*,''
  endif
 endif
endif

case (.false.)
  write(*,*) ' no input file <',trim(infile) ,'> found !! '
end select

do i=1,nat
xyz(1:3,i)=txyz(1:3,i)
iat(i)=tiat(i)
ifrez(i)=tifrez(i)
enddo
return

end subroutine


 SUBROUTINE ELEM(KEY1, NAT)
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 CHARACTER*(*) KEY1
 CHARACTER*2 ELEMNT(94),E

 DATA ELEMNT/'h ','he',                                      &
 'li','be','b ','c ','n ','o ','f ','ne',                    &
 'na','mg','al','si','p ','s ','cl','ar',                    &
 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',     &
 'zn','ga','ge','as','se','br','kr',                         &
 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',     &
 'cd','in','sn','sb','te','i ','xe',                         &
 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
 'au','hg','tl','pb','bi','po','at','rn',                    &
 'fr','ra','ac','th','pa','u ','np','pu'/

 nat=0
 e='  '
 k=1
 DO J=1,len(key1)
    if (k.gt.2)exit
    N=ICHAR(key1(J:J))
    if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
       e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
       k=k+1
    endif
    if(n.ge.ichar('a') .and. n.le.ichar('z') )then
       e(k:k)=key1(j:j)
       k=k+1
    endif
 enddo
 DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

character*2 FUNCTION ESYM(I)
CHARACTER*2 ELEMNT(94)
DATA ELEMNT/'h ','he',                                           &
  'li','be','b ','c ','n ','o ','f ','ne',                       &
  'na','mg','al','si','p ','s ','cl','ar',                       &
  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',        &
  'zn','ga','ge','as','se','br','kr',                            &
  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',        &
  'cd','in','sn','sb','te','i ','xe',                            &
  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',   &
  'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',   &
  'au','hg','tl','pb','bi','po','at','rn',                       &
  'fr','ra','ac','th','pa','u ','np','pu'/
  ESYM=ELEMNT(I)
  RETURN
END

      SUBROUTINE READL(A1,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) A1
      DIMENSION X(*)
      I=0
      IS=1
  10  I=I+1
      X(I)=READAA(A1,IS,IB,IE)
      IF(IB.GT.0 .AND. IE.GT.0) THEN
                                IS=IE
                                GOTO 10
      ENDIF
      N=I-1
      RETURN
      END


      FUNCTION READAA(A,ISTART,IEND,IEND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 READAA
      CHARACTER*(*) A
      NINE=ICHAR('9')
      IZERO=ICHAR('0')
      MINUS=ICHAR('-')
      IDOT=ICHAR('.')
      ND=ICHAR('D')
      NE=ICHAR('E')
      IBL=ICHAR(' ')
      IEND=0
      IEND2=0
      IDIG=0
      C1=0
      C2=0
      ONE=1.D0
      X = 1.D0
      NL=LEN(A)
      DO 10 J=ISTART,NL-1
         N=ICHAR(A(J:J))
         M=ICHAR(A(J+1:J+1))
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20

   10 CONTINUE
      READAA=0.D0
      RETURN
   20 CONTINUE
      IEND=J
      DO 30 I=J,NL
         N=ICHAR(A(I:I))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C1=C1*10+N-IZERO
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
            ONE=-1.D0
         ELSEIF(N.EQ.IDOT) THEN
            GOTO 40
         ELSE
            GOTO 60
         ENDIF
   30 CONTINUE
   40 CONTINUE
      IDIG=0
      DO 50 II=I+1,NL
         N=ICHAR(A(II:II))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C2=C2*10+N-IZERO
            X = X /10
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
            X=-X
         ELSE
            GOTO 60
         ENDIF
   50 CONTINUE
!C
!C PUT THE PIECES TOGETHER
!C
   60 CONTINUE
      READAA= ONE * ( C1 + C2 * X)
      DO 55 J=IEND,NL
         N=ICHAR(A(J:J))
         IEND2=J
         IF(N.EQ.IBL)RETURN
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
      RETURN

   57 C1=0.0D0
      ONE=1.0D0
      DO 31 I=J+1,NL
         N=ICHAR(A(I:I))
         IEND2=I
         IF(N.EQ.IBL)GOTO 70
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
         IF(N.EQ.MINUS)ONE=-1.0D0
   31 CONTINUE
   61 CONTINUE
   70 READAA=READAA*10**(ONE*C1)
      RETURN
      END

!******************
!* write xyz      *
!******************
subroutine wrxyz(iat,nat,xyz,infile)

!use parm
implicit none
integer i,j,k,l,nat,iat(nat)
real(8) xyz(3,*)
integer ierr
character(*) infile
character(2) esym
real(8) f

!f=0.5291770d0
f=1d0
open(unit=55,file=infile,iostat=ierr,status='replace')
if(ierr.ne.0) stop 'cannot write xopt.xyz'
write(55,'(I5)') nat
write(55,*)
!write(55,'(2F16.8)')energy,gnorm
do i=1,nat
 write(55,'(a2,5x,3(F18.14,3x))') esym(iat(i)), xyz(1,i)*f,xyz(2,i)*f,xyz(3,i)*f
enddo
close(55)
end subroutine wrxyz




subroutine DiagSM(xvar,mat,eig)
implicit none
integer i,j,k
real(8), allocatable :: aux(:)
integer info,lwork,xvar
real(8) ,intent(in) :: mat(xvar,xvar)
real(8) xx
real(8), intent(out) :: eig(xvar)

eig=0
call dsyev ('V','U',xvar,mat,xvar,eig,xx,-1,info)
lwork=int(xx)
allocate(aux(lwork))
call dsyev ('V','U',xvar,mat,xvar,eig,aux,lwork,info)
if(info/=0) print*,'Diagonalization failed !!'
end subroutine
