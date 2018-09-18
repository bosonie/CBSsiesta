program cbs

use hsx_m, only: read_hsx_file,hsx_t,read_OrbIndx_file,createHS

use FDF

use calccbs

implicit none

logical :: Ksetting, SymSetting
double precision :: EMIN, EMAX,p1,p2
type(hsx_t)  :: hw
double complex,allocatable :: S(:,:),H(:,:,:), SS(:,:), HS(:,:,:)
double complex,allocatable :: SSBIS(:,:), HSBIS(:,:,:), SSTRIS(:,:),HSTRIS(:,:,:)
character(len=10) :: pre
integer :: bins,j,l,maxsuperc(3),dir,norb,nspin,ooo
integer,allocatable :: isc(:,:)
double precision :: k(3),a1(3),a2(3),a3(3),p(2),ap(3),b1(3),b2(3)


!Reading cbs.in, mandatory call it cbs.in
CALL FDF_INIT("cbs.in","out")
Ksetting = FDF_BOOLEAN('VaringKsetting',.FALSE.)
SymSetting = FDF_BOOLEAN('SymSetting',.FALSE.)
EMIN = FDF_PHYSICAL('CBS.Emin',-1.0d10,'eV')
EMAX = FDF_PHYSICAL('CBS.Emax',1.0d10,'eV')
bins = FDF_INTEGER('CBS.bins',100)
dir = FDF_INTEGER('CBS.dir',1)
pre = FDF_STRING('SystemLabel','siesta')
p(1) = FDF_GET('CBS.2DBZk1',0d0)
p(2) = FDF_GET('CBS.2DBZk2',0d0)


!Reading .HSX file
call read_hsx_file(hw,trim(adjustl(pre))//".HSX")

!Reading .ORB_INDX file
allocate(isc(hw%no_s,3))
call read_OrbIndx_file(isc,trim(adjustl(pre))//".ORB_INDX")

!Reading .XV file
open(189,file=trim(adjustl(pre))//".XV",status="old")
read(189,*) a1(1),a1(2),a1(3)
read(189,*) a2(1),a2(2),a2(3)
read(189,*) a3(1),a3(2),a3(3)
close(189)

!Summary on the direction you calculate the complex k
!calculate the 2D reciprocal vectors! 
if (dir==1) then
    write(*,*) "You are calculating the complex band structure in the direction perpendicular"
    write(*,*) "to the plane generated by lattice vectors a2 and a3 (in next lines in bohr)"
    write(*,*) "a2=",a2
    write(*,*) "a3=",a3
    call get_2D_rec_vect(a2,a3,b1,b2,ap)
    write(*,*) "In cartesian units this direction is"
    write(*,*) ap
else if (dir==2) then
    write(*,*) "You are calculating the complex band structure in the direction perpendicular"
    write(*,*) "to the plane generated by lattice vectors a3 and a1 (in next lines in bohr)"
    write(*,*) "a3=",a3
    write(*,*) "a1=",a1
    call get_2D_rec_vect(a3,a1,b1,b2,ap)
    write(*,*) "In cartesian units this direction is"
    write(*,*) ap
elseif (dir==3) then
    write(*,*) "You are calculating the complex band structure in the direction perpendicular"
    write(*,*) "to the plane generated by lattice vectors a1 and a2 (in next lines in bohr)"
    write(*,*) "a1=",a1
    write(*,*) "a2=",a2
    call get_2D_rec_vect(a1,a2,b1,b2,ap)
    write(*,*) "In cartesian units this direction is"
    write(*,*) ap
else
    STOP "Error: you messed up with direction, allowed values are 1, 2 or 3"
endif

!Norbital and Nspin
norb=hw%no_u
nspin=hw%nspin

!Understanding how many replicas of the elementary cell
!are in the supercell in each direction
maxsuperc(:)=1
do j=1,3
  do l=1,hw%no_s
    if (isc(l,j).gt.maxsuperc(j)) then
      maxsuperc(j)=isc(l,j)
    endif
  enddo
enddo

!Last summary on the calculation properties
write(*,*) "N cell replicas in choosen direc",maxsuperc(dir), "norb",hw%no_u, "nspin", hw%nspin

if (maxsuperc(dir).gt.3) then
   write(*,*) "warning, number of replica in the choosen direction is > 3"
   write(*,*) "for the calculation we will use the 3 layers structure"
endif

k(:)=0.d0



!Fist case, I have fixed E (Emax of input) and 
!1) fix k as well and find eigenvectors for symmetries
!2) vary k
if (Ksetting.or.SymSetting) then
 
  if (SymSetting) then
        k(:)=p(1)*b1(:)+p(2)*b2(:)
        write(*,*) "analysing projections at"
        write(*,*) "E =", emax
        write(*,*) "K =", k

        allocate(S(hw%no_u,hw%no_u),H(hw%no_u,hw%no_u,hw%nspin),SS(hw%no_u,hw%no_u),HS(hw%no_u,hw%no_u,hw%nspin))
        allocate(SSBIS(hw%no_u,hw%no_u), HSBIS(hw%no_u,hw%no_u,hw%nspin))
        allocate(SSTRIS(hw%no_u,hw%no_u),HSTRIS(hw%no_u,hw%no_u,hw%nspin))

        call createHS(hw,dir,isc,k,S,H,HS,SS,HSBIS,SSBIS,HSTRIS,SSTRIS)

        if (maxsuperc(dir).eq.1) then
                call KFirstLayerInteraction(EMAX,H,S,HS,SS,p1,p2,SymSetting)
        elseif (maxsuperc(dir).eq.2) then
                call KSecondLayerInteraction(EMAX,H,S,HS,SS,HSBIS,SSBIS,p1,p2,SymSetting)
        else
                call KThirdLayerInteraction(EMAX,H,S,HS,SS,HSBIS,SSBIS,HSTRIS,SSTRIS,p1,p2,SymSetting)
        endif

        deallocate(S,H,SS,HS)
        deallocate(SSBIS,HSBIS,SSTRIS,HSTRIS)

   else
        write(*,*) "b1 and b2 are the 2D reciprocal vector:"
        write(*,*) "b1=",b1(:)
        write(*,*) "b2=",b2(:)
        write(*,*) "Varing k parallel for E =",emax
       
        p1=0.d0
        p2=0.d0
        do l=-15,15
          do j=-15,15
           p1=l*0.5/15.d0
           p2=j*0.5/15.d0
           k(:)=p1*b1(:)+p2*b2(:)
      !     write(*,*) k
       
          allocate(S(hw%no_u,hw%no_u),H(hw%no_u,hw%no_u,hw%nspin),SS(hw%no_u,hw%no_u),HS(hw%no_u,hw%no_u,hw%nspin))
          allocate(SSBIS(hw%no_u,hw%no_u), HSBIS(hw%no_u,hw%no_u,hw%nspin))
          allocate(SSTRIS(hw%no_u,hw%no_u),HSTRIS(hw%no_u,hw%no_u,hw%nspin))
       
          call createHS(hw,dir,isc,k,S,H,HS,SS,HSBIS,SSBIS,HSTRIS,SSTRIS)
       
          if (maxsuperc(dir).eq.1) then
             call KFirstLayerInteraction(EMAX,H,S,HS,SS,p1,p2,SymSetting)
          elseif (maxsuperc(dir).eq.2) then
             call KSecondLayerInteraction(EMAX,H,S,HS,SS,HSBIS,SSBIS,p1,p2,SymSetting)
          else
             call KThirdLayerInteraction(EMAX,H,S,HS,SS,HSBIS,SSBIS,HSTRIS,SSTRIS,p1,p2,SymSetting)
          endif
       
          deallocate(S,H,SS,HS)
          deallocate(SSBIS,HSBIS,SSTRIS,HSTRIS)
         enddo
        enddo

   endif

!Fixed k and change E
else

   !Define the k parallel where to calculate the CBS and print out a summary
   k(:)=p(1)*b1(:)+p(2)*b2(:)
   write(*,*) "The k parallel you choose is"
   write(*,*) p(1),"*b1 + ",p(2),"*b2"
   write(*,*) "where b1 and b2 are reciprocal lattice vectors of the 2D BZ"
   write(*,*) "b1 and b2 in the next lines:"
   write(*,*) "b1=",b1(:)
   write(*,*) "b2=",b2(:)
   write(*,*) "The k parallel you choose in cartesian units (bohr) is:"
   write(*,*) k(:)
   
   
   allocate(S(hw%no_u,hw%no_u),H(hw%no_u,hw%no_u,hw%nspin), SS(hw%no_u,hw%no_u), HS(hw%no_u,hw%no_u,hw%nspin))
   allocate(SSBIS(hw%no_u,hw%no_u), HSBIS(hw%no_u,hw%no_u,hw%nspin))
   allocate(SSTRIS(hw%no_u,hw%no_u),HSTRIS(hw%no_u,hw%no_u,hw%nspin))
   
   call createHS(hw,dir,isc,k,S,H,HS,SS,HSBIS,SSBIS,HSTRIS,SSTRIS)
  
   if (maxsuperc(dir).eq.1) then
      call FirstLayerInteraction(EMIN,EMAX,bins,H,S,HS,SS)
   elseif (maxsuperc(dir).eq.2) then
      call SecondLayerInteraction(EMIN,EMAX,bins,H,S,HS,SS,HSBIS,SSBIS)
   else
      call ThirdLayerInteraction(EMIN,EMAX,bins,H,S,HS,SS,HSBIS,SSBIS,HSTRIS,SSTRIS)
   endif
   
   deallocate(S,H,SS,HS)
   deallocate(SSBIS,HSBIS,SSTRIS,HSTRIS)

endif


deallocate(isc)

end program


subroutine get_2D_rec_vect(aa,aabis,b1,b2,ap)
!The theory behind:
!I first calculate the vector perpendicular to the aa and aabis,
!I call it ap, second I normalise it.
!Third I solve the system aa*b1=2pi,aabis*b1=0,ap*b1=0
!which is the definition of reciprocal vectors and to impose b1
!is in the same plane of aa and aabis
!Same for b2
implicit none

double precision, intent(in)  :: aa(3),aabis(3)
double precision, intent(out) :: b1(3),b2(3),ap(3)

integer :: i,info
double precision :: Matrix(3,3), OtherMat(3,3),OtherMat2(3,3),IPIV(3,3)
double precision :: PI,apmod,det,bisdet,trisdet

ap(:)=0.d0
ap(1)=aa(2)*aabis(3)-aa(3)*aabis(2)
ap(2)=aabis(1)*aa(3)-aa(1)*aabis(3)
ap(3)=aa(1)*aabis(2)-aa(2)*aabis(1)
apmod=sqrt(ap(1)**2+ap(2)**2+ap(3)**2)
ap(:)=ap(:)/apmod
Matrix(1,:)=ap(:)
Matrix(2,:)=aa(:)
Matrix(3,:)=aabis(:)
det=0.d0
bisdet=0.d0
trisdet=0.d0
call a33det(Matrix, det)
b1(:)=0d0
b2(:)=0d0
PI=4.D0*DATAN(1.D0)
do i=1,3
   OtherMat=Matrix
   OtherMat(1,i)=0
   OtherMat(2,i)=0
   OtherMat(3,i)=2*PI
   OtherMat2=Matrix
   OtherMat2(1,i)=0
   OtherMat2(2,i)=2*PI
   OtherMat2(3,i)=0
   call a33det(OtherMat, bisdet)
   call a33det(OtherMat2, trisdet)
   b1(i)=bisdet/det
   b2(i)=trisdet/det
enddo

return 

end subroutine



subroutine a33det(A,DET)

implicit none

double precision, intent(in)  :: A(3,3)
double precision, intent(out) :: DET


DET =   A(1,1)*A(2,2)*A(3,3)  &
      - A(1,1)*A(2,3)*A(3,2)  &
      - A(1,2)*A(2,1)*A(3,3)  &
      + A(1,2)*A(2,3)*A(3,1)  &
      + A(1,3)*A(2,1)*A(3,2)  &
      - A(1,3)*A(2,2)*A(3,1)

return

end subroutine a33det
