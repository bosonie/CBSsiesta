module calccbs

implicit none

double precision, parameter, private :: RyToEv = 13.605698066

contains


subroutine FirstLayerInteraction(emin,emax,nbin,HD,SD,HS,SS)

implicit none

integer :: l,i,m,j,dime,nspin,lwork,info
integer,intent(in) :: nbin
double precision,intent(in) :: emin, emax
double precision :: dummyemax,dummyemin,E,step
double complex,intent(in) :: SD(:,:),HD(:,:,:), SS(:,:), HS(:,:,:)
double complex,allocatable :: SSTRAS(:,:), HSTRAS(:,:,:)
double complex,allocatable :: BB(:,:),AA(:,:),work(:)
double complex,allocatable :: alpha(:),beta(:)
double complex :: VL(1,1)
double complex,allocatable :: VR(:,:)!,ppp(2,2),qqq(2,2),zzz(2),xxx(2)
double precision,allocatable :: rwork(:)
character(len=21) :: filename, filename2
!complex,intent(in) :: SSBIS(:,:), HSBIS(:,:,:), SSTRIS(:,:),HSTRIS(:,:,:)

dime=size(SD,dim=1)
nspin=size(HD,dim=3)

!Create conjugate transpose
allocate(SSTRAS(dime,dime),HSTRAS(dime,dime,nspin))
SSTRAS(:,:)=0.d0
HSTRAS(:,:,:)=0.d0
do m=1,dime
  do l=1,dime
    SSTRAS(m,l)=CONJG(SS(l,m))
    do i=1,nspin
     HSTRAS(m,l,i)=CONJG(HS(l,m,i))
    enddo
  enddo
enddo


allocate(BB(dime*2,dime*2),AA(dime*2,dime*2))
allocate(alpha(dime*2),beta(dime*2))
lwork=4*2*dime-1
allocate(work(lwork),rwork(8*2*dime))
allocate(VR(dime*2,dime*2))

do i=1,nspin
  write(filename,'("spin",I1,"layer1")')i
  write(filename2,'("spin",I1,"layer1purereal")')i
  open(unit=i+100,file=filename,status="unknown")
  open(unit=i+200,file=filename2,status="unknown") 

  !change energy in Ry
  dummyemax=emax/RyToEv
  dummyemin=emin/RyToEv
  step=(dummyemax-dummyemin)/nbin 
  do j=1,nbin+1
    alpha(:)=1.d0
    beta(:)=1.d0
    BB(:,:)=0.d0
    do l=1,dime
       BB(dime+l,dime+l)=1.d0
    enddo
    AA(:,:)=0.d0
    do l=1,dime
       AA(dime+l,l)=1.d0
    enddo
    E=dummyemin+(j-1)*step
    do m=1,dime
       do l=1,dime
            BB(m,l)=HSTRAS(m,l,i)-E*SSTRAS(m,l)
            AA(m,l)=E*SD(m,l)-HD(m,l,i)
            AA(m,dime+l)=E*SS(m,l)-HS(m,l,i)
       enddo
    enddo
    call ZGGEV('N', 'V', 2*dime, AA, 2*dime, BB, 2*dime, alpha, beta, VL, 1, VR, dime*2, work, lwork, rwork, info)
    if (info .ne. 0) then
         stop 'failed zggev diagonalisation'
         write(*,*) info
    end if
    do m=1,2*dime
         write(i+100,*) E*RyToEv, -real(real(log(alpha(m)/beta(m)))), real(aimag(log(alpha(m)/beta(m))))
         if (-real(real(log(alpha(m)/beta(m)))).le.0.001.and.-real(real(log(alpha(m)/beta(m)))).ge.-0.001) then
            write(i+200,*) E*RyToEv, real(aimag(log(alpha(m)/beta(m)))), -real(real(log(alpha(m)/beta(m))))
         endif
    enddo
  enddo
  close(i+100)
enddo 

deallocate(AA,BB,alpha,beta)
deallocate(work,rwork,VR)
deallocate(HSTRAS,SSTRAS)

end subroutine FirstLayerInteraction



subroutine SecondLayerInteraction(emin,emax,nbin,HD,SD,HS,SS,HSBIS,SSBIS)

implicit none

integer :: l,i,m,j,dime,nspin,lwork,info
integer,intent(in) :: nbin
double precision,intent(in) :: emin, emax
double precision :: dummyemax,dummyemin,E,step
double complex,intent(in) :: SD(:,:),HD(:,:,:), SS(:,:), HS(:,:,:)
double complex,allocatable :: SSTRAS(:,:), HSTRAS(:,:,:)
double complex,allocatable :: SSBISTRAS(:,:), HSBISTRAS(:,:,:)
double complex,allocatable :: BB(:,:),AA(:,:),work(:)
double complex,allocatable :: alpha(:),beta(:)
double complex :: VL(1,1)
double complex,allocatable :: VR(:,:)!,ppp(2,2),qqq(2,2),zzz(2),xxx(2)
double precision,allocatable :: rwork(:)
character(len=21) :: filename,filename2
double complex,intent(in) :: SSBIS(:,:), HSBIS(:,:,:)!, SSTRIS(:,:),HSTRIS(:,:,:)

dime=size(SD,dim=1)
nspin=size(HD,dim=3)

!Create conjugate transpose
allocate(SSTRAS(dime,dime),HSTRAS(dime,dime,nspin))
allocate(SSBISTRAS(dime,dime),HSBISTRAS(dime,dime,nspin))
SSTRAS(:,:)=0.d0
HSTRAS(:,:,:)=0.d0
SSBISTRAS(:,:)=0.d0
HSBISTRAS(:,:,:)=0.d0
do m=1,dime
  do l=1,dime
    SSTRAS(m,l)=CONJG(SS(l,m))
    SSBISTRAS(m,l)=CONJG(SSBIS(l,m))
    do i=1,nspin
     HSTRAS(m,l,i)=CONJG(HS(l,m,i))
     HSBISTRAS(m,l,i)=CONJG(HSBIS(l,m,i))
    enddo
  enddo
enddo


allocate(BB(dime*4,dime*4),AA(dime*4,dime*4))
allocate(alpha(dime*4),beta(dime*4))
lwork=4*4*dime-1
allocate(work(lwork),rwork(8*4*dime))
allocate(VR(dime*4,dime*4))

do i=1,nspin
  write(filename,'("spin",I1,"layer2")')i
  write(filename2,'("spin",I1,"layer2purereal")')i
  open(unit=i+100,file=filename,status="unknown")
  open(unit=i+200,file=filename2,status="unknown")
 
  !change energy in Ry
  dummyemax=emax/RyToEv
  dummyemin=emin/RyToEv
  step=(dummyemax-dummyemin)/nbin 
  do j=1,nbin+1
    alpha(:)=1.d0
    beta(:)=1.d0
    BB(:,:)=0.d0
    do l=1,3*dime
       BB(dime+l,dime+l)=1.d0
    enddo
    AA(:,:)=0.d0
    do l=1,3*dime
       AA(dime+l,l)=1.d0
    enddo
    E=dummyemin+(j-1)*step
    do m=1,dime
       do l=1,dime
            BB(m,l)=HSBISTRAS(m,l,i)-E*SSBISTRAS(m,l)
            BB(m,l+dime)=HSTRAS(m,l,i)-E*SSTRAS(m,l)
            AA(m,l+dime)=E*SD(m,l)-HD(m,l,i)
            AA(m,2*dime+l)=E*SS(m,l)-HS(m,l,i)
            AA(m,3*dime+l)=E*SSBIS(m,l)-HSBIS(m,l,i)
       enddo
    enddo
    call ZGGEV('N', 'V', 4*dime, AA, 4*dime, BB, 4*dime, alpha, beta, VL, 1, VR, 4*dime, work, lwork, rwork, info)
    if (info .ne. 0) then
         stop 'failed zggev diagonalisation'
         write(*,*) info
    end if
    do m=1,4*dime
         write(i+100,*) E*RyToEv, -real(real(log(alpha(m)/beta(m)))), real(aimag(log(alpha(m)/beta(m))))
         if (-real(real(log(alpha(m)/beta(m)))).le.0.001.and.-real(real(log(alpha(m)/beta(m)))).ge.-0.001) then
            write(i+200,*) E*RyToEv, real(aimag(log(alpha(m)/beta(m)))), -real(real(log(alpha(m)/beta(m))))
         endif
    enddo
  enddo
  close(i+100)
  close(i+200)
enddo 

deallocate(AA,BB,alpha,beta)
deallocate(work,rwork,VR)
deallocate(HSTRAS,SSTRAS,HSBISTRAS,SSBISTRAS)

end subroutine SecondLayerInteraction


subroutine ThirdLayerInteraction(emin,emax,nbin,HD,SD,HS,SS,HSBIS,SSBIS,HSTRIS,SSTRIS)

implicit none

integer :: l,i,m,j,dime,nspin,lwork,info
integer,intent(in) :: nbin
double precision,intent(in) :: emin, emax
double precision :: dummyemax,dummyemin,E,step
double complex,intent(in) :: SD(:,:),HD(:,:,:), SS(:,:), HS(:,:,:)
double complex,allocatable :: SSTRAS(:,:), HSTRAS(:,:,:)
double complex,allocatable :: SSBISTRAS(:,:), HSBISTRAS(:,:,:)
double complex,allocatable :: SSTRISTRAS(:,:), HSTRISTRAS(:,:,:)
double complex,allocatable :: BB(:,:),AA(:,:),work(:)
double complex,allocatable :: alpha(:),beta(:)
double complex :: VL(1,1)
double complex,allocatable :: VR(:,:)!,ppp(2,2),qqq(2,2),zzz(2),xxx(2)
double precision,allocatable :: rwork(:)
character(len=21) :: filename, filename2
double complex,intent(in) :: SSBIS(:,:), HSBIS(:,:,:), SSTRIS(:,:), HSTRIS(:,:,:)

dime=size(SD,dim=1)
nspin=size(HD,dim=3)

!Create conjugate transpose
allocate(SSTRAS(dime,dime),HSTRAS(dime,dime,nspin))
allocate(SSBISTRAS(dime,dime),HSBISTRAS(dime,dime,nspin))
allocate(SSTRISTRAS(dime,dime),HSTRISTRAS(dime,dime,nspin))
SSTRAS(:,:)=0.d0
HSTRAS(:,:,:)=0.d0
SSBISTRAS(:,:)=0.d0
HSBISTRAS(:,:,:)=0.d0
SSTRISTRAS(:,:)=0.d0
HSTRISTRAS(:,:,:)=0.d0
do m=1,dime
  do l=1,dime
    SSTRAS(m,l)=CONJG(SS(l,m))
    SSBISTRAS(m,l)=CONJG(SSBIS(l,m))
    SSTRISTRAS(m,l)=CONJG(SSTRIS(l,m))
    do i=1,nspin
     HSTRAS(m,l,i)=CONJG(HS(l,m,i))
     HSBISTRAS(m,l,i)=CONJG(HSBIS(l,m,i))
     HSTRISTRAS(m,l,i)=CONJG(HSTRIS(l,m,i))
    enddo
  enddo
enddo


allocate(BB(dime*6,dime*6),AA(dime*6,dime*6))
allocate(alpha(dime*6),beta(dime*6))
lwork=4*6*dime-1
allocate(work(lwork),rwork(8*6*dime))
allocate(VR(dime*6,dime*6))

do i=1,nspin
  write(filename,'("spin",I1,"layer3")')i
  write(filename2,'("spin",I1,"layer3purereal")')i
  open(unit=i+100,file=filename,status="unknown")
  open(unit=i+200,file=filename2,status="unknown")
 
  !change energy in Ry
  dummyemax=emax/RyToEv
  dummyemin=emin/RyToEv
  step=(dummyemax-dummyemin)/nbin 
  do j=1,nbin+1
    alpha(:)=1.d0
    beta(:)=1.d0
    BB(:,:)=0.d0
    do l=1,5*dime
       BB(dime+l,dime+l)=1.d0
    enddo
    AA(:,:)=0.d0
    do l=1,5*dime
       AA(dime+l,l)=1.d0
    enddo
    E=dummyemin+(j-1)*step
    do m=1,dime
       do l=1,dime
            BB(m,l)=HSTRISTRAS(m,l,i)-E*SSTRISTRAS(m,l)
            BB(m,dime+l)=HSBISTRAS(m,l,i)-E*SSBISTRAS(m,l)
            BB(m,l+2*dime)=HSTRAS(m,l,i)-E*SSTRAS(m,l)
            AA(m,l+dime*2)=E*SD(m,l)-HD(m,l,i)
            AA(m,3*dime+l)=E*SS(m,l)-HS(m,l,i)
            AA(m,4*dime+l)=E*SSBIS(m,l)-HSBIS(m,l,i)
            AA(m,5*dime+l)=E*SSTRIS(m,l)-HSTRIS(m,l,i)
       enddo
    enddo
    call ZGGEV('N', 'V', 6*dime, AA, 6*dime, BB, 6*dime, alpha, beta, VL, 1, VR, 6*dime, work, lwork, rwork, info)
    if (info .ne. 0) then
         stop 'failed zggev diagonalisation'
         write(*,*) info
    end if
    do m=1,6*dime
         write(i+100,*) E*RyToEv, -real(real(log(alpha(m)/beta(m)))), real(aimag(log(alpha(m)/beta(m))))
         if (-real(real(log(alpha(m)/beta(m)))).le.0.001.and.-real(real(log(alpha(m)/beta(m)))).ge.-0.001) then
            write(i+200,*) E*RyToEv, real(aimag(log(alpha(m)/beta(m)))), -real(real(log(alpha(m)/beta(m))))
         endif
    enddo
  enddo
  close(i+100)
  close(i+200)
enddo 

deallocate(AA,BB,alpha,beta)
deallocate(work,rwork,VR)
deallocate(HSTRAS,SSTRAS,HSBISTRAS,SSBISTRAS,HSTRISTRAS,SSTRISTRAS)

end subroutine ThirdLayerInteraction



subroutine KFirstLayerInteraction(emax,HD,SD,HS,SS,p1,p2,SymSetting)

implicit none

logical,intent(in) :: SymSetting
integer :: l,i,m,j,dime,nspin,lwork,info,ll
double precision,intent(in) :: emax,p1,p2
double precision :: dummyemax,dummyemin,E,step,posmin,reposmin
double complex,intent(in) :: SD(:,:),HD(:,:,:), SS(:,:), HS(:,:,:)
double complex,allocatable :: SSTRAS(:,:), HSTRAS(:,:,:)
double complex,allocatable :: BB(:,:),AA(:,:),work(:)
double complex,allocatable :: alpha(:),beta(:)
double complex :: VL(1,1),llproj
double complex,allocatable :: VR(:,:)!,ppp(2,2),qqq(2,2),zzz(2),xxx(2)
double precision,allocatable :: rwork(:)
character(len=12) :: filename,filename2
!complex,intent(in) :: SSBIS(:,:), HSBIS(:,:,:), SSTRIS(:,:),HSTRIS(:,:,:)


dime=size(SD,dim=1)
nspin=size(HD,dim=3)

!Create conjugate transpose
allocate(SSTRAS(dime,dime),HSTRAS(dime,dime,nspin))
SSTRAS(:,:)=0.d0
HSTRAS(:,:,:)=0.d0
do m=1,dime
  do l=1,dime
    SSTRAS(m,l)=CONJG(SS(l,m))
    do i=1,nspin
     HSTRAS(m,l,i)=CONJG(HS(l,m,i))
    enddo
  enddo
enddo


allocate(BB(dime*2,dime*2),AA(dime*2,dime*2))
allocate(alpha(dime*2),beta(dime*2))
lwork=4*2*dime-1
allocate(work(lwork),rwork(8*2*dime))
allocate(VR(dime*2,dime*2))


do i=1,nspin
  if (SymSetting) then
        write(filename,'("ProjSpin",I1,"layer1")')i
        open(unit=i+100,file=filename,Access = 'append',status="unknown")
        write(i+100,*) "#Orbitals has the same label of ORB_INDEX file"
  else
        write(filename,'("KSpin",I1,"layer1")')i
        open(unit=i+100,file=filename,Access = 'append',status="unknown")
        write(filename2,'("KMinSpin",I1,"layer1")')i
        open(unit=i+300,file=filename2,Access = 'append',status="unknown")
  endif
  E=emax/RyToEv
  alpha(:)=1.d0
  beta(:)=1.d0
  BB(:,:)=0.d0
  do l=1,dime
     BB(dime+l,dime+l)=1.d0
  enddo
  AA(:,:)=0.d0
  do l=1,dime
     AA(dime+l,l)=1.d0
  enddo
  do m=1,dime
     do l=1,dime
          BB(m,l)=HSTRAS(m,l,i)-E*SSTRAS(m,l)
          AA(m,l)=E*SD(m,l)-HD(m,l,i)
          AA(m,dime+l)=E*SS(m,l)-HS(m,l,i)
     enddo
  enddo
  call ZGGEV('N', 'V', 2*dime, AA, 2*dime, BB, 2*dime, alpha, beta, VL, 1, VR, dime*2, work, lwork, rwork, info)
  if (info .ne. 0) then
       stop 'failed zggev diagonalisation'
       write(*,*) info
  end if
  if (SymSetting) then
    do m=1,2*dime
        write(i+100,*) "eigenvalue",-real(real(log(alpha(m)/beta(m)))),real(aimag(log(alpha(m)/beta(m))))
        !do l=1,6*dime
        !   write(*,*)  abs(real(real(VR(l,m))))+abs(real(aimag(VR(l,m)))),VR(l,m)
        !enddo
        do ll=1,dime
          llproj=0.d0
          do l=1,dime
            llproj=llproj+VR(l,m)*SD(ll,l)+VR(l+1*dime,m)*SS(ll,l)
          enddo
          write(i+100,*) ll,llproj,(real(real(llproj)))**2+(real(aimag(llproj)))**2
        enddo
    enddo
  else
   posmin=300
   reposmin=0.d0
    do m=1,2*dime
       write(i+100,*) p1,p2, -real(real(log(alpha(m)/beta(m)))), real(aimag(log(alpha(m)/beta(m))))
       if (-real(real(log(alpha(m)/beta(m)))).gt.0.d0.and.-real(real(log(alpha(m)/beta(m)))).lt.posmin) then
          posmin=-real(real(log(alpha(m)/beta(m))))
          reposmin=real(aimag(log(alpha(m)/beta(m))))
       endif
    enddo
    if (p2.eq.0.5) then
     write(i+300,'(2f16.8,2f16.8)') p1,p2, posmin, reposmin
     write(i+300,*)
    else
     write(i+300,'(2f16.8,2f16.8)') p1,p2, posmin, reposmin
    endif
   close(i+300)
  endif
  close(i+100)
enddo 

deallocate(AA,BB,alpha,beta)
deallocate(work,rwork,VR)
deallocate(HSTRAS,SSTRAS)

end subroutine KFirstLayerInteraction



subroutine KSecondLayerInteraction(emax,HD,SD,HS,SS,HSBIS,SSBIS,p1,p2,SymSetting)

implicit none

logical,intent(in) :: SymSetting
integer :: l,i,m,j,dime,nspin,lwork,info,ll
double precision,intent(in) :: emax,p1,p2
double precision :: dummyemax,dummyemin,E,step,posmin,reposmin
double complex,intent(in) :: SD(:,:),HD(:,:,:), SS(:,:), HS(:,:,:)
double complex,allocatable :: SSTRAS(:,:), HSTRAS(:,:,:)
double complex,allocatable :: SSBISTRAS(:,:), HSBISTRAS(:,:,:)
double complex,allocatable :: BB(:,:),AA(:,:),work(:)
double complex,allocatable :: alpha(:),beta(:)
double complex :: VL(1,1),llproj
double complex,allocatable :: VR(:,:)!,ppp(2,2),qqq(2,2),zzz(2),xxx(2)
double precision,allocatable :: rwork(:)
character(len=12) :: filename, filename2
double complex,intent(in) :: SSBIS(:,:), HSBIS(:,:,:)!, SSTRIS(:,:),HSTRIS(:,:,:)

dime=size(SD,dim=1)
nspin=size(HD,dim=3)

!Create conjugate transpose
allocate(SSTRAS(dime,dime),HSTRAS(dime,dime,nspin))
allocate(SSBISTRAS(dime,dime),HSBISTRAS(dime,dime,nspin))
SSTRAS(:,:)=0.d0
HSTRAS(:,:,:)=0.d0
SSBISTRAS(:,:)=0.d0
HSBISTRAS(:,:,:)=0.d0
do m=1,dime
  do l=1,dime
    SSTRAS(m,l)=CONJG(SS(l,m))
    SSBISTRAS(m,l)=CONJG(SSBIS(l,m))
    do i=1,nspin
     HSTRAS(m,l,i)=CONJG(HS(l,m,i))
     HSBISTRAS(m,l,i)=CONJG(HSBIS(l,m,i))
    enddo
  enddo
enddo


allocate(BB(dime*4,dime*4),AA(dime*4,dime*4))
allocate(alpha(dime*4),beta(dime*4))
lwork=4*4*dime-1
allocate(work(lwork),rwork(8*4*dime))
allocate(VR(dime*4,dime*4))

do i=1,nspin
 if (SymSetting) then
        write(filename,'("ProjSpin",I1,"layer2")')i
        open(unit=i+100,file=filename,Access = 'append',status="unknown")
        write(i+100,*) "#Orbitals has the same label of ORB_INDEX file"
  else
        write(filename,'("KSpin",I1,"layer2")')i
        open(unit=i+100,file=filename,Access = 'append',status="unknown")
        write(filename2,'("KMinSpin",I1,"layer2")')i
        open(unit=i+300,file=filename2,Access = 'append',status="unknown")
  endif
  E=emax/RyToEv
  alpha(:)=1.d0
  beta(:)=1.d0
  BB(:,:)=0.d0
  do l=1,3*dime
     BB(dime+l,dime+l)=1.d0
  enddo
  AA(:,:)=0.d0
  do l=1,3*dime
     AA(dime+l,l)=1.d0
  enddo
  do m=1,dime
     do l=1,dime
          BB(m,l)=HSBISTRAS(m,l,i)-E*SSBISTRAS(m,l)
          BB(m,l+dime)=HSTRAS(m,l,i)-E*SSTRAS(m,l)
          AA(m,l+dime)=E*SD(m,l)-HD(m,l,i)
          AA(m,2*dime+l)=E*SS(m,l)-HS(m,l,i)
          AA(m,3*dime+l)=E*SSBIS(m,l)-HSBIS(m,l,i)
     enddo
  enddo
  call ZGGEV('N', 'V', 4*dime, AA, 4*dime, BB, 4*dime, alpha, beta, VL, 1, VR, 4*dime, work, lwork, rwork, info)
  if (info .ne. 0) then
       stop 'failed zggev diagonalisation'
       write(*,*) info
  end if
  if (SymSetting) then
    do m=1,4*dime
        write(i+100,*) "eigenvalue", -real(real(log(alpha(m)/beta(m)))),real(aimag(log(alpha(m)/beta(m))))
        !do l=1,6*dime
        !   write(*,*) abs(real(real(VR(l,m))))+abs(real(aimag(VR(l,m)))),VR(l,m)
        !enddo
        do ll=1,dime
          llproj=0.d0
          do l=1,dime
            llproj=llproj+VR(l,m)*SSTRAS(ll,l)+VR(l+1*dime,m)*SD(ll,l)+&
             +VR(l+2*dime,m)*SS(ll,l)+VR(l+3*dime,m)*SSBIS(ll,l)
          enddo
          write(i+100,*) ll,llproj,(real(real(llproj)))**2+(real(aimag(llproj)))**2
        enddo
    enddo
  else
   posmin=300
   reposmin=0.d0
    do m=1,4*dime
       write(i+100,*) p1,p2, -real(real(log(alpha(m)/beta(m)))),  real(aimag(log(alpha(m)/beta(m)))) 
       if (-real(real(log(alpha(m)/beta(m)))).gt.0.d0.and.-real(real(log(alpha(m)/beta(m)))).lt.posmin) then
          posmin=-real(real(log(alpha(m)/beta(m))))
          reposmin=real(aimag(log(alpha(m)/beta(m))))
       endif
    enddo
    if (p2.eq.0.5) then
     write(i+300,'(2f16.8,2f16.8)') p1,p2, posmin, reposmin
     write(i+300,*)
    else
     write(i+300,'(2f16.8,2f16.8)') p1,p2, posmin, reposmin
    endif
   close(i+300)
  endif
  close(i+100)
enddo 

deallocate(AA,BB,alpha,beta)
deallocate(work,rwork,VR)
deallocate(HSTRAS,SSTRAS,HSBISTRAS,SSBISTRAS)

end subroutine KSecondLayerInteraction


subroutine KThirdLayerInteraction(emax,HD,SD,HS,SS,HSBIS,SSBIS,HSTRIS,SSTRIS,p1,p2,SymSetting)

implicit none

logical,intent(in) :: SymSetting
integer :: l,i,m,j,dime,nspin,lwork,info,ll
double precision,intent(in) :: emax,p1,p2
double precision :: dummyemax,dummyemin,E,step,posmin,reposmin
double complex,intent(in) :: SD(:,:),HD(:,:,:), SS(:,:), HS(:,:,:)
double complex,allocatable :: SSTRAS(:,:), HSTRAS(:,:,:)
double complex,allocatable :: SSBISTRAS(:,:), HSBISTRAS(:,:,:)
double complex,allocatable :: SSTRISTRAS(:,:), HSTRISTRAS(:,:,:)
double complex,allocatable :: BB(:,:),AA(:,:),work(:)
double complex,allocatable :: alpha(:),beta(:)
double complex :: VL(1,1), llproj
double complex,allocatable :: VR(:,:)!,ppp(2,2),qqq(2,2),zzz(2),xxx(2)
double precision,allocatable :: rwork(:)
character(len=21) :: filename, filename2
double complex,intent(in) :: SSBIS(:,:), HSBIS(:,:,:), SSTRIS(:,:), HSTRIS(:,:,:)

dime=size(SD,dim=1)
nspin=size(HD,dim=3)

!Create conjugate transpose
allocate(SSTRAS(dime,dime),HSTRAS(dime,dime,nspin))
allocate(SSBISTRAS(dime,dime),HSBISTRAS(dime,dime,nspin))
allocate(SSTRISTRAS(dime,dime),HSTRISTRAS(dime,dime,nspin))
SSTRAS(:,:)=0.d0
HSTRAS(:,:,:)=0.d0
SSBISTRAS(:,:)=0.d0
HSBISTRAS(:,:,:)=0.d0
SSTRISTRAS(:,:)=0.d0
HSTRISTRAS(:,:,:)=0.d0
do m=1,dime
  do l=1,dime
    SSTRAS(m,l)=CONJG(SS(l,m))
    SSBISTRAS(m,l)=CONJG(SSBIS(l,m))
    SSTRISTRAS(m,l)=CONJG(SSTRIS(l,m))
    do i=1,nspin
     HSTRAS(m,l,i)=CONJG(HS(l,m,i))
     HSBISTRAS(m,l,i)=CONJG(HSBIS(l,m,i))
     HSTRISTRAS(m,l,i)=CONJG(HSTRIS(l,m,i))
    enddo
  enddo
enddo


allocate(BB(dime*6,dime*6),AA(dime*6,dime*6))
allocate(alpha(dime*6),beta(dime*6))
lwork=4*6*dime-1
allocate(work(lwork),rwork(8*6*dime))
allocate(VR(dime*6,dime*6))

do i=1,nspin
  if (SymSetting) then
        write(filename,'("ProjSpin",I1,"layer3")')i
        open(unit=i+100,file=filename,Access = 'append',status="unknown")
        write(i+100,*) "#Orbitals has the same label of ORB_INDEX file"
  else
        write(filename,'("KSpin",I1,"layer3")')i
        open(unit=i+100,file=filename,Access = 'append',status="unknown")
        write(filename2,'("KMinSpin",I1,"layer3")')i
        open(unit=i+300,file=filename2,Access = 'append',status="unknown")
  endif
  E=emax/RyToEv
  alpha(:)=1.d0
  beta(:)=1.d0
  BB(:,:)=0.d0
  do l=1,5*dime
     BB(dime+l,dime+l)=1.d0
  enddo
  AA(:,:)=0.d0
  do l=1,5*dime
     AA(dime+l,l)=1.d0
  enddo
  do m=1,dime
     do l=1,dime
          BB(m,l)=HSTRISTRAS(m,l,i)-E*SSTRISTRAS(m,l)
          BB(m,dime+l)=HSBISTRAS(m,l,i)-E*SSBISTRAS(m,l)
          BB(m,l+2*dime)=HSTRAS(m,l,i)-E*SSTRAS(m,l)
          AA(m,l+dime*2)=E*SD(m,l)-HD(m,l,i)
          AA(m,3*dime+l)=E*SS(m,l)-HS(m,l,i)
          AA(m,4*dime+l)=E*SSBIS(m,l)-HSBIS(m,l,i)
          AA(m,5*dime+l)=E*SSTRIS(m,l)-HSTRIS(m,l,i)
     enddo
  enddo
  call ZGGEV('N', 'V', 6*dime, AA, 6*dime, BB, 6*dime, alpha, beta, VL, 1, VR, 6*dime, work, lwork, rwork, info)
  if (info .ne. 0) then
       stop 'failed zggev diagonalisation'
       write(*,*) info
  end if
  if (SymSetting) then
    do m=1,6*dime
        write(i+100,*) "eigenvalue", -real(real(log(alpha(m)/beta(m)))), real(aimag(log(alpha(m)/beta(m))))
        !do l=1,6*dime
        !   write(*,*) abs(real(real(VR(l,m))))+abs(real(aimag(VR(l,m)))),VR(l,m)
        !enddo
        do ll=1,dime
          llproj=0.d0
          do l=1,dime
            llproj=llproj+VR(l,m)*SSBISTRAS(ll,l)+VR(l+dime,m)*SSTRAS(ll,l)+VR(l+2*dime,m)*SD(ll,l)+ &
             +VR(l+3*dime,m)*SS(ll,l)+VR(l+4*dime,m)*SSBIS(ll,l)+VR(l+5*dime,m)*SSTRIS(ll,l)
          enddo
          write(i+100,*) ll,llproj,(real(real(llproj)))**2+(real(aimag(llproj)))**2
        enddo
    enddo
  else
   posmin=300
   reposmin=0.d0
   do m=1,6*dime
       write(i+100,*) p1,p2, -real(real(log(alpha(m)/beta(m)))), real(aimag(log(alpha(m)/beta(m))))
       if (-real(real(log(alpha(m)/beta(m)))).gt.0.d0.and.-real(real(log(alpha(m)/beta(m)))).lt.posmin) then
          posmin=-real(real(log(alpha(m)/beta(m))))
          reposmin=real(aimag(log(alpha(m)/beta(m))))
       endif
   enddo
   if (p2.eq.0.5) then
     write(i+300,'(2f16.8,2f16.8)') p1,p2, posmin, reposmin
     write(i+300,*)
   else
     write(i+300,'(2f16.8,2f16.8)') p1,p2, posmin, reposmin
   endif
   close(i+300)
  endif
  close(i+100)
enddo 

deallocate(AA,BB,alpha,beta)
deallocate(work,rwork,VR)
deallocate(HSTRAS,SSTRAS,HSBISTRAS,SSBISTRAS,HSTRISTRAS,SSTRISTRAS)

end subroutine KThirdLayerInteraction


end module calccbs

