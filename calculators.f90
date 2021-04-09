module calculators

use readers, only: hsx_t

implicit none

double precision, parameter, private :: RyToEv = 13.605698066

type, public :: Hlayer
  double complex, pointer :: the_h(:,:,:) => null()
end type

type, public :: Slayer
  double complex, pointer :: the_s(:,:) => null()
end type

contains

subroutine CBSfixedK(hw,dir,isc,k,emin,emax,nbin,dist)
        
type(hsx_t), intent(in)  :: hw
integer, dimension(:,:), intent(in) :: isc
integer, intent(in) :: dir, nbin
double precision, intent(in) :: emin, emax, dist
double precision, intent(in) :: k(:)

type(Hlayer), allocatable :: hes(:)
type(Slayer), allocatable :: ses(:)
type(Hlayer), allocatable :: trans_hes(:)
type(Slayer), allocatable :: trans_ses(:)

complex :: i, phase
integer :: io,j,ij,jos,jo,layers,l,m,dime,nspin,lwork,info,spin,lay,integg,start_eig
integer, allocatable :: count_fail(:)
double precision :: dummyemax,dummyemin,E,step,thek
double complex,allocatable :: SSTRAS(:,:), HSTRAS(:,:,:)
double complex,allocatable :: BB(:,:),AA(:,:),work(:)
double complex,allocatable :: alpha(:),beta(:)
double complex :: VL(1,1)
double complex,allocatable :: VR(:,:)
double precision,allocatable :: rwork(:)
character(len=21) :: filename, filename2
   
   i=(0,1)

   ! First count how many layers in direction dir
   layers = 1
   do l=1,hw%no_s
     if (isc(l,dir).gt.layers) then
       layers=isc(l,dir)
     endif
   enddo

   ! Refering to theory hes(1) is HD, hes(2) is HS, hes(3) is HSBIS and so on..
   ! They are one more that the layers since also HD is included. Same for ses
   allocate(ses(layers+1))
   allocate(hes(layers+1))
   ! Allocate every matrix of ses and hes and initialize to zero
   do l=1,layers+1
     allocate(ses(l)%the_s(hw%no_u,hw%no_u),hes(l)%the_h(hw%no_u,hw%no_u,hw%nspin))
     ses(l)%the_s(1:hw%no_u,1:hw%no_u) = 0
     hes(l)%the_h(1:hw%no_u,1:hw%no_u,1:hw%nspin) = 0
   enddo

   do io = 1,hw%no_u                             ! loop on unit cell orbitals
     do j = 1,hw%numh(io)                        ! loop on connected orbitals
        ij = hw%listhptr(io)+j                   ! sparse-matrix array index
        jos = hw%listh(ij)                       ! index of connected orbital
        jo = hw%indxuo(jos)                      ! equiv. orbital in unit cell
        phase = exp(i*sum(k(:)*hw%xij(:,ij)))    !exp(i*(k(1)*hw%xij(1,ij)+k(2)*hw%xij(2,ij)))
        do l=1,layers+1
          if (isc(jos,dir).eq.-(l-1)) then
            hes(l)%the_h(io,jo,1:hw%nspin) = hes(l)%the_h(io,jo,1:hw%nspin) + phase*hw%hamilt(ij,1:hw%nspin)
            ses(l)%the_s(io,jo) = ses(l)%the_s(io,jo) + phase*hw%Sover(ij) ! overlap matrix element
          endif
        enddo
     enddo
   enddo

   dime = hw%no_u
   nspin = hw%nspin

   allocate(count_fail(nspin))

   !Create conjugate transpose, now we do not have the HD/SD
   allocate(trans_ses(layers))
   allocate(trans_hes(layers))
   ! Allocate every matrix and give value
   do lay=1,layers
     allocate(trans_ses(lay)%the_s(hw%no_u,hw%no_u),trans_hes(lay)%the_h(hw%no_u,hw%no_u,hw%nspin))
     do m=1,dime
       do l=1,dime
         trans_ses(lay)%the_s(m,l) = CONJG(ses(lay+1)%the_s(l,m))
         do spin=1,nspin
           trans_hes(lay)%the_h(m,l,spin) = CONJG(hes(lay+1)%the_h(l,m,spin))
         enddo
       enddo
     enddo
   enddo

   allocate(BB(dime*layers*2,dime*layers*2),AA(dime*layers*2,dime*layers*2))
   allocate(alpha(dime*layers*2),beta(dime*layers*2))

   lwork=4*2*layers*dime-1
   allocate(work(lwork),rwork(8*2*layers*dime))
   allocate(VR(dime*2*layers,dime*2*layers))
   
   
   do spin=1,nspin
     
     count_fail(spin) = 0

     write(filename,'("spin",I1,"layer",I1)')spin,layers
     write(filename2,'("spin",I1,"layer",I1,"purereal")')spin,layers
     open(unit=spin+100,file=filename,status="unknown")
     open(unit=spin+200,file=filename2,status="unknown") 

     !change energy in Ry
     dummyemax=emax/RyToEv
     dummyemin=emin/RyToEv
     step=(dummyemax-dummyemin)/nbin 
     do j=1,nbin+1
       alpha(:)=1.d0
       beta(:)=1.d0
       BB(:,:)=0.d0
       do l=1,(2*layers-1)*dime
          BB(dime+l,dime+l)=1.d0
       enddo
       AA(:,:)=0.d0
       do l=1,(2*layers-1)*dime
          AA(dime+l,l)=1.d0
       enddo
       E=dummyemin+(j-1)*step
       do m=1,dime
          do l=1,dime
            do integg=1,layers+1
               AA(m,(layers-2+integg)*dime+l) = E*ses(integg)%the_s(m,l)-hes(integg)%the_h(m,l,spin)
            enddo
            do integg=1,layers
               BB(m,l+dime*(integg-1)) = trans_hes(layers+1-integg)%the_h(m,l,spin) - E*trans_ses(layers+1-integg)%the_s(m,l)
            enddo
          enddo
       enddo
       start_eig = 1
       call ZGGEV('N', 'V', 2*layers*dime, AA, 2*layers*dime, BB, 2*layers*dime, &
         alpha, beta, VL, 1, VR, dime*2*layers, work, lwork, rwork, info)
       if (info .ne. 0) then
         !write(*,*) info
         start_eig = info+1
         count_fail(spin) = count_fail(spin)+1
         !stop 'failed zggev diagonalisation'
         !cycle
       end if
       do m=start_eig,2*layers*dime
         thek=-real(real(log(alpha(m)/beta(m))))/dist
         if (thek.gt.10000000) cycle
         if (thek.ne.thek) cycle
         write(spin+100,*) E*RyToEv, thek, real(aimag(log(alpha(m)/beta(m))))/dist
         if (thek.le.0.001.and.thek.ge.-0.001) then
           write(spin+200,*) E*RyToEv, real(aimag(log(alpha(m)/beta(m))))/dist, thek
         endif
       enddo
     enddo
     close(spin+100)
   enddo 

do jo=1,nspin
  write(*,*) "For spin", jo, ",", count_fail(jo), "bins out of", nbin, "could not produce all the eigenvalues."
enddo
if (any(count_fail>0)) then
  write(*,*) "Problem due to a failure of LAPACK routine ZGGEV, &
          probably linked to numerical instability of Schur decomposition within LAPACK."
endif 

deallocate(AA,BB,alpha,beta)
deallocate(work,rwork,VR)
deallocate(trans_ses,trans_hes,hes,ses,count_fail)

end subroutine CBSfixedK


end module calculators

