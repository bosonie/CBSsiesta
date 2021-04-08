module readers

implicit none

integer, parameter, private :: dp = selected_real_kind(14,100)
integer, parameter, private :: sp = selected_real_kind(6,30)

public  :: read_hsx_file
public  :: read_OrbIndx_file

! Set derived type hsx_t to hold info of HSX file, containing:
!   nspecies                : number of chemical species
!   na_u                    : number of atoms in unit cell
!   no_u                    : number of orbitals in unit cell
!   no_s                    : number of orbitals in supercell
!   nspin                   : number of spin components
!   nh                      : dimension of arrays hamilt, Sover, and xij
!   gamma                   : was this a gamma-only calculation?
!   has_xij                 : does the file contain xij vectors?
!   no(nspecies)            : number of atomic orbitals of each species
!   nquant(nspecies,naoatx) : principal quantum number of each atomic orbital
!                             with naoatx=max(no)
!   lquant(nspecies,naoatx) : ang. momentum number of each atomic orbital
!   zeta(nspecies,naoatx)   : zeta-index of each atomic orbital
!   iaorb(no_u)             : atom to which each orbital belongs
!   iphorb(no_u)            : index of each orbital within its atom
!   label(nspecies)         : atomic label (symbol) of each species
!   numh(no_u)              : num of nonzero elements in each row of hamiltonian
!   listhptr(no_u)          : row-start index in sparse-matrix arrays
!   listh(nh)               : orbital index of nonzero matrix elements
!   indxuo(no_s)            : index of equivalent orbital in first unit cell
!   hamilt(nh,nspin)        : hamiltonian matrix elements in sparse format
!   Sover(nh)               : overlap matrix elements in sparse format
!   xij(3,nh)               : vector between each pair of connected orbitals
!   isa(na_u)               : species index of each atom
!   zval(nspecies)          : atomic number of each species
! To transform from sparse to full format, for a given k point:
!   S(1:no_u,1:no_u) = 0                       ! full complex overlap matrix
!   H(1:no_u,1:no_u,1:nspin) = 0               ! full complex hamiltonian
!   do io = 1,no_u                             ! loop on unit cell orbitals
!     do j = 1,numh(io)                        ! loop on connected orbitals
!        ij = listhptr(io)+j                   ! sparse-matrix array index
!        jos = listh(ij)                       ! index of connected orbital
!        jo = indxuo(jos)                      ! equiv. orbital in unit cell
!        phase = exp(i*sum(k(:)*xij(:,ij)))    ! phase factor between orbs.
!        H(io,jo,1:nspin) = H(io,jo,1:nspin) + ! hamiltonian matrix element
!                           phase*hamilt(ij,1:spin)
!        S(io,jo) = S(io,jo) + phase*Sover(ij) ! overlap matrix element
!     enddo
!   enddo
! Notice that io,jo are within unit cell, and jos is within supercell

type, public :: hsx_t
  integer :: nspecies
  integer :: na_u
  integer :: no_u
  integer :: no_s
  integer :: nspin
  integer :: nh
  logical :: gamma
  logical :: has_xij = .false.
  integer, pointer :: no(:) => null()
  integer, pointer :: nquant(:,:) => null()
  integer, pointer :: lquant(:,:) => null()
  integer, pointer :: zeta(:,:) => null()
  integer, pointer :: iaorb(:) => null()
  integer, pointer :: iphorb(:) => null()
  character(len=20), pointer :: label(:) => null()
  integer, pointer  :: numh(:) => null() 
  integer, pointer  :: listhptr(:) => null()
  integer, pointer  :: listh(:) => null()  
  integer, pointer  :: indxuo(:) => null()
  real(dp), pointer :: hamilt(:,:) => null()
  real(dp), pointer :: Sover(:) => null()
  real(dp), pointer :: xij(:,:) => null()
  integer, pointer  :: isa(:) => null()
  real(dp), pointer :: zval(:) => null()
  real(dp)          :: qtot, temp           !  fossils
end type

private

CONTAINS


subroutine read_OrbIndx_file(isc,fname)
!integer, intent(in) :: tot !Total number of orbitals, meaning no_u
!integer, intent(out)  :: isc(tot,3)
integer, dimension(:,:), intent(out) :: isc
character(len=*), intent(in) :: fname

  integer :: hsx_u, tot, nn, totorb, io, ia, iuo, l, dumm
  character(len=10) :: a,b 
  double precision :: du

  tot = size(isc,dim=1)
  nn  = size(isc,dim=2)
  call get_unit_number(hsx_u)
  !print *, "Using unit: ", hsx_u

  open(hsx_u,file=trim(fname),status="old")
  read(hsx_u,*) dumm, totorb
  if (totorb.ne.tot) then
     print *, "Something wrong, total number of orbitals doesn't match"
     STOP
  endif
 ! allocate(isc(tot,3))

  read(hsx_u,*)
  read(hsx_u,*)
  do l=1,totorb
    read(hsx_u,*) io,ia,dumm,b,dumm,dumm,dumm,dumm,dumm,a,a,du,isc(l,:),iuo
  enddo
  close(hsx_u)

end subroutine read_OrbIndx_file


!--------------------------------------

subroutine read_hsx_file(hsx,fname)
type(hsx_t), intent(out)  :: hsx
character(len=*), intent(in) :: fname

!
! Reads HSX file "fname" and stores the info in the hsx data structure
! (Real arrays are stored in double precision)

  integer, allocatable  :: ibuff(:)
  real(sp), allocatable  :: hbuff(:)
  real(sp), allocatable  :: buff3(:,:)

  integer numx, ind, no_u, nnz, na_u, nspecies, nspin, nh, i
  integer :: im, is, hsx_u, ia, io, iostat, k, naoatx, no_s
  logical  :: debug = .false.

  call get_unit_number(hsx_u)
  !print *, "Using unit: ", hsx_u

  open(hsx_u,file=trim(fname),status='old',form='unformatted')

!read n of orbitals in unit-cell, supercell, n of spin components, dimension of H S ..
  read(hsx_u,iostat=iostat) hsx%no_u, hsx%no_s, hsx%nspin, hsx%nh
  if (iostat /= 0) STOP "nnao, no_s..."

  no_u = hsx%no_u
  no_s = hsx%no_s

  read(hsx_u,iostat=iostat) hsx%gamma !only gamma?
  if (iostat /= 0) STOP "gamma"
  IF (DEBUG) PRINT *, "GAMMA=", hsx%gamma
  if (.not. hsx%gamma) then
     allocate(hsx%indxuo(no_s))
     read(hsx_u) (hsx%indxuo(i),i=1,hsx%no_s)
  else
     allocate(hsx%indxuo(hsx%no_u))
     do i=1,hsx%no_u
        hsx%indxuo(i) = i
     enddo
  endif

  nh  = hsx%nh
  nspin = hsx%nspin
  !print *, "nh: ", nh
  allocate (hsx%numh(no_u), hsx%listhptr(no_u), hsx%listh(nh))

       allocate (hsx%xij(3,nh),stat=iostat)
       allocate (hsx%hamilt(nh,nspin),stat=iostat)
       allocate (hsx%Sover(nh),stat=iostat)

  read(hsx_u,iostat=iostat) (hsx%numh(io), io=1,no_u)      
  if (iostat /= 0) STOP "numh"

  numx = maxval(hsx%numh(1:no_u))
  allocate(ibuff(numx), hbuff(numx), buff3(3,numx))

  nnz = sum(hsx%numh(1:hsx%no_u))
   if (nnz > nh) STOP "nh overflow in HS"
  ! Create listhptr 
  hsx%listhptr(1)=0
  do io=2,hsx%no_u
     hsx%listhptr(io)=hsx%listhptr(io-1)+hsx%numh(io-1)
  enddo

  do io=1,hsx%no_u
     read(hsx_u,iostat=iostat) (ibuff(im), im=1,hsx%numh(io))
     if (iostat /= 0) STOP "listh"
     do im=1,hsx%numh(io)
        hsx%listh(hsx%listhptr(io)+im) = ibuff(im)
     enddo
  enddo

  do is=1,hsx%nspin
     do io=1,hsx%no_u
        read(hsx_u,iostat=iostat) (hbuff(im), im=1,hsx%numh(io))
        if (iostat /= 0) STOP "Hamilt"
        do im=1,hsx%numh(io)
           hsx%hamilt(hsx%listhptr(io)+im,is) = hbuff(im)
           if (debug) print *, "Hamilt ", io, im, hbuff(im)
        enddo
     enddo
  enddo
  !
  !       Read overlap matrix
  !
  do io=1,hsx%no_u
     read(hsx_u,iostat=iostat) (hbuff(im), im=1,hsx%numh(io))
     if (iostat /= 0) STOP "Overlap matrix read error"
     do im=1,hsx%numh(io)
        hsx%Sover(hsx%listhptr(io)+im) = hbuff(im)
        if (debug) print *, "S ", io, im, hbuff(im)
     enddo
  enddo

  read(hsx_u,iostat=iostat) hsx%qtot, hsx%temp           ! fossils
  if (iostat /= 0) STOP "Qtot, temp, read error"

  !
  !        Always read xijk
  !
  do io=1,hsx%no_u
     read(hsx_u,iostat=iostat) ((buff3(k,im), k=1,3), im=1,hsx%numh(io))
     if (iostat /= 0) STOP "xij(k) read error"
     do im=1,hsx%numh(io)
        ind = hsx%listhptr(io)+im
        if (debug) print *, "xijk ", buff3(:,im)
        hsx%xij(1:3,ind) = buff3(1:3,im)
     enddo
  enddo
  hsx%has_xij = .true.

  !
  !        Read auxiliary info
  !
  read(hsx_u) hsx%nspecies
  nspecies = hsx%nspecies
  !print *, "nspecies: ", nspecies
  allocate(hsx%label(nspecies), hsx%zval(nspecies), hsx%no(nspecies))
  read(hsx_u) (hsx%label(is),hsx%zval(is),hsx%no(is), is=1,nspecies)
  naoatx = maxval(hsx%no(1:nspecies))
  allocate (hsx%nquant(nspecies,naoatx), hsx%lquant(nspecies,naoatx), &
       hsx%zeta(nspecies,naoatx))
  do is=1, nspecies
     do io=1, hsx%no(is)
        read(hsx_u) hsx%nquant(is,io), hsx%lquant(is,io), hsx%zeta(is,io)
     enddo
  enddo
  read(hsx_u) hsx%na_u
  na_u = hsx%na_u
  allocate(hsx%isa(na_u))
  allocate(hsx%iaorb(no_u), hsx%iphorb(no_u))
  read(hsx_u) (hsx%isa(ia), ia=1,na_u)
  read(hsx_u) (hsx%iaorb(io), hsx%iphorb(io), io=1,no_u)

  close(hsx_u)
  deallocate(ibuff, hbuff, buff3)

end subroutine read_hsx_file

!--------------------------------------------------------------

subroutine get_unit_number(lun)
integer, intent(out) :: lun

logical :: used
integer :: iostat

do lun= 10, 99
   inquire(unit=lun, opened=used, iostat=iostat)
   if (iostat .ne. 0) used = .true.
   if (.not. used) return
enddo
STOP "Cannot get unit"
end subroutine get_unit_number



end module readers
