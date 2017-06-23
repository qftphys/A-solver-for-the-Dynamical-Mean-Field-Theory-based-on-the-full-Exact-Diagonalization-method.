!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN
  USE SF_CONSTANTS,only:zero
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_SETUP
  implicit none
  private

  !>build sparse hamiltonian of the sector
  public  :: ed_buildH_c
  !
  !
  !> Related auxiliary routines:
  public  :: setup_Hv_sector
  public  :: delete_Hv_sector

  integer                      :: Hsector=0
  logical                      :: Hstatus=.false.
  type(sector_map)             :: H,Hup,Hdw


contains



  subroutine setup_Hv_sector(isector)
    integer                   :: isector
    Hsector=isector
    Hstatus=.true.
    call build_sector(isector,H)
  end subroutine setup_Hv_sector


  subroutine delete_Hv_sector()
    call delete_sector(Hsector,H)
    Hsector=0
    Hstatus=.false.
  end subroutine delete_Hv_sector


  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_c(Hmat)
    complex(8),dimension(:,:)              :: Hmat
    integer                                :: isector
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: dim,dimUp,dimDw
    integer                                :: i,iup,idw
    integer                                :: m,mup,mdw
    integer                                :: ishift,ishift_up,ishift_dw
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw
    !
    if(.not.Hstatus)stop "ed_buildH_c ERROR: Hsector NOT set"
    isector=Hsector
    !
    dim=getdim(isector)
    !
    call assert_shape(Hmat,[dim,dim],"ed_buildH_c","Hmat")
    !
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
          enddo
       enddo
    enddo
    !
    Hmat=zero
    !
    !-----------------------------------------------!
    states: do i=1,Dim
       m = H%map(i)
       impi = i
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "FED_HAMILTONIAN/Himp.f90"
       !
       !LOCAL INTERACTION
       include "FED_HAMILTONIAN/Hint.f90"
       !
       !BATH HAMILTONIAN
       include "FED_HAMILTONIAN/Hbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "FED_HAMILTONIAN/Himp_bath.f90"
       !
       !
    enddo states
    !-----------------------------------------------!
    !
  end subroutine ed_buildH_c


end MODULE ED_HAMILTONIAN
