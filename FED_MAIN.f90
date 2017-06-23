module ED_MAIN
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG
  USE SF_LINALG
  USE SF_IOTOOLS, only: str
  USE SF_TIMER,only: start_timer,stop_timer
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  !
  !>INIT ED SOLVER
  !
  interface ed_init_solver
     module procedure :: ed_init_solver_single
     module procedure :: ed_init_solver_lattice
#ifdef _MPI
     module procedure :: ed_init_solver_lattice_mpi
#endif
  end interface ed_init_solver
  !>
  public :: ed_init_solver


  !
  !> ED SOLVER
  !
  interface ed_solve
     module procedure :: ed_solve_single
     module procedure :: ed_solve_lattice
#ifdef _MPI
     module procedure :: ed_solve_lattice_mpi
#endif
  end interface ed_solve
  !>
  public :: ed_solve



  real(8),dimension(:),allocatable                   :: wr,wm
  character(len=64)                                  :: suffix



contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_single(bath,Hloc)
    real(8),dimension(:),intent(inout) :: bath
    complex(8),intent(in)              :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    logical                            :: MPI_MASTER=.true.
    integer                            :: MPI_ERR
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    !Init bath:
    call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath(dmft_bath)
    call init_dmft_bath(dmft_bath)
    call get_dmft_bath(dmft_bath,bath)
    !
    if(isetup)call setup_pointers_normal
    call deallocate_dmft_bath(dmft_bath)
    isetup=.false.
    !
  end subroutine ed_init_solver_single


  !+-----------------------------------------------------------------------------+!
  !                           INEQUVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_lattice(bath,Hloc)
    real(8),dimension(:,:)         :: bath ![Nlat][:]
    complex(8),intent(in)          :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer                        :: ilat,Nineq,Nsect
    character(len=5)               :: tmp_suffix
    !
    Nineq = size(bath,1)
    do ilat=1,Nineq             !all nodes check the bath, u never know...
       ed_file_suffix="_site"//str(ilat,Npad=4)
       call ed_init_solver_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
    end do
    ed_file_suffix=""
  end subroutine ed_init_solver_lattice

  !+-----------------------------------------------------------------------------+!
  !                           INEQUVALENT SITES - MPI                             !
  !+-----------------------------------------------------------------------------+!
#ifdef _MPI
  subroutine ed_init_solver_lattice_mpi(MpiComm,bath,Hloc)
    integer                        :: MpiComm
    real(8),dimension(:,:)         :: bath ![Nlat][:]
    complex(8),intent(in)          :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer                        :: ilat,Nineq,Nsect
    logical                        :: check_dim
    character(len=5)               :: tmp_suffix
    integer                        :: MPI_ERR
    !
    Nineq = size(bath,1)
    do ilat=1,Nineq             !all nodes check the bath, u never know...
       ed_file_suffix="_site"//str(ilat,Npad=4)
       call ed_init_solver_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
    end do
    call MPI_Barrier(MpiComm,MPI_ERR)
    ed_file_suffix=""
    !
  end subroutine ed_init_solver_lattice_mpi
#endif




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single(bath,Hloc)
    real(8),dimension(:),intent(in) :: bath
    complex(8),optional,intent(in)  :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                         :: check
    !
    if(present(Hloc))call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    call write_dmft_bath(dmft_bath,LOGfile)
    call save_dmft_bath(dmft_bath,used=.true.)
    !
    !
    call setup_eigenspace()
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()         !find target states by digonalization of Hamiltonian
    call observables_impurity()         !obtain impurity observables as thermal averages.  
    call buildgf_impurity()             !build the one-particle impurity Green's functions  & Self-energy
    if(chiflag)call buildchi_impurity() !build the local susceptibilities (spin [todo charge])
    call local_energy_impurity()        !obtain the local energy of the effective impurity problem.
    !
    call delete_eigenspace()
    !    
    call deallocate_dmft_bath(dmft_bath)
    !
  end subroutine ed_solve_single




  !+-----------------------------------------------------------------------------+!
  !                           INEQUVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_lattice(bath,Hloc,Uloc_ii,Ust_ii,Jh_ii)
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    complex(8)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nsites
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(mii))deallocate(mii)
    if(allocated(eii))deallocate(eii)
    if(allocated(ddii))deallocate(ddii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(mii(Nsites,Norb))
    allocate(eii(Nsites,4))
    allocate(ddii(Nsites,4))
    !
    !Allocate the self-energies global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    allocate(Smatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Allocate the imp GF global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Gmatsii))deallocate(Gmatsii)
    if(allocated(Grealii))deallocate(Grealii)
    allocate(Gmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Grealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Check the dimensions of the bath are ok:
    do ilat=1,Nsites
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smatsii  = zero 
    Srealii  = zero 
    Gmatsii  = zero 
    Grealii  = zero 
    nii      = 0d0  
    dii      = 0d0  
    mii      = 0d0  
    eii      = 0d0  
    ddii     = 0d0  
    !
    call start_timer
    !
    do ilat = 1, Nsites
       write(LOGfile,'(A)')" solves site: "//str(ilat,Npad=4)
       !
       ed_file_suffix="_site"//str(ilat,Npad=4)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
       !
       Smatsii(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Srealii(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       Gmatsii(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
       Grealii(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
       nii(ilat,1:Norb)       = ed_dens(1:Norb)
       dii(ilat,1:Norb)       = ed_docc(1:Norb)
       mii(ilat,1:Norb)       = ed_dens_up(1:Norb)-ed_dens_dw(1:Norb)
       eii(ilat,:)            = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       ddii(ilat,:)           = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
    enddo
    !
    call stop_timer(LOGfile)
    !
    ed_file_suffix=""
    !
  end subroutine ed_solve_lattice



  !+-----------------------------------------------------------------------------+!
  !                           INEQUVALENT SITES - MPI                             !
  !+-----------------------------------------------------------------------------+!
#ifdef _MPI
  subroutine ed_solve_lattice_mpi(MpiComm,bath,Hloc,Uloc_ii,Ust_ii,Jh_ii)
    integer          :: MpiComm
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    complex(8)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    !MPI  auxiliary vars
    complex(8)       :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Gmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Greal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)          :: nii_tmp(size(bath,1),Norb)
    real(8)          :: dii_tmp(size(bath,1),Norb)
    real(8)          :: mii_tmp(size(bath,1),Norb)
    real(8)          :: eii_tmp(size(bath,1),4)
    real(8)          :: ddii_tmp(size(bath,1),4)
    !
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nsites
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !
    integer          :: MPI_ID=0
    integer          :: MPI_SIZE=1
    logical          :: MPI_MASTER=.true.
    !
    integer          :: mpi_err 
    !
    MPI_ID     = get_Rank_MPI(MpiComm)
    MPI_SIZE   = get_Size_MPI(MpiComm)
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(mii))deallocate(mii)
    if(allocated(eii))deallocate(eii)
    if(allocated(ddii))deallocate(ddii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(mii(Nsites,Norb))
    allocate(eii(Nsites,4))
    allocate(ddii(Nsites,4))
    !
    !Allocate the self-energies global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    allocate(Smatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Allocate the imp GF global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Gmatsii))deallocate(Gmatsii)
    if(allocated(Grealii))deallocate(Grealii)
    allocate(Gmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Grealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !
    !Check the dimensions of the bath are ok:
    do ilat=1+MPI_ID,Nsites,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smats_tmp  = zero
    Sreal_tmp  = zero
    Gmats_tmp  = zero
    Greal_tmp  = zero
    nii_tmp    = 0d0
    dii_tmp    = 0d0
    mii_tmp    = 0d0
    eii_tmp    = 0d0
    ddii_tmp   = 0d0
    !
    call start_timer
    !
    do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
       write(LOGfile,"(A)")str(MPI_ID)//" solves site: "//str(ilat,Npad=4)
       !
       ed_file_suffix="_site"//str(ilat,Npad=4)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
       !
       Smats_tmp(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Sreal_tmp(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       Gmats_tmp(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
       Greal_tmp(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)       = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)       = ed_docc(1:Norb)
       mii_tmp(ilat,1:Norb)       = ed_dens_up(1:Norb)-ed_dens_dw(1:Norb)
       eii_tmp(ilat,:)            = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       ddii_tmp(ilat,:)           = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
    enddo
    !
    call MPI_Barrier(MpiComm,MPI_ERR)
    !
    call stop_timer(LOGfile)
    !
    ed_file_suffix=""
    !
    Smatsii  = zero 
    Srealii  = zero 
    Gmatsii  = zero 
    Grealii  = zero 
    nii      = 0d0  
    dii      = 0d0  
    mii      = 0d0  
    eii      = 0d0  
    ddii     = 0d0  
    call MPI_ALLREDUCE(Smats_tmp,Smatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Sreal_tmp,Srealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Gmats_tmp,Gmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Greal_tmp,Grealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(mii_tmp,mii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(ddii_tmp,ddii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    !
  end subroutine ed_solve_lattice_mpi
#endif


end module ED_MAIN
