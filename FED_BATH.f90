MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  private



  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  !explicit symmetries:
  interface break_symmetry_bath
     module procedure break_symmetry_bath_site
     module procedure break_symmetry_bath_lattice
  end interface break_symmetry_bath
  !
  interface spin_symmetrize_bath
     module procedure spin_symmetrize_bath_site
     module procedure spin_symmetrize_bath_lattice
  end interface spin_symmetrize_bath
  !
  interface ph_symmetrize_bath
     module procedure ph_symmetrize_bath_site
     module procedure ph_symmetrize_bath_lattice
  end interface ph_symmetrize_bath
  !
  interface ph_trans_bath
     module procedure ph_trans_bath_site
     module procedure ph_trans_bath_lattice
  end interface ph_trans_bath
  !
  interface enforce_normal_bath
     module procedure enforce_normal_bath_site
     module procedure enforce_normal_bath_lattice
  end interface enforce_normal_bath


  public :: get_bath_dimension
  public :: check_bath_dimension
  public :: break_symmetry_bath
  public :: spin_symmetrize_bath
  public :: ph_symmetrize_bath
  public :: ph_trans_bath
  public :: enforce_normal_bath



  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  !PUBLIC   = transparent to the final user
  !INTERNAL = opaque to the user but available for internal use in the code.
  !
  !DMFT BATH procedures:
  public :: allocate_dmft_bath               !INTERNAL (for effective_bath)
  public :: deallocate_dmft_bath             !INTERNAL (for effective_bath)
  public :: init_dmft_bath                   !INTERNAL (for effective_bath)
  public :: write_dmft_bath                  !INTERNAL (for effective_bath)
  public :: save_dmft_bath                   !INTERNAL (for effective_bath)
  public :: set_dmft_bath                    !INTERNAL (for effective_bath)
  public :: get_dmft_bath                    !INTERNAL (for effective_bath)






contains


  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  !
  ! Get size of each dimension of the component array. 
  ! The Result is an rank 1 integer array Ndim with dimension:
  ! 3 for get_component_size_bath
  ! 2 for get_spin_component_size_bath & get_orb_component_size_bath
  ! 1 for get_spin_orb_component_size_bath
  !+-------------------------------------------------------------------+
  function get_bath_dimension(Hloc_nn,ispin_) result(bath_size)
    complex(8),optional,intent(in) :: Hloc_nn(:,:,:,:)
    integer,optional               :: ispin_
    integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo
    complex(8),allocatable         :: Hloc(:,:,:,:)
    !
    select case(bath_type)
    case default
       !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
       bath_size = Norb*Nbath + Norb*Nbath
       if(.not.present(ispin_))bath_size=Nspin*bath_size
    case('hybrid')
       !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
       bath_size = Nbath + Norb*Nbath
       if(.not.present(ispin_))bath_size=Nspin*bath_size
    end select
  end function get_bath_dimension



  !##################################################################
  !
  !     USER BATH PREDEFINED SYMMETRIES:
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : given a bath array apply a specific transformation or 
  ! impose a given symmetry:
  ! - break spin symmetry by applying a symmetry breaking field
  ! - given a bath array set both spin components to have 
  !    the same bath, i.e. impose non-magnetic solution
  ! - given a bath array enforces the particle-hole symmetry 
  !    by setting the positive energies in modulo identical to the negative
  !    ones.
  ! - given a bath enforce normal (i.e. non superconducting) solution
  ! - given a dmft bath pull/push the components W^{ss'}_\a(l) of the Hybridization 
  !    matrix
  !+-------------------------------------------------------------------+
  subroutine break_symmetry_bath_site(bath_,field,sign,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    real(8)                :: field
    real(8)                :: sign
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine break_symmetry_bath_site
  subroutine break_symmetry_bath_lattice(bath_,field,sign,save)
    real(8),dimension(:,:) :: bath_
    real(8)                :: field
    real(8)                :: sign
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix="_site"//reg(txtfy(ilat,Npad=4))
       call break_symmetry_bath_site(bath_(ilat,:),field,sign,save_)
    enddo
    ed_file_suffix=""
  end subroutine break_symmetry_bath_lattice

     !---------------------------------------------------------!

  subroutine spin_symmetrize_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    if(Nspin==1)then
       write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    !
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    !
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine spin_symmetrize_bath_site
  subroutine spin_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix="_site"//reg(txtfy(ilat,Npad=4))
       call spin_symmetrize_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine spin_symmetrize_bath_lattice

  !---------------------------------------------------------!

  subroutine ph_symmetrize_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    if(mod(Nbath,2)==0)then
       do i=1,Nbath/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
       enddo
    else
       do i=1,(Nbath-1)/2
          dmft_bath_%e(:,:,Nbath+1-i)=-dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,Nbath+1-i)= dmft_bath_%v(:,:,i)
       enddo
       dmft_bath_%e(:,:,(Nbath-1)/2+1)=0.d0
    endif
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine ph_symmetrize_bath_site
  subroutine ph_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix="_site"//reg(txtfy(ilat,Npad=4))
       call ph_symmetrize_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine ph_symmetrize_bath_lattice

  !---------------------------------------------------------!

  subroutine ph_trans_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    type(effective_bath)   :: tmp_dmft_bath
    integer                :: i
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call allocate_dmft_bath(tmp_dmft_bath)
    call set_dmft_bath(bath_,dmft_bath_)
    if(Nbath==1)return
    do i=1,Nbath
       select case(Norb)
       case default
          ! do nothing
          dmft_bath_%e(:,:,i)= dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)= dmft_bath_%v(:,:,i)
       case(1)
          dmft_bath_%e(:,:,i)= -dmft_bath_%e(:,:,i)
          dmft_bath_%v(:,:,i)=  dmft_bath_%v(:,:,i)
       case(2)
          tmp_dmft_bath%e(:,1,i) = -dmft_bath_%e(:,2,i)
          tmp_dmft_bath%e(:,2,i) = -dmft_bath_%e(:,1,i)
          dmft_bath_%e(:,:,i)    = tmp_dmft_bath%e(:,:,i)
          tmp_dmft_bath%v(:,1,i) = dmft_bath_%v(:,2,i)
          tmp_dmft_bath%v(:,2,i) = dmft_bath_%v(:,1,i)
          dmft_bath_%v(:,:,i)    = tmp_dmft_bath%v(:,:,i)
       end select
    end do
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine ph_trans_bath_site
  subroutine ph_trans_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix="_site"//reg(txtfy(ilat,Npad=4))
       call ph_trans_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine ph_trans_bath_lattice

  !---------------------------------------------------------!

  subroutine enforce_normal_bath_site(bath_,save)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    logical,optional       :: save
    logical                :: save_
    save_=.true.;if(present(save))save_=save
    call allocate_dmft_bath(dmft_bath_)
    call set_dmft_bath(bath_,dmft_bath_)
    if(save_)call save_dmft_bath(dmft_bath_)
    call get_dmft_bath(dmft_bath_,bath_)
    call deallocate_dmft_bath(dmft_bath_)
  end subroutine enforce_normal_bath_site
  subroutine enforce_normal_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix="_site"//reg(txtfy(ilat,Npad=4))
       call enforce_normal_bath_site(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine enforce_normal_bath_lattice





  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_dmft_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(dmft_bath_%status)call deallocate_dmft_bath(dmft_bath_)
    !
    select case(bath_type)
    case default
       !
       allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
       allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
       !
    case('hybrid')
       !
       allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
       allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
       !
    end select
    dmft_bath_%status=.true.
  end subroutine allocate_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_dmft_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(allocated(dmft_bath_%e))   deallocate(dmft_bath_%e)
    if(allocated(dmft_bath_%v))   deallocate(dmft_bath_%v)
    dmft_bath_%status=.false.
  end subroutine deallocate_dmft_bath




  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_dmft_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    integer              :: i,unit,flen,Nh
    integer              :: io,jo,iorb,ispin,jorb,jspin
    logical              :: IOfile
    real(8)              :: de,noise_tot
    character(len=21)    :: space
    if(.not.dmft_bath_%status)stop "init_dmft_bath error: bath not allocated"
    !Get energies:
    dmft_bath_%e(:,:,1)    =-hwband
    dmft_bath_%e(:,:,Nbath)= hwband
    Nh=Nbath/2
    if(mod(Nbath,2)==0.and.Nbath>=4)then
       de=hwband/max(Nh-1,1)
       dmft_bath_%e(:,:,Nh)  = -1.d-3
       dmft_bath_%e(:,:,Nh+1)=  1.d-3
       do i=2,Nh-1
          dmft_bath_%e(:,:,i)   =-hwband + (i-1)*de
          dmft_bath_%e(:,:,Nbath-i+1)= hwband - (i-1)*de
       enddo
    elseif(mod(Nbath,2)/=0.and.Nbath>=3)then
       de=hwband/Nh
       dmft_bath_%e(:,:,Nh+1)= 0.0d0
       do i=2,Nh
          dmft_bath_%e(:,:,i)        =-hwband + (i-1)*de
          dmft_bath_%e(:,:,Nbath-i+1)= hwband - (i-1)*de
       enddo
    endif
    !Get spin-keep yhbridizations
    do i=1,Nbath
       dmft_bath_%v(:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath)))
    enddo
    !
    !Read from file if exist:
    !
    inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")'Reading bath from file'//trim(Hfile)//trim(ed_file_suffix)//".restart"
       unit = free_unit()
       flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
       !
       open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
       !
       select case(bath_type)
       case default
          !
          read(unit,*)
          do i=1,min(flen,Nbath)
             read(unit,*)((&
                  dmft_bath_%e(ispin,iorb,i),&
                  dmft_bath_%v(ispin,iorb,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
          !
       case ('hybrid')
          read(unit,*)
          !
          do i=1,min(flen,Nbath)
             read(unit,*)(&
                  dmft_bath_%e(ispin,1,i),&
                  (&
                  dmft_bath_%v(ispin,iorb,i),&
                  iorb=1,Norb),&
                  ispin=1,Nspin)
          enddo
          !
       end select
       close(unit)
    endif
  end subroutine init_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_dmft_bath(dmft_bath_,unit)
    type(effective_bath) :: dmft_bath_
    integer,optional     :: unit
    integer              :: unit_
    integer              :: i
    integer              :: io,jo,iorb,ispin
    complex(8)           :: hybr_aux
    complex(8)           :: hrep_aux(Nspin*Norb,Nspin*Norb)
    unit_=LOGfile;if(present(unit))unit_=unit
    if(.not.dmft_bath_%status)stop "write_dmft_bath error: bath not allocated"
    select case(bath_type)
    case default
       !
       write(unit_,"(90(A21,1X))")&
            ((&
            "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
            "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
            iorb=1,Norb),ispin=1,Nspin)
       do i=1,Nbath
          write(unit_,"(90(F21.12,1X))")((&
               dmft_bath_%e(ispin,iorb,i),&
               dmft_bath_%v(ispin,iorb,i),&
               iorb=1,Norb),ispin=1,Nspin)
       enddo
       !
    case('hybrid')
       !
       write(unit_,"(90(A21,1X))")(&
            "#Ek_s"//reg(txtfy(ispin)),&
            ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),&
            ispin=1,Nspin)
       do i=1,Nbath
          write(unit_,"(90(F21.12,1X))")(&
               dmft_bath_%e(ispin,1,i),&
               (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
               ispin=1,Nspin)
       enddo
       !
    end select
  end subroutine write_dmft_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : save the bath to a given file using the write bath
  ! procedure and formatting: 
  !+-------------------------------------------------------------------+
  subroutine save_dmft_bath(dmft_bath_,file,used)
    type(effective_bath)      :: dmft_bath_
    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_
    character(len=16)         :: extension
    integer                   :: unit_
    if(.not.dmft_bath_%status)stop "save_dmft_bath error: bath is not allocated"
    used_=.false.;if(present(used))used_=used
    extension=".restart";if(used_)extension=".used"
    file_=reg(reg(Hfile)//reg(ed_file_suffix)//reg(extension))
    if(present(file))file_=reg(file)
    unit_=free_unit()
    open(unit_,file=reg(file_))
    call write_dmft_bath(dmft_bath_,unit_)
    close(unit_)
  end subroutine save_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_dmft_bath(bath_,dmft_bath_)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath
    logical                :: check
    complex(8)             :: hrep_aux(Nspin*Norb,Nspin*Norb)
    complex(8)             :: U(Nspin*Norb,Nspin*Norb)
    complex(8)             :: Udag(Nspin*Norb,Nspin*Norb)
    real(8)                :: element_R,element_I,eps_k,lambda_k
    if(.not.dmft_bath_%status)stop "set_dmft_bath error: bath not allocated"
    check = check_bath_dimension(bath_)
    if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       !
       stride = 0
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Nbath
                io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                dmft_bath_%e(ispin,iorb,i) = bath_(io)
             enddo
          enddo
       enddo
       stride = Nspin*Norb*Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Nbath
                io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                dmft_bath_%v(ispin,iorb,i) = bath_(io)
             enddo
          enddo
       enddo
       !
       !
    case ('hybrid')
       !
       stride = 0
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             dmft_bath_%e(ispin,1,i) = bath_(io)
          enddo
       enddo
       stride = Nspin*Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Nbath
                io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                dmft_bath_%v(ispin,iorb,i) = bath_(io)
             enddo
          enddo
       enddo
       !
       !
    end select
  end subroutine set_dmft_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  subroutine get_dmft_bath(dmft_bath_,bath_)
    type(effective_bath)   :: dmft_bath_
    real(8),dimension(:)   :: bath_
    complex(8)             :: hrep_aux(Nspin*Norb,Nspin*Norb)
    complex(8)             :: U(Nspin*Norb,Nspin*Norb)
    complex(8)             :: Udag(Nspin*Norb,Nspin*Norb)
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath
    logical                :: check
    if(.not.dmft_bath_%status)stop "get_dmft_bath error: bath not allocated"
    check=check_bath_dimension(bath_)
    if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       !
       stride = 0
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Nbath
                io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                bath_(io) = dmft_bath_%e(ispin,iorb,i) 
             enddo
          enddo
       enddo
       stride = Nspin*Norb*Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Nbath
                io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                bath_(io) = dmft_bath_%v(ispin,iorb,i)
             enddo
          enddo
       enddo
       !
    case ('hybrid')
       !
       stride = 0
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             bath_(io) =  dmft_bath_%e(ispin,1,i)
          enddo
       enddo
       stride = Nspin*Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,Nbath
                io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                bath_(io) =  dmft_bath_%v(ispin,iorb,i)
             enddo
          enddo
       enddo
       !
    end select
  end subroutine get_dmft_bath







  !##################################################################
  !
  !     USER BATH CHECKS:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_,Hloc_nn) result(bool)
    real(8),dimension(:)           :: bath_
    integer                        :: Ntrue
    logical                        :: bool
    complex(8),optional,intent(in) :: Hloc_nn(:,:,:,:)
    if (present(Hloc_nn))then
       Ntrue = get_bath_dimension(Hloc_nn)
    else
       Ntrue = get_bath_dimension()
    endif
    bool  = ( size(bath_) == Ntrue )
  end function check_bath_dimension






END MODULE ED_BATH
