!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!NOTE: in the MPI implementation we may require all the nodes to 
!evaluate the GF, this is safer, simpler and works for both Lanc &
!Ed. For Lanc we can indeed assign the contribution from each state 
!to different node and accumulate the result at the end.
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GREENS_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy,splot
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_LINALG,  only: inv,inv_sym,inv_her,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  implicit none
  private 



  public :: buildGf_impurity

  public :: buildChi_impurity



  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer                :: state_vec
  complex(8),dimension(:),pointer             :: state_cvec
  real(8)                                     :: state_e

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable            :: wm,tau,wr,vm

  !Auxiliary functions GF
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:) :: impDeltamats,impDeltareal
  complex(8),allocatable,dimension(:,:,:,:,:) :: invimpG0mats,invimpG0real
  complex(8),allocatable,dimension(:,:,:,:,:) :: invimpGmats,invimpGreal



contains



  subroutine buildgf_impurity()
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*dble(2*arange(1,Lmats)-1)
    wr     = linspace(wini,wfin,Lreal)
    !
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    if(.not.allocated(impDeltamats)) allocate(impDeltamats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpG0mats)) allocate(invimpG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpGmats))  allocate( invimpGmats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impDeltareal)) allocate(impDeltareal(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpG0real)) allocate(invimpG0real(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpGreal))  allocate( invimpGreal(Nspin,Nspin,Norb,Norb,Lreal))
    impDeltamats=zero
    invimpGmats=zero
    invimpG0mats=zero
    impDeltareal=zero
    invimpGreal=zero
    invimpG0real=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    call get_sigma_normal()
    !
    if(ed_print_Sigma)call ed_print_impSigma()
    if(ed_print_G)call ed_print_impG()
    if(ed_print_G0)call ed_print_impG0()
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    if(allocated(invimpG0mats))deallocate(invimpG0mats)
    if(allocated(invimpGmats))deallocate(invimpGmats)
    if(allocated(impDeltamats))deallocate(impDeltamats)
    if(allocated(invimpG0real))deallocate(invimpG0real)
    if(allocated(invimpGreal))deallocate(invimpGreal)
    if(allocated(impDeltareal))deallocate(impDeltareal)
  end subroutine buildgf_impurity





  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Green's functions
  !+------------------------------------------------------------------+
  subroutine build_gf_normal()
    integer :: iorb,jorb,ispin,i
    logical :: MaskBool
    !
    !
    !NORMAL: (default)
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_s"//str(ispin)
          call full_build_gf_normal_c_mix(iorb,ispin)
       enddo
    enddo
    !
  end subroutine build_gf_normal











  !+------------------------------------------------------------------+
  !PURPOSE  : DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine full_build_gf_normal_c(iorb,ispin)
    real(8)          :: cdgmat,cc
    integer          :: iorb,ispin,isite,istate
    integer          :: idim,isector
    integer          :: jdim,jsector
    integer          :: ib(Nlevels)
    integer          :: m,i,j,r,k,ll
    real(8)          :: sgn
    real(8)          :: Ei,Ej,matcdg
    real(8)          :: expterm,peso,de,w0,it,chij1
    complex(8)       :: iw
    type(sector_map) :: HI,HJ
    isite=impIndex(iorb,ispin)
    !
    call start_timer
    !
    do isector=1,Nsectors
       !
       jsector=getCDGsector(ispin,isector)
       if(jsector==0)cycle
       !
       call eta(isector,Nsectors)
       idim=getdim(isector)     !i-th sector dimension
       jdim=getdim(jsector)     !j-th sector dimension
       call build_sector(isector,HI)
       call build_sector(jsector,HJ)
       do i=1,idim          !loop over the states in the i-th sect.
          do j=1,jdim       !loop over the states in the j-th sect.
             cdgmat=0.d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm < cutoff)cycle
             !
             do ll=1,idim              !loop over the component of |j> (IN state!)
                m = HI%map(ll)
                ib = bdecomp(m,2*Ns)
                if(ib(isite) == 0)then
                   call cdg(isite,m,k,cc)
                   r = binary_search(HJ%map,k)
                   cdgmat=cdgmat+espace(jsector)%M(r,j)*cc*espace(isector)%M(ll,i)
                endif
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(jsector)%e(j)
             de=Ej-Ei
             peso=expterm/zeta_function
             matcdg=peso*cdgmat**2
             !
             do m=1,Lmats
                iw=xi*wm(m)
                impGmats(ispin,ispin,iorb,iorb,m)=impGmats(ispin,ispin,iorb,iorb,m)+matcdg/(iw+de)
             enddo
             !
             do m=1,Lreal
                w0=wr(m);iw=cmplx(w0,eps)
                impGreal(ispin,ispin,iorb,iorb,m)=impGreal(ispin,ispin,iorb,iorb,m)+matcdg/(iw+de)
             enddo
             !
          enddo
       enddo
       deallocate(HI%map,HJ%map)
    enddo
    call stop_progress
  end subroutine full_build_gf_normal_c



  !+------------------------------------------------------------------+
  !PURPOSE  : DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine full_build_gf_normal_c_mix(iorb,ispin)
    real(8)          :: spectral_weight,op_mat(2),sgn_cdg,sgn_c
    integer          :: iorb,ispin,isite,istate
    integer          :: idim,isector
    integer          :: jdim,jsector
    integer          :: ib(Nlevels)
    integer          :: li,rj
    integer          :: m,i,j,r,k,p
    real(8)          :: sgn
    real(8)          :: Ei,Ej
    real(8)          :: expterm,peso,de,w0,it,chij1
    complex(8)       :: iw
    type(sector_map) :: HI,HJ
    !
    isite=impIndex(iorb,ispin)
    !
    call start_timer
    !
    do isector=1,Nsectors
       jsector=getCDGsector(ispin,isector)
       if(jsector==0)cycle
       call eta(isector,Nsectors)
       idim=getdim(isector)     !i-th sector dimension
       jdim=getdim(jsector)     !j-th sector dimension
       call build_sector(isector,HI)
       call build_sector(jsector,HJ)
       do i=1,idim          !loop over the states in the i-th sect.
          do j=1,jdim       !loop over the states in the j-th sect.
             op_mat=0.d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm < cutoff)cycle
             !
             do li=1,idim              !loop over the component of |I> (IN state!)
                m = HI%map(li)
                !
                ib = bdecomp(m,2*Ns)
                if(ib(isite) == 1)cycle
                call cdg(isite,m,k,sgn_cdg)
                rj = binary_search(HJ%map,k)
                !
                ib = bdecomp(k,2*Ns)
                if(ib(isite) == 0)cycle
                call c(isite,k,p,sgn_c)
                if(p/=m)cycle
                !
                op_mat(1)=op_mat(1) + conjg(espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
                op_mat(2)=op_mat(2) + conjg(espace(isector)%M(li,i))*sgn_c*espace(jsector)%M(rj,j)
                !
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(jsector)%e(j)
             de=Ej-Ei
             peso=expterm/zeta_function
             spectral_weight=peso*product(op_mat)
             !
             do m=1,Lmats
                iw=xi*wm(m)
                impGmats(ispin,ispin,iorb,iorb,m)=impGmats(ispin,ispin,iorb,iorb,m)+spectral_weight/(iw+de)
             enddo
             !
             do m=1,Lreal
                w0=wr(m);iw=cmplx(w0,eps)
                impGreal(ispin,ispin,iorb,iorb,m)=impGreal(ispin,ispin,iorb,iorb,m)+spectral_weight/(iw+de)
             enddo
             !
          enddo
       enddo
       deallocate(HI%map,HJ%map)
    enddo
    call stop_progress
  end subroutine full_build_gf_normal_c_mix





  !+------------------------------------------------------------------+
  !                    SELF-ENERGY FUNCTIONS 
  !+------------------------------------------------------------------+
  subroutine get_sigma_normal
    integer                                           :: i,ispin,iorb
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
    complex(8),dimension(Norb,Norb)                   :: invGimp
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    !
    !Get G0^-1
    invG0mats(:,:,:,:,:) = invg0_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real(:,:,:,:,:) = invg0_bath_real(dcmplx(wr(:),eps),dmft_bath)
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       !
       !Get Gimp^-1
       do ispin=1,Nspin
          do iorb=1,Norb
             invGmats(ispin,ispin,iorb,iorb,:) = one/impGmats(ispin,ispin,iorb,iorb,:)
             invGreal(ispin,ispin,iorb,iorb,:) = one/impGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             impSmats(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
             impSreal(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !
    case ("hybrid")   !Diagonal in spin only. Full Orbital structure
       !
       !Get Gimp^-1
       do ispin=1,Nspin
          do i=1,Lmats
             invGimp = impGmats(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGmats(ispin,ispin,:,:,i)=invGimp
          enddo
          !
          do i=1,Lreal
             invGimp = impGreal(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGreal(ispin,ispin,:,:,i)=invGimp
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          impSmats(ispin,ispin,:,:,:) = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
          !
          impSreal(ispin,ispin,:,:,:) = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       enddo
       !
    end select
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_real(dcmplx(wr(:),eps),dmft_bath)
    !!
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine get_sigma_normal






  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildChi_impurity()
    integer :: i
    !
    call allocate_grids
    !
    !
    !BUILD SPIN SUSCEPTIBILITY
    if(.not.allocated(spinChi_tau)) stop "buildChi_impurity: spinChi_tau not allocated"
    if(.not.allocated(spinChi_w))  stop "buildChi_impurity: spinChi_w not allocated"
    if(.not.allocated(spinChi_iv)) stop "buildChi_impurity: spinChi_iv not allocated"
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    call build_chi_spin()
    !PRINTING:
    call ed_print_impChi()
    !
    call deallocate_grids
  end subroutine buildChi_impurity




  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Spin Susceptibility 
  !+------------------------------------------------------------------+
  subroutine build_chi_spin()
    integer :: iorb
    write(LOGfile,"(A)")"Get impurity spin Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_spin_l"//reg(txtfy(iorb))
       call full_ed_build_spinChi_c(iorb)
    enddo
    ! if(Norb>1)then
    !    write(LOGfile,"(A)")"Get Chi_spin_tot"
    !    call lanc_ed_build_spinChi_tot_c()
    ! endif
    spinChi_tau = SpinChi_tau/zeta_function
    spinChi_w   = spinChi_w/zeta_function
    spinChi_iv  = spinChi_iv/zeta_function
  end subroutine build_chi_spin





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine full_ed_build_spinChi_c(iorb)
    integer                    :: iorb
    real(8)                    :: cdgmat(Norb),chij(Norb),chitot,spin,spintot
    integer,dimension(Nlevels) :: ib
    integer                    :: i,j,k,r,ll,m,in,is,ispin,jorb,isector
    integer                    :: idim,ia
    real(8)                    :: Ei,Ej,cc,peso(Norb),pesotot
    real(8)                    :: expterm,de,w0,it
    complex(8)                 :: iw 
    type(sector_map)           :: HI    !map of the Sector S to Hilbert space H
    !
    !Spin susceptibility \X(tau). |<i|S_z|j>|^2
    call start_timer()
    do isector=1,Nsectors !loop over <i| total particle number
       call eta(isector,Nsectors,LOGfile)
       idim=getdim(isector)
       call build_sector(isector,HI)
       do i=1,idim 
          do j=1,idim
             chij=0.d0
             chitot=0.d0
             expterm=exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idim 
                ia=HI%map(ll)
                ib = bdecomp(ia,2*Ns)
                spintot=0.d0
                spin=dble(ib(iorb))-dble(ib(iorb+Ns)) !nup - ndw
                chij(iorb)=chij(iorb)+espace(isector)%M(ll,i)*spin*espace(isector)%M(ll,j)
                ! spintot=spintot+dble(ib(iorb))-dble(ib(iorb+Ns))
                ! chitot=chitot+espace(isector)%M(ll,i)*spintot*espace(isector)%M(ll,j)
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ei-Ej
             peso=chij/zeta_function
             pesotot=chitot/zeta_function
             !
             !Matsubara (bosonic) frequency
             if(de>cutoff)spinChi_iv(iorb,0)=spinChi_iv(iorb,0)-peso(iorb)*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/de
             do m=1,Lmats
                iw=xi*vm(m)
                spinChi_iv(iorb,m)=spinChi_iv(iorb,m)+peso(iorb)*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
             enddo
             !
             !Real-frequency
             do m=1,Lreal
                w0=wr(m);iw=cmplx(w0,eps,8)
                !Retarded = Commutator = response function
                spinChi_w(iorb,m)=spinChi_w(iorb,m)+peso(iorb)*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
             enddo
             !
             !Imaginary time:
             do m=0,Ltau 
                it=tau(m)
                spinChi_tau(iorb,m)=spinChi_tau(iorb,m) + exp(-it*Ei)*exp(-(beta-it)*Ej)*peso(iorb)
             enddo
             !
          enddo
       enddo
    enddo
    call stop_progress
  end subroutine full_ed_build_spinChi_c




  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2.d0*dble(i)
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
  end subroutine deallocate_grids




end MODULE ED_GREENS_FUNCTIONS
