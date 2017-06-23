!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_BATH
  USE ED_AUX_FUNX

  implicit none
  private
  !
  public :: observables_impurity

  public :: local_energy_impurity


  logical,save                       :: iolegend=.true.
  real(8),dimension(:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable   :: docc
  real(8),dimension(:),allocatable   :: magz
  real(8),dimension(:,:),allocatable :: sz2,n2
  real(8),dimensioN(:,:),allocatable :: zimp,simp
  real(8)                            :: s2tot
  real(8)                            :: Egs
  !



contains 





  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_impurity()
    integer,dimension(Nlevels)      :: ib
    integer                         :: i,j
    integer                         :: izero,istate
    integer                         :: isector,jsector
    integer                         :: idim,jdim
    integer                         :: isz,jsz
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite,ibath
    integer                         :: numstates
    integer                         :: r,m,k
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: boltzman_weight
    real(8)                         :: state_weight
    real(8)                         :: weight
    real(8)                         :: Ei
    real(8)                         :: norm
    real(8),dimension(Norb)         :: nup,ndw,Sz,nt
    type(sector_map)                :: H,HJ
    complex(8),allocatable          :: vvinit(:)
    complex(8),dimension(:),pointer :: evec
    !
    !
    !LOCAL OBSERVABLES:
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
    !
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    !
    do isector=1,Nsectors
       idim    = getdim(isector)
       call build_sector(isector,H)
       !
       do istate=1,idim
          Ei=espace(isector)%e(istate)
          boltzman_weight=exp(-beta*Ei)/zeta_function
          if(boltzman_weight < cutoff)cycle
          !
          evec => espace(isector)%M(:,istate)
          !
          do i=1,idim
             m=H%map(i)
             ib = bdecomp(m,2*Ns)
             !
             state_weight=conjg(evec(i))*evec(i)
             weight = boltzman_weight*state_weight
             !
             !Get operators:
             do iorb=1,Norb
                nup(iorb)= dble(ib(iorb))
                ndw(iorb)= dble(ib(iorb+Ns))
                sz(iorb) = (nup(iorb) - ndw(iorb))/2.d0
                nt(iorb) =  nup(iorb) + ndw(iorb)
             enddo
             !
             !Evaluate averages of observables:
             do iorb=1,Norb
                dens(iorb)     = dens(iorb)      +  nt(iorb)*weight
                dens_up(iorb)  = dens_up(iorb)   +  nup(iorb)*weight
                dens_dw(iorb)  = dens_dw(iorb)   +  ndw(iorb)*weight
                docc(iorb)     = docc(iorb)      +  nup(iorb)*ndw(iorb)*weight
                magz(iorb)     = magz(iorb)      +  (nup(iorb)-ndw(iorb))*weight
                sz2(iorb,iorb) = sz2(iorb,iorb)  +  (sz(iorb)*sz(iorb))*weight
                n2(iorb,iorb)  = n2(iorb,iorb)   +  (nt(iorb)*nt(iorb))*weight
                do jorb=iorb+1,Norb
                   sz2(iorb,jorb) = sz2(iorb,jorb)  +  (sz(iorb)*sz(jorb))*weight
                   sz2(jorb,iorb) = sz2(jorb,iorb)  +  (sz(jorb)*sz(iorb))*weight
                   n2(iorb,jorb)  = n2(iorb,jorb)   +  (nt(iorb)*nt(jorb))*weight
                   n2(jorb,iorb)  = n2(jorb,iorb)   +  (nt(jorb)*nt(iorb))*weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*weight
          enddo
       enddo
       deallocate(H%map)
    enddo


    call get_szr
    if(iolegend)call write_legend
    call write_observables()

    write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
    write(LOGfile,"(A,10f18.12,A)")"docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
    if(Nspin==2)write(LOGfile,"(A,10f18.12,A)") "mag "//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)

    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
    enddo
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2,n2)
    deallocate(simp,zimp)    
  end subroutine observables_impurity




  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_impurity()
    integer,dimension(Nlevels)      :: ib
    integer                         :: i,j
    integer                         :: izero,istate
    integer                         :: isector
    integer                         :: idim
    integer                         :: iorb,jorb,ispin
    integer                         :: numstates
    integer                         :: m,k1,k2,k3,k4
    real(8)                         :: sg1,sg2,sg3,sg4
    real(8)                         :: Egs
    real(8)                         :: Ei
    real(8)                         :: boltzman_weight
    real(8)                         :: state_weight
    real(8)                         :: weight
    real(8)                         :: norm
    real(8),dimension(Norb)         :: nup,ndw
    real(8),dimension(Nspin,Norb)   :: eloc
    complex(8),dimension(:),pointer :: evec
    type(sector_map)                :: H
    logical                         :: Jcondition
    !
    !
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !
    do isector=1,Nsectors
       idim    = getdim(isector)
       call build_sector(isector,H)
       !
       do istate=1,idim
          Ei=espace(isector)%e(istate)
          boltzman_weight=exp(-beta*Ei)/zeta_function
          if(boltzman_weight < cutoff)cycle
          !
          evec => espace(isector)%M(:,istate)
          !
          do i=1,idim
             m=H%map(i)
             ib = bdecomp(m,2*Ns)
             !
             state_weight = conjg(evec(i))*evec(i)
             weight = boltzman_weight*state_weight
             !
             !Get operators:
             do iorb=1,Norb
                nup(iorb)= dble(ib(iorb))
                ndw(iorb)= dble(ib(iorb+Ns))
             enddo
             !
             !start evaluating the Tr(H_loc) to estimate potential energy
             !
             !LOCAL ENERGY
             ed_Eknot = ed_Eknot + dot_product(eloc(1,:),nup)*weight + dot_product(eloc(Nspin,:),ndw)*weight
             !==> HYBRIDIZATION TERMS I: same or different orbitals, same spins.
             do iorb=1,Norb
                do jorb=1,Norb
                   !SPIN UP
                   if((ib(iorb)==0).AND.(ib(jorb)==1))then
                      call c(jorb,m,k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      j=binary_search(H%map,k2)
                      if(j==0)cycle
                      ed_Eknot = ed_Eknot + impHloc(1,1,iorb,jorb)*sg1*sg2*evec(i)*conjg(evec(j))*boltzman_weight
                   endif
                   !SPIN DW
                   if((ib(iorb+Ns)==0).AND.(ib(jorb+Ns)==1))then
                      call c(jorb+Ns,m,k1,sg1)
                      call cdg(iorb+Ns,k1,k2,sg2)
                      j=binary_search(H%map,k2)
                      if(j==0)cycle
                      ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*evec(i)*conjg(evec(j))*boltzman_weight
                   endif
                enddo
             enddo
             !
             !
             !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*weight
             do iorb=1,Norb
                ed_Epot = ed_Epot + Uloc(iorb)*nup(iorb)*ndw(iorb)*weight
             enddo
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*weight
                      ed_Dust = ed_Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*weight
                   enddo
                enddo
             endif
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*weight
                      ed_Dund = ed_Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*weight
                   enddo
                enddo
             endif
             !
             !SPIN-EXCHANGE (S-E) TERMS
             !S-E: Jh *( c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up )  (i.ne.j) 
             if(Norb>1.AND.Jhflag)then
                do iorb=1,Norb
                   do jorb=1,Norb
                      Jcondition=((iorb/=jorb).AND.&
                           (ib(jorb)==1)      .AND.&
                           (ib(iorb+Ns)==1)   .AND.&
                           (ib(jorb+Ns)==0)   .AND.&
                           (ib(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,m,k1,sg1)
                         call c(iorb+Ns,k1,k2,sg2)
                         call cdg(jorb+Ns,k2,k3,sg3)
                         call cdg(iorb,k3,k4,sg4)
                         j=binary_search(H%map,k4)
                         if(j==0)cycle
                         ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                         ed_Dse  = ed_Dse  + sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                      endif
                   enddo
                enddo
             endif
             !
             !PAIR-HOPPING (P-H) TERMS
             !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
             !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
             if(Norb>1.AND.Jhflag)then
                do iorb=1,Norb
                   do jorb=1,Norb
                      Jcondition=((iorb/=jorb).AND.&
                           (ib(jorb)==1)      .AND.&
                           (ib(jorb+Ns)==1)   .AND.&
                           (ib(iorb+Ns)==0)   .AND.&
                           (ib(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,m,k1,sg1)
                         call c(jorb+Ns,k1,k2,sg2)
                         call cdg(iorb+Ns,k2,k3,sg3)
                         call cdg(iorb,k3,k4,sg4)
                         j=binary_search(H%map,k4)
                         if(j==0)cycle
                         ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                         ed_Dph  = ed_Dph  + sg1*sg2*sg3*sg4*evec(i)*conjg(evec(j))*boltzman_weight
                      endif
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                do iorb=1,Norb
                   ed_Ehartree=ed_Ehartree - 0.5d0*uloc(iorb)*(nup(iorb)+ndw(iorb))*weight + 0.25d0*uloc(iorb)*weight
                enddo
                if(Norb>1)then
                   do iorb=1,Norb
                      do jorb=iorb+1,Norb
                         ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*weight + 0.25d0*Ust*weight
                         ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*weight + 0.25d0*(Ust-Jh)*weight
                      enddo
                   enddo
                endif
             endif
          enddo
       enddo
       deallocate(H%map)
    enddo
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose==3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
       write(LOGfile,"(A,10f18.12)")"Dse     =",ed_Dse
       write(LOGfile,"(A,10f18.12)")"Dph     =",ed_Dph
    endif
    call write_energy_info()
    call write_energy()
    !
    !
  end subroutine local_energy_impurity



  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: ispin,iorb
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3d0*pi/beta
    do ispin=1,Nspin
       do iorb=1,Norb
          simp(iorb,ispin) = dimag(impSmats(ispin,ispin,iorb,iorb,1)) - &
               wm1*(dimag(impSmats(ispin,ispin,iorb,iorb,2))-dimag(impSmats(ispin,ispin,iorb,iorb,1)))/(wm2-wm1)
          zimp(iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,ispin,iorb,iorb,1))/wm1 ))
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_info.ed")
    write(unit,"(A1,90(A10,6X))")"#",&
         (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(5*Norb+1))//"s2",&
         reg(txtfy(5*Norb+2))//"egs",&
         ((reg(txtfy(5*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((5+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((5+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         ((reg(txtfy((6+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend

  subroutine write_energy_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>",&
         reg(txtfy(7))//"<Dse>",&
         reg(txtfy(8))//"<Dph>"
    close(unit)
  end subroutine write_energy_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_all"//reg(ed_file_suffix)//".ed",position='append')
    write(unit,"(90(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         (magz(iorb),iorb=1,Norb),&
         s2tot,egs,&
         ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)    
    !
    unit = free_unit()
    open(unit,file="parameters_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
    close(unit)
    !
    unit = free_unit()
    open(unit,file="observables_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90(F15.9,1X))")&
         (dens(iorb),iorb=1,Norb),&
         (docc(iorb),iorb=1,Norb),&
         (dens_up(iorb),iorb=1,Norb),&
         (dens_dw(iorb),iorb=1,Norb),&
         (magz(iorb),iorb=1,Norb),&
         s2tot,egs,&
         ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         
  end subroutine write_observables

  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy



end MODULE ED_OBSERVABLES







! !
! !SUPERCONDUCTING ORDER PARAMETER
! if(ed_mode=="superc")then
!    do ispin=1,Nspin
!       do iorb=1,Norb
!          do isector=1,Nsectors
!             idim    = getdim(isector)
!             call build_sector(isector,H)
!             do istate=1,idim
!                Ei=espace(isector)%e(istate)
!                boltzman_weight=exp(-beta*Ei)/zeta_function
!                if(boltzman_weight < cutoff)cycle
!                !
!                evec => espace(isector)%M(:,istate)
!                !
!                !GET <(C_UP + CDG_DW)(CDG_UP + C_DW)> = <N_UP> + < 1 - N_DW> + 2*<PHI>
!                isz = getsz(isector)
!                if(isz<Ns)then
!                   jsz     = isz+1
!                   jsector = getsector(jsz,1)
!                   jdim    = getdim(jsector)
!                   allocate(vvinit(jdim))
!                   vvinit=zero
!                   call build_sector(jsector,HJ)
!                   do i=1,idim
!                      m=H%map(i)
!                      ib = bdecomp(m,2*Ns)
!                      if(ib(iorb)==0)then
!                         call cdg(iorb,m,r,sgn)
!                         j=binary_search(HJ%map,r)
!                         vvinit(j) = sgn*evec(i)
!                      endif
!                   enddo
!                   do i=1,idim
!                      m=H%map(i)
!                      ib = bdecomp(m,2*Ns)
!                      if(ib(iorb+Ns)==1)then
!                         call c(iorb+Ns,m,r,sgn)
!                         j=binary_search(HJ%map,r)
!                         vvinit(j) = vvinit(j) + sgn*evec(i)
!                      endif
!                   enddo
!                   deallocate(HJ%map)
!                   phisc(iorb) = phisc(iorb) + dot_product(vvinit,vvinit)*boltzman_weight
!                   deallocate(vvinit)
!                endif
!                if(associated(evec)) nullify(evec)
!             enddo
!             deallocate(H%map)
!             phisc(iorb) = 0.5d0*(phisc(iorb) - dens_up(iorb) - (1.d0-dens_dw(iorb)))
!          enddo
!       enddo
!    enddo
! end if


! !IMPURITY DENSITY MATRIX
! if(allocated(imp_density_matrix)) deallocate(imp_density_matrix)
! allocate(imp_density_matrix(Nspin,Nspin,Norb,Norb))
! imp_density_matrix=zero
! do isector=1,Nsectors
!    idim    = getdim(isector)
!    call build_sector(isector,H)
!    do istate=1,idim
!       Ei=espace(isector)%e(istate)
!       boltzman_weight=exp(-beta*Ei)/zeta_function
!       if(boltzman_weight < cutoff)cycle
!       !
!       evec => espace(isector)%M(:,istate)
!       !
!       !
!       !Diagonal densities
!       do ispin=1,Nspin
!          do iorb=1,Norb
!             isite=impIndex(iorb,ispin)
!             do i=1,idim
!                m=H%map(i)
!                ib = bdecomp(m,2*Ns)
!                imp_density_matrix(ispin,ispin,iorb,iorb) = imp_density_matrix(ispin,ispin,iorb,iorb) + &
!                     ib(isite)*conjg(evec(i))*evec(i)*boltzman_weight
!             enddo
!          enddo
!       enddo
!       !off-diagonal
!       do ispin=1,Nspin
!          do jspin=1,Nspin
!             do iorb=1,Norb
!                do jorb=1,Norb
!                   if((ed_mode=="normal").and.(ispin/=jspin))cycle
!                   if((bath_type=="normal").and.(iorb/=jorb))cycle
!                   if((.not.dmft_bath%mask(ispin,jspin,iorb,jorb,1)).and.(.not.dmft_bath%mask(ispin,jspin,iorb,jorb,2)))cycle
!                   isite=impIndex(iorb,ispin)
!                   jsite=impIndex(jorb,jspin)
!                   do i=1,idim
!                      m=H%map(i)
!                      ib = bdecomp(m,2*Ns)
!                      if((ib(jsite)==1).and.(ib(isite)==0))then
!                         call c(jsite,i,r,sgn1)
!                         call cdg(isite,r,k,sgn2)
!                         j=binary_search(H%map,k)
!                         imp_density_matrix(ispin,jspin,iorb,jorb) = imp_density_matrix(ispin,jspin,iorb,jorb) + &
!                              sgn1*evec(i)*sgn2*conjg(evec(j))*boltzman_weight
!                      endif
!                   enddo
!                enddo
!             enddo
!          enddo
!       enddo
!    enddo
!    deallocate(H%map)
! enddo



