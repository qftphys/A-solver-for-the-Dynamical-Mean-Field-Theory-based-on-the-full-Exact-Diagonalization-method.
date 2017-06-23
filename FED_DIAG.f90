!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module ED_DIAG
  USE SF_CONSTANTS
  USE SF_LINALG, only: eigh
  USE SF_TIMER,  only: start_timer,stop_timer,eta
  USE SF_IOTOOLS, only:reg,free_unit
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_SETUP
  USE ED_HAMILTONIAN
  !
  implicit none
  private


  public :: diagonalize_impurity


contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine diagonalize_impurity
    integer                     :: nup,ndw,isector,dim
    integer                     :: sz,nt
    integer                     :: i,j,unit
    real(8),dimension(Nsectors) :: e0 
    real(8)                     :: egs
    logical                     :: Tflag
    !
    e0=1000.d0
    write(LOGfile,"(A)")"Diagonalize impurity H:"
    call start_timer()
    !
    !
    sector: do isector=1,Nsectors
       !
       Dim      = getdim(isector)
       !
       if(ed_verbose==3)then
          nup  = getnup(isector)
          ndw  = getndw(isector)
          write(LOGfile,"(A,I4,A6,I2,A6,I2,A6,I15)")"Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector)
       elseif(ed_verbose==1.OR.ed_verbose==2)then
          call eta(isector,Nsectors,LOGfile)
       endif
       !
       call setup_Hv_sector(isector)
       call ed_buildH_c(espace(isector)%M)
       call delete_Hv_sector()
       call eigh(espace(isector)%M,espace(isector)%e,'V','U')
       if(dim==1)espace(isector)%M=one
       !
       e0(isector)=minval(espace(isector)%e)
       !
    enddo sector
    !
    call stop_timer(LOGfile)
    !
    !Get the ground state energy and rescale energies
    egs=minval(e0)
    forall(isector=1:Nsectors)espace(isector)%e = espace(isector)%e - egs
    !
    !Get the partition function Z
    zeta_function=0.d0;zeta_function=0.d0
    do isector=1,Nsectors
       dim=getdim(isector)
       do i=1,dim
          zeta_function=zeta_function+exp(-beta*espace(isector)%e(i))
       enddo
    enddo
    !
    write(LOGfile,"(A)")"DIAG resume:"
    open(free_unit(unit),file='egs'//reg(ed_file_suffix)//".ed",position='append')
    do isector=1,Nsectors
       if(e0(isector)/=0d0)cycle
       nup  = getnup(isector)
       ndw  = getndw(isector)
       write(LOGfile,"(A,F20.12,2I4)")'Egs =',e0(isector),nup,ndw
       write(unit,"(F20.12,2I4)")e0(isector),nup,ndw
    enddo
    write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    close(unit)
    return
  end subroutine diagonalize_impurity


end MODULE ED_DIAG









