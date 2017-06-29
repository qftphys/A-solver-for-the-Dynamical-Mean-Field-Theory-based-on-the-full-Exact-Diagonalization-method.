MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
  implicit none
  private

  public :: init_ed_structure
  public :: setup_eigenspace
  public :: delete_eigenspace
  public :: setup_pointers_normal
  public :: build_sector
  public :: delete_sector
  public :: bdecomp
  public :: c,cdg
  public :: binary_search



contains





  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure()
    logical                                           :: control
    real(8),dimension(Nspin,Nspin,Norb,Norb)          :: reHloc         !local hamiltonian, real part 
    real(8),dimension(Nspin,Nspin,Norb,Norb)          :: imHloc         !local hamiltonian, imag part
    integer                                           :: i,dim_sector_max(2),iorb,jorb,ispin,jspin
    integer                                           :: maxtwoJz,inJz,dimJz
    integer                                           :: isector,in,shift
    logical                                           :: MPI_MASTER=.true.
    !
    !
    ! Norb    = # of impurity orbitals
    ! Nbath   = # of bath levels (depending on bath_type)
    ! Ns      = # of levels (per spin)
    ! Nlevels = 2*Ns = Total # of levels (counting spin degeneracy 2) 
    select case(bath_type)
    case default
       Ns = (Nbath+1)*Norb
    case ('hybrid')
       Ns = Nbath+Norb
    end select
    Nlevels  = 2*Ns
    !
    Nsectors = (Ns+1)*(Ns+1) !nup=0:Ns;ndw=0:Ns
    !
    dim_sector_max=0
    dim_sector_max(1)=get_normal_sector_dimension(Ns/2)
    dim_sector_max(2)=get_normal_sector_dimension(Ns-Ns/2)
    !
    write(LOGfile,"(A)")"Summary:"
    write(LOGfile,"(A)")"--------------------------------------------"
    write(LOGfile,"(A,I15)")'# of levels/spin      = ',Ns
    write(LOGfile,"(A,I15)")'Total size            = ',Nlevels
    write(LOGfile,"(A,I15)")'# of impurities       = ',Norb
    write(LOGfile,"(A,I15)")'# of bath/impurity    = ',Nbath
    write(LOGfile,"(A,I15)")'# of Bath levels/spin = ',Ns-Norb
    write(LOGfile,"(A,2I15)")'Largest Sector(s)    = ',dim_sector_max
    write(LOGfile,"(A,I15)")'Number of sectors     = ',Nsectors
    write(LOGfile,"(A)")"--------------------------------------------"
    !
    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    reHloc = 0d0 ; imHloc = 0d0
    !
    !Search and read impHloc
    inquire(file=trim(HLOCfile),exist=control)
    if(control)then
       write(LOGfile,*)"Reading impHloc from file: "//reg(HLOCfile)
       open(50,file=trim(HLOCfile),status='old')
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       close(50)
    else
       write(LOGfile,*)"impHloc file not found."
       write(LOGfile,*)"impHloc should be defined elsewhere..."
       call sleep(2)
    endif
    impHloc = dcmplx(reHloc,imHloc)
    if(control)then
       write(LOGfile,"(A)")"H_local:"
       call print_Hloc(impHloc)
    endif
    !
    !
    !Allocate indexing arrays
    allocate(impIndex(Norb,2));impIndex=0
    allocate(getDim(Nsectors));getDim=0
    !
    allocate(getDimUp(Nsectors),getDimDw(Nsectors));getDimUp=0;getDimDw=0
    allocate(getNup(Nsectors),getNdw(Nsectors));getNup=0;getNdw=0
    !
    !
    allocate(getSector(0:Ns,0:Ns))
    getSector=0
    !
    allocate(getCsector(2,Nsectors));getCsector=0
    allocate(getCDGsector(2,Nsectors));getCDGsector=0
    !
    !
    allocate(getBathStride(Norb,Nbath));getBathStride=0
    !
    !
    !CHECKS:
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "ED ERROR: Nspin > 2 is currently not supported"
    if(Norb>2)stop "ED ERROR: Norb > 2 is currently not supported"
    ! if(bath_type=="hybrid")stop "ED ERROR: bath_type==hybrid is currently not supported"
    !
    !
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
    endif
    if(Nspin>1.AND.ed_twin.eqv..true.)then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
       call sleep(1)
    end if
    !
    !
    !
    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    impSmats=zero
    impSreal=zero
    !
    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nspin,Nspin,Norb,Norb,Lreal))
    impGmats=zero
    impGreal=zero
    !
    allocate(impG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impG0real(Nspin,Nspin,Norb,Norb,Lreal))
    impG0mats=zero
    impG0real=zero
    !
    !
    !allocate observables
    allocate(ed_dens(Norb),ed_docc(Norb),ed_dens_up(Norb),ed_dens_dw(Norb))
    ed_dens=0d0
    ed_docc=0d0
    ed_dens_up=0d0
    ed_dens_dw=0d0
    !
    if(chiflag)then
       allocate(spinChi_tau(Norb+1,0:Ltau))
       allocate(spinChi_w(Norb+1,Lreal))
       allocate(spinChi_iv(Norb+1,0:Lmats))
    endif
    !
  end subroutine init_ed_structure




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_eigenspace
    integer :: isector,dim,jsector
    if(allocated(espace)) deallocate(espace)
    allocate(espace(1:Nsectors))
    do isector=1,Nsectors
       dim=getdim(isector);if(dim==0)stop "setup_eigenspace: dim==0!"
       allocate(espace(isector)%e(dim))
       allocate(espace(isector)%M(dim,dim))
    enddo
  end subroutine setup_eigenspace


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine delete_eigenspace
    integer :: isector
    if(allocated(espace))then
       do isector=1,size(espace)
          deallocate(espace(isector)%e)
          deallocate(espace(isector)%M)
       end do
       deallocate(espace)
    endif
  end subroutine delete_eigenspace



  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  !
  !NORMAL CASE
  !
  subroutine setup_pointers_normal
    integer                                           :: i,in,dim,isector,jsector
    integer                                           :: nup,ndw,jup,jdw,iorb
    integer                                           :: unit,status,istate
    logical                                           :: IOfile
    integer                                           :: anint
    real(8)                                           :: adouble
    integer                                           :: list_len
    integer,dimension(:),allocatable                  :: list_sector
    isector=0
    do nup=0,Ns
       do ndw=0,Ns
          isector=isector+1
          getSector(nup,ndw)=isector
          getNup(isector)=nup
          getNdw(isector)=ndw
          dim = get_normal_sector_dimension(nup,ndw)
          getDim(isector)=dim
          getDimUp(isector)=get_normal_sector_dimension(nup)
          getDimDw(isector)=get_normal_sector_dimension(ndw)
       enddo
    enddo
    !
    !
    !
    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo
    !
    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)       = Norb + i
       enddo
    end select
    !
    getCsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup-1;jdw=ndw;if(jup < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw-1;if(jdw < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup+1;jdw=ndw;if(jup > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_normal







  !+------------------------------------------------------------------+
  !PURPOSE  : return the dimension of a sector
  !+------------------------------------------------------------------+
  !NORMAL
  function get_normal_sector_dimension(nup,ndw) result(dim)
    integer :: nup
    integer,optional :: ndw
    integer :: dim,dimup,dimdw
    if(present(ndw))then
       dimup = binomial(Ns,nup)    !this ensures better evaluation of the dimension
       dimdw = binomial(Ns,ndw)    !as it avoids large numbers
    else
       dimup = binomial(Ns,nup)
       dimdw = 1
    endif
    dim=dimup*dimdw
  end function get_normal_sector_dimension




  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the 
  !states i\in Hilbert_space from the states count in H_sector.
  !|ImpUP,BathUP>|ImpDW,BathDW >
  !+------------------------------------------------------------------+
  subroutine build_sector(isector,Hup)
    integer                                      :: isector
    type(sector_map)                             :: Hup
    integer                                      :: nup,ndw,sz,nt,twoJz
    integer                                      :: nup_,ndw_,sz_,nt_
    integer                                      :: twoSz_,twoLz_
    integer                                      :: i,ibath,iorb
    integer                                      :: iup,idw
    integer                                      :: dim
    integer                                      :: ivec(Ns),jvec(Ns)
    nup = getNup(isector)
    ndw = getNdw(isector)
    dim = getDim(isector)
    call map_allocate(Hup,dim)
    dim=0
    do idw=0,2**Ns-1
       jvec  = bdecomp(idw,Ns)
       ndw_  = sum(jvec)
       if(ndw_ /= ndw)cycle
       do iup=0,2**Ns-1
          ivec  = bdecomp(iup,Ns)
          nup_  = sum(ivec)
          if(nup_ /= nup)cycle
          dim      = dim+1
          Hup%map(dim) = iup + idw*2**Ns
       enddo
    enddo
  end subroutine build_sector





  subroutine delete_sector(isector,Hup)!,Hdw)
    integer                   :: isector
    type(sector_map)          :: Hup
    ! type(sector_map),optional :: Hdw
    call map_deallocate(Hup)
  end subroutine delete_sector




  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg















  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  function binomial(n1,n2) result(nchoos)
    real(8) :: xh
    integer :: n1,n2,i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search








end MODULE ED_SETUP
