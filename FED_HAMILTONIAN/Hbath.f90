  !
  !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
  htmp=zero
  do iorb=1,size(dmft_bath%e,2)
     do kp=1,Nbath
        alfa=getBathStride(iorb,kp)
        htmp =htmp + dmft_bath%e(1,iorb,kp)*ib(alfa)        !UP
        htmp =htmp + dmft_bath%e(Nspin,iorb,kp)*ib(alfa+Ns) !DW
     enddo
  enddo
  !
  Hmat(i,i) = Hmat(i,i) + htmp
  !


