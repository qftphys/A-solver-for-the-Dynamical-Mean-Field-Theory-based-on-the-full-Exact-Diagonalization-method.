  !Diagonal Elements, i.e. local part
  htmp = zero
  htmp = htmp - xmu*(sum(nup)+sum(ndw))
  !
  do iorb=1,Norb
     htmp = htmp + impHloc(1,1,iorb,iorb)*nup(iorb)
     htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
  enddo
  !
  Hmat(i,i)=Hmat(i,i) + htmp
  !
  !
  !Off-diagonal elements, i.e. non-local part
  !1. same spin:
  do iorb=1,Norb
     do jorb=1,Norb
        !this loop considers only the orbital off-diagonal terms
        !because iorb=jorb can not have simultaneously
        !occupation 0 and 1, as required by this if Jcondition:
        !UP
        Jcondition = &
             (impHloc(1,1,iorb,jorb)/=zero) .AND. &
             (ib(jorb)==1)                  .AND. &
             (ib(iorb)==0)
        if (Jcondition) then
           call c(jorb,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           j = binary_search(H%map,k2)
           htmp = impHloc(1,1,iorb,jorb)*sg1*sg2
           !
           Hmat(i,j) = Hmat(i,j) + htmp
           !
        endif
        !DW
        Jcondition = &
             (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
             (ib(jorb+Ns)==1)                       .AND. &
             (ib(iorb+Ns)==0)
        if (Jcondition) then
           call c(jorb+Ns,m,k1,sg1)
           call cdg(iorb+Ns,k1,k2,sg2)
           j = binary_search(H%map,k2)
           htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
           !
           Hmat(i,j) = Hmat(i,j) + htmp
           !
        endif
     enddo
  enddo
  !
