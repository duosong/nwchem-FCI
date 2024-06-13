      subroutine fci_string(rtdb,
     &    d_v2,k_v2_offset,
     &    d_t1,k_t1_offset,
     &    d_t2,k_t2_offset,
     &    nos,nvs,noas,nobs,nvas,nvbs,
     &    nstot,
     &    h,v,repulsion)
c
c  Only for RHF
c
c all varaible below have to be defined in the tce_energy.f
c
c please do not freeze the occupied orbitals for now
c
c nstot - total number of spinorbitals nstot=nos+nvs
c nos   - number of occupied spinorbitals
c nvs   - number of virtual spinorbitals
c noas  - number of occupied alpha spinorbitals
c nobs  - number of occupied beta spinorbitals
c nvas  - number of virtual  alhoa spinorbitals
c nvbs  - number of virtual beta spinorbitals
c
c spinorbital convention
c
c | noas | nobs | nvas | nvbs |
c
c ene_orb contains orbital energies
c
c COMMENTS: 
c iext1 & iext2 - redundant
c
c
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_v2    ! GA handle for v2
      integer k_v2_offset
      integer d_t1,d_t2
      integer k_t1_offset,k_t2_offset
      integer nos,nvs ! # of occupied/virtual spinorbitals
      integer noas,nobs,nvas,nvbs ! # of occupied/virtual alphas betas
      integer nstot
cc      double precision ene_orb(nos+nvs)  ! a copy of dbl_mb(k_sorted)
      integer size
      integer rtdb
      integer l_aux1,k_aux1 ! for local memory allocator: loc. mem. buffer 1
      integer i,j,k,l,m,n         ! auxiliary indices
      integer ia,ib,ic,id,ie      !auxiliary indices
      integer iaux
c      integer nact              !number of active virtual orbitals
c h created here from "orbital" matrix horb
      double precision h(nos+nvs,nos+nvs)
c horb valid only for RHF case
      double precision horb((nos+nvs)/2,(nos+nvs)/2)
      double precision v(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
c FCI dimension
      integer qnumber
      integer irrep_g
      integer dim_fci
      integer alpha_str
      integer l_str,k_str
      integer k_matrix,l_matrix
      logical sym_test
c
c
c
      double precision xxx,yyy,zzz
c
      double precision repulsion
c
      logical nodezero
      logical oprint_qa
c
c
      nodezero = (ga_nodeid().eq.0)
      oprint_qa=util_print('ducc_qa', print_high)
c
c
c
c
      do i=1,(nos+nvs)/2
      do j=1,(nos+nvs)/2
        horb(i,j)=0.0d0
      enddo
      enddo
c
      do i=1,nos+nvs
      do j=1,nos+nvs
        h(i,j)=0.0d0
      enddo
      enddo
c
      call kinetic_hcore_1(rtdb,horb,h,nos,nvs,noas,nobs,nvas,nvbs)
c
c     
c forming v matrix
c
      do i=1,nos+nvs
      do j=1,nos+nvs
      do k=1,nos+nvs
      do l=1,nos+nvs
        v(i,j,k,l)=0.0d0 
      enddo
      enddo
      enddo
      enddo
c
      call mapping_v2_m(rtdb,d_v2,k_v2_offset,v,nos,nvs)
c
c
c
c  you have everything here for FCI: v & h are ready
c
c  ATTENTION:   Works only for the RHF case
c
c
c
c  Matrix dimensions - simplest case Ms symmetry binom(n,k)*binom(n,k) 
c  n=noas+nvas / k=noas
c
c  
        irrep_g = 0
        dim_fci = 0
        alpha_str = 0
        sym_test = .true.
c
      if (.not.MA_PUSH_GET(mt_int,(noas+nvas),'a_string',l_str,
     1 k_str)) CALL ERRQUIT('a_string',0,MA_ERR)
c
        qnumber = 0
        do m = 0, noas+nvas-1
         qnumber = qnumber + 2**m
        enddo 
c
        do i=1,qnumber   !main loop 
c
         do m=1,noas+nvas
          int_mb(k_str+m-1) = 0
         enddo
c         
         k=i
c
         do j = noas+nvas, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str+j-1) = 0
            else
             int_mb(k_str+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of alpha particle check
         m=0
         do j=1,noas+nvas
c additional symmetry checks possible here
          m = m+int_mb(k_str+j-1)
         enddo
         if(m.eq.noas) alpha_str = alpha_str + 1
c
        enddo !main loop
c
c
c        Determination of dim_fci with symmetry
c        returns dim_fci (with symmetry inclusion)
         call dim_det(noas,nobs,nvas,nvbs,nos,nvs,nstot,
     1        irrep_g,alpha_str,dim_fci)
c
c
      if (.not.MA_POP_STACK(l_str))
     1   call errquit("l_str",0,MA_ERR)
c
      if(nodezero) then
        write(6,*)'alpha strings: ',alpha_str
        write(6,*)'dim_fci: ',dim_fci
        call util_flush
      endif
c
      if (.not.MA_PUSH_GET(mt_dbl,dim_fci*dim_fci,'fcim',
     1      l_matrix,k_matrix))
     1      call errquit('fcim matrix',1,MA_ERR)
c
c 
c initialization
c
       do i=1,dim_fci*dim_fci
         dbl_mb(k_matrix+i-1) = 0.0d0
       enddo
c
c
c
c *** debug ***
       if(nodezero) then
        write(6,*)'before fci_matrix'
        call util_flush(6)
       endif
c *************
c
       call fci_matrix(rtdb,
     &      d_t1,k_t1_offset,d_t2,k_t2_offset,
     &      dim_fci,alpha_str,dbl_mb(k_matrix),h,v,
     &      nos,nvs,noas,nobs,nvas,nvbs,nstot,repulsion,irrep_g)
c
c
c *** debug ***
       if(nodezero) then
        write(6,*)'after fci_matrix'
        call util_flush(6)
       endif
c *************
c
c
      if(.not.MA_POP_STACK(l_matrix))
     &      call errquit('rrst3',4,MA_ERR)
c
c
c
c
c
      return 
      end
c
c
c
c
c ------------- linearized approach RRST(1)
c
c
c
c
c
       subroutine fci_matrix(rtdb,
     &      d_t1,k_t1_offset,d_t2,k_t2_offset, 
     &      dim_fci,alpha_str,matrix,h,v,
     &      nos,nvs,noas,nobs,nvas,nvbs,nstot,repulsion,irrep_g)
c
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c auxiliary integers
        integer i,j,k,l,m,n,iaux
        integer p,q,r,s,t,u,w,x,y,z
        integer ipl1,ipl2
c
        integer rtdb
        integer d_t1,d_t2
        integer k_t1_offset,k_t2_offset
        integer nstot
        integer dim_fci
        integer alpha_str
        double precision h(nstot,nstot)
        double precision v(nstot,nstot,nstot,nstot)
        double precision matrix(dim_fci,dim_fci)
        double precision emt(dim_fci,dim_fci)
        double precision ept(dim_fci,dim_fci)
        double precision trial(dim_fci)
        double precision vec1(dim_fci),vec2(dim_fci),vec3(dim_fci)
        double precision xxx,yyy,zzz
        double precision repulsion
        double precision my_rhf
c *** debug for H4 ****
        double precision m1(dim_fci,dim_fci),m2(dim_fci,dim_fci)
        double precision xa,xb,xc,xd,aaa,bbb,ccc,xroot1,xroot2
        double precision xxdel
c *********************
c
        integer nos,nvs,noas,nobs,nvas,nvbs
        integer irrep_g
        integer maxm
c
        integer qnumber
        integer k_str,l_str
        integer half_str(alpha_str,noas+nvas)
        integer strings(dim_fci,nstot)
        integer string1(nstot)
        integer string2(nstot)
        integer string(nstot)
        integer string_diff(nstot)
        integer ind_set(2*noas)
        integer dim_fci_test
        integer iphase
        integer isum
c active
        integer dim_act
        integer actoi(6)
        integer actspin(nstot)
c
        logical nodezero
        logical sym_test
c  eigensolver
        double precision wr(dim_fci),wi(dim_fci)
        double precision vl,vr,work(3*dim_fci)
        integer info
c
c
        sym_test=.true.
        nodezero = (ga_nodeid().eq.0)
c
c
c
        my_rhf = 0.0d0
        do i=1,noas+nobs
         my_rhf = my_rhf + h(i,i)
        enddo
        do i=1,noas+nobs
        do j=1,noas+nobs
         my_rhf = my_rhf + (0.50d0)*v(i,j,i,j)
        enddo
        enddo
c
        my_rhf = my_rhf + repulsion
c
c
        if(nodezero) then 
         write(6,*)'nos: ',nos
         write(6,*)'nvs: ',nvs
         write(6,*)'noas: ',noas
         write(6,*)'nobs: ',nobs
         write(6,*)'nvas: ',nvas
         write(6,*)'nvbs" ',nvbs
         write(6,*)'nstot = ',nstot
         write(6,*)'dim_fci = ',dim_fci
         write(6,*)'alpha_str = ',alpha_str
         write(6,*)'repulsion = ',repulsion
         write(6,*)'my_rhf = ',my_rhf
         call util_flush(6)
        endif
c
        do i=1,noas+nvas
        do j=1,alpha_str
         half_str(j,i) = 0
        enddo
        enddo
c
        do i=1,nstot
        do j=1,dim_fci
         strings(j,i) = 0
        enddo
        enddo
c
        do i=1,dim_fci
        do j=1,dim_fci
         matrix(j,i)=0.0d0
        enddo
        enddo
c
c
c
        sym_test = .true.
        dim_fci_test = 0
        alpha_str_test = 0
c
      if (.not.MA_PUSH_GET(mt_int,(noas+nvas),'a_string',l_str,
     1 k_str)) CALL ERRQUIT('a_string',0,MA_ERR)
c ----- half_str starts here --------------------
        qnumber = 0
        do m = 0, noas+nvas-1
         qnumber = qnumber + 2**m
        enddo
c
        do i=1,qnumber   !main loop 
c
         do m=1,noas+nvas
          int_mb(k_str+m-1) = 0
         enddo
c         
         k=i
c
         do j = noas+nvas, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str+j-1) = 0
            else
             int_mb(k_str+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of alpha particle check
         m=0
         do j=1,noas+nvas
c additional symmetry checks possible here
          m = m+int_mb(k_str+j-1)
         enddo
         if(m.eq.noas) then
c additional symmetry checks possible here
           alpha_str_test = alpha_str_test + 1
           do n=1,noas+nvas
            half_str(alpha_str_test,n) = 
     1        int_mb(k_str+n-1)    
           enddo
         endif
c
       enddo !main loop
c
      if (.not.MA_POP_STACK(l_str))
     1   call errquit("l_str",0,MA_ERR)
c
      if(alpha_str_test.ne.alpha_str) 
     1   call errquit("dim_fci_test",0,MA_ERR)

c   *** debug ***
c       if(nodezero) then 
c        do k=1,alpha_str
c         write(6,220) half_str(k,1),half_str(k,2),half_str(k,3),
c     1                half_str(k,4)
c        enddo
c       endif
c 220   format(4i5)
c   *************

c  ------- done with half_str -----------------------
c
c
c  ------- total number of strings (only valid fro RHF ref.)
c  
      m = 0
      do i=1,alpha_str
      do j=1,alpha_str
c 
       do l=1,nstot
        string(l)=0
       enddo
c
       do l=1,noas
        string(l)=half_str(i,l)
       enddo
       do l=1,nvas
       string(l+nos)=half_str(i,l+noas)
       enddo
       do l=1,noas
        string(l+noas)=half_str(j,l)
       enddo
       do l=1,nvas
        string(l+nos+nvas)=half_str(j,l+noas)
       enddo
c
         do l=1,2*noas
          ind_set(l)=0
         enddo
         k=1
         do l=1,nstot
          if(string(l).eq.1) then
           ind_set(k)=l
           k=k+1
          endif
         enddo
c
         irrep1=0
         do k=1,2*noas
          iaux = ind_set(k)
          irrep1 = ieor(irrep1,int_mb(k_irs_sorted+iaux-1))
         enddo
         if(irrep1.eq.irrep_g) then
          m = m + 1
          do k=1,nstot
           strings(m,k)=string(k)
          enddo
         endif
c
      enddo !j
      enddo !i
c
c
c *** debug ****
c      if(nodezero) then 
c       write(6,*)'full strings'
c       do i=1,dim_fci
c        write(6,224)(strings(i,j),j=1,nstot)
c       enddo
c       call util_flush(6)
c      endif
c 224  format(8i5)
c **************

 
c
c ------- done with strings list --------------------
c
c upper half only
c
      do i=1,dim_fci  !main loop1
      do j=i,dim_fci  !main loop2
c
       idiff=0
       do k=1,nstot
        string1(k)=strings(i,k)
        string2(k)=strings(j,k)
        string_diff(k)=string1(k)-string2(k)
        idiff=idiff+abs(string_diff(k))
       enddo
c
       idiff=idiff/2
c ---------------------
       if(idiff.gt.2) then
        matrix(i,j)=0.0d0
        matrix(j,i)=0.0d0
       endif
c ---------------------
       if(idiff.eq.2) then
        p=0
        q=0
        r=0
        s=0
c       v^{pq}_{rs} a_p^+ a_q^+ a_s a_r (p<q & r<s)
        do k=1,nstot !k ----
         if(string_diff(k).eq.1) then
          if(p.eq.0) then
           p=k
          else  
           q=k
          endif
         endif
         if(string_diff(k).eq.-1) then
          if(r.eq.0) then
           r=k
          else
           s=k
          endif
         endif
        enddo !k ---
c
        iphase = 1
c r part 
        isum=0
        do k=1,r-1
         isum=isum+string2(k) 
        enddo 
        string2(r)=0
c s part
        do k=1,s-1
         isum=isum+string2(k) 
        enddo
        string2(s)=0
c q part
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=1
c p part
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
        iphase=(-1)**isum
c      
        matrix(i,j)=matrix(i,j)+dble(iphase)*v(p,q,r,s) 
        matrix(j,i)=matrix(j,i)+dble(iphase)*v(p,q,r,s)
c check
        do k=1,nstot
         if((string1(k)-string2(k)).ne.0) 
     1   call errquit("string_1_2_test",0,MA_ERR)
        enddo  
c
       endif !idiff=2
c ------------------------------------------------
       if(idiff.eq.1) then
        p=0
        q=0
        do k=1,nstot !k ----
         if(string_diff(k).eq.1) p=k
         if(string_diff(k).eq.-1) q=k
        enddo
c
         k=1
         do l=1,nstot
          if(string2(l).eq.1) then
           ind_set(k)=l
           k=k+1
          endif
         enddo
c
        iphase = 1
c q part 
        isum=0
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=0
c p part        
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
         iphase=(-1)**isum
c
         xxx = 0.0d0
         xxx = xxx + h(p,q)*dble(iphase)
         do m=1,2*noas 
          t=ind_set(m)
          xxx = xxx + v(p,t,q,t)*dble(iphase) 
         enddo
c
         matrix(i,j) = matrix(i,j) + xxx
         matrix(j,i) = matrix(j,i) + xxx
c 
         endif !idiff=1
c ------------------------------------------------
       if(idiff.eq.0) then
         k=1
         do l=1,nstot
          if(string1(l).eq.1) then
           ind_set(k)=l
           k=k+1
          endif
         enddo
c *** debug ***
      if(k-1.ne.2*noas)
     1   call errquit("2noas_test",0,MA_ERR)
c *************
         do k=1,2*noas
          p=ind_set(k)
          matrix(i,i) = matrix(i,i) + h(p,p)
          do l=1,2*noas
           q=ind_set(l) 
           matrix(i,i) = matrix(i,i) + (0.50d0)*v(p,q,p,q)
          enddo
         enddo
         matrix(i,i) = matrix(i,i) + repulsion
c *** debug ***
c      if(nodezero) then 
c       write(6,*)'diagonal elem. ',i,matrix(i,i)
c       call util_flush(6)
c      endif
c *************
       endif !idiff=0
c ------------------------------------------------
c
      enddo           !main loop2
      enddo           !main loop1
c
c  *** debug ****
c      if(nodezero) then
c       write(6,*)' matrix diag '
c       do i=1,dim_fci
c        write(6,223) matrix(i,i),i,i
c       enddo
c       call util_flush(6)
c      endif 
c 223  format(f16.6,2x,2i5)
c ***************
c
c  Moments evaluation
c
        maxm=20
        do m=1,dim_fci
         trial(m) = 0.0d0
        enddo
        do i=1,dim_fci
         n=0
         do j=1,2*noas
          n = n + strings(i,j)
         enddo
         if(n.eq.2*noas) then
          trial(i)=1.0d0
          go to 30
         endif
        enddo
c
 30     continue
c
        call moments(matrix,trial,dim_fci,maxm)
c form now iact aact
        iact = 1
        aact = 2
c
        call  texponents(rtdb,ept,emt,
     1     d_t1,k_t1_offset,d_t2,k_t2_offset,
     1     nos,nvs,noas,nobs,nvas,nvbs,nstot,
     1     dim_fci,strings,v,iact,aact)
c e^{-T} H e^{T}
c        do i=1,dim_fci
c         vec1(i) = 0.0d0
c         vec2(i) = 0.0d0
c         vec3(i) = 0.0d0
c        enddo
cc
c        do i=1,dim_fci
c        do j=1,dim_fci
c         vec1(i) = vec1(i) + ept(i,j)*trial(j)
c        enddo
c        enddo
cc  
c        do i=1,dim_fci
c        do j=1,dim_fci
c         vec2(i) = vec2(i) + matrix(i,j)*vec1(j)
c        enddo
c        enddo
cc
c        do i=1,dim_fci
c        do j=1,dim_fci
c         vec3(i) = vec3(i) + emt(i,j)*vec2(j)
c        enddo
c        enddo
cc
c        if(nodezero) then
c         write(6,*)'--- moments of the CC equations ---'
c         do i=1,dim_fci
cc debug for H4 ***
c          write(6,110) vec3(i),(strings(i,j),j=1,8)
c         enddo
cc ****************
c         call util_flush(6)
c        endif
c 110    format(f20.10)
c e^{-T_ext}He^{T_ext}
        do i=1,dim_fci
        do j=1,dim_fci
         m1(i,j)=0.0d0
         m2(i,j)=0.0d0
        enddo
        enddo
        do i=1,dim_fci
        do j=1,dim_fci
        do k=1,dim_fci
         m1(i,k) = m1(i,k) + matrix(i,j)*ept(j,k)
        enddo
        enddo
        enddo
        do i=1,dim_fci
        do j=1,dim_fci
        do k=1,dim_fci
         m2(i,k) = m2(i,k) + emt(i,j)*m1(j,k)
        enddo
        enddo
        enddo
c *** debug for H4 ****
        xa=m2(1,1)
        xb=m2(1,15)
        xc=m2(15,1)
        xd=m2(15,15)
        aaa=1.0d0
        bbb=(-1.0d0)*(xa+xd)
        ccc=(-1.0d0)*xb*xc+xa*xd
        xxdel=bbb*bbb-(4.0d0)*aaa*ccc
        xroot1=((-1.0d0)*bbb-dsqrt(xxdel))/(2.0d0*aaa)
        xroot2=((-1.0d0)*bbb+dsqrt(xxdel))/(2.0d0*aaa)
        if(nodezero) then
         write(6,*)'SES-CC 2 & 3 active orbitals'
         write(6,*)'xa ',xa
         write(6,*)'xb ',xb
         write(6,*)'xc ',xc
         write(6,*)'xd ',xd
         write(6,*)'xroot1: ',xroot1
         write(6,*)'xroot2: ',xroot2
         xxx=(xroot2-xroot1)*27.211390d0
         write(6,*)'excit. ene. : ',xxx,' eV'
         call util_flush(6)
        endif

c  ********************

c
c diagonalization
c

      call util_dgeev('n','n',dim_fci,matrix,dim_fci,
     1  wr,wi,vl,1,vr,1,
     2  work,3*dim_fci,info)
c
c  *** debug ****
c      if(nodezero) then
c       write(6,*)' wr       wi '
c       do i=1,dim_fci
c        write(6,222) wr(i),wi(i)
c       enddo
c       call util_flush(6)
c      endif 
c 222  format(2f16.6)
c ***************
c
      xxx=wr(1)
      do i=1,dim_fci
       if(wr(i).lt.xxx) xxx=wr(i)
      enddo
c
      if(nodezero) then
      write(6,*)'----------'
      write(6,*)'lowest FCI energy: ',xxx
      write(6,*)'----------'
      call util_flush(6)
      endif
c
c      if(nodezero) then
c      write(6,*)'----------'
c      write(6,*)'Excitation energies (Hartree) (eV) '
c       do i=1,dim_fci
c        write(6,303) wr(i)-xxx,(wr(i)-xxx)*27.211390d0
c       enddo
c      write(6,*)'----------'
c      call util_flush(6)
c      endif
c
 303  format(f15.6,4x,f15.6)
c
        return 
        end
c
c
c
c
c
c    
c
      subroutine dim_det(noas,nobs,nvas,nvbs,nos,nvs,nstot,
     1           irrep_g,alpha_str,dim_fci)
c
c
c nstot - total number of spinorbitals nstot=nos+nvs
c nos   - number of occupied spinorbitals
c nvs   - number of virtual spinorbitals
c noas  - number of occupied alpha spinorbitals
c nobs  - number of occupied beta spinorbitals
c nvas  - number of virtual  alhoa spinorbitals
c nvbs  - number of virtual beta spinorbitals
c
c spinorbital convention
c
c | noas | nobs | nvas | nvbs |
c
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
        integer nos,nvs,noas,nobs,nvas,nvbs,nstot
        integer i,j,k,l,m,n,iaux
        integer qnumber
        integer alpha_str,alpha_str_test
        integer dim_fci
        integer half_str(alpha_str,noas+nvas)
        integer k_str,l_str
        integer ind_set(2*noas)
        integer string(nstot)
        integer irrep1,irrep_g
c
        dim_fci = 0
        alpha_str_test = 0
c
      if (.not.MA_PUSH_GET(mt_int,(noas+nvas),'a_string',l_str,
     1 k_str)) CALL ERRQUIT('a_string',0,MA_ERR)
c ----- half_str starts here --------------------
        qnumber = 0
        do m = 0, noas+nvas-1
         qnumber = qnumber + 2**m
        enddo
c
        do i=1,qnumber   !main loop 
c
         do m=1,noas+nvas
          int_mb(k_str+m-1) = 0
         enddo
c         
         k=i
c
         do j = noas+nvas, 1, -1 !bin. dec.
            l=k/(2**(j-1))
            if(l.eq.0) then
             int_mb(k_str+j-1) = 0
            else
             int_mb(k_str+j-1) = 1
             k=k-2**(j-1)
            endif
         enddo !bin. dec.
c numer of alpha particle check
         m=0
         do j=1,noas+nvas
c additional symmetry checks possible here
          m = m+int_mb(k_str+j-1)
         enddo
         if(m.eq.noas) then
c additional symmetry checks possible here
           alpha_str_test = alpha_str_test + 1
           do n=1,noas+nvas
            half_str(alpha_str_test,n) =
     1        int_mb(k_str+n-1)
           enddo
         endif
c
       enddo !main loop
c
      if (.not.MA_POP_STACK(l_str))
     1   call errquit("l_str",0,MA_ERR)
c
      if(alpha_str_test.ne.alpha_str)
     1   call errquit("dim_fci_test",0,MA_ERR)
c



      dim_fci = 0
      do i=1,alpha_str
      do j=1,alpha_str
c 
       do l=1,nstot
        string(l)=0
       enddo
c
       do l=1,noas
        string(l)=half_str(i,l)
       enddo
       do l=1,nvas
       string(l+nos)=half_str(i,l+noas)
       enddo
       do l=1,noas
        string(l+noas)=half_str(j,l)
       enddo
       do l=1,nvas
        string(l+nos+nvas)=half_str(j,l+noas)
       enddo
c
         do l=1,2*noas
          ind_set(l)=0
         enddo
         k=1
         do l=1,nstot
          if(string(l).eq.1) then
           ind_set(k)=l
           k=k+1
          endif
         enddo
c
         irrep1=0
         do k=1,2*noas
          iaux = ind_set(k)  
          irrep1 = ieor(irrep1,int_mb(k_irs_sorted+iaux-1))
         enddo
         if(irrep1.eq.irrep_g) dim_fci = dim_fci+1
      enddo !j
      enddo !i
c
        return
        end
c
c
c
c
      subroutine moments(h,x,dim_fci,maxm)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      double precision h(dim_fci,dim_fci),x(dim_fci)
      integer dim_fci,maxm
      double precision y(dim_fci),z(dim_fci)
      double precision xxx,yyy,zzz
      double precision vmom(maxm)
      integer i,j,k,l,m,n
c
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
c
      do i=1,maxm
       vmom(i) = 0.0d0
      enddo 
c
      do i=1,dim_fci
       y(i) = x(i)
      enddo
c
      do i=1,maxm
c
       do m=1,dim_fci
        z(m) = 0.0d0
       enddo
       do m=1,dim_fci
       do n=1,dim_fci
        z(m) = z(m) + h(m,n)*y(n)
       enddo
       enddo
       do m=1,dim_fci
        vmom(i) = vmom(i) + x(m)*z(m)
       enddo
c 
       do m=1,dim_fci
        y(m)=z(m)
       enddo
c
      enddo !maxm      
c
       if(nodezero) then
        write(6,*)'-------- Moments --------'
         do m=1,maxm
          write(6,10) m,vmom(m)
         enddo 
        write(6,*)'-------------------------'
        call util_flush(6)
       endif 
c
 10    format('K(',i3,')',4x,f16.6)
c
      return 
      end
c
c
c
c in texponents added #include "tce_amps.fh"
c v - only for testing purposes: once tested get rid of it.
c
      subroutine  texponents(rtdb,ept,emt,
     1     d_t1,k_t1_offset,d_t2,k_t2_offset,
     1     nos,nvs,noas,nobs,nvas,nvbs,nstot,
     1     dim_fci,strings,v,iact,aact) 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
        integer rtdb
        integer d_t1,d_t2
        integer k_t1_offset,k_t2_offset
        integer dim_fci
        double precision ept(dim_fci,dim_fci)
        double precision emt(dim_fci,dim_fci)
        double precision tm(dim_fci,dim_fci)
        double precision m1(dim_fci,dim_fci)
        double precision m2(dim_fci,dim_fci) 
        double precision v(nstot,nstot,nstot,nstot)
        integer nos,nvs,noas,nobs,nvas,nvbs,nstot
        integer strings(dim_fci,nstot)
        integer iact !active occupied orbital
        integer aact !active virtual orbital
        integer nfact
c
c t2 in a nice representation : be careful here (you may want to reindex
c virtuals
        double precision t1(nos,nos+1:nos+nvs)
        double precision t2(nos,nos,nos+1:nos+nvs,nos+1:nos+nvs)        
        double precision ft1(nstot,nstot)
        double precision ft2(nstot,nstot,nstot,nstot)
        double precision eccsd
        integer inda(nstot)
c
        integer string1(nstot)
        integer string2(nstot)
        integer string(nstot)
        integer string_diff(nstot)
c
        integer p,q,r,s,t,u,w,x,y,z
        integer i,j,k,l,m,n,iaux
        integer ia,ib,ic,id
        integer idiff,isum,iphase
        double precision xxx,yyy,zzz
c
        logical nodezero
c
c
         nodezero = (ga_nodeid().eq.0)
c
c
       
      do i=1,nos
      do ia=nos+1,nos+nvs
        t1(i,ia)=0.0d0
      enddo
      enddo
c
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        t2(i,j,ia,ib)=0.0d0
      enddo
      enddo
      enddo
      enddo
c 
c
      call mapping_t1(rtdb,d_t1,k_t1_offset,t1,nos,nvs)
c
c
      call mapping_t2(rtdb,d_t2,k_t2_offset,t2,nos,nvs)
c
c t1 -> ft1 / t2 -> ft2
c
c changing index ordering in ft1 and ft2 consistent with
c creation/annihilation operators ordering
c
      do i=1,nstot
      do j=1,nstot
       ft1(i,j) = 0.0d0
      enddo
      enddo
c
      do i=1,nstot
      do j=1,nstot
      do k=1,nstot
      do l=1,nstot
       ft2(i,j,k,l) = 0.0d0
      enddo
      enddo
      enddo
      enddo
c
      do i=1,nos
      do ia=nos+1,nos+nvs
        ft1(ia,i) = t1(i,ia)
      enddo
      enddo
c
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        ft2(ia,ib,i,j) =  t2(i,j,ia,ib)
      enddo
      enddo
      enddo
      enddo
c
c check point: calculating CCSD correlation energy 
c 
      eccsd=0.0d0
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        eccsd=eccsd+0.250d0*v(i,j,ia,ib)*t2(i,j,ia,ib)
        eccsd=eccsd+0.50d0*v(i,j,ia,ib)*t1(i,ia)*t1(j,ib)
      enddo
      enddo
      enddo
      enddo
c 
      if(nodezero) then
        write(6,*)'From  CCSD corr. ene.',eccsd
        call util_flush(6)
      endif
c
c    inda =1 for external indices
c
      do i=1,nstot
       inda(i)=1
      enddo
c attention: works only with the the c1 symmetry
c            fix is for any type of symetry

c *** debug for H4 ***
      iact = 1
      aact = 3
c ********************
 
      inda(iact) = 0
      inda(iact+nvas) = 0
      inda(nos+(aact-noas)) = 0
      inda(nos+nvas+(aact-noas)) = 0

c *** debug for H4 ***
      if(nodezero) then
      write(6,*)'inda vector - length nstot'
      write(6,112)(inda(j),j=1,8)
      call util_flush(6)
      endif
 112  format(8i5)
c ********************

  
c
c done with indact: external Ts operators
c
      do i=1,nos
      do ia=nos+1,nos+nvs
        iaux = inda(i) + inda(ia)
        if(iaux.eq.0) then 
         write(6,*)'canc.',ia,i,ft1(ia,i)
         ft1(ia,i) = 0.0d0
        endif
c *** debug for H4 ***
        if(nodezero) then
        if(iaux.eq.0) write(6,*)'ext.Singles i ia ',i,ia
        call util_flush(6)
        endif
c ********************
      enddo
      enddo      
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
       iaux = inda(i) + inda(j) + inda(ia) + inda(ib)
       if(iaux.eq.0) then
          write(6,*)'canc.',ia,ib,i,j,ft2(ia,ib,i,j)
          ft2(ia,ib,i,j) =  0.0d0
       endif
c *** debug for H4 ***
        if(nodezero) then
            if(iaux.eq.0) write(6,*)'ext.Doubles ij iaib ',i,j,ia,ib
        call util_flush(6)
        endif
c ********************
      enddo
      enddo
      enddo
      enddo
c
c Formation of the T_ext  matrix
c
      do i=1,dim_fci  
      do j=1,dim_fci  
       tm(i,j) = 0.0d0
      enddo
      enddo
c
      do i=1,dim_fci  !main loop1
      do j=1,dim_fci  !main loop2
c
       idiff=0
       do k=1,nstot
        string1(k)=strings(i,k)
        string2(k)=strings(j,k)
        string_diff(k)=string1(k)-string2(k)
        idiff=idiff+abs(string_diff(k))
       enddo
c
       idiff=idiff/2
c ---------------------
       if(idiff.gt.2) then
        tm(i,j)=0.0d0
       endif
c ---------------------
       if(idiff.eq.2) then
        p=0
        q=0
        r=0
        s=0
c       t^{ab}_{ij} a_p^+ a_q^+ a_s a_r (p<q & r<s)
        do k=1,nstot !k ----
         if(string_diff(k).eq.1) then
          if(p.eq.0) then
           p=k
          else
           q=k
          endif
         endif
         if(string_diff(k).eq.-1) then
          if(r.eq.0) then
           r=k
          else
           s=k
          endif
         endif
        enddo !k ---
c
        iphase = 1
c r part 
        isum=0
        do k=1,r-1
         isum=isum+string2(k)
        enddo
        string2(r)=0
c s part
        do k=1,s-1
         isum=isum+string2(k)
        enddo
        string2(s)=0
c q part
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=1
c p part
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
        iphase=(-1)**isum
c      
        tm(i,j)=tm(i,j)+dble(iphase)*ft2(p,q,r,s)
c
c check
        do k=1,nstot
         if((string1(k)-string2(k)).ne.0)
     1   call errquit("string_1_2_test",0,MA_ERR)
        enddo
c
       endif !idiff=2
c ------------------------------------------------
       if(idiff.eq.1) then
        p=0
        q=0
        do k=1,nstot !k ----
         if(string_diff(k).eq.1) p=k
         if(string_diff(k).eq.-1) q=k
        enddo
c
        iphase = 1
c q part 
        isum=0
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=0
c p part        
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
         iphase=(-1)**isum
c
c q part 
        isum=0
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=0
c p part        
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
         iphase=(-1)**isum
c
         tm(i,j) = tm(i,j) + ft1(p,q)*dble(iphase)
c 
         endif !idiff=1
c ------------------------------------------------
      enddo           !main loop2
      enddo           !main loop1
c calculating ept,emt
      do i=1,dim_fci
      do j=1,dim_fci
       if(i.eq.j) then
         ept(i,i) = 1.0d0
         emt(i,i) = 1.0d0
       else
         ept(i,j) = 0.0d0
         ept(i,j) = 0.0d0
       endif
      enddo
      enddo
c
      do i=1,dim_fci
      do j=1,dim_fci
       m1(i,j) = tm(i,j)
       m2(i,j) = 0.0d0
      enddo
      enddo  
c
      do i=1,dim_fci
      do j=1,dim_fci
       ept(i,j) = ept(i,j) + m1(i,j)
       emt(i,j) = emt(i,j) - m1(i,j)
      enddo
      enddo
c
c
      nfact=1
c
      do k=1,2*noas-1
c
      nfact = nfact*(k+1)
c
       do l=1,dim_fci
        do m=1,dim_fci
         do n=1,dim_fci
          m2(l,m)=m2(l,m)+tm(l,n)*m1(n,m)
         enddo
        enddo
       enddo
c
       xxx = ((1.0d0)/(dble(nfact)))
       yyy = dble((-1)**(k+1))
c
       do m=1,dim_fci
       do n=1,dim_fci
        ept(m,n) = ept(m,n) + xxx*m2(m,n)
        emt(m,n) = emt(m,n) + xxx*yyy*m2(m,n)
       enddo
       enddo
c
       do l=1,dim_fci
       do m=1,dim_fci
        m1(l,m) = m2(l,m)
        m2(l,m) = 0.0d0
       enddo 
       enddo
c
      enddo !k 2*noa
c
c Test
c number of electrons in the system = 2*noas
c tm^{2*noas+1} should be equal to zero: tm is a nilpotent matrix
c
c      do i=1,dim_fci
c      do j=1,dim_fci
c       m1(i,j) = tm(i,j)
c       m2(i,j) = 0.0d0
c      enddo
c      enddo  
cc
c      do k=1,2*noas
cc
c       do l=1,dim_fci
c        do m=1,dim_fci
c         do n=1,dim_fci
c          m2(l,m)=m2(l,m)+tm(l,n)*m1(n,m)
c         enddo
c        enddo
c       enddo
cc
c       do l=1,dim_fci
c       do m=1,dim_fci
c        m1(l,m) = m2(l,m)
c        m2(l,m) = 0.0d0
c       enddo 
c       enddo
cc
c      enddo !k 2*noas
cc 
c      xxx = 0.0d0
cc
c      do i=1,dim_fci
c      do j=1,dim_fci
c       if(dabs(m1(i,j)).gt.xxx) xxx=dabs(m1(i,j))
c      enddo
c      enddo
cc
c      if(nodezero) then 
c        write(6,*)'---- max tm^(2noas+1) ---'
c        write(6,*)'max : ',xxx
c        write(6,*)'-------------------------'
c        call util_flush(6)
c      endif
cc
      return 
      end
c
c
c
c
c
c
c
c actoi - from fci_matrix
c actspin  - from fci_matrix 
c
      subroutine  act_space_size(rtdb,
     1     nos,nvs,noas,nobs,nvas,nvbs,nstot,
     1     strings,actoi,actspin,dim_fci,dim_act)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
        integer rtdb
        integer dim_fci, dim_act
        integer nos,nvs,noas,nobs,nvas,nvbs,nstot
        integer strings(dim_fci,nstot)
        integer string_aux(nstot)
        integer string_ref(nstot)
        integer string_diff(nstot)
        integer iact1,iact2,iact3,iact4,iact5,iact6 !max.6 act.orb.
        integer actoi(6)
        integer actspin(nstot) ! 0-active 1-inactive
        integer nacto !total number of active orbitals
c
        integer p,q,r,s,t,u,w,x,y,z
        integer i,j,k,l,m,n,iaux
        integer idiff,isum,iphase
c
        logical nodezero
c
        nodezero = (ga_nodeid().eq.0)
c
        do i=1,dim_fci
         n=0
         do j=1,noas+nobs
          n = n + strings(i,j)
         enddo
         if(n.eq.noas+nobs) then
          do m=1,nstot
           string_ref(m)=strings(i,m)
           go to 30
          enddo 
         endif
        enddo
c
 30     continue
c
        if(nodezero) then 
         write(6,*)'ref. string'
         do i=1,nstot
          write(6,100) string_ref(i)
         enddo
        endif
c
        do i=1,6
         actoi(i) = 0
        enddo
c
        iact1 = 0
        iact2 = 0 
        iact3 = 0
        iact4 = 0
        iact5 = 0
        iact6 = 0
        nacto = 0
c
        if (.not.rtdb_get(rtdb,'tce:act1',mt_int,1,iact1)) then
     1    call errquit('fci:iact1',0,RTDB_ERR)
        else
          if(iact1.ne.0) then
           nacto = nacto + 1 
           actoi(nacto) = iact1
          endif
        end if
c
        if (.not.rtdb_get(rtdb,'tce:act2',mt_int,1,iact2)) then
     1    call errquit('fci:iact2',0,RTDB_ERR)
        else
          if(iact2.ne.0) then
           nacto = nacto + 1
           actoi(nacto) = iact2
          endif
        end if
c
        if (.not.rtdb_get(rtdb,'tce:act3',mt_int,1,iact3)) then
     1    call errquit('fci:iact3',0,RTDB_ERR)
        else
          if(iact3.ne.0) then
           nacto = nacto + 1
           actoi(nacto) = iact3
          endif
        end if
c
        if (.not.rtdb_get(rtdb,'tce:act4',mt_int,1,iact4)) then
     1    call errquit('fci:iact4',0,RTDB_ERR)
        else
          if(iact4.ne.0) then
           nacto = nacto + 1
           actoi(nacto) = iact4
          endif
        end if
c
        if (.not.rtdb_get(rtdb,'tce:act5',mt_int,1,iact5)) then
     1    call errquit('fci:iact5',0,RTDB_ERR)
        else
          if(iact5.ne.0) then
           nacto = nacto + 1
           actoi(nacto) = iact5
          endif
        end if
c
        if (.not.rtdb_get(rtdb,'tce:act6',mt_int,1,iact6)) then
     1    call errquit('fci:iact6',0,RTDB_ERR)
        else
          if(iact6.ne.0) then
           nacto = nacto + 1
           actoi(nacto) = iact6
          endif
        end if
c
        if(nodezero) then 
         write(6,*)'actoi vector'
         do i=1,6
          write(6,100) actoi(i)
         enddo
         call util_flush(6)
        endif
c actspin
        do i=1,nstot
         actspin(i) = 1
        enddo
c 
        do i=1,nacto
         m = actoi(i)
         if(m.le.noas+nobs) then
          actspin(m) = 0
          actspin(nvas+m) = 0
         else
          actspin(nos+(m-noas)) = 0
          actspin(nos+nvas+(m-nobs)) = 0
        enddo
c dim_act 
        dim_act = 0
c
      do i=1,dim_fci  !main loop1
c
       do k=1,nstot
        string_aux(k)=strings(i,k)
        string_diff(k)=string_aux(k)-string_ref(k)
       enddo
c
       m = 0
       do k=1,nstot
         if(string_diff(k).ne.0) then
          m = m + actspin(k)
         endif
       enddo
c
       if(m.eq.0) dim_act = dim_act + 1
c
      enddo !main loop1
c
c
 100    format(i5)
c
        return 
        end
c
c
c
c
c        
c
c
c
c
c
c matrix === e^{-Text} H e^{Text}
c actoi - from fci_matrix
c actspin  - from fci_matrix 
c
      subroutine  act_matrix(rtdb,
     1     nos,nvs,noas,nobs,nvas,nvbs,nstot,
     1     strings,actspin,matrix,matrix_a,dim_fci,dim_act)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
        integer rtdb
        integer dim_fci,dim_act,dim_act_test
        integer nos,nvs,noas,nobs,nvas,nvbs,nstot
        integer strings(dim_fci,nstot)
cc        integer strings_act(dim_act,nstot) !not needed
        integer maps(dim_act) ! rows in full matrix
        integer string_aux(nstot)
        integer string_ref(nstot)
        integer string_diff(nstot)
        integer actspin(nstot) ! 0-active 1-inactive
c
        integer p,q,r,s,t,u,w,x,y,z
        integer i,j,k,l,m,n,iaux
        integer idiff,isum,iphase
        double precision xxx,yyy,zzz
c
        double precision matrix(dim_fci,dim_fci) !! e^-Text H e^Text
        double precision matrix_a(dim_act,dim_act)
c  eigensolver
        double precision wr(dim_act),wi(dim_act)
        double precision vl,vr,work(3*dim_act)
        integer info
c
        logical nodezero
c
        nodezero = (ga_nodeid().eq.0)
c
        do i=1,dim_act
         maps(i) = 0
        enddo
c
c dim_act 
        dim_act_test = 0
c
      do i=1,dim_fci  !main loop1
c
       do k=1,nstot
        string_aux(k)=strings(i,k)
        string_diff(k)=string_aux(k)-string_ref(k)
       enddo
c
       m = 0
       do k=1,nstot
         if(string_diff(k).ne.0) then
          m = m + actspin(k)
         endif
       enddo
c
       if(m.eq.0) then
        dim_act_test = dim_act_test + 1
        maps(dim_act_test) = i
       endif
c
      enddo !main loop1
c
      if(dim_act_test.ne.dim_act)
     1 call errquit("dim_act_test",0,MA_ERR)
c maping matrix --> matrix_a
      do i=1,dim_act
      do j=1,dim_act
       k=maps(i)
       l=maps(j)
       matrix_a(i,j) = matrix(k,l)
      enddo
      enddo
c
      call util_dgeev('n','n',dim_act,matrix_a,dim_act,
     1  wr,wi,vl,1,vr,1,
     2  work,3*dim_act,info)
c
c  *** debug ****
      if(nodezero) then
       write(6,*)' wr       wi '
       do i=1,dim_fci
        write(6,222) wr(i),wi(i)
       enddo
       call util_flush(6)
      endif 
 222  format(2f16.6)
c ***************
c
      xxx=wr(1)
      do i=1,dim_act
       if(wr(i).lt.xxx) xxx=wr(i)
      enddo
c
      if(nodezero) then
      write(6,*)'----------'
      write(6,*)'lowest active FCI energy: ',xxx
      write(6,*)'----------'
      call util_flush(6)
      endif
c
      return
      end



c
c
c
c
c
c
c
c
      subroutine  t_ext_exp(rtdb,ept,emt,
     1     d_t1,k_t1_offset,d_t2,k_t2_offset,
     1     nos,nvs,noas,nobs,nvas,nvbs,nstot,
     1     dim_fci,strings,v,iact,aact) 
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c 
c t1/t2 - contain CCSD amplitudes
c ft1 and ft2 - contain the external components only
c
        integer rtdb
        integer d_t1,d_t2
        integer k_t1_offset,k_t2_offset
        integer dim_fci
        double precision ept(dim_fci,dim_fci)
        double precision emt(dim_fci,dim_fci)
        double precision tm(dim_fci,dim_fci)
        double precision m1(dim_fci,dim_fci)
        double precision m2(dim_fci,dim_fci) 
        double precision v(nstot,nstot,nstot,nstot)
        integer nos,nvs,noas,nobs,nvas,nvbs,nstot
        integer strings(dim_fci,nstot)
        integer iact !active occupied orbital
        integer aact !active virtual orbital
        integer nfact
c
c t2 in a nice representation : be careful here (you may want to reindex
c virtuals
        double precision t1(nos,nos+1:nos+nvs)
        double precision t2(nos,nos,nos+1:nos+nvs,nos+1:nos+nvs)        
        double precision ft1(nstot,nstot)
        double precision ft2(nstot,nstot,nstot,nstot)
        double precision eccsd
        integer inda(nstot)
c
        integer string1(nstot)
        integer string2(nstot)
        integer string(nstot)
        integer string_diff(nstot)
c
        integer p,q,r,s,t,u,w,x,y,z
        integer i,j,k,l,m,n,iaux
        integer ia,ib,ic,id
        integer idiff,isum,iphase
        double precision xxx,yyy,zzz
c
        logical nodezero
c
c
         nodezero = (ga_nodeid().eq.0)
c
c
       
      do i=1,nos
      do ia=nos+1,nos+nvs
        t1(i,ia)=0.0d0
      enddo
      enddo
c
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        t2(i,j,ia,ib)=0.0d0
      enddo
      enddo
      enddo
      enddo
c 
c
      call mapping_t1(rtdb,d_t1,k_t1_offset,t1,nos,nvs)
c
c
      call mapping_t2(rtdb,d_t2,k_t2_offset,t2,nos,nvs)
c
c t1 -> ft1 / t2 -> ft2
c
c changing index ordering in ft1 and ft2 consistent with
c creation/annihilation operators ordering
c
      do i=1,nstot
      do j=1,nstot
       ft1(i,j) = 0.0d0
      enddo
      enddo
c
      do i=1,nstot
      do j=1,nstot
      do k=1,nstot
      do l=1,nstot
       ft2(i,j,k,l) = 0.0d0
      enddo
      enddo
      enddo
      enddo
c
      do i=1,nos
      do ia=nos+1,nos+nvs
        ft1(ia,i) = t1(i,ia)
      enddo
      enddo
c
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        ft2(ia,ib,i,j) =  t2(i,j,ia,ib)
      enddo
      enddo
      enddo
      enddo
c
c check point: calculating CCSD correlation energy 
c 
      eccsd=0.0d0
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        eccsd=eccsd+0.250d0*v(i,j,ia,ib)*t2(i,j,ia,ib)
        eccsd=eccsd+0.50d0*v(i,j,ia,ib)*t1(i,ia)*t1(j,ib)
      enddo
      enddo
      enddo
      enddo
c 
      if(nodezero) then
        write(6,*)'From  CCSD corr. ene.',eccsd
        call util_flush(6)
      endif
c
c    inda =1 for external indices
c
      do i=1,nstot
       inda(i)=1
      enddo
c attention: works only with the the c1 symmetry
c            fix is for any type of symetry

c *** debug for H4 ***
      iact = 1
      aact = 3
c ********************
 
      inda(iact) = 0
      inda(iact+nvas) = 0
      inda(nos+(aact-noas)) = 0
      inda(nos+nvas+(aact-noas)) = 0

c *** debug for H4 ***
      if(nodezero) then
      write(6,*)'inda vector - length nstot'
      write(6,112)(inda(j),j=1,8)
      call util_flush(6)
      endif
 112  format(8i5)
c ********************

  
c
c done with indact: external Ts operators
c
      do i=1,nos
      do ia=nos+1,nos+nvs
        iaux = inda(i) + inda(ia)
        if(iaux.eq.0) then 
         write(6,*)'canc.',ia,i,ft1(ia,i)
         ft1(ia,i) = 0.0d0
        endif
c *** debug for H4 ***
        if(nodezero) then
        if(iaux.eq.0) write(6,*)'ext.Singles i ia ',i,ia
        call util_flush(6)
        endif
c ********************
      enddo
      enddo      
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
       iaux = inda(i) + inda(j) + inda(ia) + inda(ib)
       if(iaux.eq.0) then
          write(6,*)'canc.',ia,ib,i,j,ft2(ia,ib,i,j)
          ft2(ia,ib,i,j) =  0.0d0
       endif
c *** debug for H4 ***
        if(nodezero) then
            if(iaux.eq.0) write(6,*)'ext.Doubles ij iaib ',i,j,ia,ib
        call util_flush(6)
        endif
c ********************
      enddo
      enddo
      enddo
      enddo
c
c Formation of the T_ext  matrix
c
      do i=1,dim_fci  
      do j=1,dim_fci  
       tm(i,j) = 0.0d0
      enddo
      enddo
c
      do i=1,dim_fci  !main loop1
      do j=1,dim_fci  !main loop2
c
       idiff=0
       do k=1,nstot
        string1(k)=strings(i,k)
        string2(k)=strings(j,k)
        string_diff(k)=string1(k)-string2(k)
        idiff=idiff+abs(string_diff(k))
       enddo
c
       idiff=idiff/2
c ---------------------
       if(idiff.gt.2) then
        tm(i,j)=0.0d0
       endif
c ---------------------
       if(idiff.eq.2) then
        p=0
        q=0
        r=0
        s=0
c       t^{ab}_{ij} a_p^+ a_q^+ a_s a_r (p<q & r<s)
        do k=1,nstot !k ----
         if(string_diff(k).eq.1) then
          if(p.eq.0) then
           p=k
          else
           q=k
          endif
         endif
         if(string_diff(k).eq.-1) then
          if(r.eq.0) then
           r=k
          else
           s=k
          endif
         endif
        enddo !k ---
c
        iphase = 1
c r part 
        isum=0
        do k=1,r-1
         isum=isum+string2(k)
        enddo
        string2(r)=0
c s part
        do k=1,s-1
         isum=isum+string2(k)
        enddo
        string2(s)=0
c q part
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=1
c p part
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
        iphase=(-1)**isum
c      
        tm(i,j)=tm(i,j)+dble(iphase)*ft2(p,q,r,s)
c
c check
        do k=1,nstot
         if((string1(k)-string2(k)).ne.0)
     1   call errquit("string_1_2_test",0,MA_ERR)
        enddo
c
       endif !idiff=2
c ------------------------------------------------
       if(idiff.eq.1) then
        p=0
        q=0
        do k=1,nstot !k ----
         if(string_diff(k).eq.1) p=k
         if(string_diff(k).eq.-1) q=k
        enddo
c
        iphase = 1
c q part 
        isum=0
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=0
c p part        
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
         iphase=(-1)**isum
c
c q part 
        isum=0
        do k=1,q-1
         isum=isum+string2(k)
        enddo
        string2(q)=0
c p part        
        do k=1,p-1
         isum=isum+string2(k)
        enddo
        string2(p)=1
c
         iphase=(-1)**isum
c
         tm(i,j) = tm(i,j) + ft1(p,q)*dble(iphase)
c 
         endif !idiff=1
c ------------------------------------------------
      enddo           !main loop2
      enddo           !main loop1
c calculating ept,emt
      do i=1,dim_fci
      do j=1,dim_fci
       if(i.eq.j) then
         ept(i,i) = 1.0d0
         emt(i,i) = 1.0d0
       else
         ept(i,j) = 0.0d0
         ept(i,j) = 0.0d0
       endif
      enddo
      enddo
c
      do i=1,dim_fci
      do j=1,dim_fci
       m1(i,j) = tm(i,j)
       m2(i,j) = 0.0d0
      enddo
      enddo  
c
      do i=1,dim_fci
      do j=1,dim_fci
       ept(i,j) = ept(i,j) + m1(i,j)
       emt(i,j) = emt(i,j) - m1(i,j)
      enddo
      enddo
c
c
      nfact=1
c
      do k=1,2*noas-1
c
      nfact = nfact*(k+1)
c
       do l=1,dim_fci
        do m=1,dim_fci
         do n=1,dim_fci
          m2(l,m)=m2(l,m)+tm(l,n)*m1(n,m)
         enddo
        enddo
       enddo
c
       xxx = ((1.0d0)/(dble(nfact)))
       yyy = dble((-1)**(k+1))
c
       do m=1,dim_fci
       do n=1,dim_fci
        ept(m,n) = ept(m,n) + xxx*m2(m,n)
        emt(m,n) = emt(m,n) + xxx*yyy*m2(m,n)
       enddo
       enddo
c
       do l=1,dim_fci
       do m=1,dim_fci
        m1(l,m) = m2(l,m)
        m2(l,m) = 0.0d0
       enddo 
       enddo
c
      enddo !k 2*noa
c
c Test
c number of electrons in the system = 2*noas
c tm^{2*noas+1} should be equal to zero: tm is a nilpotent matrix
c
c      do i=1,dim_fci
c      do j=1,dim_fci
c       m1(i,j) = tm(i,j)
c       m2(i,j) = 0.0d0
c      enddo
c      enddo  
cc
c      do k=1,2*noas
cc
c       do l=1,dim_fci
c        do m=1,dim_fci
c         do n=1,dim_fci
c          m2(l,m)=m2(l,m)+tm(l,n)*m1(n,m)
c         enddo
c        enddo
c       enddo
cc
c       do l=1,dim_fci
c       do m=1,dim_fci
c        m1(l,m) = m2(l,m)
c        m2(l,m) = 0.0d0
c       enddo 
c       enddo
cc
c      enddo !k 2*noas
cc 
c      xxx = 0.0d0
cc
c      do i=1,dim_fci
c      do j=1,dim_fci
c       if(dabs(m1(i,j)).gt.xxx) xxx=dabs(m1(i,j))
c      enddo
c      enddo
cc
c      if(nodezero) then 
c        write(6,*)'---- max tm^(2noas+1) ---'
c        write(6,*)'max : ',xxx
c        write(6,*)'-------------------------'
c        call util_flush(6)
c      endif
cc
      return 
      end
c
c
c
c
c
