      subroutine rrst(rtdb,
     &    d_v2,k_v2_offset,nos,nvs,noas,nobs,nvas,nvbs,nactv,
     &    nstot,
     &    h,v,xtol)
c
c  Only for RHF
c
c all varaible below have to be defined in the tce_energy.f
c the nact should be defined in the input set tce:nactv
c
c All occupied orbitals are treated as active
c nactv defines the number of active virtual orbitals
c please do not freeze the occupied orbitals for now
c
c nactv - number of active virtual orbitals
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
      integer nos,nvs ! # of occupied/virtual spinorbitals
      integer noas,nobs,nvas,nvbs ! # of occupied/virtual alphas betas
      integer nactv,nstot
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
c RRST 
      double precision v1(nos+nvs,nos+nvs)
      double precision h0(nos+nvs,nos+nvs)
      integer iext(nos+nvs)  ! iext(i)=0/1  (0 i=active) (1 i=ext)
      integer vnext   !number of virtual externa orb.
      integer iext1(nstot,nstot)
      integer iext2(nstot,nstot,nstot,nstot)
      integer ispin(nstot)
      integer isize1,isize2,isize3,isize4
      integer isize1_e1,isize1_e2
      integer isize2_e1,isize2_e2,isize2_e3,isize2_e4
      integer lin_size 
c
c      integer ext3(nstot,nstot,nstot,nstot,nstot,nstot)
c v2 <===> v(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
      double precision a1_1o(nos+nvs,nos+nvs)
      double precision a2_1o(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
      double precision a1_2o(nos+nvs,nos+nvs)
      double precision a2_2o(nstot,nstot,nstot,nstot)
c      double precision a3_2o(nstot,nstot,nstot,nstot,nstot,nstot)
c
c
c t2 in a nice representation : be careful here (you may want to reindex
c virtuals
c      double precision t1(nos,nos+1:nos+nvs)
c      double precision t2(nos,nos,nos+1:nos+nvs,nos+1:nos+nvs)
c transformed matrices
ccx      double precision ht(nos+nvs,nos+nvs)
!      double precision FOCKt(nos+nvs,nos+nvs)
ccx      double precision vt(nos+nvs,nos+nvs,nos+nvs,nos+nvs)  
c correlation energies
      double precision eccsd
      double precision xxx,yyy,zzz
c RRST norms and linear equation settings
      double precision xtol
      double precision h1_0ext,h1_1ext,h1_2ext
      double precision h1_0ext_diag,h1_2ext_diag
      double precision h0_0ext,h0_1ext,h0_2ext
      double precision v1_0ext,v1_1ext,v1_2ext
      double precision v2_0ext,v2_1ext,v2_2ext,v2_3ext,v2_4ext
      double precision v2_0ext_diag,v2_2ext_diag,v2_4ext_diag
      double precision oene(nstot)
      integer l_amp1, k_amp1
      integer l_amp2, k_amp2
      integer l_matrix,k_matrix
      integer l_free,k_free
      integer amp1n(nstot,nstot)
      integer amp2n(nstot,nstot,nstot,nstot)
      integer ind1,ind2
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
        iext(i)=0
      do j=1,nos+nvs
        h(i,j)=0.0d0
        h0(i,j)=0.0d0
        v1(i,j)=0.0d0
        a1_1o(i,j)=0.0d0
        a1_2o(i,j)=0.0d0
        iext1(i,j)=0
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
        a2_1o(i,j,k,l)=0.0d0
        a2_2o(i,j,k,l)=0.0d0 
        iext2(i,j,k,l)=0
      enddo
      enddo
      enddo
      enddo
c
      call mapping_v2_m(rtdb,d_v2,k_v2_offset,v,nos,nvs)
c
c
c
c  you have everything here for RRST: v & h are ready
c
c  ATTENTION:   Works only for the RHF case
c
      vnext=nvas-nactv !number of virtual external orbitals
c
      do i=1,nvas-nactv
       iext(noas+nobs+nactv+i)=1
       iext(noas+nobs+nvas+nactv+i)=1
      enddo
c
c
c *** debug ***
      if(nodezero) then
        write(6,*)'xtol:  ',xtol
        write(6,*)'nstot: ',nstot
        write(6,*)'nactv: ',nactv
        write(6,*)'vnext: ',vnext
        write(6,*)'------- IEXT -----'
        do i=1,nstot
         write(6,200) iext(i),i
        enddo
        write(6,*)'------- IEXT -----'
        call util_flush(6)
      endif !nodezero
c *************
c
c
      do i=1,nstot
      do j=1,nstot
       iext1(i,j)=iext(i)+iext(j)
       amp1n(i,j)=0
      enddo
      enddo
c
      do i=1,nstot
      do j=1,nstot
      do k=1,nstot
      do l=1,nstot
       iext2(i,j,k,l)=iext(i)+iext(j)+iext(k)+iext(l)
       amp2n(i,j,k,l)=0
      enddo
      enddo
      enddo
      enddo
c
c *** debug ***
c      if(nodezero) then
c        write(6,*)'------- iext1 -----'
c        do i=1,nstot
c        do j=1,nstot
c          write(6,202) i,j,iext1(i,j)
c        enddo
c        enddo
c        write(6,*)'------- iext1 -----'
c        call util_flush(6)
 202    format(3i4)
c      endif !nodezero
c *************
c
c
c
c ispin
c
      do i=1,noas
       ispin(i)=1
      enddo
      do i=noas+1,noas+nobs
       ispin(i)=2
      enddo
      do i=nos+1,nos+nvas
       ispin(i)=1
      enddo
      do i=nos+nvas+1,nstot
       ispin(i)=2
      enddo
c *** debug ***
c      if(nodezero) then
c        write(6,*)'nstot: ',nstot
c        write(6,*)'------- ISPIN -----'
c        do i=1,nstot
c         write(6,200) ispin(i),i
c        enddo
c        write(6,*)'------- ISPIN -----'
c        call util_flush(6)
 200    format(2i7)
c      endif !nodezero
c *************
c
c
c
c constructing h0 & v1
c
      do i=1,nstot
       oene(i)=0.0d0
      enddo
c
      do i=1,nstot
      do j=1,nstot
       h0(i,j)=h(i,j)
      enddo
      enddo
c
      do i=1,nstot
      do j=1,nstot
        do m=1,nos
         h0(i,j)=h0(i,j)+v(i,m,j,m)
         v1(i,j)=v1(i,j)-v(i,m,j,m)
        enddo
      enddo
      enddo
c
      do i=1,nstot
       oene(i)=h0(i,i)
      enddo
c
c *** debug ***
      if(nodezero) then
        write(6,*)'------- h matrix -----'
        do i=1,nstot
        do j=1,nstot
         xxx=dabs(h(i,j))
         if(xxx.gt.xtol) then
          write(6,201) i,j,h(i,j)
         endif
        enddo
        enddo
      endif
c        write(6,*)'------- h matrix -----'
c        call util_flush(6)
c      endif !nodezero
c      if(nodezero) then 
c       write(6,*)'---- oene vector ------'
c       do i=1,nstot
c        write(6,230) i,oene(i)
c       enddo
c       call util_flush(6)
c      endif
c 230  format(2x,i5,f16.6)
c      if(nodezero) then
c        write(6,*)'------- h0 matrix -----'
c        do i=1,nstot
c        do j=1,nstot
c         xxx=dabs(h0(i,j))
c         if(xxx.gt.1.0d-9) then
c          write(6,201) i,j,h0(i,j)
c         endif
c        enddo
c        enddo
c        write(6,*)'------- h0 matrix -----'
c        call util_flush(6)
 201    format(2i4,3x,f14.6)
c      endif !nodezero
c *************
c
c h0 &  v1 done here
c
c  NORMS 
c
      h0_0ext=0.0d0
      h0_1ext=0.0d0
      h0_2ext=0.0d0
c
      h1_0ext=0.0d0
      h1_1ext=0.0d0
      h1_2ext=0.0d0
c
      h1_0ext_diag=0.0d0
      h1_2ext_diag=0.0d0
c
      v1_0ext=0.0d0
      v1_1ext=0.0d0
      v1_2ext=0.0d0
c
      v2_0ext=0.0d0
      v2_1ext=0.0d0
      v2_2ext=0.0d0
      v2_3ext=0.0d0
      v2_4ext=0.0d0
c
      v2_0ext_diag=0.0d0
      v2_2ext_diag=0.0d0
      v2_4ext_diag=0.0d0
c
      do i=1,nstot
      do j=1,nstot
        xxx=h0(i,j)
        yyy=v1(i,j)
        zzz=h(i,j)
        iaux=iext(i)+iext(j)
        if(iaux.eq.0) then
         h0_0ext=h0_0ext+xxx*xxx
         v1_0ext=v1_0ext+yyy*yyy
         h1_0ext=h1_0ext+zzz*zzz
         if(i.eq.j) then 
          h1_0ext_diag=h1_0ext_diag+zzz*zzz
         endif 
        else if(iaux.eq.1) then
         h0_1ext=h0_1ext+xxx*xxx
         v1_1ext=v1_1ext+yyy*yyy
         h1_1ext=h1_1ext+zzz*zzz
        else
         h0_2ext=h0_2ext+xxx*xxx
         v1_2ext=v1_2ext+yyy*yyy
         h1_2ext=h1_2ext+zzz*zzz
         if(i.eq.j) then
          h1_2ext_diag=h1_2ext_diag+zzz*zzz
         endif
        endif
      enddo
      enddo
c
      do i=1,nstot
      do j=i+1,nstot
      do k=1,nstot
      do l=k+1,nstot
       iaux=iext(i)+iext(j)+iext(k)+iext(l)
       xxx=v(i,j,k,l)
       if(iaux.eq.0) then
          v2_0ext=v2_0ext+xxx*xxx
           if((i.eq.k).and.(j.eq.l)) then
            v2_0ext_diag=v2_0ext_diag+xxx*xxx
           endif
       else if(iaux.eq.1) then 
          v2_1ext=v2_1ext+xxx*xxx
       else if(iaux.eq.2) then
          v2_2ext=v2_2ext+xxx*xxx
           if((i.eq.k).and.(j.eq.l)) then
            v2_2ext_diag=v2_2ext_diag+xxx*xxx
           endif
       else if(iaux.eq.3) then
          v2_3ext=v2_3ext+xxx*xxx
       else 
          v2_4ext=v2_4ext+xxx*xxx
           if((i.eq.k).and.(j.eq.l)) then
            v2_4ext_diag=v2_4ext_diag+xxx*xxx
           endif
       endif
      enddo
      enddo
      enddo
      enddo
c
      h1_0ext=dsqrt(h1_0ext)
      h1_1ext=dsqrt(h1_1ext)
      h1_2ext=dsqrt(h1_2ext)
c
      h1_0ext_diag=dsqrt(h1_0ext_diag)
      h1_2ext_diag=dsqrt(h1_2ext_diag)

c
      h0_0ext=dsqrt(h0_0ext)
      h0_1ext=dsqrt(h0_1ext)
      h0_2ext=dsqrt(h0_2ext)
c
      v1_0ext=dsqrt(v1_0ext)
      v1_1ext=dsqrt(v1_1ext)
      v1_2ext=dsqrt(v1_2ext)
c
      v2_0ext=dsqrt(v2_0ext)
      v2_1ext=dsqrt(v2_1ext)
      v2_2ext=dsqrt(v2_2ext)
      v2_3ext=dsqrt(v2_3ext)
      v2_4ext=dsqrt(v2_4ext)
c
      v2_0ext_diag=dsqrt(v2_0ext_diag)
      v2_2ext_diag=dsqrt(v2_2ext_diag)
      v2_4ext_diag=dsqrt(v2_4ext_diag)
c
      if(nodezero) then 
c
         write(6,97) vnext,h1_0ext
         write(6,98) vnext,h1_1ext
         write(6,99) vnext,h1_2ext
c
         write(6,100) vnext,h0_0ext
         write(6,101) vnext,h0_1ext
         write(6,102) vnext,h0_2ext
c
         write(6,103) vnext,v1_0ext
         write(6,104) vnext,v1_1ext
         write(6,105) vnext,v1_2ext
c
         write(6,106) vnext,v2_0ext
         write(6,107) vnext,v2_1ext
         write(6,108) vnext,v2_2ext
         write(6,109) vnext,v2_3ext
         write(6,110) vnext,v2_4ext
c
         write(6,111) vnext,h1_0ext_diag
         write(6,112) vnext,h1_2ext_diag
         write(6,113) vnext,v2_0ext_diag
         write(6,114) vnext,v2_2ext_diag
         write(6,115) vnext,v2_4ext_diag
      endif
c
c
c
c
c
c
      isize1 = 0
      isize2 = 0
      isize3 = 0
      isize4 = 0
c
      isize1_e1 = 0
      isize1_e2 = 0
c
      isize2_e1 = 0
      isize2_e2 = 0
      isize2_e3 = 0
      isize2_e4 = 0
c
c number of irreducible A1 and A2
c
c selection of A1
c xxx check won 
      do i=1,nstot
      do j=1,nstot
cf        xxx=oene(i)-oene(j)
        zzz=dabs(h(i,j))
        iaux=iext(i)+iext(j)
        if(i.gt.j) then 
        if((iaux.gt.0).and.(zzz.gt.xtol)) then
cf        if(dabs(xxx).gt.xtol) then
         isize1=isize1+1
         amp1n(i,j)=isize1
         amp1n(j,i)=-isize1
         if(iaux.eq.1) isize1_e1=isize1_e1+1
         if(iaux.eq.2) isize1_e2=isize1_e2+1
cf        endif
        endif
        endif
      enddo
      enddo
c selection of A2
      do i=1,nstot
      do j=i+1,nstot
      do k=1,nstot
      do l=k+1,nstot
       iaux=iext(i)+iext(j)+iext(k)+iext(l)
cf       xxx=oene(i)+oene(j)-oene(k)-oene(l)
       zzz=dabs(v(i,j,k,l))
       ind1=i+((j-1)*(j-2))/2
       ind2=k+((l-1)*(l-2))/2
       if(ind2.lt.ind1) then
ccc error here !!!!!
ccc should be i.ne or j.ne.l (equivalent to ind2<ind1 
ccc restricting domain (domain restriction)i
cf       if((i.ne.k).and.(j.ne.l)) then 
       if((iaux.gt.0).and.(zzz.gt.xtol)) then 
cf       if(dabs(xxx).gt.xtol) then
        isize2=isize2+1
        amp2n(i,j,k,l)=isize2+isize1
        amp2n(i,j,l,k)=-(isize2+isize1)
        amp2n(j,i,k,l)=-(isize2+isize1)
        amp2n(j,i,l,k)=(isize2+isize1)
c
        amp2n(k,l,i,j)=-(isize2+isize1)
        amp2n(k,l,j,i)=(isize2+isize1)
        amp2n(l,k,i,j)=(isize2+isize1)
        amp2n(l,k,j,i)=-(isize2+isize1)
c flipped terms - automaticallt included in the same loop 
c (with different numbers) - IMPOSE ANTISYMMETRY !!!!
c fix this problem
        if(iaux.eq.1) isize2_e1=isize2_e1+1
        if(iaux.eq.2) isize2_e2=isize2_e2+1
        if(iaux.eq.3) isize2_e3=isize2_e3+1
        if(iaux.eq.4) isize2_e4=isize2_e4+1
       endif
cf       endif
       endif
cf       endif
      enddo
      enddo
      enddo
      enddo
c
c
      if(nodezero) then 
       write(6,116) isize1
       write(6,118) isize1_e1
       write(6,119) isize1_e2
       write(6,*)'     '
       write(6,117) isize2
       write(6,120) isize2_e1
       write(6,121) isize2_e2
       write(6,122) isize2_e3
       write(6,123) isize2_e4
      endif
c
c
c
c
c  LINEAR EQUATIONS SD
c
      lin_size = isize1+isize2 
      if(nodezero) then
        write(6,131) lin_size
        call util_flush
      endif
c
      if (.not.MA_PUSH_GET(mt_int,2*isize1,'rrst1',l_amp1,
     1 k_amp1)) CALL ERRQUIT('rrst1',0,MA_ERR)
c
      if (.not.MA_PUSH_GET(mt_int,4*isize2,'rrst1a',l_amp2,
     1 k_amp2)) CALL ERRQUIT('rrst1a',0,MA_ERR)
c
      if (.not.MA_PUSH_GET(mt_dbl,lin_size*lin_size,'rrst2',
     1      l_matrix,k_matrix))
     1      call errquit('rrst2',1,MA_ERR)
c
      if (.not.MA_PUSH_GET(mt_dbl,lin_size,'rrst2a',
     1      l_free,k_free))
     1      call errquit('rrst2a',1,MA_ERR)
c 
c initialization
c
       do i=1,2*isize1
        int_mb(k_amp1+i-1) = 0
       enddo
       do i=1,4*isize2
        int_mb(k_amp2+i-1) = 0
       enddo
       do i=1,lin_size
        dbl_mb(k_free+i-1) = 0.0d0
       enddo
       do i=1,lin_size*lin_size
         dbl_mb(k_matrix+i-1) = 0.0d0
       enddo
c
c
c
c *** debug ***
       if(nodezero) then
        write(6,*)'before grad. desc.'
        call util_flush(6)
       endif
c *************
c
       call grad_desc(rtdb,
     &      int_mb(k_amp1),int_mb(k_amp2),isize1,isize2,
     &      lin_size,iext,dbl_mb(k_matrix),dbl_mb(k_free),h,v,
     &      amp1n,amp2n,nstot,xtol,h0,v1,oene,
     &      nos,nvs,noas,nobs,nvas,nvbs)
c
c
c *** debug ***
       if(nodezero) then
        write(6,*)'after grad. desc.'
        call util_flush(6)
       endif
c *************
c
c
      if(.not.MA_POP_STACK(l_free))
     &      call errquit('rrst3a',4,MA_ERR)
c
      if(.not.MA_POP_STACK(l_matrix))
     &      call errquit('rrst3',4,MA_ERR)
c
        if (.not.MA_POP_STACK(l_amp2))
     1    call errquit("rrst4a",0,MA_ERR)
c
        if (.not.MA_POP_STACK(l_amp1))
     1    call errquit("rrst4",0,MA_ERR)
c
c
c
  97  format('# of external orbitals: ',i4,2x,'h1_0ext norm: ',f10.5)
  98  format('# of external orbitals: ',i4,2x,'h1_1ext norm: ',f10.5)
  99  format('# of external orbitals: ',i4,2x,'h1_2ext norm: ',f10.5)
c
 100  format('# of external orbitals: ',i4,2x,'h0_0ext norm: ',f10.5)
 101  format('# of external orbitals: ',i4,2x,'h0_1ext norm: ',f10.5)
 102  format('# of external orbitals: ',i4,2x,'h0_2ext norm: ',f10.5)
c
 103  format('# of external orbitals: ',i4,2x,'v1_0ext norm: ',f10.5)
 104  format('# of external orbitals: ',i4,2x,'v1_1ext norm: ',f10.5)
 105  format('# of external orbitals: ',i4,2x,'v1_2ext norm: ',f10.5)
c
 106  format('# of external orbitals: ',i4,2x,'v2_0ext norm: ',f10.5)
 107  format('# of external orbitals: ',i4,2x,'v2_1ext norm: ',f10.5)
 108  format('# of external orbitals: ',i4,2x,'v2_2ext norm: ',f10.5)
 109  format('# of external orbitals: ',i4,2x,'v2_3ext norm: ',f10.5)
 110  format('# of external orbitals: ',i4,2x,'v2_4ext norm: ',f10.5)
c
 111  format('# of ext. orb.: ',i4,2x,'h1_0ext_diag norm: ',f10.5)
 112  format('# of ext. orb.: ',i4,2x,'h1_2ext_diag norm: ',f10.5)
 113  format('# of ext. orb.: ',i4,2x,'v2_0ext_diag norm: ',f10.5)
 114  format('# of ext. orb.: ',i4,2x,'v2_2ext_diag norm: ',f10.5)
 115  format('# of ext. orb.: ',i4,2x,'v2_4ext_diag norm: ',f10.5)
c
 116  format('# of A1 amplitudes: ',2x,i5)
 117  format('# of A2 amplitudes: ',2x,i5)
c
 118  format('# of A1(1ext) amplitudes: ',2x,i5)
 119  format('# of A1(2ext) amplitudes: ',2x,i5)
c
 120  format('# of A2(1ext) amplitudes: ',2x,i5)
 121  format('# of A2(2ext) amplitudes: ',2x,i5)
 122  format('# of A2(3ext) amplitudes: ',2x,i5)
 123  format('# of A2(4ext) amplitudes: ',2x,i5)
c
 130  format(f14.7)
 131  format('size of the linear problem lin_size: ',1x,i7)
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
       subroutine grad_desc(rtdb,amp1_list,amp2_list,isize1,isize2,
     &      lin_size,iext,matrix,free,h,v,
     &      amp1n,amp2n,nstot,xtol,h0,v1,oene,
     &      nos,nvs,noas,nobs,nvas,nvbs)
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
c
c auxiliary integers doubles
c
        integer i,j,k,l,m,n
c
        integer rtdb
        integer isize1,isize2,lin_size
        integer nstot
        double precision h(nstot,nstot)
        double precision h0(nstot,nstot)
        double precision v1(nstot,nstot)
        double precision oene(nstot)
        double precision v(nstot,nstot,nstot,nstot)
        double precision matrix(lin_size,lin_size)
        double precision free(lin_size)
        double precision xxx,yyy,zzz
        double precision xtol
        integer iext(nstot)
        integer amp1n(nstot,nstot)
        integer amp2n(nstot,nstot,nstot,nstot)
        integer amp1_list(isize1,2)
        integer amp2_list(isize2,4)
        integer ind1,ind2
        integer nos,nvs,noas,nobs,nvas,nvbs
c gauge transformation
        double precision c(nstot,nstot)
c gradient algorithms
        double precision cnew(nstot,nstot)
        double precision cold(nstot,nstot)
        double precision cwork(nstot,nstot)
        double precision grad(nstot,nstot)
        double precision cmax
        double precision penalty 
c penalty is a part of sol_norm evaluaton
c it is also printed independently
c
        double precision sol_norm,sol_norm_old,sol_norm_new
        double precision xgamma
        double precision delta_h
c
        integer maxiter_gd
c
        logical nodezero
c
        nodezero = (ga_nodeid().eq.0)
c initial settings
        xgamma=0.50d0
        maxiter_gd=100000000
        delta_h=0.000010d0
        sol_norm=0.d0
        sol_norm_old=0.0d0
        sol_norm_new=0.0d0
        if(nodezero) then
         write(6,*)'xgamma = ',xgamma
         write(6,*)'maxiter_gd = ',maxiter_gd
         write(6,*)'delta_h = ',delta_h
         call util_flush(6)
        endif 
c initial guess
         do m=1,nstot
         do n=1,nstot
          cold(m,n)=0.0d0
         enddo
         enddo
c
         do m=1,nstot
         do n=m+1,nstot
          xxx=dble(m+n)/(2.0d0)
          yyy=dsin(xxx)
          cold(m,n)=yyy
          cold(n,m)=-yyy
         enddo
         enddo
c
        cmax=0.0d0
        do m=1,nstot
        do n=m+1,nstot
         if(dabs(cold(m,n)).gt.cmax) cmax=dabs(cold(m,n))
        enddo
        enddo
c
        call  rrst_1(rtdb,amp1_list,amp2_list,isize1,isize2,
     &      lin_size,iext,matrix,free,h,v,
     &      amp1n,amp2n,nstot,xtol,h0,v1,oene,
     &    nos,nvs,noas,nobs,nvas,nvbs,cold,sol_norm_old,penalty)
c
        if(nodezero) then
         write(6,*)'Initial guess sol_norm: ',sol_norm_old
         write(6,*)'Cmax : ',cmax
        endif
c
        do i=1,maxiter_gd  !main iteration loop 
c
c gradient calculations
c
         do m=1,nstot
         do n=1,nstot
          grad(m,n)=0.0d0
         enddo
         enddo
c
        xxx=0.0d0
        call  rrst_1(rtdb,amp1_list,amp2_list,isize1,isize2,
     &      lin_size,iext,matrix,free,h,v,
     &      amp1n,amp2n,nstot,xtol,h0,v1,oene,
     &    nos,nvs,noas,nobs,nvas,nvbs,cold,xxx,penalty)
c
         do m=1,nstot
         do n=m+1,nstot
c setting cwork for m,n entry
            do k=1,nstot
            do l=1,nstot
             cwork(k,l)=cold(k,l)
            enddo !l
            enddo !k
            cwork(m,n)=cwork(m,n)+delta_h
            cwork(n,m)=cwork(n,m)-delta_h
c
        yyy=0.0d0
        call  rrst_1(rtdb,amp1_list,amp2_list,isize1,isize2,
     &      lin_size,iext,matrix,free,h,v,
     &      amp1n,amp2n,nstot,xtol,h0,v1,oene,
     &    nos,nvs,noas,nobs,nvas,nvbs,cwork,yyy,penalty)
c
        grad(m,n)=(yyy-xxx)/(delta_h)
        grad(n,m)=(-1.0d0)*grad(m,n)
c
         enddo !n>m
         enddo !m
c
c cnew update
c
         do m=1,nstot
         do n=1,nstot
          cnew(m,n)=0.0d0
         enddo
         enddo
c
         do m=1,nstot
         do n=m+1,nstot
          cnew(m,n)=cold(m,n)-xgamma*grad(m,n)
          cnew(n,m)=(-1.0d0)*cnew(m,n)
         enddo
         enddo
        sol_norm_new=0.0d0
        call  rrst_1(rtdb,amp1_list,amp2_list,isize1,isize2,
     &      lin_size,iext,matrix,free,h,v,
     &      amp1n,amp2n,nstot,xtol,h0,v1,oene,
     &    nos,nvs,noas,nobs,nvas,nvbs,cnew,sol_norm_new,penalty)
c
c checks here
        cmax=0.0d0
        do m=1,nstot
        do n=m+1,nstot
         if(dabs(cnew(m,n)).gt.cmax) cmax=dabs(cnew(m,n))
        enddo
        enddo
c
        yyy=dabs(sol_norm_new-sol_norm_old)
        if(nodezero) then
         if(yyy.lt.1.0d-6) then
          write(6,100) 
          write(6,101) i,yyy,sol_norm_new
          write(6,104) penalty
          write(6,*)' Converged min. norm C matrix '
           do m=1,nstot
           do n=m+1,nstot
            write(6,102) cnew(m,n),m,n
           enddo
           enddo
           call util_flush(6)
           go to 10
         else 
          if(mod(i,10000).eq.0) then
           write(6,101) i,yyy,sol_norm_new
           write(6,103) cmax
          endif
          call util_flush(6)
         endif
        endif !nodezero
c new -> old
        sol_norm_old=sol_norm_new
c
        do m=1,nstot
        do n=1,nstot
         cold(m,n)=cnew(m,n) 
        enddo
        enddo       
c
        enddo !main iteration loop
c
 10     continue
c
 100    format('WOW !!! Convergence !!! WOW')
 101    format('Iter: ',i6,2x,'diff: ',f15.7,2x,'sol.norm: ',f15.7)
 102    format(f15.7,3x,2(i5))
 103    format('Max C amp. ',2x,f15.4)
 104    format('value of penalty function: ',2x,f15.7)
c
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
c
c
c
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
       subroutine rrst_1(rtdb,amp1_list,amp2_list,isize1,isize2,
     &      lin_size,iext,matrix,free,h,v,
     &      amp1n,amp2n,nstot,xtol,h0,v1,oene,
     &    nos,nvs,noas,nobs,nvas,nvbs,c,sol_norm,penalty)
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
        integer isize1,isize2,lin_size
        integer nstot
        double precision h(nstot,nstot)
        double precision h0(nstot,nstot)
        double precision v1(nstot,nstot)
        double precision oene(nstot)
        double precision v(nstot,nstot,nstot,nstot)
        double precision matrix(lin_size,lin_size)
        double precision mp1s(isize1)
        double precision mp1d(isize2)
        double precision mp2s(isize1)
        double precision mp2d(isize2)
        double precision des(isize1),ded(isize2)
        double precision aux1(isize1),aux2(isize2)
        double precision free(lin_size)
        double precision xxx,yyy,zzz
        double precision xtol
        integer iext(nstot)
        integer amp1n(nstot,nstot)
        integer amp2n(nstot,nstot,nstot,nstot)
        integer amp1_list(isize1,2)
        integer amp2_list(isize2,4)
        integer ind1,ind2
        integer nos,nvs,noas,nobs,nvas,nvbs
c gauge transformation
        double precision c(nstot,nstot)
        double precision sol_norm
        double precision penalty
c linear solver
        integer iwork(lin_size)
        integer info
c NR solver
        integer iter
        double precision xvec(lin_size)
        double precision deltax(lin_size)
        double precision fun(lin_size)
        double precision jac(lin_size,lin_size)
c Gauge matrix C - parameters (breaking spin symmetry)
        double precision alphax,betax
        double precision x1max,x2max,xnorm1,xnorm2,xnorm
c
        double precision  xsg
        logical nodezero
c
c
        nodezero = (ga_nodeid().eq.0)
c
c
c        if(nodezero) then 
c         write(6,*)'nos: ',nos
c         write(6,*)'nvs: ',nvs
c         write(6,*)'noas: ',noas
c         write(6,*)'nobs: ',nobs
c         write(6,*)'nvas: ',nvas
c         write(6,*)'nvbs" ',nvbs
c         call util_flush(6)
c        endif
c
c
c       if(nodezero) then
c        write(6,*)'rrst_1 xtol: ',xtol
c        call util_flush(6)
c       endif
c
c selection of A1 (irreducible) --- MP1
c
      m=0
      do i=1,nstot
      do j=1,nstot
cf      xxx=oene(i)-oene(j)
        zzz=dabs(h(i,j))
        iaux=iext(i)+iext(j)
        if(i.gt.j) then
        if((iaux.gt.0).and.(zzz.gt.xtol)) then
cf        if(dabs(xxx).gt.xtol) then
         m=m+1
         amp1_list(m,1)=i
         amp1_list(m,2)=j
        endif
        endif
cf        endif
      enddo
      enddo
c
      if(m.ne.isize1) then
         write(6,*)'isize1= ',isize1
         write(6,*)'m = ',m
         call errquit("m-isize1",0,MA_ERR)
      endif
c
c
c selection of A2 (irreducible)  --- MP1
c
      m=0
c
      do i=1,nstot
      do j=i+1,nstot
      do k=1,nstot
      do l=k+1,nstot
       iaux=iext(i)+iext(j)+iext(k)+iext(l)
cf       xxx=oene(i)+oene(j)-oene(k)-oene(l)
       zzz=dabs(v(i,j,k,l))
       ind1=i+((j-1)*(j-2))/2
       ind2=k+((l-1)*(l-2))/2
       if(ind2.lt.ind1) then
c error domain restriction
c condition repetition
cf       if((i.ne.k).and.(j.ne.l)) then
       if((iaux.gt.0).and.(zzz.gt.xtol)) then
cf       if(dabs(xxx).gt.xtol) then
        m=m+1
        amp2_list(m,1)=i
        amp2_list(m,2)=j
        amp2_list(m,3)=k
        amp2_list(m,4)=l
cf       endif
       endif
cf       endif
       endif
      enddo
      enddo
      enddo
      enddo
c
      if(m.ne.isize2) then
         write(6,*)'isize2= ',isize2
         write(6,*)'m = ',m
         call errquit("m-isize2",0,MA_ERR)
      endif
c
 203  format('AMP1[1st-order]^{',i3,'}_{',i3,'}: ',f12.6)
 204  format('AMP2[1st-order]^{',2i3,'}_{',2i3,'}: ',f12.6)
c
c
c
c
c
 205  format(f14.6,2x,i4)
c
c
 201    format(2i4,3x,f14.6)
c
c
c
c
c
c
c
c 
c    ATTENTION: linear matrix is assumed to be here!!!
c    Be carefull! ADD removed part of the code!!
c    free vectors seems to be OK, but check it.
c
c
c
c
c
c linear equations reinstated (make sure these are linear equations for
c A amplitudes!!!! (at the end matrix = - matrix ) 
c
c
c
       do i=1,lin_size
       do j=1,lin_size
         matrix(i,j)=0.0d0
       enddo
       enddo
c
c
c
c ----------------------------------------------------------------
c
c
c singles equations (no contributions from doubles)
      do i=1,isize1
        p=amp1_list(i,1)
        q=amp1_list(i,2)
c S1
        do s=1,nstot
         ipl1=abs(amp1n(s,q))
         xxx=xsg(amp1n(s,q))
         if(ipl1.ne.0) then
          matrix(i,ipl1)=matrix(i,ipl1)+h(p,s)*xxx
         endif
        enddo
c S2
        do s=1,nstot
         ipl1=abs(amp1n(p,s))
         xxx=xsg(amp1n(p,s))
         if(ipl1.ne.0) then
          matrix(i,ipl1)=matrix(i,ipl1)-h(s,q)*xxx
         endif
        enddo
      enddo
c
c doubles equations (contributions from one- and two-body amplit.)
      do i=1,isize2
       p=amp2_list(i,1)
       q=amp2_list(i,2)
       r=amp2_list(i,3)
       s=amp2_list(i,4)
c D1
       do u=1,nstot
        ipl1=abs(amp1n(u,r))
        xxx=xsg(amp1n(u,r))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+v(p,q,u,s)*xxx
        endif
       enddo
c D2
       do u=1,nstot
        ipl1=abs(amp1n(u,s))
        xxx=xsg(amp1n(u,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-v(p,q,u,r)*xxx
        endif
       enddo
c D3
       do t=1,nstot
        ipl1=abs(amp1n(p,t))
        xxx=xsg(amp1n(p,t))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-v(t,q,r,s)*xxx
        endif
       enddo
c D4
       do t=1,nstot
        ipl1=abs(amp1n(q,t))
        xxx=xsg(amp1n(q,t))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+v(t,p,r,s)*xxx
        endif
       enddo
c D5
       do t=1,nstot
        ipl1=abs(amp2n(t,q,r,s))
        xxx=xsg(amp2n(t,q,r,s))
        if(ipl1.ne.0)then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+h(p,t)*xxx
        endif
       enddo
c D6
       do t=1,nstot
        ipl1=abs(amp2n(t,p,r,s))
        xxx=xsg(amp2n(t,p,r,s))
        if(ipl1.ne.0)then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-h(q,t)*xxx
        endif
       enddo
c D7
       do t=1,nstot
        ipl1=abs(amp2n(p,q,t,s))
        xxx=xsg(amp2n(p,q,t,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-h(t,r)*xxx
        endif
       enddo
c D8
       do t=1,nstot
        ipl1=abs(amp2n(p,q,t,r))
        xxx=xsg(amp2n(p,q,t,r))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+h(t,s)*xxx
        endif
       enddo
c D9
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(t,u,r,s))
        xxx=xsg(amp2n(t,u,r,s))
        if(ipl1.ne.0) then
          matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+
     &    (0.50d0)*v(p,q,t,u)*xxx
        endif
       enddo
       enddo
c D10
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(p,q,t,u))
        xxx=xsg(amp2n(p,q,t,u))
        if(ipl1.ne.0) then
          matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-
     &    (0.50d0)*v(t,u,r,s)*xxx
        endif
       enddo
       enddo
c
      enddo !isize2
c
c sign inversion (A <--> B) Gauge Transformation
c
      do i=1,isize1+isize2
      do j=1,isize1+isize2
        matrix(i,j)=(-1.0d0)*matrix(i,j)
      enddo
      enddo
c
c
c
c
c  end of linear B equations
c ----------------------------------------------------------------
c
c      
c
c    GAUGE TRANSFORMATION BASED FROMULATION
c 
c  
c
c
c
c    zeroing  
c
        do i=1,lin_size
         free(i)=0.0d0
        enddo
c
c  no zeroing C matrix here!!!!!
c
c  Gauge trans. matrix C defined through grad_desc
c
c
c
c Defning free vector (with sign already reversed)
c sign change by [H,C]
      do i=1,isize1
        p=amp1_list(i,1)
        q=amp1_list(i,2)
        free(i)= (-1.0d0)*h(p,q)
        do t=1,nstot
         free(i)=free(i)+h(p,t)*c(t,q)
         free(i)=free(i)-c(p,t)*h(t,q)
c [H,C] sign change OK Gauge Transformation
        enddo 
      enddo
      do i=1,isize2
        p=amp2_list(i,1)
        q=amp2_list(i,2)
        r=amp2_list(i,3)
        s=amp2_list(i,4)
        free(i+isize1)=(-1.0d0)*v(p,q,r,s)
c free(i+isize1) below corrected July 25th
        do t=1,nstot
         free(i+isize1)=free(i+isize1)+v(p,q,t,s)*c(t,r)
         free(i+isize1)=free(i+isize1)-v(p,q,t,r)*c(t,s)
         free(i+isize1)=free(i+isize1)-c(p,t)*v(t,q,r,s)
         free(i+isize1)=free(i+isize1)+c(q,t)*v(t,p,r,s)
c [H,C] sign change OK Gauge Transformation
        enddo
      enddo
c Matrix update: singles
      do i=1,isize1 !singles
        p=amp1_list(i,1)
        q=amp1_list(i,2)
c G4
        do t=1,nstot
        do u=1,nstot
         ipl1=abs(amp1n(u,q))
         xxx=xsg(amp1n(u,q))
         if(ipl1.ne.0) then
          matrix(i,ipl1)=matrix(i,ipl1)+(0.50d0)*h(p,t)*c(t,u)*xxx
         endif
        enddo
        enddo
c G8
        do t=1,nstot
        do u=1,nstot
         ipl1=abs(amp1n(u,q))
         xxx=xsg(amp1n(u,q))
         if(ipl1.ne.0) then
          matrix(i,ipl1)=matrix(i,ipl1)-(0.50d0)*h(t,u)*c(p,t)*xxx
         endif
        enddo
        enddo
c G12
        do t=1,nstot
        do u=1,nstot
         ipl1=abs(amp1n(p,t))
         xxx=xsg(amp1n(p,t))
         if(ipl1.ne.0) then
          matrix(i,ipl1)=matrix(i,ipl1)-(0.50d0)*h(t,u)*c(u,q)*xxx
         endif
        enddo
        enddo
c G17
        do t=1,nstot
        do u=1,nstot
         ipl1=abs(amp1n(p,t))
         xxx=xsg(amp1n(p,t))
         if(ipl1.ne.0) then
          matrix(i,ipl1)=matrix(i,ipl1)+(0.50d0)*h(u,q)*c(t,u)*xxx
         endif
        enddo
        enddo
c
      enddo !singles
c Matrix update: doubles
      do i=1,isize2  !doubles
       p=amp2_list(i,1)
       q=amp2_list(i,2)
       r=amp2_list(i,3)
       s=amp2_list(i,4)
c G1a 
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,r))
        xxx=xsg(amp1n(u,r))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(p,q,t,s)*c(t,u)*xxx
        endif
       enddo
       enddo
c G1b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,s))
        xxx=xsg(amp1n(u,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         v(p,q,t,r)*c(t,u)*xxx
        endif
       enddo
       enddo
c G2a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,r))
        xxx=xsg(amp1n(u,r))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         v(p,q,t,u)*c(t,s)*xxx
        endif
       enddo
       enddo
c G2b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,s))
        xxx=xsg(amp1n(u,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(p,q,t,u)*c(t,r)*xxx
        endif
       enddo
       enddo
c G3
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(u,w,r,s))
        xxx=xsg(amp2n(u,w,r,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(p,q,t,w)*c(t,u)*xxx
        endif
       enddo
       enddo
       enddo
c G5a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(u,q,r,s))
        xxx=xsg(amp2n(u,q,r,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         h(p,t)*c(t,u)*xxx
        endif
       enddo
       enddo
c G5b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(u,p,r,s))
        xxx=xsg(amp2n(u,p,r,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         h(q,t)*c(t,u)*xxx
        endif
       enddo
       enddo
c G6a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,r))
        xxx=xsg(amp1n(u,r))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(t,q,u,s)*c(p,t)*xxx
        endif
       enddo
       enddo
c G6b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,s))
        xxx=xsg(amp1n(u,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         v(t,q,u,r)*c(p,t)*xxx
        endif
       enddo
       enddo
c G6c
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,r))
        xxx=xsg(amp1n(u,r))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         v(t,p,u,s)*c(q,t)*xxx
        endif
       enddo
       enddo
c G6d
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,s))
        xxx=xsg(amp1n(u,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(t,p,u,r)*c(q,t)*xxx
        endif
       enddo
       enddo
c G7a
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(u,w,r,s))
        xxx=xsg(amp2n(u,w,r,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.250d0)*
     &                         v(t,q,u,w)*c(p,t)*xxx
        endif
       enddo
       enddo
       enddo
c G7b
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(u,w,r,s))
        xxx=xsg(amp2n(u,w,r,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.250d0)*
     &                         v(t,p,u,w)*c(q,t)*xxx
        endif
       enddo
       enddo
       enddo
c G9a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(u,q,r,s))
        xxx=xsg(amp2n(u,q,r,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         h(t,u)*c(p,t)*xxx
        endif
       enddo
       enddo
c G9b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(u,p,r,s))
        xxx=xsg(amp2n(u,p,r,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         h(t,u)*c(q,t)*xxx
        endif
       enddo
       enddo
c G10a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(p,t))
        xxx=xsg(amp1n(p,t))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         v(t,q,u,s)*c(u,r)*xxx
        endif
       enddo
       enddo
c G10b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(p,t))
        xxx=xsg(amp1n(p,t))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(t,q,u,r)*c(u,s)*xxx
        endif
       enddo
       enddo
c G10c
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(q,t))
        xxx=xsg(amp1n(q,t))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(t,p,u,s)*c(u,r)*xxx
        endif
       enddo
       enddo
c G10d
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(q,t))
        xxx=xsg(amp1n(q,t))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         v(t,p,u,r)*c(u,s)*xxx
        endif
       enddo
       enddo
c G11a
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(p,q,t,u))
        xxx=xsg(amp2n(p,q,t,u))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.250d0)*
     &                         v(t,u,w,s)*c(w,r)*xxx
        endif
       enddo
       enddo
       enddo
c G11b
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(p,q,t,u))
        xxx=xsg(amp2n(p,q,t,u))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.250d0)*
     &                         v(t,u,w,r)*c(w,s)*xxx
        endif
       enddo
       enddo
       enddo
c G13a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(p,q,t,s))
        xxx=xsg(amp2n(p,q,t,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         h(t,u)*c(u,r)*xxx
        endif
       enddo
       enddo
c G13b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(p,q,t,r))
        xxx=xsg(amp2n(p,q,t,r))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         h(t,u)*c(u,s)*xxx
        endif
       enddo
       enddo
c G14a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(p,t))
        xxx=xsg(amp1n(p,t))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(u,q,r,s)*c(t,u)*xxx
        endif
       enddo
       enddo
c G14b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(q,t))
        xxx=xsg(amp1n(q,t))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         v(u,p,r,s)*c(t,u)*xxx
        endif
       enddo
       enddo
c G15a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(p,u))
        xxx=xsg(amp1n(p,u))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         v(t,u,r,s)*c(q,t)*xxx
        endif
       enddo
       enddo
c G15b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(q,u))
        xxx=xsg(amp1n(q,u))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(t,u,r,s)*c(p,t)*xxx
        endif
       enddo
       enddo
c G16
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(p,q,t,w))
        xxx=xsg(amp2n(p,q,t,w))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         v(u,w,r,s)*c(t,u)*xxx
        endif
       enddo
       enddo
       enddo
c G18a
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(p,q,t,s))
        xxx=xsg(amp2n(p,q,t,s))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)+(0.50d0)*
     &                         h(u,r)*c(t,u)*xxx
        endif
       enddo
       enddo
c G18b
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp2n(p,q,t,r))
        xxx=xsg(amp2n(p,q,t,r))
        if(ipl1.ne.0) then
         matrix(i+isize1,ipl1)=matrix(i+isize1,ipl1)-(0.50d0)*
     &                         h(u,s)*c(t,u)*xxx
        endif
       enddo
       enddo
c
      enddo !doubles
c
c
c
c
c       if(nodezero) then
c        write(6,*)'    '
c        write(6,*)'Solving Linearized Gauge Transformation Equations'
c        write(6,*)'    '
c        call util_flush(6)
c       endif
c
c
c Linearized Gauge approximation
c
        call dgesv(lin_size,1,matrix,lin_size,iwork,free,lin_size,info)
        if (info .ne. 0)
     1    call errquit('rrts: lin. solver failed',info,
     2    UNKNOWN_ERR)
c
c        if(nodezero) then ! nodezero
c         do i=1,isize1
c          write(6,10) amp1_list(i,1),amp1_list(i,2),free(i)
c         enddo
c         do i=1,isize2
c          write(6,11) amp2_list(i,1),amp2_list(i,2),amp2_list(i,3),
c     &                amp2_list(i,4),free(i+isize1)
c         enddo
c         call util_flush(6)
c        endif !nodezero
c
c
c Solution Analysis
c
        sol_norm=0.0d0
        penalty=0.0d0
c
        xxx=0.0d0
        do i=1,nstot
        do j=i+1,nstot
         if(dabs(c(i,j)).gt.xxx) xxx=dabs(c(i,j))
        enddo
        enddo
c
        if(xxx.lt.1.0d0) then
         penalty = 0.0d0
        else
c         yyy=xxx*xxx*xxx*xxx*xxx*xxx*xxx*xxx*xxx*xxx
c         zzz=yyy*yyy
         penalty=xxx*xxx*xxx-1.0d0
c         penalty = dexp(zzz-1.0d0)-1.0d0
        endif
c
        do i=1,lin_size
          sol_norm=sol_norm+free(i)*free(i)
        enddo
        sol_norm=dsqrt(sol_norm)
        sol_norm=sol_norm+penalty
c
         x1max=0.0d0
         do i=1,isize1
          xxx=dabs(free(i))
          if(xxx.gt.x1max) x1max=xxx
         enddo
c
         x2max=0.0d0
         do i=1,isize2
          j=i+isize1
          xxx=dabs(free(j))
          if(xxx.gt.x2max) x2max=xxx
         enddo
c
c         if(nodezero) then
c           write(6,*)'-----'
c           write(6,500) x1max
c           write(6,501) x2max
c           write(6,*)'-----'
c          call util_flush(6)
c         endif
c
 498     format(' alphax : ',f16.6)
 499     format(' betax  : ',f16.6)
 500     format(' max. B1 amp. = ',f16.6)
 501     format(' max. B2 amp. = ',f16.6)
 502     format(' norm B1 sol. = ',f16.6)
 503     format(' norm B2 sol. = ',f16.6)
 504     format(' tot. norm of sol. = ',f16.6)
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
c
c
c
c
c
c
  9     format('SD Quadratic approximation: Newton-Raphson solver')
 10     format(2i5,2x,f14.6)
 11     format(4i5,2x,f14.6)
 12     format('NR iterative process has converged')
 13     format('iteration: ',i5,2x,'norm of residual vector: ',f16.10)
 14     format('Converged SD amplitudes')
c
       return 
       end
c
c
c
c
c
c ------------------ Newton-Raphson infrastructure
c
c
c
      subroutine nrquantum(xvec,fun,jac,v,h,amp1n,amp2n,
     &           amp1_list,amp2_list,nstot,isize1,isize2,lin_size)
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
        integer i,j,k,l,m,n,iaux,ii
        integer p,q,r,s,t,u,w,x,y,z
        integer ipl1,ipl2
c
        integer isize1,isize2,lin_size
        integer nstot
        double precision h(nstot,nstot)
        double precision v(nstot,nstot,nstot,nstot)
        double precision matrix(lin_size,lin_size)
        double precision xvec(lin_size)
        double precision fun(lin_size)
        double precision jac(lin_size,lin_size)
        double precision xx,yy,zz
        integer amp1n(nstot,nstot)
        integer amp2n(nstot,nstot,nstot,nstot)
        integer amp1_list(isize1,2)
        integer amp2_list(isize2,4)
c
        double precision  xsg
        logical nodezero
c
        nodezero = (ga_nodeid().eq.0)
c
c zeroing function & jacobian
        do i=1,lin_size
         fun(i) = 0.0d0
         do j=1,lin_size
          jac(i,j) = 0.0d0
         enddo
        enddo
c free terms
      do i=1,isize1
        p=amp1_list(i,1)
        q=amp1_list(i,2)
        fun(i)= fun(i) + h(p,q)
      enddo
c
      do i=1,isize2
        p=amp2_list(i,1)
        q=amp2_list(i,2)
        r=amp2_list(i,3)
        s=amp2_list(i,4)
        fun(i+isize1) = fun(i+isize1) + v(p,q,r,s)
      enddo
c linear part
      do i=1,lin_size
      do j=1,lin_size
       fun(i) = fun(i) + matrix(i,j)*xvec(j)
       jac(i,j) = jac(i,j) + matrix(i,j)
      enddo
      enddo
c *** debug ***
       if(nodezero) then
          write(6,*)'in nrquantum step1'
          call util_flush(6)
       endif
c *************
c quadratic contributions to fun & jac
c Qsingles ----------------------------------------
      do i=1,isize1
        p=amp1_list(i,1)
        q=amp1_list(i,2)
c QS1 (1h)
        do t=1,nstot
        do u=1,nstot
         ipl1=abs(amp1n(t,u))
         xx=xsg(amp1n(t,u))
         ipl2=abs(amp1n(u,q))
         yy=xsg(amp1n(u,q))
         if((ipl1.ne.0).and.(ipl2.ne.0)) then
          fun(i)=fun(i)+(0.50d0)*h(p,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
          jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(p,t)*xvec(ipl2)*xx*yy
          jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(p,t)*xvec(ipl1)*xx*yy
         endif
        enddo
        enddo
c QS2+QS3 (6h/13h) - change in the weight factor
        do t=1,nstot
        do u=1,nstot
         ipl1=abs(amp1n(p,t))
         xx=xsg(amp1n(p,t))
         ipl2=abs(amp1n(u,q))
         yy=xsg(amp1n(u,q))
         if((ipl1.ne.0).and.(ipl2.ne.0)) then
          fun(i)=fun(i)-(1.0d0)*h(t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
          jac(i,ipl1)=jac(i,ipl1)-(1.0d0)*h(t,u)*xvec(ipl2)*xx*yy      
          jac(i,ipl2)=jac(i,ipl2)-(1.00d0)*h(t,u)*xvec(ipl1)*xx*yy
         endif
        enddo
        enddo
c QS4 (20h)
        do t=1,nstot
        do u=1,nstot
         ipl1=abs(amp1n(p,t))
         xx=xsg(amp1n(p,t))
         ipl2=abs(amp1n(t,u))
         yy=xsg(amp1n(t,u))
         if((ipl1.ne.0).and.(ipl2.ne.0)) then
          fun(i)=fun(i)+(0.50d0)*h(u,q)*xvec(ipl1)*xvec(ipl2)*xx*yy
          jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(u,q)*xvec(ipl2)*xx*yy
          jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(u,q)*xvec(ipl1)*xx*yy
         endif
        enddo
        enddo
      enddo !isize1
c *** debug ***
       if(nodezero) then
          write(6,*)'in nrquantum step2'
          call util_flush(6)
       endif
c *************
c Qdoubles  !-------------------------------------------------
      do ii=1,isize2 !ii - local scale
       p=amp2_list(ii,1)
       q=amp2_list(ii,2)
       r=amp2_list(ii,3)
       s=amp2_list(ii,4)
       i=ii+isize1  !i - global scale
c QD1 (2h)
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(t,u))
        xx=xsg(amp1n(t,u))
        ipl2=abs(amp2n(u,q,r,s))
        yy=xsg(amp2n(u,q,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*h(p,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(p,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(p,t)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(t,u))
        xx=xsg(amp1n(t,u))
        ipl2=abs(amp2n(u,p,r,s))
        yy=xsg(amp2n(u,p,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*h(q,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*h(q,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*h(q,t)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c QD2 (3h)
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,r))
        xx=xsg(amp1n(u,r))
        ipl2=abs(amp2n(t,q,u,s))
        yy=xsg(amp2n(t,q,u,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*h(p,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(p,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(p,t)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(u,s))
        xx=xsg(amp1n(u,s))
        ipl2=abs(amp2n(t,q,u,r))
        yy=xsg(amp2n(t,q,u,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*h(p,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*h(p,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*h(p,t)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(u,r))
        xx=xsg(amp1n(u,r))
        ipl2=abs(amp2n(t,q,u,s))
        yy=xsg(amp2n(t,q,u,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*h(q,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*h(q,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*h(q,t)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(u,s))
        xx=xsg(amp1n(u,s))
        ipl2=abs(amp2n(t,p,u,r))
        yy=xsg(amp2n(t,p,u,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*h(q,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(q,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(q,t)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c QD3 (4h)
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(t,q,u,w))
        xx=xsg(amp2n(t,q,u,w))
        ipl2=abs(amp2n(u,w,r,s))
        yy=xsg(amp2n(u,w,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.25d0)*h(p,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.25d0)*h(p,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.25d0)*h(p,t)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp2n(t,p,u,w))
        xx=xsg(amp2n(t,p,u,w))
        ipl2=abs(amp2n(u,w,r,s))
        yy=xsg(amp2n(u,w,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.25d0)*h(q,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.25d0)*h(q,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.25d0)*h(q,t)*xvec(ipl1)*xx*yy
        endif 
       enddo
       enddo
       enddo
c QD4 (7h) + QD9 (15h)
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(p,t))
        xx=xsg(amp1n(p,t))
        ipl2=abs(amp2n(u,q,r,s))
        yy=xsg(amp2n(u,q,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(1.00d0)*h(t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(1.00d0)*h(t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(1.00d0)*h(t,u)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,t))
        xx=xsg(amp1n(q,t))
        ipl2=abs(amp2n(u,p,r,s))
        yy=xsg(amp2n(u,p,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(1.00d0)*h(t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(1.00d0)*h(t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(1.00d0)*h(t,u)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c QD5 (8h) + QD8 (14h)
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,r))
        xx=xsg(amp1n(u,r))
        ipl2=abs(amp2n(p,q,t,s))
        yy=xsg(amp2n(p,q,t,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(1.00d0)*h(t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(1.00d0)*h(t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(1.00d0)*h(t,u)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(u,s))
        xx=xsg(amp1n(u,s))
        ipl2=abs(amp2n(p,q,t,r))
        yy=xsg(amp2n(p,q,t,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(1.00d0)*h(t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(1.00d0)*h(t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(1.00d0)*h(t,u)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c QD6 (9h)
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(u,r))
        xx=xsg(amp1n(u,r))
        ipl2=abs(amp2n(p,q,t,u))
        yy=xsg(amp2n(p,q,t,u))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*h(t,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(t,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(t,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(u,s))
        xx=xsg(amp1n(u,s))
        ipl2=abs(amp2n(p,q,t,u))
        yy=xsg(amp2n(p,q,t,u))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*h(t,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*h(t,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*h(t,r)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c QD7 (10h) + QD11 (17h)
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(p,q,t,w))
        xx=xsg(amp2n(p,q,t,w))
        ipl2=abs(amp2n(u,w,r,s))
        yy=xsg(amp2n(u,w,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(1.00d0)*h(t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(1.00d0)*h(t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(1.00d0)*h(t,u)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c QD10 (16h)
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(p,u))
        xx=xsg(amp1n(p,u))
        ipl2=abs(amp2n(t,u,r,s))
        yy=xsg(amp2n(t,u,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*h(q,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(q,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(q,t)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,u))
        xx=xsg(amp1n(q,u))
        ipl2=abs(amp2n(t,u,r,s))
        yy=xsg(amp2n(t,u,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*h(p,t)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*h(p,t)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*h(p,t)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c QD12 (21h)
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(t,u))
        xx=xsg(amp1n(t,u))
        ipl2=abs(amp2n(p,q,t,s))
        yy=xsg(amp2n(p,q,t,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*h(u,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(u,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(u,r)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(t,u))
        xx=xsg(amp1n(t,u))
        ipl2=abs(amp2n(p,q,t,r))
        yy=xsg(amp2n(p,q,t,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*h(u,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*h(u,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*h(u,s)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c QD13 (22h)
c QD12 (21h)
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(p,t))
        xx=xsg(amp1n(p,t))
        ipl2=abs(amp2n(t,q,u,s))
        yy=xsg(amp2n(t,q,u,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*h(u,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(u,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(u,r)*xvec(ipl1)*xx*yy
        endif
c 
        ipl1=abs(amp1n(p,t))
        xx=xsg(amp1n(p,t))
        ipl2=abs(amp2n(t,q,u,r))
        yy=xsg(amp2n(t,q,u,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*h(u,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*h(u,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*h(u,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,t))
        xx=xsg(amp1n(q,t))
        ipl2=abs(amp2n(t,p,u,s))
        yy=xsg(amp2n(t,p,u,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*h(u,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*h(u,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*h(u,r)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,t))
        xx=xsg(amp1n(q,t))
        ipl2=abs(amp2n(t,p,u,r))
        yy=xsg(amp2n(t,p,u,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*h(u,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*h(u,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*h(u,s)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c Q14 (23h)
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp2n(p,q,t,u))
        xx=xsg(amp2n(p,q,t,u))
        ipl2=abs(amp2n(t,u,w,s))
        yy=xsg(amp2n(t,u,w,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.25d0)*h(w,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.25d0)*h(w,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.25d0)*h(w,r)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp2n(p,q,t,u))
        xx=xsg(amp2n(p,q,t,u))
        ipl2=abs(amp2n(t,u,w,r))
        yy=xsg(amp2n(t,u,w,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.25d0)*h(w,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.25d0)*h(w,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.25d0)*h(w,s)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q15 1v
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(t,s))
        xx=xsg(amp1n(t,s))
        ipl2=abs(amp1n(u,r))
        yy=xsg(amp1n(u,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*v(p,q,t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*v(p,q,t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*v(p,q,t,u)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(t,r))
        xx=xsg(amp1n(t,r))
        ipl2=abs(amp1n(u,s))
        yy=xsg(amp1n(u,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*v(p,q,t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*v(p,q,t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*v(p,q,t,u)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c Q16 2v
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(t,u))
        xx=xsg(amp1n(t,u))
        ipl2=abs(amp1n(u,r))
        yy=xsg(amp1n(u,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*v(p,q,t,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*v(p,q,t,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*v(p,q,t,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(t,u))
        xx=xsg(amp1n(t,u))
        ipl2=abs(amp1n(u,s))
        yy=xsg(amp1n(u,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*v(p,q,t,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*v(p,q,t,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*v(p,q,t,r)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c Q17 3v
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp1n(w,r))
        xx=xsg(amp1n(w,r))
        ipl2=abs(amp2n(t,u,w,s))
        yy=xsg(amp2n(t,u,w,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.25d0)*v(p,q,t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.25d0)*v(p,q,t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.25d0)*v(p,q,t,u)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(w,s))
        xx=xsg(amp1n(w,s))
        ipl2=abs(amp2n(t,u,w,r))
        yy=xsg(amp2n(t,u,w,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.25d0)*v(p,q,t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.25d0)*v(p,q,t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.25d0)*v(p,q,t,u)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q18 4v
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp1n(t,u))
        xx=xsg(amp1n(t,u))
        ipl2=abs(amp2n(u,w,r,s))
        yy=xsg(amp2n(u,w,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*v(p,q,t,w)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*v(p,q,t,w)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*v(p,q,t,w)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q19 5v
       do t=1,nstot
       do u=1,nstot
       do x=1,nstot
       do y=1,nstot
        ipl1=abs(amp2n(t,u,x,y))
        xx=xsg(amp2n(t,u,x,y))
        ipl2=abs(amp2n(x,y,r,s))
        yy=xsg(amp2n(x,y,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.1250d0)*v(p,q,t,u)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.1250d0)*v(p,q,t,u)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.1250d0)*v(p,q,t,u)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
       enddo
c Q20 6v + Q24 10v
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp1n(p,t))
        xx=xsg(amp1n(p,t))
        ipl2=abs(amp1n(u,r))
        yy=xsg(amp1n(u,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(1.00d0)*v(t,q,u,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(1.00d0)*v(t,q,u,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(1.00d0)*v(t,q,u,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(p,t))
        xx=xsg(amp1n(p,t))
        ipl2=abs(amp1n(u,s))
        yy=xsg(amp1n(u,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(1.00d0)*v(t,q,u,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(1.00d0)*v(t,q,u,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(1.00d0)*v(t,q,u,r)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,t))
        xx=xsg(amp1n(q,t))
        ipl2=abs(amp1n(u,r))
        yy=xsg(amp1n(u,r))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(1.00d0)*v(t,p,u,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(1.00d0)*v(t,p,u,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(1.00d0)*v(t,p,u,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,t))
        xx=xsg(amp1n(q,t))
        ipl2=abs(amp1n(u,s))
        yy=xsg(amp1n(u,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(1.00d0)*v(t,p,u,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(1.00d0)*v(t,p,u,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(1.00d0)*v(t,p,u,r)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q21 7v + Q25 11v
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp1n(p,t))
        xx=xsg(amp1n(p,t))
        ipl2=abs(amp2n(u,w,r,s))
        yy=xsg(amp2n(u,w,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*v(t,q,u,w)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*v(t,q,u,w)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*v(t,q,u,w)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,t))
        xx=xsg(amp1n(q,t))
        ipl2=abs(amp2n(u,w,r,s))
        yy=xsg(amp2n(u,w,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*v(t,p,u,w)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*v(t,p,u,w)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*v(t,p,u,w)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q22 8v + Q26 12v
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp1n(w,r))
        xx=xsg(amp1n(w,r))
        ipl2=abs(amp2n(p,q,t,u))
        yy=xsg(amp2n(p,q,t,u))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*v(t,u,w,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*v(t,u,w,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*v(t,u,w,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(w,s))
        xx=xsg(amp1n(w,s))
        ipl2=abs(amp2n(p,q,t,u))
        yy=xsg(amp2n(p,q,t,u))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*v(t,u,w,r)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*v(t,u,w,r)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*v(t,u,w,r)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q23 9v + Q27 13v
       do t=1,nstot
       do u=1,nstot
       do x=1,nstot
       do y=1,nstot
        ipl1=abs(amp2n(p,q,t,u))
        xx=xsg(amp2n(p,q,t,u))
        ipl2=abs(amp2n(x,y,r,s))
        yy=xsg(amp2n(x,y,r,s))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.250d0)*v(t,u,x,y)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.250d0)*v(t,u,x,y)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.250d0)*v(t,u,x,y)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
       enddo
c Q28 14v
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp1n(p,t))
        xx=xsg(amp1n(p,t))
        ipl2=abs(amp1n(t,u))
        yy=xsg(amp1n(t,u))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*v(u,q,r,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*v(u,q,r,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*v(u,q,r,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,t))
        xx=xsg(amp1n(q,t))
        ipl2=abs(amp1n(t,u))
        yy=xsg(amp1n(t,u))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*v(u,p,r,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*v(u,p,r,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*v(u,p,r,s)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q29 15v
       do t=1,nstot
       do u=1,nstot
        ipl1=abs(amp1n(p,u))
        xx=xsg(amp1n(p,u))
        ipl2=abs(amp1n(q,t))
        yy=xsg(amp1n(q,t))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*v(t,u,r,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*v(t,u,r,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*v(t,u,r,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,u))
        xx=xsg(amp1n(q,u))
        ipl2=abs(amp1n(p,t))
        yy=xsg(amp1n(p,t))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.50d0)*v(t,u,r,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.50d0)*v(t,u,r,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.50d0)*v(t,u,r,s)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
c Q30 16v
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp1n(t,w))
        xx=xsg(amp1n(t,w))
        ipl2=abs(amp2n(p,q,t,u))
        yy=xsg(amp2n(p,q,t,u))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.50d0)*v(u,w,r,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.50d0)*v(u,w,r,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.50d0)*v(u,w,r,s)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q31 17v
       do t=1,nstot
       do u=1,nstot
       do w=1,nstot
        ipl1=abs(amp1n(p,t))
        xx=xsg(amp1n(p,t))
        ipl2=abs(amp2n(t,q,u,w))
        yy=xsg(amp2n(t,q,u,w))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.250d0)*v(u,w,r,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.250d0)*v(u,w,r,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.250d0)*v(u,w,r,s)*xvec(ipl1)*xx*yy
        endif
c
        ipl1=abs(amp1n(q,t))
        xx=xsg(amp1n(q,t))
        ipl2=abs(amp2n(t,p,u,w))
        yy=xsg(amp2n(t,p,u,w))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)-(0.250d0)*v(u,w,r,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)-(0.250d0)*v(u,w,r,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)-(0.250d0)*v(u,w,r,s)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
c Q32 18v
       do t=1,nstot
       do u=1,nstot
       do x=1,nstot
       do y=1,nstot
        ipl1=abs(amp2n(p,q,t,u))
        xx=xsg(amp2n(p,q,t,u))
        ipl2=abs(amp2n(t,u,x,y))
        yy=xsg(amp2n(t,u,x,y))
        if((ipl1.ne.0).and.(ipl2.ne.0))  then
        fun(i)=fun(i)+(0.1250d0)*v(x,y,r,s)*xvec(ipl1)*xvec(ipl2)*xx*yy
        jac(i,ipl1)=jac(i,ipl1)+(0.1250d0)*v(x,y,r,s)*xvec(ipl2)*xx*yy
        jac(i,ipl2)=jac(i,ipl2)+(0.1250d0)*v(x,y,r,s)*xvec(ipl1)*xx*yy
        endif
       enddo
       enddo
       enddo
       enddo
c
      enddo !isize2
c
c *** debug ***
       if(nodezero) then
          write(6,*)'in nrquantum step3'
          call util_flush(6)
       endif
c *************
c 
c
      return 
      end

c
c
c
c
c
c ------------------
c
c
c
       double precision function xsg(i)
       implicit none
       integer i
       xsg=0.0d0
       if(i.ge.0) then 
         xsg=1.0d0
       else 
         xsg=-1.0d0
       endif
       return 
       end
