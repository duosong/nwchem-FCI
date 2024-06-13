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
