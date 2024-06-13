c
      subroutine kinetic_hcore(rtdb,spust)
C     $Id: tce_mo2e_zones_4a_disk_ga_N5.F 25705 2014-06-02 18:58:56Z
C     kowalski $
      IMPLICIT NONE
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb
      integer g_hcore
      integer nbf2
      integer k_core,l_core
      integer iha,ihb,ipa,ipb
      integer isigma1
      integer l_movecs_orb,k_movecs_orb
      integer l_aux,k_aux
      integer dim_mol_a,dim_mol_b
      integer l_final,k_final
      integer i,j,k
      integer uhf_a_num, uhf_b_num
      double precision xxx
cc      logical oskel
      logical spust
      logical nodezero
      oskel=.false.
      nodezero=(ga_nodeid().eq.0)
c
      if (.not. ga_create(mt_dbl,nbf,nbf, 'AO ovl',
     &        0, 0, g_hcore))
     &        call errquit('Error creating ga_hcore',0,GA_ERR)
c
      call ga_zero(g_hcore)
c 
      call int_1e_ga(ao_bas_han,ao_bas_han,g_hcore,'kinetic',oskel)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_hcore,'potential',oskel)
c create MA l_core
      nbf2=nbf*nbf
      if (.not.ma_push_get(mt_dbl,nbf2,'h_core',l_core,k_core))
     1    call errquit('tce_core: MA problem',0,MA_ERR)
      call dfill(nbf2, 0.0d0, dbl_mb(k_core), 1)
c h_hcore ==> l_core ! attention - the last parameter (leading
c             dimension of local memory buffer)
        call ga_get(g_hcore,1,nbf,1,nbf,dbl_mb(k_core),nbf)
c 
c alpha orbitals only
c
        iha = nocc(1)-nfc(1)                 !number of alpha holes
        ihb = nocc(ipol)-nfc(ipol)           !number of beta holes
        ipa = nmo(1)-nocc(1)-nfv(1)          !number of alpha particles
        ipb = nmo(ipol)-nocc(ipol)-nfv(ipol) !number of beta particles
c
c
        dim_mol_a=iha+ipa    !block dimension for alpha
        dim_mol_b-ihb+ipb    !block dimension for beta
c
c
c a-a block
c
c
      if (.not.ma_push_get(mt_dbl,nbf*(iha+ipa)
     1  ,"sorted MO coeffs",
     2  l_movecs_orb,k_movecs_orb))
     3  call errquit('tce_mo2e_zone: MA problem 1',0,
     2    BASIS_ERR)
      call dfill(nbf*(iha+ipa),0.0d0, dbl_mb(k_movecs_orb), 1)
      do i=1,iha
      do isigma1=1,nbf
       dbl_mb(k_movecs_orb+(i-1)*nbf+isigma1-1)=
     & dbl_mb(k_movecs_sorted+(i-1)*nbf+isigma1-1)
      enddo
      enddo
      do i=iha+1,iha+ipa
      do isigma1=1,nbf
       dbl_mb(k_movecs_orb+(i-1)*nbf+isigma1-1)=
     & dbl_mb(k_movecs_sorted+(i+ihb-1)*nbf+isigma1-1)
      enddo
      enddo
c
      if (.not.ma_push_get(mt_dbl,nbf*dim_mol_a,'l_aux',l_aux,k_aux))
     1    call errquit('l_aux_core: MA problem',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,dim_mol_a*dim_mol_a,'l_final',
     1    l_final,k_final))
     1    call errquit('tce_core: MA problem',0,MA_ERR)
c   l_aux (nbf,dim_mol_a)= k_core(nbf,nbf) * k_movecs_orb(nbf,dim_mol_a)
       call dgemm('N','N',nbf,dim_mol_a,nbf,1.0d0,dbl_mb(k_core),nbf,
     1            dbl_mb(k_movecs_orb),nbf,0.0d0,dbl_mb(k_aux),nbf)
c   l_final= k_movecs_orb(nbf,dim_mol_a)^T * l_aux(nbf,dim_mol_a)
c   l_final= k_movecs_orb_T(dim_mol_a,nbf) * l_aux(nbf,dim_mol_a)
       call dgemm('T','N',dim_mol_a,dim_mol_a,nbf,1.0d0,
     1            dbl_mb(k_movecs_orb),nbf,dbl_mb(k_aux),
     1            nbf,0.0d0,dbl_mb(k_final),dim_mol_a)
c
       if(nodezero) then !nodezero`
        if(spust) then  !spust
        k=0
        do j=1,dim_mol_a
        do i=1,dim_mol_a
         k=k+1
         xxx=dabs(dbl_mb(k_final+k-1))
         if(xxx.gt.1.0d-10) then
          if(spust) then
c yaml one_electron_integrals
           write(6,100) uhf_a_num(i,iha,ihb,ipa,ipb),
     &                  uhf_a_num(j,iha,ihb,ipa.ipb),dbl_mb(k_final+k-1)
          endif
         end if 
        enddo
        enddo
        call util_flush(6)
        end if   !spust
       end if            !nodezero
c delete k_final
      if (.not.ma_pop_stack(l_final))
     1  call errquit('tcc_hcore_final: MA problem',15,MA_ERR)
c delete k_aux
      if (.not.ma_pop_stack(l_aux))
     1  call errquit('tcc_hcore_aux: MA problem',15,MA_ERR)
c delete k_movecs_orb
      if (.not.ma_pop_stack(l_movecs_orb))
     1  call errquit('tcc_hcore: MA problem',15,MA_ERR)
c
c
c
c b-b blcok
c
c
c
      if (.not.ma_push_get(mt_dbl,nbf*(ihb+ipb)
     1  ,"sorted MO coeffs",
     2  l_movecs_orb,k_movecs_orb))
     3  call errquit('tce_mo2e_zone: MA problem 1',0,
     2    BASIS_ERR)
      call dfill(nbf*(ihb+ipb),0.0d0, dbl_mb(k_movecs_orb), 1)
      do i=1,ihb
      do isigma1=1,nbf
       dbl_mb(k_movecs_orb+(i-1)*nbf+isigma1-1)=
     & dbl_mb(k_movecs_sorted+(i-1+iha)*nbf+isigma1-1)
      enddo
      enddo
      do i=ihb+1,iha+ipb
      do isigma1=1,nbf
       dbl_mb(k_movecs_orb+(i-1)*nbf+isigma1-1)=
     & dbl_mb(k_movecs_sorted+(i+iha+ipa-1)*nbf+isigma1-1)
      enddo
      enddo
c
      if (.not.ma_push_get(mt_dbl,nbf*dim_mol_b,'l_aux',l_aux,k_aux))
     1    call errquit('l_aux_core: MA problem',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,dim_mol_b*dim_mol_b,'l_final',
     1    l_final,k_final))
     1    call errquit('tce_core: MA problem',0,MA_ERR)
c   l_aux (nbf,dim_mol_b)= k_core(nbf,nbf) * k_movecs_orb(nbf,dim_mol_b)
       call dgemm('N','N',nbf,dim_mol_b,nbf,1.0d0,dbl_mb(k_core),nbf,
     1            dbl_mb(k_movecs_orb),nbf,0.0d0,dbl_mb(k_aux),nbf)
c   l_final= k_movecs_orb(nbf,dim_mol_b)^T * l_aux(nbf,dim_mol_b)
c   l_final= k_movecs_orb_T(dim_mol_b,nbf) * l_aux(nbf,dim_mol_b)
       call dgemm('T','N',dim_mol_b,dim_mol_b,nbf,1.0d0,
     1            dbl_mb(k_movecs_orb),nbf,dbl_mb(k_aux),
     1            nbf,0.0d0,dbl_mb(k_final),dim_mol_b)
c
       if(nodezero) then !nodezero`
        if(spust) then  !spust
        k=0
        do j=1,dim_mol_b
        do i=1,dim_mol_b
         k=k+1
         xxx=dabs(dbl_mb(k_final+k-1))
         if(xxx.gt.1.0d-10) then
          if(spust) then
c yaml one_electron_integrals
           write(6,100) uhf_b_num(i,iha,ihb,ipa,ipb),
     &                  uhf_b_num(j,iha,ihb,ipa.ipb),dbl_mb(k_final+k-1)
          endif
         end if
        enddo
        enddo
        call util_flush(6)
        end if   !spust
       end if            !nodezero
c delete k_final
      if (.not.ma_pop_stack(l_final))
     1  call errquit('tcc_hcore_final: MA problem',15,MA_ERR)
c delete k_aux
      if (.not.ma_pop_stack(l_aux))
     1  call errquit('tcc_hcore_aux: MA problem',15,MA_ERR)
c delete k_movecs_orb
      if (.not.ma_pop_stack(l_movecs_orb))
     1  call errquit('tcc_hcore: MA problem',15,MA_ERR)
c
c
c
c delete MA l_core
      if (.not.ma_pop_stack(l_core))
     1     call errquit('l_core',15,MA_ERR)
c
        if (.not.ga_destroy(g_hcore))
     1    call errquit('deletefile: hcore GA problem',0,GA_ERR)
c
 100  format(2i5,3x,f18.10)
c
      call ga_sync()
      return 
      end
c
c
c
      integer function uhf_a_num(i,iha,ihb,ipa,ipb)
      implicit none
      integer i,iha,ihb,ipa,opb
      uhf_a_num=0
      if(i.le.iha) then
        uhf_a_num=i
      else
        uhf_a_num=i+ihb
      end if
      return 
      end
c
      integer function uhf_b_num(i,iha,ihb,ipa,ipb)
      implicit none
      integer i,iha,ihb,ipa,opb
      uhf_b_num=0
      if(i.le.ihb) then
        uhf_b_num=i+iha
      else
        uhf_b_num=i+iha+ipa
      end if
      return
      end


