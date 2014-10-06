      real    snd_tick

cif T3D_PVM
c     real    t_str, t_end
c     data    snd_tick/6.6e-9/
c#elif ASA_I860
       real    t_str, t_end
       real    ttime
       data    snd_tick/1./
c#else
c     real*4  t_str, t_end
c     real*4  rdum, etime, time_used(2)
c     data    snd_tick/1./
c#endif

c     IN SUBROUTINE sltdma
      real    t_fact, t_per_p1, t_ord_p1, t_com1_p1, t_com2_p1
      common  /c_sltdma_1/
     &        t_fact, t_per_p1, t_ord_p1, t_com1_p1, t_com2_p1

c    IN SUBROUTINE prolong_p1
      integer  nt_prol_p1
      parameter (nt_prol_p1=5)
      real     t_prol_p1(nt_prol_p1+1,l_max)
      common   /c_prolong_p1_1/
     &         t_prol_p1

c    IN SUBROUTINE restrict_p1
      integer  nt_rstr_p1
      parameter (nt_rstr_p1=5)
      real     t_rstr_p1(nt_rstr_p1+1,l_max)
      common   /c_restrict_p1_1/
     &         t_rstr_p1
