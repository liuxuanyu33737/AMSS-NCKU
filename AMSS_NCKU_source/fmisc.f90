

#include "macrodef.fh"
#if ENABLE_ABE_PROFILING
!=========================================================================
! Level-6 profiling support for f_global_interp() (the Cell-center
! global_interp below).  Pure diagnostics - no algorithmic or numerical
! change.  The module accumulates internal phase timings of every
! global_interp() call while abe_l6_gi_arm(1) is in force.  MPatch.C arms it
! only around the Level-5-marked Patch::Interp_Points() calls (surf_MassPAng
! cgh path), so all other f_global_interp() callers keep near-zero overhead.
! The accumulators are fetched/reduced/printed from C++ (abe_l5_interp_report
! in MPatch.C) through the bind(C) entry points below.
!=========================================================================
module abe_l6_gi_prof
  use iso_c_binding
  implicit none

  integer, parameter :: GI_NPHASE = 4
  integer, parameter :: GI_INDEX  = 1   ! index/bounds prep + box clamp + coord
  integer, parameter :: GI_DECIDE = 2   ! decide3d(): stencil fetch (sym/SoA)
  integer, parameter :: GI_POLIN  = 3   ! polin3(): 3-D polynomial interpolation
  integer, parameter :: GI_TOTAL  = 4   ! whole global_interp call

  logical(c_bool) :: l6_armed = .false.
  logical         :: l6_init  = .false.
  real(c_double)  :: l6_cnt(GI_NPHASE) = 0.0_c_double
  real(c_double)  :: l6_tot(GI_NPHASE) = 0.0_c_double
  real(c_double)  :: l6_mn(GI_NPHASE)  = huge(1.0_c_double)
  real(c_double)  :: l6_mx(GI_NPHASE)  = 0.0_c_double

contains

  subroutine abe_l6_gi_arm(on) bind(C, name="abe_l6_gi_arm")
    integer(c_int), value :: on
    if (on /= 0) then
      l6_armed = .true.
    else
      l6_armed = .false.
    end if
  end subroutine abe_l6_gi_arm

  subroutine abe_l6_gi_reset() bind(C, name="abe_l6_gi_reset")
    l6_cnt = 0.0_c_double
    l6_tot = 0.0_c_double
    l6_mn  = huge(1.0_c_double)
    l6_mx  = 0.0_c_double
    l6_init = .false.
  end subroutine abe_l6_gi_reset

  subroutine abe_l6_gi_get(ph, cnt, tot, mn, mx) bind(C, name="abe_l6_gi_get")
    integer(c_int), value :: ph
    real(c_double) :: cnt, tot, mn, mx
    cnt = 0.0_c_double
    tot = 0.0_c_double
    mn  = huge(1.0_c_double)
    mx  = 0.0_c_double
    if (ph >= 1 .and. ph <= GI_NPHASE) then
      cnt = l6_cnt(ph)
      tot = l6_tot(ph)
      mn  = l6_mn(ph)
      mx  = l6_mx(ph)
    end if
  end subroutine abe_l6_gi_get

  subroutine abe_l6_gi_rec(ph, dt)
    integer, intent(in)        :: ph
    real(c_double), intent(in) :: dt
    if (.not. l6_init) then
      call abe_l6_gi_reset()
      l6_init = .true.
    end if
    if (ph >= 1 .and. ph <= GI_NPHASE) then
      l6_cnt(ph) = l6_cnt(ph) + 1.0_c_double
      l6_tot(ph) = l6_tot(ph) + dt
      if (dt < l6_mn(ph)) l6_mn(ph) = dt
      if (dt > l6_mx(ph)) l6_mx(ph) = dt
    end if
  end subroutine abe_l6_gi_rec

  function abe_l6_wtime() result(t)
    real(c_double) :: t
    interface
      function gi_mpi_wtime() bind(C, name="MPI_Wtime") result(wt)
        import :: c_double
        real(c_double) :: wt
      end function gi_mpi_wtime
    end interface
    t = gi_mpi_wtime()
  end function abe_l6_wtime

end module abe_l6_gi_prof
!=========================================================================
! Level-7 profiling support for polin3() (the 3-D polynomial-interpolation
! kernel shared by global_interp()/interp_2()/global_interp_ss()).  Pure
! diagnostics - no algorithmic or numerical change.
!
! polin3() follows zyx order: it first interpolates along x3a (the third
! axis of ya, i.e. z for the Cell-center global_interp() callers) at m*n
! grid lines, then along x2a (second axis / y) at m lines, and finally a
! single polint along x1a (first axis / x).  This module accumulates the
! per-call wall time of each block (plus the trivial entry/setup segment)
! while abe_l7_p3_arm(1) is in force.  MPatch.C arms/disarms it together
! with the Level-6 f_global_interp() profile (abe_l6_gi_arm), i.e. only
! around the surf_MassPAng-marked Interp_Points() calls, so unmarked
! polin3() callers keep near-zero overhead.
!
! polin3() and polint() are not modified algorithmically.  The polint()
! invocation count per polin3() call is n*m + m + 1 (here m = n = ordn),
! accumulated in l7_polint so the report can separate call overhead from
! the polint() body cost.
!=========================================================================
module abe_l7_p3_prof
  use iso_c_binding
  implicit none

  integer, parameter :: P3_NPHASE = 5
  integer, parameter :: P3_PREP  = 1   ! entry/setup before the first polint
  integer, parameter :: P3_ZDIR  = 2   ! x3a-direction block (m*n polints)
  integer, parameter :: P3_YDIR  = 3   ! x2a-direction block (m polints)
  integer, parameter :: P3_XDIR  = 4   ! x1a-direction final polint (1)
  integer, parameter :: P3_TOTAL = 5   ! whole polin3() call

  logical(c_bool) :: l7_armed = .false.
  logical         :: l7_init  = .false.
  real(c_double)  :: l7_cnt(P3_NPHASE) = 0.0_c_double
  real(c_double)  :: l7_tot(P3_NPHASE) = 0.0_c_double
  real(c_double)  :: l7_mn(P3_NPHASE)  = huge(1.0_c_double)
  real(c_double)  :: l7_mx(P3_NPHASE)  = 0.0_c_double
  real(c_double)  :: l7_polint         = 0.0_c_double   ! polint() calls

contains

  subroutine abe_l7_p3_arm(on) bind(C, name="abe_l7_p3_arm")
    integer(c_int), value :: on
    if (on /= 0) then
      l7_armed = .true.
    else
      l7_armed = .false.
    end if
  end subroutine abe_l7_p3_arm

  subroutine abe_l7_p3_reset() bind(C, name="abe_l7_p3_reset")
    l7_cnt = 0.0_c_double
    l7_tot = 0.0_c_double
    l7_mn  = huge(1.0_c_double)
    l7_mx  = 0.0_c_double
    l7_polint = 0.0_c_double
    l7_init = .false.
  end subroutine abe_l7_p3_reset

  subroutine abe_l7_p3_get(ph, cnt, tot, mn, mx) bind(C, name="abe_l7_p3_get")
    integer(c_int), value :: ph
    real(c_double) :: cnt, tot, mn, mx
    cnt = 0.0_c_double
    tot = 0.0_c_double
    mn  = huge(1.0_c_double)
    mx  = 0.0_c_double
    if (ph >= 1 .and. ph <= P3_NPHASE) then
      cnt = l7_cnt(ph)
      tot = l7_tot(ph)
      mn  = l7_mn(ph)
      mx  = l7_mx(ph)
    end if
  end subroutine abe_l7_p3_get

  subroutine abe_l7_p3_get_polint(cnt) bind(C, name="abe_l7_p3_get_polint")
    real(c_double) :: cnt
    cnt = l7_polint
  end subroutine abe_l7_p3_get_polint

  subroutine p3_rec(ph, dt)
    integer, intent(in)        :: ph
    real(c_double), intent(in) :: dt
    if (.not. l7_init) then
      call abe_l7_p3_reset()
      l7_init = .true.
    end if
    if (ph >= 1 .and. ph <= P3_NPHASE) then
      l7_cnt(ph) = l7_cnt(ph) + 1.0_c_double
      l7_tot(ph) = l7_tot(ph) + dt
      if (dt < l7_mn(ph)) l7_mn(ph) = dt
      if (dt > l7_mx(ph)) l7_mx(ph) = dt
    end if
  end subroutine p3_rec

  function p3_wtime() result(t)
    real(c_double) :: t
    interface
      function p3_mpi_wtime() bind(C, name="MPI_Wtime") result(wt)
        import :: c_double
        real(c_double) :: wt
      end function p3_mpi_wtime
    end interface
    t = p3_mpi_wtime()
  end function p3_wtime

end module abe_l7_p3_prof
#endif


#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
!---------------------------------------------------------------------------------------------------
! copy a point of data into data target for vertext center code
!---------------------------------------------------------------------------------------------------
  subroutine pointcopy(wei,llbout,uubout,ext_out,data_out,xx,yy,zz,dv)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out
  real*8,dimension(3) :: llbout,uubout
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,intent(in) :: xx,yy,zz,dv

  real*8,dimension(3) :: ho
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::pointcopy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  if(any(ext_out == 1))then
       write(*,*)"fmisc.f90::pointcopy: meets iolated points for out data"
       write(*,*) llbout,uubout
       stop
  else
    ho = (uubout-llbout)/(ext_out-1)
  endif    
  i = idint((xx-llbout(1))/ho(1)+0.4)+1
  j = idint((yy-llbout(2))/ho(2)+0.4)+1
  k = idint((zz-llbout(3))/ho(3)+0.4)+1

  if(i<1 .or. i>ext_out(1) .or. &
     j<1 .or. j>ext_out(2) .or. &
     k<1 .or. k>ext_out(3) )then
     write(*,*)"i,j,k = ",i,j,k
     write(*,*)"ext = ",ext_out
     stop
  endif
  if(dabs(llbout(1)+(i-1)*ho(1)-xx)>ho(1)/2 .or. &
     dabs(llbout(2)+(j-1)*ho(2)-yy)>ho(2)/2 .or. &
     dabs(llbout(3)+(k-1)*ho(3)-zz)>ho(3)/2 )then    
     write(*,*)"fmisc.f90::pointcopy: llbout = ",llbout
     write(*,*)"fmisc.f90::pointcopy: ho = ",ho
     write(*,*)"fmisc.f90::pointcopy: x,y,z = ",llbout(1)+(i-1)*ho(1),llbout(2)+(j-1)*ho(2),llbout(3)+(k-1)*ho(3)
     write(*,*)"fmisc.f90::pointcopy: point = ",xx,yy,zz
     stop
   endif

  data_out(i,j,k)=dv

  return

  end subroutine pointcopy
!---------------------------------------------------------------------------------------------------
! copy a part of data from data source, for vertex center code
!---------------------------------------------------------------------------------------------------
  subroutine copy(wei,llbout,uubout,ext_out,data_out,llbin,uubin,ext_in,data_in,lcopy,ucopy)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out,ext_in
  real*8,dimension(3),intent(in) :: lcopy,ucopy
  real*8,dimension(3) :: llbout,uubout,llbin,uubin
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,dimension(ext_in(1),ext_in(2),ext_in(3)),intent(in)::data_in

  real*8,dimension(3) :: ho,hi
  integer,dimension(3) :: illo,iuuo,illi,iuui

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  if(any(ext_out == 1))then
    if(any(ext_in == 1))then
       write(*,*)"fmisc.f90::copy: meets iolated points for both in and out data"
       write(*,*) llbin,uubin
       write(*,*) llbout,uubout
       stop
     else
       hi = (uubin-llbin)/(ext_in-1)
       ho = hi
    endif
  else
    ho = (uubout-llbout)/(ext_out-1)
    if(any(ext_in == 1))then
      hi = ho
    else
      hi = (uubin-llbin)/(ext_in-1)
      if(any(abs(hi-ho) > min(hi,ho)/2))then
         write(*,*)"fmisc.f90::copy: meets copy reqest for different numerical grid"
         write(*,*)hi,ho
         stop
      endif
    endif
  endif    
  illo = idint((lcopy-llbout)/ho+0.4)+1
  iuuo = ext_out - idint((uubout-ucopy)/ho+0.4)
  illi = idint((lcopy-llbin)/hi+0.4)+1
  iuui = ext_in - idint((uubin-ucopy)/hi+0.4)

  if(any(llbout-lcopy>ho/2) .or. any(ucopy-uubout>ho/2))then
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     write(*,*)"fmisc.f90::copy: ho = ",ho
     write(*,*)llbout-lcopy,ucopy-uubout
     stop
  elseif(any(llbin -lcopy>hi/2) .or. any(ucopy-uubin >hi/2))then
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
  elseif(any(illo<1) .or. any(illi<1) .or. any(illo-iuuo>0) .or. any(illi-iuui>0) .or. &
         any(iuui-ext_in>0) .or. any(iuuo-ext_out>0))then
     write(*,*)"fmisc.f90::copy: illi = ",illi
     write(*,*)"fmisc.f90::copy: iuui = ",iuui
     write(*,*)"fmisc.f90::copy: illo = ",illo
     write(*,*)"fmisc.f90::copy: iuuo = ",iuuo     
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
   endif

  data_out(illo(1):iuuo(1),illo(2):iuuo(2),illo(3):iuuo(3))=data_in(illi(1):iuui(1),illi(2):iuui(2),illi(3):iuui(3))

  return

  end subroutine copy
!-----------------------------------------------------------------------------------------------------------------
! three dimensional interpolation for vertex center grid structure  
  subroutine global_interp(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0
  logical::decide3d

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1
       
  cmin = 1
  cmax = ex
  if(Symmetry == OCTANT  .and.dabs(X(1))<dX) cmin(1) = -ORDN/2+2
  if(Symmetry == OCTANT  .and.dabs(Y(1))<dY) cmin(2) = -ORDN/2+2
  if(Symmetry /= NO_SYMM .and.dabs(Z(1))<dZ) cmin(3) = -ORDN/2+2
  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 if(cxB(1)>0)then
  cx(1) = (x1 - X(cxB(1)))/dX
 else
  cx(1) = (x1 + X(2-cxB(1)))/dX
 endif
 if(cxB(2)>0)then
  cx(2) = (y1 - Y(cxB(2)))/dY
 else
  cx(2) = (y1 + Y(2-cxB(2)))/dY
 endif
 if(cxB(3)>0)then
  cx(3) = (z1 - Z(cxB(3)))/dZ
 else
  cx(3) = (z1 + Z(2-cxB(3)))/dZ
 endif

  if(decide3d(ex,f,f,cxB,cxT,SoA,ya,ORDN,Symmetry))then
     write(*,*)"global_interp position: ",x1,y1,z1
     write(*,*)"data range: ",X(1),X(ex(1)),Y(1),Y(ex(2)),Z(1),Z(ex(3))
     stop
  endif
  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp

!----------------------------------------------------------------
! decide which 3d data to be used does not surport PI-Symmetry yet 
!----------------------------------------------------------------
  function decide3d(ex,f,fpi,cxB,cxT,SoA,ya,ORDN,Symmetry)  result(gont)
  implicit none

  integer,                                 intent(in) :: ORDN,Symmetry
  integer,dimension(1:3)                 , intent(in) :: ex,cxB,cxT
  real*8, dimension(1:3)                 , intent(in) :: SoA
  real*8, dimension(ex(1),ex(2),ex(3))   , intent(in) :: f,fpi
  real*8, dimension(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3)), intent(out):: ya
  logical::gont

  integer,dimension(1:3) :: fmin1,fmin2,fmax1,fmax2
  integer::i,j,k,m

  gont=.false.
  do m=1,3
! check cxB and cxT are NaN or not  
    if(.not.(iabs(cxB(m)).ge.0)) gont=.true.
    if(.not.(iabs(cxT(m)).ge.0)) gont=.true.
    fmin1(m) = max(1,cxB(m))
    fmax1(m) = cxT(m)
    fmin2(m) = cxB(m)
    fmax2(m) = min(0,cxT(m))
    if((fmin1(m).le.fmax1(m)).and.(  fmin1(m)<1.or.  fmax1(m)>ex(m)))gont=.true.
    if((fmin2(m).le.fmax2(m)).and.(2-fmax2(m)<1.or.2-fmin2(m)>ex(m)))gont=.true.
  enddo
!sanity check
  if(gont)then
          write(*,*)"error in decide3d"
          write(*,*)((fmin1.le.fmax1).and.(  fmin1<1.or.  fmax1>ex))
          write(*,*)((fmin2.le.fmax2).and.(2-fmax2<1.or.2-fmin2>ex))
          write(*,*)"cxB, cxT and data shape:"
          write(*,*)cxB,cxT,ex
          write(*,*)"resulted fmin1, fmax1 and fmin2, fmax2:"
          write(*,*)fmin1,fmax1,fmin2,fmax2
  else

  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,k)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,j,k)*SoA(1)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,2-j,k)*SoA(2)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,2-j,k)*SoA(1)*SoA(2)
        enddo
     enddo
  enddo
           
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,2-k)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,j,2-k)*SoA(1)*SoA(3)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,2-j,2-k)*SoA(2)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(2-i,2-j,2-k)*SoA(1)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

  endif

  end function decide3d

!---------------------------------------------------------------------------------------
subroutine symmetry_bd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3)),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+2)*SoA(3)
   enddo

end subroutine symmetry_bd

subroutine symmetry_tbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,-ord+1:extc(3)+ord),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-1-i,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+2)*SoA(3)
      funcc(:,:,extc(3)+1+i) = funcc(:,:,extc(3)-1-i)*SoA(3)
   enddo

end subroutine symmetry_tbd

subroutine symmetry_stbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, dimension(2), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+2,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-1-i,1:extc(3))*SoA(2)
   enddo

end subroutine symmetry_stbd

subroutine symmetry_sntbd(ord,extc,func,funcc,SoA,actd)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord,actd
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
  if(actd==0)then
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+2,1:extc(2),1:extc(3))*SoA
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-1-i,1:extc(2),1:extc(3))*SoA
   enddo
  elseif(actd==1)then
   do i=0,ord-1
      funcc(1:extc(1),-i,1:extc(3)) = funcc(1:extc(1),i+2,1:extc(3))*SoA
      funcc(1:extc(1),extc(2)+1+i,1:extc(3)) = funcc(1:extc(1),extc(2)-1-i,1:extc(3))*SoA
   enddo
  else
    write(*,*)"symmetry_sntbd: not recognized actd = ",actd
  endif

end subroutine symmetry_sntbd


subroutine d2dump(wei,llb,uub,ext,data_in,data_out,gord,SoA)
  implicit none
  integer,             intent(in) :: wei,gord
  integer,dimension(3),intent(in) :: ext
  real*8, dimension(3),intent(in) :: SoA
  real*8, dimension(3) :: llb,uub
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in)   ::data_in
  real*8, dimension(ext(1),ext(2)),       intent(inout)::data_out

  real*8  :: dZ
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  dZ = (uub(3)-llb(3))/(ext(3)-1)
  k = idint((0-llb(3))/dZ+0.4)+1
     
  if(k < 1)then
      write(*,*) "d2dump: something must be wrong"
      return
  endif

  data_out(i,j) = data_in(i,j,k)

end subroutine d2dump

#else
#ifdef Cell
!subroutine interp_2 support cell center only
!-----------------------------------------------------------------------------
!
! Interpolate function f using weights Delx, Dely and Delz
!
!-----------------------------------------------------------------------------

  subroutine interp_2(ex,f,f_int,il,iu,jl,ju,kl,ku,Dx,Dy,Dz,&
                                     ordn,SoA,symmetry)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  integer,                             intent(in) :: il,iu,jl,ju,kl,ku,ordn
  real*8,                              intent(in) :: Dx,Dy,Dz,SoA(3)

!~~~~~~> Other parameters:

  integer :: j,imin,jmin,kmin
  real*8, dimension(1:ordn) :: x1a
  real*8, dimension(1:ordn,1:ordn,1:ordn) :: ya
  real*8, parameter :: ONE=1.d0
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: ddy,symX,symY,symZ

  symX = SoA(1)
  symY = SoA(2)
  symZ = SoA(3)
 
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  ya(2:ordn,2:ordn,2:ordn) = f(il+1:iu,jl+1:ju,kl+1:ku)
   
  if( il < imin .and. symmetry < OCTANT ) then
      write(*,*) 'Error in interp_2!!!'
      stop
  endif
  if( il < imin ) then
    ya(1,2:ordn,2:ordn) = f(imin,jl+1:ju,kl+1:ku)* symX
  else
    ya(1,2:ordn,2:ordn) = f(il  ,jl+1:ju,kl+1:ku)
  endif
  
  if( jl < jmin .and. symmetry < OCTANT ) then
     write(*,*) 'Error in interp_2!!!'
     stop
  endif
   
  if( jl < jmin ) then
    ya(2:ordn,1,2:ordn) = f(il+1:iu,jmin,kl+1:ku)* symY
  else
    ya(2:ordn,1,2:ordn) = f(il+1:iu,jl,kl+1:ku)
  endif

  if( kl < kmin .and. symmetry < EQUATORIAL ) then
    write(*,*)  'Error in interp_2!!!'
    stop
  endif

  if( kl < kmin ) then
   ya(2:ordn,2:ordn,1) = f(il+1:iu,jl+1:ju,kmin)* symZ
  else
   ya(2:ordn,2:ordn,1) = f(il+1:iu,jl+1:ju,kl  )
  endif

  if( il < imin .and. jl < jmin ) then
   ya(1,1,2:ordn) = f(imin,jmin,kl+1:ku)* symX * symY
  else if( il >= imin .and. jl < jmin ) then
   ya(1,1,2:ordn) = f(il,jmin,kl+1:ku)* symY
  else if( il < imin .and. jl >= jmin ) then
   ya(1,1,2:ordn) = f(imin,jl,kl+1:ku)* symX
  else 
   ya(1,1,2:ordn) = f(il,jl,kl+1:ku)
  endif
   
  if( il < imin .and. kl < kmin ) then
   ya(1,2:ordn,1) = f(imin,jl+1:ju,kmin)* symX * symZ
  else if( il >= imin .and. kl < kmin ) then
   ya(1,2:ordn,1) = f(il,jl+1:ju,kmin)* symZ
  else if( il < imin .and. kl >= kmin ) then
   ya(1,2:ordn,1) = f(imin,jl+1:ju,kl)* symX
  else
   ya(1,2:ordn,1) = f(il,jl+1:ju,kl)
  endif

  if( jl < jmin .and. kl < kmin ) then
    ya(2:ordn,1,1) = f(il+1:iu,jmin,kmin)* symY * symZ
  else if( jl >= jmin .and. kl < kmin ) then
    ya(2:ordn,1,1) =   f(il+1:iu,jl,kmin)* symZ
  else if( jl < jmin .and. kl >= kmin ) then
    ya(2:ordn,1,1) =   f(il+1:iu,jmin,kl)* symY
  else
   ya(2:ordn,1,1) = f(il+1:iu,jl,kl)
  endif

  if( il < imin ) then
   if( jl < jmin .and. kl < kmin) then
    ya(1,1,1) = f(imin,jmin,kmin)* symX * symY * symZ
   else if( jl >= jmin .and. kl < kmin ) then
    ya(1,1,1) = f(imin,jl,kmin)* symX * symZ
   else if( jl < jmin .and. kl >= kmin ) then
    ya(1,1,1) = f(imin,jmin,kl)* symX * symY
   else
    ya(1,1,1) = f(imin,jl,kl)* symX
   endif
  else
   if( jl < jmin .and. kl < kmin) then
    ya(1,1,1) = f(il,jmin,kmin)* symY * symZ
   else if( jl >= jmin .and. kl < kmin ) then
    ya(1,1,1) = f(il,jl,kmin)* symZ
   else if( jl < jmin .and. kl >= kmin ) then
    ya(1,1,1) = f(il,jmin,kl)* symY
   else
    ya(1,1,1) = f(il,jl,kl)
   endif
  endif
  
  call polin3(x1a,x1a,x1a,ya,Dx,Dy,Dz,f_int,ddy,ordn)

  if(.not.(dabs(f_int).ge.0))then
    write(*,*)"find nan in interp_2:",f_int,"inputs are:"
!    write(*,*)ya
!    write(*,*)"-----------------------------------------"
!    write(*,*)f(il:iu,jl:ju,kl:ku)
    write(*,*)Dx,Dy,Dz,symx,symy,symz,ordn
    write(*,*)il,iu,jl,ju,kl,ku,ex,symmetry
  endif

  return

  end subroutine interp_2
!---------------------------------------------------------------------------------------------------
! copy a point of data into data target for vertext center code
!---------------------------------------------------------------------------------------------------
  subroutine pointcopy(wei,llbout,uubout,ext_out,data_out,xx,yy,zz,dv)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out
  real*8,dimension(3) :: llbout,uubout
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,intent(in) :: xx,yy,zz,dv

  real*8,dimension(3) :: ho
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::pointcopy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  ho = (uubout-llbout)/ext_out
  i = idint((xx-llbout(1))/ho(1)+0.4)+1
  j = idint((yy-llbout(2))/ho(2)+0.4)+1
  k = idint((zz-llbout(3))/ho(3)+0.4)+1

  if(i<1 .or. i>ext_out(1) .or. &
     j<1 .or. j>ext_out(2) .or. &
     k<1 .or. k>ext_out(3) )then
     write(*,*)"i,j,k = ",i,j,k
     write(*,*)"ext = ",ext_out
     stop
  endif
  if(dabs(llbout(1)+(i-0.5)*ho(1)-xx)>ho(1)/2 .or. &
     dabs(llbout(2)+(j-0.5)*ho(2)-yy)>ho(2)/2 .or. &
     dabs(llbout(3)+(k-0.5)*ho(3)-zz)>ho(3)/2 )then    
     write(*,*)"fmisc.f90::pointcopy: llbout = ",llbout
     write(*,*)"fmisc.f90::pointcopy: ho = ",ho
     write(*,*)"fmisc.f90::pointcopy: x,y,z = ",llbout(1)+(i-0.5)*ho(1),llbout(2)+(j-0.5)*ho(2),llbout(3)+(k-0.5)*ho(3)
     write(*,*)"fmisc.f90::pointcopy: point = ",xx,yy,zz
     stop
   endif

  data_out(i,j,k)=dv

  return

  end subroutine pointcopy
!---------------------------------------------------------------------------------------------------
! copy a part of data from data source, for cell center code
!---------------------------------------------------------------------------------------------------
  subroutine copy(wei,llbout,uubout,ext_out,data_out,llbin,uubin,ext_in,data_in,lcopy,ucopy)
  implicit none
  integer,intent(in) :: wei
  integer,dimension(3),intent(in) ::ext_out,ext_in
  real*8,dimension(3),intent(in) :: lcopy,ucopy
  real*8,dimension(3) :: llbout,uubout,llbin,uubin
  real*8,dimension(ext_out(1),ext_out(2),ext_out(3)),intent(inout)::data_out
  real*8,dimension(ext_in(1),ext_in(2),ext_in(3)),intent(in)::data_in

  real*8,dimension(3) :: ho,hi
  integer,dimension(3) :: illo,iuuo,illi,iuui

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

!!!   
  ho = (uubout-llbout)/ext_out
  hi = (uubin-llbin)/ext_in
  illo = idint((lcopy-llbout)/ho+0.4)+1
  iuuo = ext_out - idint((uubout-ucopy)/ho+0.4)
  illi = idint((lcopy-llbin)/hi+0.4)+1
  iuui = ext_in - idint((uubin-ucopy)/hi+0.4)

  if(any(llbout-lcopy>ho/2) .or. any(ucopy-uubout>ho/2))then
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     write(*,*)"fmisc.f90::copy: ho = ",ho
     write(*,*)llbout-lcopy,ucopy-uubout
     stop
  elseif(any(llbin -lcopy>hi/2) .or. any(ucopy-uubin >hi/2))then
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
  elseif(any(illo<1) .or. any(illi<1) .or. any(illo-iuuo>0) .or. any(illi-iuui>0) .or. &
         any(iuui-ext_in>0) .or. any(iuuo-ext_out>0))then
     write(*,*)"fmisc.f90::copy: illi = ",illi
     write(*,*)"fmisc.f90::copy: iuui = ",iuui
     write(*,*)"fmisc.f90::copy: illo = ",illo
     write(*,*)"fmisc.f90::copy: iuuo = ",iuuo     
     write(*,*)"fmisc.f90::copy: llbout = ",llbout
     write(*,*)"fmisc.f90::copy: uubout = ",uubout
     write(*,*)"fmisc.f90::copy: llbin = ",llbin
     write(*,*)"fmisc.f90::copy: uubin = ",uubin
     write(*,*)"fmisc.f90::copy: llbcopy = ",lcopy
     write(*,*)"fmisc.f90::copy: uubcopy = ",ucopy
     stop
   endif

  data_out(illo(1):iuuo(1),illo(2):iuuo(2),illo(3):iuuo(3))=data_in(illi(1):iuui(1),illi(2):iuui(2),illi(3):iuui(3))

  return

  end subroutine copy
!--------------------------------------------------------------------------
! three dimensional interpolation for cell center grid structure  
  subroutine global_interp(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry)
#if ENABLE_ABE_PROFILING
  use abe_l6_gi_prof
#endif
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0
  logical::decide3d
#if ENABLE_ABE_PROFILING
  real*8 :: l6_t0, l6_t1, l6_t2, l6_t3

  if(l6_armed) l6_t0 = abe_l6_wtime()
#endif
  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1
       
  cmin = 1
  cmax = ex
  if(Symmetry == OCTANT  .and.dabs(X(1))<dX) cmin(1) = -ORDN/2+1
  if(Symmetry == OCTANT  .and.dabs(Y(1))<dY) cmin(2) = -ORDN/2+1
  if(Symmetry /= NO_SYMM .and.dabs(Z(1))<dZ) cmin(3) = -ORDN/2+1

  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 if(cxB(1)>0)then
  cx(1) = (x1 - X(cxB(1)))/dX
 else
  cx(1) = (x1 + X(1-cxB(1)))/dX
 endif
 if(cxB(2)>0)then
  cx(2) = (y1 - Y(cxB(2)))/dY
 else
  cx(2) = (y1 + Y(1-cxB(2)))/dY
 endif
 if(cxB(3)>0)then
  cx(3) = (z1 - Z(cxB(3)))/dZ
 else
  cx(3) = (z1 + Z(1-cxB(3)))/dZ
 endif
#if ENABLE_ABE_PROFILING
  if(l6_armed) l6_t1 = abe_l6_wtime()
#endif

  if(decide3d(ex,f,f,cxB,cxT,SoA,ya,ORDN,Symmetry))then
     write(*,*)"global_interp position: ",x1,y1,z1
     write(*,*)"data range: ",X(1),X(ex(1)),Y(1),Y(ex(2)),Z(1),Z(ex(3))
     stop
  endif
#if ENABLE_ABE_PROFILING
  if(l6_armed) l6_t2 = abe_l6_wtime()
#endif
  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)
#if ENABLE_ABE_PROFILING
  if(l6_armed)then
    l6_t3 = abe_l6_wtime()
    call abe_l6_gi_rec(GI_INDEX, l6_t1 - l6_t0)
    call abe_l6_gi_rec(GI_DECIDE, l6_t2 - l6_t1)
    call abe_l6_gi_rec(GI_POLIN,  l6_t3 - l6_t2)
    call abe_l6_gi_rec(GI_TOTAL,  l6_t3 - l6_t0)
  endif
#endif

  return

  end subroutine global_interp

!-----------------------------------------------------------------------------
! Cell-only fixed-order interpolation of 17 surf_MassPAng fields at one point.
! Geometry, source indices, and weights are shared; SoA signs are per field.
!-----------------------------------------------------------------------------
  subroutine global_interp_batch17(ex,X,Y,Z, &
       f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17, &
       yout,x1,y1,z1,ORDN,soa,Symmetry,status)
  implicit none

  integer, intent(in) :: ex(3),ORDN,Symmetry
  real*8, intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, intent(in) :: f1(ex(1),ex(2),ex(3)),f2(ex(1),ex(2),ex(3))
  real*8, intent(in) :: f3(ex(1),ex(2),ex(3)),f4(ex(1),ex(2),ex(3))
  real*8, intent(in) :: f5(ex(1),ex(2),ex(3)),f6(ex(1),ex(2),ex(3))
  real*8, intent(in) :: f7(ex(1),ex(2),ex(3)),f8(ex(1),ex(2),ex(3))
  real*8, intent(in) :: f9(ex(1),ex(2),ex(3)),f10(ex(1),ex(2),ex(3))
  real*8, intent(in) :: f11(ex(1),ex(2),ex(3)),f12(ex(1),ex(2),ex(3))
  real*8, intent(in) :: f13(ex(1),ex(2),ex(3)),f14(ex(1),ex(2),ex(3))
  real*8, intent(in) :: f15(ex(1),ex(2),ex(3)),f16(ex(1),ex(2),ex(3))
  real*8, intent(in) :: f17(ex(1),ex(2),ex(3))
  real*8, intent(out) :: yout(17)
  real*8, intent(in) :: x1,y1,z1,soa(3,17)
  integer, intent(out) :: status

  integer :: i,j,k,m,q
  integer :: cxB(3),cxT(3),cxI(3),cmin(3),cmax(3)
  integer :: src_i(6,6,6),src_j(6,6,6),src_k(6,6,6)
  integer :: reflect_mask(6,6,6),qi,qj,qk,mask
  integer, parameter :: NO_SYMM=0,EQUATORIAL=1,OCTANT=2
  real*8 :: dX,dY,dZ,cx(3),wx(6),wy(6),wz(6)
  real*8 :: ya(6,6,6),ztmp(6,6),ytmp(6)

  status=1
  yout=0.d0
  if(ORDN.ne.6) return
  if(minval(ex).lt.2) return

  dX=X(2)-X(1)
  dY=Y(2)-Y(1)
  dZ=Z(2)-Z(1)
  cxI(1)=idint((x1-X(1))/dX+0.4d0)+1
  cxI(2)=idint((y1-Y(1))/dY+0.4d0)+1
  cxI(3)=idint((z1-Z(1))/dZ+0.4d0)+1
  cxB=cxI-ORDN/2+1
  cxT=cxB+ORDN-1

  cmin=1
  cmax=ex
  if(Symmetry.eq.OCTANT .and. dabs(X(1)).lt.dX) cmin(1)=-ORDN/2+1
  if(Symmetry.eq.OCTANT .and. dabs(Y(1)).lt.dY) cmin(2)=-ORDN/2+1
  if(Symmetry.ne.NO_SYMM .and. dabs(Z(1)).lt.dZ) cmin(3)=-ORDN/2+1
  do m=1,3
    if(cxB(m).lt.cmin(m))then
      cxB(m)=cmin(m)
      cxT(m)=cxB(m)+ORDN-1
    endif
    if(cxT(m).gt.cmax(m))then
      cxT(m)=cmax(m)
      cxB(m)=cxT(m)+1-ORDN
    endif
  enddo

  if(cxB(1).gt.0)then
    cx(1)=(x1-X(cxB(1)))/dX
  else
    cx(1)=(x1+X(1-cxB(1)))/dX
  endif
  if(cxB(2).gt.0)then
    cx(2)=(y1-Y(cxB(2)))/dY
  else
    cx(2)=(y1+Y(1-cxB(2)))/dY
  endif
  if(cxB(3).gt.0)then
    cx(3)=(z1-Z(cxB(3)))/dZ
  else
    cx(3)=(z1+Z(1-cxB(3)))/dZ
  endif

  do k=1,6
    qk=cxB(3)+k-1
    do j=1,6
      qj=cxB(2)+j-1
      do i=1,6
        qi=cxB(1)+i-1
        src_i(i,j,k)=merge(qi,1-qi,qi.gt.0)
        src_j(i,j,k)=merge(qj,1-qj,qj.gt.0)
        src_k(i,j,k)=merge(qk,1-qk,qk.gt.0)
        mask=0
        if(qi.le.0) mask=ibset(mask,0)
        if(qj.le.0) mask=ibset(mask,1)
        if(qk.le.0) mask=ibset(mask,2)
        reflect_mask(i,j,k)=mask
        if(src_i(i,j,k).lt.1 .or. src_i(i,j,k).gt.ex(1) .or. &
           src_j(i,j,k).lt.1 .or. src_j(i,j,k).gt.ex(2) .or. &
           src_k(i,j,k).lt.1 .or. src_k(i,j,k).gt.ex(3)) return
      enddo
    enddo
  enddo

  do i=1,6
    wx(i)=1.d0
    wy(i)=1.d0
    wz(i)=1.d0
    do q=1,6
      if(q.ne.i)then
        wx(i)=wx(i)*(cx(1)-dble(q-1))/dble(i-q)
        wy(i)=wy(i)*(cx(2)-dble(q-1))/dble(i-q)
        wz(i)=wz(i)*(cx(3)-dble(q-1))/dble(i-q)
      endif
    enddo
  enddo

  call eval_field(f1, soa(:,1), yout(1))
  call eval_field(f2, soa(:,2), yout(2))
  call eval_field(f3, soa(:,3), yout(3))
  call eval_field(f4, soa(:,4), yout(4))
  call eval_field(f5, soa(:,5), yout(5))
  call eval_field(f6, soa(:,6), yout(6))
  call eval_field(f7, soa(:,7), yout(7))
  call eval_field(f8, soa(:,8), yout(8))
  call eval_field(f9, soa(:,9), yout(9))
  call eval_field(f10,soa(:,10),yout(10))
  call eval_field(f11,soa(:,11),yout(11))
  call eval_field(f12,soa(:,12),yout(12))
  call eval_field(f13,soa(:,13),yout(13))
  call eval_field(f14,soa(:,14),yout(14))
  call eval_field(f15,soa(:,15),yout(15))
  call eval_field(f16,soa(:,16),yout(16))
  call eval_field(f17,soa(:,17),yout(17))
  status=0
  return

  contains

  subroutine eval_field(f,soav,value)
  implicit none
  real*8, intent(in) :: f(ex(1),ex(2),ex(3)),soav(3)
  real*8, intent(out) :: value
  real*8 :: signv
  integer :: ii,jj,kk,local_mask

  do kk=1,6
    do jj=1,6
      do ii=1,6
        signv=1.d0
        local_mask=reflect_mask(ii,jj,kk)
        if(btest(local_mask,0)) signv=signv*soav(1)
        if(btest(local_mask,1)) signv=signv*soav(2)
        if(btest(local_mask,2)) signv=signv*soav(3)
        ya(ii,jj,kk)=f(src_i(ii,jj,kk),src_j(ii,jj,kk), &
                          src_k(ii,jj,kk))*signv
      enddo
    enddo
  enddo
  do jj=1,6
    do ii=1,6
      ztmp(ii,jj)=0.d0
      do kk=1,6
        ztmp(ii,jj)=ztmp(ii,jj)+wz(kk)*ya(ii,jj,kk)
      enddo
    enddo
  enddo
  do ii=1,6
    ytmp(ii)=0.d0
    do jj=1,6
      ytmp(ii)=ytmp(ii)+wy(jj)*ztmp(ii,jj)
    enddo
  enddo
  value=0.d0
  do ii=1,6
    value=value+wx(ii)*ytmp(ii)
  enddo
  end subroutine eval_field

  end subroutine global_interp_batch17
!----------------------------------------------------------------
! decide which 3d data to be used does not surport PI-Symmetry yet 
!----------------------------------------------------------------
  function decide3d(ex,f,fpi,cxB,cxT,SoA,ya,ORDN,Symmetry)  result(gont)
  implicit none

  integer,                                 intent(in) :: ORDN,Symmetry
  integer,dimension(1:3)                 , intent(in) :: ex,cxB,cxT
  real*8, dimension(1:3)                 , intent(in) :: SoA
  real*8, dimension(ex(1),ex(2),ex(3))   , intent(in) :: f,fpi
  real*8, dimension(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3)), intent(out):: ya
  logical::gont

  integer,dimension(1:3) :: fmin1,fmin2,fmax1,fmax2
  integer::i,j,k,m

  gont=.false.
  do m=1,3
! check cxB and cxT are NaN or not  
    if(.not.(iabs(cxB(m)).ge.0)) gont=.true.
    if(.not.(iabs(cxT(m)).ge.0)) gont=.true.
    fmin1(m) = max(1,cxB(m))
    fmax1(m) = cxT(m)
    fmin2(m) = cxB(m)
    fmax2(m) = min(0,cxT(m))
    if((fmin1(m).le.fmax1(m)).and.(  fmin1(m)<1.or.  fmax1(m)>ex(m)))gont=.true.
    if((fmin2(m).le.fmax2(m)).and.(1-fmax2(m)<1.or.1-fmin2(m)>ex(m)))gont=.true.
  enddo
!sanity check
  if(gont)then
          write(*,*)"error in decide3d"
          write(*,*)((fmin1.le.fmax1).and.(  fmin1<1.or.  fmax1>ex))
          write(*,*)((fmin2.le.fmax2).and.(1-fmax2<1.or.1-fmin2>ex))
          write(*,*)"cxB, cxT and data shape:"
          write(*,*)cxB,cxT,ex
          write(*,*)"resulted fmin1, fmax1 and fmin2, fmax2:"
          write(*,*)fmin1,fmax1,fmin2,fmax2
  else

  do k=fmin1(3),fmax1(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,k)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,j,k)*SoA(1)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,1-j,k)*SoA(2)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,1-j,k)*SoA(1)*SoA(2)
        enddo
     enddo
  enddo
           
  do k=fmin2(3),fmax2(3)
     do j=fmin1(2),fmax1(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,j,1-k)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,j,1-k)*SoA(1)*SoA(3)
        enddo
     enddo
     do j=fmin2(2),fmax2(2)
        do i=fmin1(1),fmax1(1)
           ya(i,j,k) = f(i,1-j,1-k)*SoA(2)*SoA(3)
        enddo
        do i=fmin2(1),fmax2(1)
           ya(i,j,k) = f(1-i,1-j,1-k)*SoA(1)*SoA(2)*SoA(3)
        enddo
     enddo
  enddo

  endif

  end function decide3d

!---------------------------------------------------------------------------------------
subroutine symmetry_bd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3)),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3)
   enddo

end subroutine symmetry_bd

subroutine symmetry_tbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,-ord+1:extc(3)+ord),intent(out):: funcc
  real*8, dimension(1:3), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-i,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3)
      funcc(:,:,extc(3)+1+i) = funcc(:,:,extc(3)-i)*SoA(3)
   enddo

end subroutine symmetry_tbd

subroutine symmetry_stbd(ord,extc,func,funcc,SoA)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, dimension(2), intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
      funcc(:,extc(2)+1+i,1:extc(3)) = funcc(:,extc(2)-i,1:extc(3))*SoA(2)
   enddo

end subroutine symmetry_stbd

subroutine symmetry_sntbd(ord,extc,func,funcc,SoA,actd)
  implicit none

!~~~~~~> input arguments
  integer,intent(in) :: ord,actd
  integer,dimension(3),   intent(in) :: extc
  real*8, dimension(extc(1),extc(2),extc(3)),intent(in ):: func
  real*8, dimension(-ord+1:extc(1)+ord,-ord+1:extc(2)+ord,extc(3)),intent(out):: funcc
  real*8, intent(in) :: SoA

  integer::i

  funcc = 0.d0
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
  if(actd==0)then
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA
      funcc(extc(1)+1+i,1:extc(2),1:extc(3)) = funcc(extc(1)-i,1:extc(2),1:extc(3))*SoA
   enddo
  elseif(actd==1)then
   do i=0,ord-1
      funcc(1:extc(1),-i,1:extc(3)) = funcc(1:extc(1),i+1,1:extc(3))*SoA
      funcc(1:extc(1),extc(2)+1+i,1:extc(3)) = funcc(1:extc(1),extc(2)-i,1:extc(3))*SoA
   enddo
  else
    write(*,*)"symmetry_sntbd: not recognized actd = ",actd
  endif

end subroutine symmetry_sntbd

subroutine d2dump(wei,llb,uub,ext,data_in,data_out,gord,SoA)
  implicit none
  integer,intent(in) :: wei,gord
  integer,dimension(3),intent(in) ::ext
  real*8,dimension(3),intent(in) :: SoA
  real*8,dimension(3) :: llb,uub
  real*8,dimension(ext(1),ext(2),ext(3)),intent(in)::data_in
  real*8,dimension(ext(1),ext(2)),intent(inout)::data_out

  real*8 :: dZ
  integer :: i,j,k

!sanity check
  if(wei.ne.3)then
     write(*,*)"fmisc.f90::copy: this routine only surport 3 dimension"
     write(*,*)"dim = ",wei
     stop
  endif

  dZ = (uub(3)-llb(3))/ext(3)
  k = idint((0-llb(3))/dZ+0.4)+1

  select case (gord)
  case (2)
     if(k > 2)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5625d0*(data_in(i,j,k)+data_in(i,j,k-1))-0.0625d0*(data_in(i,j,k+1)+data_in(i,j,k-2))
         enddo
         enddo
     else if(k == 1)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5625d0*(data_in(i,j,k)+SoA(3)*data_in(i,j,k))-0.0625d0*(data_in(i,j,k+1)+SoA(3)*data_in(i,j,k+1))
         enddo
         enddo
     else
      write(*,*) "d2dump: something must be wrong, k = ",k
      return
     endif
  case (3)
     if(k > 3)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5859375d0*(data_in(i,j,k)+data_in(i,j,k-1))   &
                           -0.9765625d-1*(data_in(i,j,k+1)+data_in(i,j,k-2)) &
                           +0.1171875d-1*(data_in(i,j,k+2)+data_in(i,j,k-3))
         enddo
         enddo
     else if(k == 1)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5859375d0*(data_in(i,j,k)+SoA(3)*data_in(i,j,k))   &
                           -0.9765625d-1*(data_in(i,j,k+1)+SoA(3)*data_in(i,j,k+1)) &
                           +0.1171875d-1*(data_in(i,j,k+2)+SoA(3)*data_in(i,j,k+2))
         enddo
         enddo
     else
      write(*,*) "d2dump: something must be wrong, k = ",k
      return
     endif
  case (4)
     if(k > 4)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5981445312d0*(data_in(i,j,k)+data_in(i,j,k-1))   &
                            -0.1196289063d0*(data_in(i,j,k+1)+data_in(i,j,k-2)) &
                           +0.2392578125d-1*(data_in(i,j,k+2)+data_in(i,j,k-3)) &
                           -0.2441406250d-2*(data_in(i,j,k+3)+data_in(i,j,k-4))
         enddo
         enddo
     else if(k == 1)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.5981445312d0*(data_in(i,j,k)+SoA(3)*data_in(i,j,k))   &
                            -0.1196289063d0*(data_in(i,j,k+1)+SoA(3)*data_in(i,j,k+1)) &
                           +0.2392578125d-1*(data_in(i,j,k+2)+SoA(3)*data_in(i,j,k+2)) &
                           -0.2441406250d-2*(data_in(i,j,k+3)+SoA(3)*data_in(i,j,k+3))
         enddo
         enddo
     else
      write(*,*) "d2dump: something must be wrong, k = ",k
      return
     endif
  case (5)
     if(k > 5)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.6056213378d0*(data_in(i,j,k)+data_in(i,j,k-1))   &
                            -0.1345825196d0*(data_in(i,j,k+1)+data_in(i,j,k-2)) &
                           +0.3460693359d-1*(data_in(i,j,k+2)+data_in(i,j,k-3)) &
                           -0.6179809571d-2*(data_in(i,j,k+3)+data_in(i,j,k-4)) &
                           +0.5340576171d-3*(data_in(i,j,k+4)+data_in(i,j,k-5))
         enddo
         enddo
     else if(k == 1)then
         do i=1,ext(1)
         do j=1,ext(2)
             data_out(i,j) = 0.6056213378d0*(data_in(i,j,k)+SoA(3)*data_in(i,j,k))   &
                            -0.1345825196d0*(data_in(i,j,k+1)+SoA(3)*data_in(i,j,k+1)) &
                           +0.3460693359d-1*(data_in(i,j,k+2)+SoA(3)*data_in(i,j,k+2)) &
                           -0.6179809571d-2*(data_in(i,j,k+3)+SoA(3)*data_in(i,j,k+3)) &
                           +0.5340576171d-3*(data_in(i,j,k+4)+SoA(3)*data_in(i,j,k+4))
         enddo
         enddo
     else
      write(*,*) "d2dump: something must be wrong, k = ",k
      return
     endif
  case default
    write(*,*) "d2dump: not recognized ord = ",gord
    return
  end select

end subroutine d2dump

#else
#error Not define Vertex nor Cell
#endif  
#endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! common code for cell and vertex
!------------------------------------------------------------------------------
! Fixed-size six-point Lagrangian polynomial interpolation.
! This is the hot path used by prolongation when ghost_width == 3.
! Keep polint below as the general-order implementation used by polin2/polin3.
!------------------------------------------------------------------------------

  subroutine polint6_fast(xa,ya,x,y,dy)

  implicit none

!~~~~~~> Input Parameter:
  real*8, dimension(6), intent(in) :: xa,ya
  real*8, intent(in) :: x
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameter:

  integer :: ns
  real*8 :: c1,c2,c3,c4,c5,c6
  real*8 :: d1,d2,d3,d4,d5,d6
  real*8 :: h1,h2,h3,h4,h5,h6
  real*8 :: den,dif,dift

!~~~~~~>

  c1=ya(1); c2=ya(2); c3=ya(3)
  c4=ya(4); c5=ya(5); c6=ya(6)
  d1=ya(1); d2=ya(2); d3=ya(3)
  d4=ya(4); d5=ya(5); d6=ya(6)
  h1=xa(1)-x; h2=xa(2)-x; h3=xa(3)-x
  h4=xa(4)-x; h5=xa(5)-x; h6=xa(6)-x

  ns=1
  dif=dabs(x-xa(1))
  dift=dabs(x-xa(2)); if(dift < dif)then; ns=2; dif=dift; endif
  dift=dabs(x-xa(3)); if(dift < dif)then; ns=3; dif=dift; endif
  dift=dabs(x-xa(4)); if(dift < dif)then; ns=4; dif=dift; endif
  dift=dabs(x-xa(5)); if(dift < dif)then; ns=5; dif=dift; endif
  dift=dabs(x-xa(6)); if(dift < dif)then; ns=6; endif

  y=ya(ns)
  ns=ns-1

! First Neville column
  den=h1-h2; if(den == 0.d0) goto 900
  den=(c2-d1)/den; d1=h2*den; c1=h1*den
  den=h2-h3; if(den == 0.d0) goto 900
  den=(c3-d2)/den; d2=h3*den; c2=h2*den
  den=h3-h4; if(den == 0.d0) goto 900
  den=(c4-d3)/den; d3=h4*den; c3=h3*den
  den=h4-h5; if(den == 0.d0) goto 900
  den=(c5-d4)/den; d4=h5*den; c4=h4*den
  den=h5-h6; if(den == 0.d0) goto 900
  den=(c6-d5)/den; d5=h6*den; c5=h5*den
  if(2*ns < 5)then
    if(ns == 0) dy=c1
    if(ns == 1) dy=c2
    if(ns == 2) dy=c3
  else
    if(ns == 3) dy=d3
    if(ns == 4) dy=d4
    if(ns == 5) dy=d5
    ns=ns-1
  endif
  y=y+dy

! Second Neville column
  den=h1-h3; if(den == 0.d0) goto 900
  den=(c2-d1)/den; d1=h3*den; c1=h1*den
  den=h2-h4; if(den == 0.d0) goto 900
  den=(c3-d2)/den; d2=h4*den; c2=h2*den
  den=h3-h5; if(den == 0.d0) goto 900
  den=(c4-d3)/den; d3=h5*den; c3=h3*den
  den=h4-h6; if(den == 0.d0) goto 900
  den=(c5-d4)/den; d4=h6*den; c4=h4*den
  if(2*ns < 4)then
    if(ns == 0) dy=c1
    if(ns == 1) dy=c2
  else
    if(ns == 2) dy=d2
    if(ns == 3) dy=d3
    if(ns == 4) dy=d4
    ns=ns-1
  endif
  y=y+dy

! Third Neville column
  den=h1-h4; if(den == 0.d0) goto 900
  den=(c2-d1)/den; d1=h4*den; c1=h1*den
  den=h2-h5; if(den == 0.d0) goto 900
  den=(c3-d2)/den; d2=h5*den; c2=h2*den
  den=h3-h6; if(den == 0.d0) goto 900
  den=(c4-d3)/den; d3=h6*den; c3=h3*den
  if(2*ns < 3)then
    if(ns == 0) dy=c1
    if(ns == 1) dy=c2
  else
    if(ns == 2) dy=d2
    if(ns == 3) dy=d3
    ns=ns-1
  endif
  y=y+dy

! Fourth Neville column
  den=h1-h5; if(den == 0.d0) goto 900
  den=(c2-d1)/den; d1=h5*den; c1=h1*den
  den=h2-h6; if(den == 0.d0) goto 900
  den=(c3-d2)/den; d2=h6*den; c2=h2*den
  if(2*ns < 2)then
    if(ns == 0) dy=c1
  else
    if(ns == 1) dy=d1
    if(ns == 2) dy=d2
    ns=ns-1
  endif
  y=y+dy

! Fifth Neville column
  den=h1-h6; if(den == 0.d0) goto 900
  den=(c2-d1)/den; d1=h6*den; c1=h1*den
  if(2*ns < 1)then
    dy=c1
  else
    dy=d1
  endif
  y=y+dy

  return

900 continue
  write(*,*) 'failure in polint6_fast for point',x
  write(*,*) 'with input points: ',xa
  stop

  end subroutine polint6_fast
!------------------------------------------------------------------------------
! Lagrangian polynomial interpolation
!------------------------------------------------------------------------------

  subroutine polint(xa,ya,x,y,dy,ordn)

  implicit none

!~~~~~~> Input Parameter:
  integer,intent(in) :: ordn
  real*8, dimension(ordn), intent(in) :: xa,ya
  real*8, intent(in) :: x
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameter:

  integer :: m,n,ns
  real*8, dimension(ordn) :: c,d,den,ho
  real*8 :: dif,dift

!~~~~~~>

  n=ordn
  m=ordn

  if (n == 6) then
  c(1:6)=ya(1:6)
  d(1:6)=ya(1:6)
  ho(1:6)=xa(1:6)-x

  ns=1
  dif=abs(x-xa(1))
  dift=abs(x-xa(2))
  if(dift < dif) then
    ns=2
    dif=dift
  end if
  dift=abs(x-xa(3))
  if(dift < dif) then
    ns=3
    dif=dift
  end if
  dift=abs(x-xa(4))
  if(dift < dif) then
    ns=4
    dif=dift
  end if
  dift=abs(x-xa(5))
  if(dift < dif) then
    ns=5
    dif=dift
  end if
  dift=abs(x-xa(6))
  if(dift < dif) then
    ns=6
  end if

  y=ya(ns)
  ns=ns-1

  den(1:5)=ho(1:5)-ho(2:6)
  if (any(den(1:5) == 0.0))then
    write(*,*) 'failure in polint for point',x
    write(*,*) 'with input points: ',xa
    stop
  endif
  den(1:5)=(c(2:6)-d(1:5))/den(1:5)
  d(1:5)=ho(2:6)*den(1:5)
  c(1:5)=ho(1:5)*den(1:5)
  if (2*ns < 5) then
    dy=c(ns+1)
  else
    dy=d(ns)
    ns=ns-1
  end if
  y=y+dy

  den(1:4)=ho(1:4)-ho(3:6)
  if (any(den(1:4) == 0.0))then
    write(*,*) 'failure in polint for point',x
    write(*,*) 'with input points: ',xa
    stop
  endif
  den(1:4)=(c(2:5)-d(1:4))/den(1:4)
  d(1:4)=ho(3:6)*den(1:4)
  c(1:4)=ho(1:4)*den(1:4)
  if (2*ns < 4) then
    dy=c(ns+1)
  else
    dy=d(ns)
    ns=ns-1
  end if
  y=y+dy

  den(1:3)=ho(1:3)-ho(4:6)
  if (any(den(1:3) == 0.0))then
    write(*,*) 'failure in polint for point',x
    write(*,*) 'with input points: ',xa
    stop
  endif
  den(1:3)=(c(2:4)-d(1:3))/den(1:3)
  d(1:3)=ho(4:6)*den(1:3)
  c(1:3)=ho(1:3)*den(1:3)
  if (2*ns < 3) then
    dy=c(ns+1)
  else
    dy=d(ns)
    ns=ns-1
  end if
  y=y+dy

  den(1:2)=ho(1:2)-ho(5:6)
  if (any(den(1:2) == 0.0))then
    write(*,*) 'failure in polint for point',x
    write(*,*) 'with input points: ',xa
    stop
  endif
  den(1:2)=(c(2:3)-d(1:2))/den(1:2)
  d(1:2)=ho(5:6)*den(1:2)
  c(1:2)=ho(1:2)*den(1:2)
  if (2*ns < 2) then
    dy=c(ns+1)
  else
    dy=d(ns)
    ns=ns-1
  end if
  y=y+dy

  den(1)=ho(1)-ho(6)
  if (den(1) == 0.0)then
    write(*,*) 'failure in polint for point',x
    write(*,*) 'with input points: ',xa
    stop
  endif
  den(1)=(c(2)-d(1))/den(1)
  d(1)=ho(6)*den(1)
  c(1)=ho(1)*den(1)
  if (2*ns < 1) then
    dy=c(ns+1)
  else
    dy=d(ns)
    ns=ns-1
  end if
  y=y+dy

  else
  c=ya
  d=ya
  ho=xa-x

  ns=1
  dif=abs(x-xa(1))
  do m=1,n
   dift=abs(x-xa(m))
   if(dift < dif) then
    ns=m
    dif=dift
   end if
  end do

  y=ya(ns)
  ns=ns-1
  do m=1,n-1
    den(1:n-m)=ho(1:n-m)-ho(1+m:n)
    if (any(den(1:n-m) == 0.0))then
      write(*,*) 'failure in polint for point',x
      write(*,*) 'with input points: ',xa
      stop
    endif
    den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
    d(1:n-m)=ho(1+m:n)*den(1:n-m)
    c(1:n-m)=ho(1:n-m)*den(1:n-m)
    if (2*ns < n-m) then
      dy=c(ns+1)
    else
      dy=d(ns)
      ns=ns-1
    end if
    y=y+dy
  end do
  end if

  return

  end subroutine polint
!------------------------------------------------------------------------------
!
! interpolation in 2 dimensions, follow yx order
!
!------------------------------------------------------------------------------
  subroutine polin2(x1a,x2a,ya,x1,x2,y,dy,ordn)

  implicit none

!~~~~~~> Input parameters:
  integer,intent(in) :: ordn
  real*8, dimension(1:ordn), intent(in) :: x1a,x2a
  real*8, dimension(1:ordn,1:ordn), intent(in) :: ya
  real*8, intent(in) :: x1,x2
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameters:

  integer  :: i,m
  real*8, dimension(ordn) :: ymtmp
  real*8, dimension(ordn) :: yntmp

  m=size(x1a)
  
  do i=1,m

    yntmp=ya(i,:)
    call polint(x2a,yntmp,x2,ymtmp(i),dy,ordn)

  end do

  call polint(x1a,ymtmp,x1,y,dy,ordn)

  return

  end subroutine polin2
!------------------------------------------------------------------------------
!
! interpolation in 3 dimensions, follow zyx order
!
!------------------------------------------------------------------------------
  subroutine polin3(x1a,x2a,x3a,ya,x1,x2,x3,y,dy,ordn)
#if ENABLE_ABE_PROFILING
  use abe_l7_p3_prof
#endif

  implicit none

!~~~~~~> Input parameters:
  integer,intent(in) :: ordn
  real*8, dimension(1:ordn), intent(in) :: x1a,x2a,x3a
  real*8, dimension(1:ordn,1:ordn,1:ordn), intent(in) :: ya
  real*8, intent(in) :: x1,x2,x3
  real*8, intent(out) :: y,dy

!~~~~~~> Other parameters:

  integer  :: i,j,m,n
  real*8, dimension(ordn,ordn) :: yatmp
  real*8, dimension(ordn) :: ymtmp
  real*8, dimension(ordn) :: yntmp
  real*8, dimension(ordn) :: yqtmp

#if ENABLE_ABE_PROFILING
! Level-7 profiling locals (pure diagnostics, no numerical effect; module
! abe_l7_p3_prof above).  Timing is gated by the l7_armed flag that MPatch.C
! raises only around the surf_MassPAng-marked Interp_Points() calls (same
! window that arms the Level-6 f_global_interp() profile), so unarmed
! polin3() callers pay one logical test per call.
  real*8  :: p3_t0,p3_tb,p3_te,p3_tz,p3_ty,p3_tx
  logical :: p3_on

  p3_on = l7_armed
  if(p3_on) p3_t0 = p3_wtime()
#endif

  m=size(x1a)
  n=size(x2a)
#if ENABLE_ABE_PROFILING
  if(p3_on)then
    call p3_rec(P3_PREP,p3_wtime()-p3_t0)
    p3_tz = 0.d0
    p3_ty = 0.d0
  endif
#endif

  do i=1,m
#if ENABLE_ABE_PROFILING
   if(p3_on) p3_tb = p3_wtime()
#endif
   do j=1,n

    yqtmp=ya(i,j,:)
    call polint(x3a,yqtmp,x3,yatmp(i,j),dy,ordn)

   end do
#if ENABLE_ABE_PROFILING
   if(p3_on) p3_tz = p3_tz + (p3_wtime()-p3_tb)
   if(p3_on) p3_tb = p3_wtime()
#endif

    yntmp=yatmp(i,:)
    call polint(x2a,yntmp,x2,ymtmp(i),dy,ordn)

#if ENABLE_ABE_PROFILING
   if(p3_on) p3_ty = p3_ty + (p3_wtime()-p3_tb)
#endif
  end do

#if ENABLE_ABE_PROFILING
  if(p3_on) p3_te = p3_wtime()
#endif
  call polint(x1a,ymtmp,x1,y,dy,ordn)
#if ENABLE_ABE_PROFILING
  if(p3_on)then
    p3_tx = p3_wtime()-p3_te
    call p3_rec(P3_TOTAL,p3_wtime()-p3_t0)
    call p3_rec(P3_ZDIR,p3_tz)
    call p3_rec(P3_YDIR,p3_ty)
    call p3_rec(P3_XDIR,p3_tx)
    l7_polint = l7_polint + n*m + m + 1
  endif
#endif

  return

  end subroutine polin3
!--------------------------------------------------------------------------------------
! calculate L2norm  
  subroutine l2normhelper(ex, X, Y, Z,xmin,ymin,zmin,xmax,ymax,zmax,&
                          f,f_out,gw)

  implicit none
!~~~~~~> Input parameters:
  integer,intent(in ):: ex(1:3)
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3)),xmin,ymin,zmin,xmax,ymax,zmax
  integer,intent(in)::gw
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8, intent(out) :: f_out
!~~~~~~> Other variables:

  real*8, parameter :: ZEO = 0.D0
  real*8            :: dX, dY, dZ
  integer::imin,jmin,kmin
  integer::imax,jmax,kmax
  integer::i,j,k

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

! for ghost zone
   imin = gw+1
   jmin = gw+1
   kmin = gw+1

   imax = ex(1) - gw
   jmax = ex(2) - gw
   kmax = ex(3) - gw

!for patch boundary (i.e., not ghost boundary)

if(dabs(X(ex(1))-xmax) < dX) imax = ex(1)
if(dabs(Y(ex(2))-ymax) < dY) jmax = ex(2)
if(dabs(Z(ex(3))-zmax) < dZ) kmax = ex(3)
if(dabs(X(1)-xmin) < dX) imin = 1
if(dabs(Y(1)-ymin) < dY) jmin = 1
if(dabs(Z(1)-zmin) < dZ) kmin = 1

f_out = sum(f(imin:imax,jmin:jmax,kmin:kmax)*f(imin:imax,jmin:jmax,kmin:kmax))

f_out = f_out*dX*dY*dZ

  return

  end subroutine l2normhelper
!--------------------------------------------------------------------------------------
! calculate L2norm especially for shell Blocks
  subroutine l2normhelper_sh(ex, X, Y, Z,xmin,ymin,zmin,xmax,ymax,zmax,&
                          f,f_out,gw,ogw,Symmetry)

  implicit none
!~~~~~~> Input parameters:
  integer,intent(in ):: ex(1:3),Symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3)),xmin,ymin,zmin,xmax,ymax,zmax
  integer,intent(in)::gw,ogw
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8, intent(out) :: f_out
!~~~~~~> Other variables:

  real*8, parameter :: ZEO = 0.D0
  real*8            :: dX, dY, dZ
  integer::imin,jmin,kmin
  integer::imax,jmax,kmax
  integer::i,j,k

  real*8 :: PIo4

  PIo4 = dacos(-1.d0)/4.d0

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

! for ghost zone
   imin = gw+1
   jmin = gw+1
   kmin = gw+1

   imax = ex(1) - gw
   jmax = ex(2) - gw
   kmax = ex(3) - gw

!for patch boundary (i.e., not ghost boundary)

if(dabs(X(ex(1))-xmax) < dX)then
  if(X(ex(1))-PIo4 > dX)then
    imax = ex(1)-ogw             ! for overlap zone
  else
    imax = ex(1)
  endif
endif
if(dabs(Y(ex(2))-ymax) < dY)then
  if(Y(ex(2))-PIo4 > dY)then
    jmax = ex(2)-ogw             ! for overlap zone
  else
    jmax = ex(2)
  endif
endif
if(dabs(Z(ex(3))-zmax) < dZ) kmax = ex(3)

if(dabs(X(1)-xmin) < dX)then
  if(X(1)+PIo4 < dX)then
    imin = 1+ogw             ! for overlap zone
  else
    imin = 1
  endif
endif
if(dabs(Y(1)-ymin) < dY)then
  if(Y(1)+PIo4 < dY)then
    jmin = 1+ogw             ! for overlap zone
  else
    jmin = 1
  endif
endif
if(dabs(Z(1)-zmin) < dZ) kmin = 1

!for Symmetry ghost points
if(Symmetry==1)then
  if(dabs(ymin+gw*dY)<dY.and.Y(1)<0.d0) jmin = gw+1
  if(dabs(ymax-gw*dY)<dY.and.Y(ex(2))>0.d0) jmax = ex(2)-gw
endif
if(Symmetry==2)then
  if(dabs(xmin+gw*dX)<dX.and.X(1)<0.d0) imin = gw+1
  if(dabs(ymin+gw*dY)<dY.and.Y(1)<0.d0) jmin = gw+1
endif

f_out = sum(f(imin:imax,jmin:jmax,kmin:kmax)*f(imin:imax,jmin:jmax,kmin:kmax))

f_out = f_out*dX*dY*dZ

  return

  end subroutine l2normhelper_sh
!--------------------------------------------------------------------------------------
! calculate L2norm especially for shell Blocks
! use root mean sqrt method
  subroutine l2normhelper_sh_rms(ex, X, Y, Z,xmin,ymin,zmin,xmax,ymax,zmax,&
                          f,f_out,gw,ogw,Symmetry,Nout)

  implicit none
!~~~~~~> Input parameters:
  integer,intent(in ):: ex(1:3),Symmetry
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3)),xmin,ymin,zmin,xmax,ymax,zmax
  integer,intent(in)::gw,ogw
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8, intent(out) :: f_out
  integer,intent(out) :: Nout
!~~~~~~> Other variables:

  real*8, parameter :: ZEO = 0.D0
  real*8            :: dX, dY, dZ
  integer::imin,jmin,kmin
  integer::imax,jmax,kmax
  integer::i,j,k

  real*8 :: PIo4

  PIo4 = dacos(-1.d0)/4.d0

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

! for ghost zone
   imin = gw+1
   jmin = gw+1
   kmin = gw+1

   imax = ex(1) - gw
   jmax = ex(2) - gw
   kmax = ex(3) - gw

!for patch boundary (i.e., not ghost boundary)

if(dabs(X(ex(1))-xmax) < dX)then
  if(X(ex(1))-PIo4 > dX)then
    imax = ex(1)-ogw             ! for overlap zone
  else
    imax = ex(1)
  endif
endif
if(dabs(Y(ex(2))-ymax) < dY)then
  if(Y(ex(2))-PIo4 > dY)then
    jmax = ex(2)-ogw             ! for overlap zone
  else
    jmax = ex(2)
  endif
endif
if(dabs(Z(ex(3))-zmax) < dZ) kmax = ex(3)

if(dabs(X(1)-xmin) < dX)then
  if(X(1)+PIo4 < dX)then
    imin = 1+ogw             ! for overlap zone
  else
    imin = 1
  endif
endif
if(dabs(Y(1)-ymin) < dY)then
  if(Y(1)+PIo4 < dY)then
    jmin = 1+ogw             ! for overlap zone
  else
    jmin = 1
  endif
endif
if(dabs(Z(1)-zmin) < dZ) kmin = 1

!for Symmetry ghost points
if(Symmetry==1)then
  if(dabs(ymin+gw*dY)<dY.and.Y(1)<0.d0) jmin = gw+1
  if(dabs(ymax-gw*dY)<dY.and.Y(ex(2))>0.d0) jmax = ex(2)-gw
endif
if(Symmetry==2)then
  if(dabs(xmin+gw*dX)<dX.and.X(1)<0.d0) imin = gw+1
  if(dabs(ymin+gw*dY)<dY.and.Y(1)<0.d0) jmin = gw+1
endif

f_out = sum(f(imin:imax,jmin:jmax,kmin:kmax)*f(imin:imax,jmin:jmax,kmin:kmax))

f_out = f_out

Nout = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)

  return

  end subroutine l2normhelper_sh_rms
! locating the position of NaN
  subroutine ScaneNaN(ext,X,Y,Z,f)
  implicit none
  integer,dimension(3),intent(in) :: ext
  real*8,dimension(ext(1)) :: X
  real*8,dimension(ext(2)) :: Y
  real*8,dimension(ext(3)) :: Z
  real*8,dimension(ext(1),ext(2),ext(3)) :: f

  integer :: i,j,k

  do k=1,ext(3)
  do j=1,ext(2)
  do i=1,ext(1)
     if(abs(f(i,j,k)) .ne. abs(f(i,j,k))) write(*,*)X(i),Y(j),Z(k),f(i,j,k)
  enddo
  enddo
  enddo

  end subroutine ScaneNaN
! fortran version writefile
  subroutine fwritefile(time,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,NN,filename,data_out)
  implicit none

  real*8,intent(in) :: time,xmin,xmax,ymin,ymax,zmin,zmax
  integer,intent(in) :: nx,ny,nz,NN
  real*8,dimension(nx,ny,nz),intent(in) :: data_out
  Character(Len=NN) :: filename

 
!  Open( 12 , File = filename,form='BINARY', access="SEQUENTIAL",status="replace",action='Write')
  Open( 12 , File = filename,form='UNFORMATTED', access="DIRECT",status="replace",action='Write')
  Write( 12 ) time,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,data_out

  Close( 12 )

  end subroutine fwritefile
!--------------------------------------------------------------------------
!
! average for interface 
!
!--------------------------------------------------------------------------
  subroutine average(ext,f1,f2,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout

  real*8,parameter::HLF=0.5d0

  fout = HLF*(f1+f2)

  return

  end subroutine average
!-----------------------------------------------------------------------------  
  subroutine average3(ext,f1,f2,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------p  | t
! f2 ----------                |
! 2 points, 1st order interpolation
! 1   2
! f2  f1
! *---*--> t
!    ^
! f=3/4*f_1 + 1/4*f_2

  real*8,parameter::C1=0.75d0,C2=0.25d0

  fout = C1*f1+C2*f2

  return

  end subroutine average3
!-----------------------------------------------------------------------------  
  subroutine average2(ext,f1,f2,f3,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------   |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!       ^
! f=3/8*f_1 + 3/4*f_2 - 1/8*f_3

  real*8,parameter::C1=3.d0/8.d0,C2=3.d0/4.d0,C3=-1.d0/8.d0

  fout = C1*f1+C2*f2+C3*f3

  return

  end subroutine average2
!-----------------------------------------------------------------------------  
  subroutine average2p(ext,f1,f2,f3,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------p  |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!        ^
! f=21/32*f_1 + 7/16*f_2 - 3/32*f_3

  real*8,parameter::C1=5.d0/3.2d1,C2=1.5d1/1.6d1,C3=-3.d0/3.2d1

  fout = C1*f1+C2*f2+C3*f3

  return

  end subroutine average2p
!-----------------------------------------------------------------------------  
  subroutine average2m(ext,f1,f2,f3,fout)
  implicit none
  integer,dimension(3),   intent(in) :: ext
  real*8, dimension(ext(1),ext(2),ext(3)),intent(in):: f1,f2,f3
  real*8, dimension(ext(1),ext(2),ext(3)),intent(out):: fout
! f1 ----------                ^
!                fout ------m  |
! f2 ----------                | t
!                              |
! f3 ----------                |
! 3 points, 2nd order interpolation
! 1   2   3
! f3  f2  f1
! *---*---*--> t
!      ^
! f=5/32*f_1 + 15/16*f_2 - 3/32*f_3

  real*8,parameter::C1=5.d0/3.2d1,C2=1.5d1/1.6d1,C3=-3.d0/3.2d1

  fout = C1*f1+C2*f2+C3*f3

  return

  end subroutine average2m
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine lowerboundset(ex,chi0,TINNY)
  implicit none

!~~~~~~% Input parameters:

  integer ,intent(in):: ex(1:3)
  real*8  ,intent(in):: TINNY
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) ::chi0

  where(chi0 < TINNY) chi0 = TINNY

  return   

  end subroutine lowerboundset
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
  subroutine global_interpind(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,inds,coef,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: inds
  real*8, dimension(3*ORDN),           intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,-ORDN+1:ex(3)+ORDN) :: fh
  integer :: m
  integer,dimension(3) :: cxB,cxT
  real*8, dimension(ORDN,ORDN,ORDN) :: ya
  real*8, dimension(ORDN,ORDN) :: tmp2
  real*8, dimension(ORDN) :: tmp1
  real*8, dimension(3) :: SoAh

! +1 because c++ gives 0 for first point
  cxB = inds+1  
  cxT = cxB + ORDN - 1
      
  if(all(cxB>0).and.all(cxT<ex+1))then
     ya=f(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))
  elseif(any(cxB<-ORDN+1).or.any(cxT>ex+ORDN))then
     write(*,*)"error in global_interpind, cxB = ",cxB
     write(*,*)"                           cxT = ",cxT 
     write(*,*)"                           ext = ",ex
     stop
  else
     if(sst==-1)then
        SoAh = SoA
        if(any(cxT>ex)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==0.or.sst==1)then
        SoAh = SoA
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==2.or.sst==3)then
        SoAh(1) = SoA(2)
        SoAh(2) = SoA(3)
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst
     elseif(sst==4.or.sst==5)then
        SoAh(1) = SoA(1)
        SoAh(2) = SoA(3)
        SoAh(3) = 0
        if(cxB(3)<1.or.cxT(3)>ex(3)) write(*,*)"error global_interpind sst =",sst,cxB(3),cxT(3)
     endif
     call symmetry_tbd(ORDN,ex,f,fh,SoAh)
     ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))
  endif 

  tmp2=0
  do m=1,ORDN
    tmp2 = tmp2 + coef(2*ORDN+m)*ya(:,:,m)
  enddo

  tmp1=0
  do m=1,ORDN
    tmp1 = tmp1 + coef(ORDN+m)*tmp2(:,m)
  enddo

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*tmp1(m)
  enddo

  return

  end subroutine global_interpind
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
! special for shell to shell
  subroutine global_interpind2d(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,inds,coef,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: inds
  real*8, dimension(2*ORDN),           intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,ex(3)) :: fh
  integer :: m
  integer,dimension(2) :: cxB,cxT
  real*8, dimension(ORDN,ORDN) :: ya
  real*8, dimension(ORDN) :: tmp1
  real*8, dimension(2) :: SoAh

! +1 because c++ gives 0 for first point
  cxB = inds(1:2)+1  
  cxT = cxB + ORDN - 1
      
  if(all(cxB>0).and.all(cxT<ex(1:2)+1))then
     ya=f(cxB(1):cxT(1),cxB(2):cxT(2),inds(3))
  elseif(any(cxB<-ORDN+1).or.any(cxT>ex(1:2)+ORDN))then
     write(*,*)"error in global_interpind2d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(1:2)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind2d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(1:2)
     elseif(sst==2.or.sst==3)then
        SoAh(1) = SoA(2)
        SoAh(2) = SoA(3)
     elseif(sst==4.or.sst==5)then
        SoAh(1) = SoA(1)
        SoAh(2) = SoA(3)
     endif
     call symmetry_stbd(ORDN,ex,f,fh,SoAh)
     ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),inds(3))
  endif 

  tmp1=0
  do m=1,ORDN
    tmp1 = tmp1 + coef(ORDN+m)*ya(:,m)
  enddo

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*tmp1(m)
  enddo

  return

  end subroutine global_interpind2d
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!global interpolation with given index and coeffients
! special for shell to shell
! dumyd refer to source
  subroutine global_interpind1d(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,indsi,coef,sst,dumyd)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3),symmetry,ORDN,sst,dumyd
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA
  integer,dimension(3),                intent(in) :: indsi
  real*8, dimension(ORDN),             intent(in) :: coef

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,ex(3)) :: fh
  integer :: m
  integer :: cxB,cxT
  real*8, dimension(ORDN) :: ya
  real*8 :: SoAh
  integer,dimension(3) :: inds

! +1 because c++ gives 0 for first point
  inds = indsi + 1
  cxB = inds(1)  
  cxT = cxB + ORDN - 1
     
! active is rho  
  if(dumyd==1)then

  if(cxB>0.and.cxT<ex(1)+1)then
     ya=f(cxB:cxT,inds(2),inds(3))
  elseif(cxB<-ORDN+1.or.cxT>ex(1)+ORDN)then
     write(*,*)"error in global_interpind1d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(1)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind1d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(1)
     elseif(sst==2.or.sst==3)then
        SoAh = SoA(2)
     elseif(sst==4.or.sst==5)then
        SoAh = SoA(1)
     endif
     call symmetry_sntbd(ORDN,ex,f,fh,SoAh,1-dumyd)
     ya=fh(cxB:cxT,inds(2),inds(3))
  endif 

! active is sigma
  elseif(dumyd==0)then

  if(cxB>0.and.cxT<ex(2)+1)then
     ya=f(inds(2),cxB:cxT,inds(3))
  elseif(cxB<-ORDN+1.or.cxT>ex(2)+ORDN)then
     write(*,*)"error in global_interpind1d, cxB = ",cxB
     write(*,*)"                             cxT = ",cxT 
     write(*,*)"                             ext = ",ex(2)
     stop
  else
     if(sst==-1)then
       write(*,*)"error in global_interpind1d, sst = ",sst
       stop
     elseif(sst==0.or.sst==1)then
        SoAh = SoA(2)
     elseif(sst==2.or.sst==3)then
        SoAh = SoA(3)
     elseif(sst==4.or.sst==5)then
        SoAh = SoA(3)
     endif
     call symmetry_sntbd(ORDN,ex,f,fh,SoAh,1-dumyd)
     ya=fh(inds(2),cxB:cxT,inds(3))
  endif 

  else
          write(*,*)"error in global_interpind1d, not recognized dumyd = ",dumyd
  endif

  f_int=0
  do m=1,ORDN
    f_int = f_int + coef(m)*ya(m)
  enddo

  return

  end subroutine global_interpind1d
!-----------------------------------------------------------------------------------------------------------------
! three dimensional interpolation for both vertex and cell center grid structure  
! for distinguishing shell and Cartesian
  subroutine global_interp_ss(ex,X,Y,Z,f,f_int,x1,y1,z1,ORDN,SoA,symmetry,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3), symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1,z1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,-ORDN+1:ex(3)+ORDN) :: fh
  real*8, dimension(3) :: SoAh
  integer :: j,m,imin,jmin,kmin
  integer,dimension(3) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(3) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,dZ,ddy
  real*8, parameter :: ONE=1.d0

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)
  dZ = Z(kmin+1)-Z(kmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1
  cxI(3) = idint((z1-Z(1))/dZ+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1

  cmin = 1
  cmax = ex

  if(sst==-1)then
     SoAh = SoA
     cmin = -ORDN+1
  elseif(sst==0.or.sst==1)then
     SoAh = SoA
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==2.or.sst==3)then
     SoAh(1) = SoA(2)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==4.or.sst==5)then
     SoAh(1) = SoA(1)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  endif
  do m =1,3
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 cx(1) = (x1 - X(1))/dX-cxB(1)+1
 cx(2) = (y1 - Y(1))/dY-cxB(2)+1
 cx(3) = (z1 - Z(1))/dZ-cxB(3)+1

  call symmetry_tbd(ORDN,ex,f,fh,SoAh)
  ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),cxB(3):cxT(3))

  call polin3(x1a,x1a,x1a,ya,cx(1),cx(2),cx(3),f_int,ddy,ORDN)

  return

  end subroutine global_interp_ss
!-----------------------------------------------------------------------------------------------------------------
! two dimensional interpolation for both vertex and cell center grid structure  
! for distinguishing shell and Cartesian
  subroutine global_interp_ss_2d(ex,X,Y,indZ,f,f_int,x1,y1,ORDN,SoA,symmetry,sst)
  implicit none

!~~~~~~> Input parameters:

  integer,                             intent(in) :: ex(1:3),indZ,symmetry,ORDN,sst
  real*8,intent(in) :: X(ex(1)),Y(ex(2))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in) :: f
  real*8,                              intent(out):: f_int
  real*8,                              intent(in) :: x1,y1
  real*8, dimension(3),                intent(in) :: SoA

!~~~~~~> Other parameters:

  real*8, dimension(-ORDN+1:ex(1)+ORDN,-ORDN+1:ex(2)+ORDN,-ORDN+1:ex(3)+ORDN) :: fh
  real*8, dimension(3) :: SoAh
  integer :: j,m,imin,jmin,kmin
  integer,dimension(2) :: cxB,cxT,cxI,cmin,cmax
  real*8,dimension(2) :: cx
  real*8, dimension(1:ORDN) :: x1a
  real*8, dimension(1:ORDN,1:ORDN) :: ya
  integer, parameter :: NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2
  real*8 :: dX,dY,ddy
  real*8, parameter :: ONE=1.d0

! sanity check
  if(indZ < 1 .or. indZ > ex(3))then
      write(*,*)"error in global_interp_ss_2d, ext = ",ex(3),"ind = ",indZ
      return
  endif

  imin = lbound(f,1)
  jmin = lbound(f,2)
  kmin = lbound(f,3)

  dX = X(imin+1)-X(imin)
  dY = Y(jmin+1)-Y(jmin)

  forall( j = 1:ordn ) x1a(j) = ( j - 1 )* ONE

  cxI(1) = idint((x1-X(1))/dX+0.4)+1
  cxI(2) = idint((y1-Y(1))/dY+0.4)+1

  cxB = cxI - ORDN/2+1
  cxT = cxB + ORDN - 1

  cmin = 1
  cmax = ex(1:2)

  if(sst==-1)then
     SoAh = SoA
     cmin = -ORDN+1
  elseif(sst==0.or.sst==1)then
     SoAh = SoA
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==2.or.sst==3)then
     SoAh(1) = SoA(2)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  elseif(sst==4.or.sst==5)then
     SoAh(1) = SoA(1)
     SoAh(2) = SoA(3)
     SoAh(3) = 0
     cmin(1:2) = -ORDN+1
     cmax(1:2) = ex(1:2)+ORDN
  endif
  do m =1,2
   if(cxB(m) < cmin(m))then
      cxB(m) = cmin(m)
      cxT(m) = cxB(m) + ORDN - 1
   endif
   if(cxT(m) > cmax(m))then
      cxT(m) = cmax(m)
      cxB(m) = cxT(m) + 1 - ORDN
   endif
 enddo
 cx(1) = (x1 - X(1))/dX-cxB(1)+1
 cx(2) = (y1 - Y(1))/dY-cxB(2)+1

  call symmetry_tbd(ORDN,ex,f,fh,SoAh)
  ya=fh(cxB(1):cxT(1),cxB(2):cxT(2),indZ)

  call polin2(x1a,x1a,ya,cx(1),cx(2),f_int,ddy,ORDN)

  return

  end subroutine global_interp_ss_2d
!------------------------------------------
!fortran version of Wigner d function
!Eq.(42) of PRD 77, 024027 (2008)
!we consider only theta in [0,pi]
!------------------------------------------
  function fWigner_d_function(l,m,s,costheta) result(gont)
  implicit none
  integer,intent(in) :: l,m,s
  real*8,intent(in) :: costheta

  real*8 :: gont

  integer :: t,C1,C2
  real*8 :: ffact,vv,sinht,cosht

  C1=max(0,m-s)
  C2=min(l+m,l-s)
  vv=0
  sinht=dsqrt((1.d0-costheta)/2.d0)
  cosht=dsqrt((1.d0+costheta)/2.d0);
  if(C1/2*2==C1)then
    do t=C1,C2,2
       vv=vv+cosht**(2*l+m-s-2*t)*sinht**(2*t+s-m)/(ffact(l+m-t)*ffact(l-s-t)*ffact(t)*ffact(t+s-m))
    enddo
    do t=C1+1,C2,2
       vv=vv-cosht**(2*l+m-s-2*t)*sinht**(2*t+s-m)/(ffact(l+m-t)*ffact(l-s-t)*ffact(t)*ffact(t+s-m))
    enddo
  else
    do t=C1,C2,2
       vv=vv-cosht**(2*l+m-s-2*t)*sinht**(2*t+s-m)/(ffact(l+m-t)*ffact(l-s-t)*ffact(t)*ffact(t+s-m))
    enddo
    do t=C1+1,C2,2
       vv=vv+cosht**(2*l+m-s-2*t)*sinht**(2*t+s-m)/(ffact(l+m-t)*ffact(l-s-t)*ffact(t)*ffact(t+s-m))
    enddo
  endif
  
  gont = vv*dsqrt(ffact(l+m)*ffact(l-m)*ffact(l+s)*ffact(l-s))

  return

  end function fWigner_d_function
!----------------------------------
  function ffact(N) result(gont)
  implicit none
  integer,intent(in) :: N

  real*8 :: gont

  integer :: i

! sanity check
  if(N < 0)then
     write(*,*) "ffact: error input for factorial"
     return
  endif

  gont = 1.d0
  do i=1,N
     gont = gont*i
  enddo

  return

  end function ffact
!---------------------------  
!Eq.(41) of PRD 77, 024027 (2008)
!----------------------------------
  function Yslm(s,l,m,the,phi) result(gont)
  implicit none
  integer,intent(in) :: s,l,m
  real*8,intent(in) :: the,phi

  double complex :: gont

  real*8 :: fWigner_d_function,PI,rp

  PI = dacos(-1.d0)

  rp = fWigner_d_function(l,m,s,dcos(the))
  rp = rp*dsqrt((2*l+1.d0)/4.d0/PI)
  if(s/2*2.ne.s) rp = -rp

  gont = dcmplx(dcos(m*phi),dsin(m*phi))

  gont = rp*gont

  return

  end function Yslm
!------------------------------------------------------------------------------------  
subroutine set_value(ext,data_out,rr)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(out) :: data_out
  REAL*8, intent(in) :: rr

  data_out = rr

  return
end subroutine set_value
subroutine add_value(ext,data_out,rr)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(inout) :: data_out
  REAL*8, intent(in) :: rr

  data_out = data_out + rr

  return
end subroutine add_value
! copy array2 to array1  
subroutine array_copy(ext,data1,data2)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(out) :: data1
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(in) :: data2

  data1 = data2

  return
  end subroutine array_copy
! add array2 to array1  
subroutine array_add(ext,data1,data2)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(inout) :: data1
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(in) :: data2

  data1 = data1 + data2

  return
  end subroutine array_add
! subtract array2 from array1  
subroutine array_subtract(ext,data1,data2)

  IMPLICIT NONE

  integer, intent(in) :: ext(3)
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(inout) :: data1
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(in) :: data2

  data1 = data1 - data2

  return
  end subroutine array_subtract
! find out the maximum
subroutine find_maximum(ext,X,Y,Z,fun,val,pos,llb,uub)

  implicit none

  integer,intent(in) :: ext(3),llb(3),uub(3)
  real*8 :: X(ext(1)),Y(ext(2)),Z(ext(3))
  REAL*8, DIMENSION(ext(1),ext(2),ext(3)), intent(in) :: fun
  real*8,intent(out) :: val,pos(3)

  integer :: i,j,k,ii,jj,kk
  real*8 :: tmp

  tmp = 0.d0

  ii=1
  jj=1
  kk=1

  do k=llb(3)+1,ext(3)-uub(3)
  do j=llb(2)+1,ext(2)-uub(2)
  do i=llb(1)+1,ext(1)-uub(1)
      if(dabs(fun(i,j,k)) > tmp)then
         tmp = dabs(fun(i,j,k))
         ii = i
         jj = j
         kk = k
      endif
  enddo
  enddo
  enddo

  pos(1) = X(ii)
  pos(2) = Y(jj)
  pos(3) = Z(kk)
  val = tmp

  return

end subroutine

#ifndef VERIFY_V7_SYMMETRY_NOZERO
#define VERIFY_V7_SYMMETRY_NOZERO 0
#endif

! V7-13-2: derivative-only preparation. The public symmetry_bd is unchanged.
subroutine symmetry_bd_deriv_fast(ord,extc,func,funcc,SoA)
#if VERIFY_V7_SYMMETRY_NOZERO
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
#endif
  implicit none
  integer,intent(in) :: ord,extc(3)
  real*8,intent(in) :: func(extc(1),extc(2),extc(3)),SoA(3)
  real*8,intent(out) :: funcc(-ord+1:extc(1),-ord+1:extc(2),-ord+1:extc(3))
#if defined(Cell) && !defined(Vertex) && ghost_width == 3
  integer :: i
#if VERIFY_V7_SYMMETRY_NOZERO
  integer :: ix,iy,iz,sx,sy,sz
  integer*8 :: mismatch_count
  real*8 :: old_value,new_value,max_abs,max_rel,err
#endif
  if (ord == 2 .and. all(extc >= 2)) then
#if VERIFY_V7_SYMMETRY_NOZERO
    ! Poison the existing buffer; a missing write cannot inherit a valid zero.
    funcc=ieee_value(0.d0,ieee_quiet_nan)
#endif
    ! Coverage, in order (P=1:n, L=-1:0): PPP, LPP, (L+P)LP,
    ! (L+P)(L+P)L. Every reflection reads already initialized elements.
  funcc(1:extc(1),1:extc(2),1:extc(3)) = func
   do i=0,ord-1
      funcc(-i,1:extc(2),1:extc(3)) = funcc(i+1,1:extc(2),1:extc(3))*SoA(1)
   enddo
   do i=0,ord-1
      funcc(:,-i,1:extc(3)) = funcc(:,i+1,1:extc(3))*SoA(2)
   enddo
   do i=0,ord-1
      funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3)
   enddo


#if VERIFY_V7_SYMMETRY_NOZERO
    ! Stream the old copy/reflection result using scalars, including the
    ! original x then y then z multiplication order. No reference grid.
    mismatch_count=0
    max_abs=0.d0
    max_rel=0.d0
    do iz=-1,extc(3)
    do iy=-1,extc(2)
    do ix=-1,extc(1)
      sx=ix
      sy=iy
      sz=iz
      if(ix<=0) sx=1-ix
      if(iy<=0) sy=1-iy
      if(iz<=0) sz=1-iz
      old_value=func(sx,sy,sz)
      if(ix<=0) old_value=old_value*SoA(1)
      if(iy<=0) old_value=old_value*SoA(2)
      if(iz<=0) old_value=old_value*SoA(3)
      new_value=funcc(ix,iy,iz)
      if(transfer(old_value,0_8)/=transfer(new_value,0_8)) then
        mismatch_count=mismatch_count+1
        if(ieee_is_finite(old_value).and.ieee_is_finite(new_value)) then
          err=abs(new_value-old_value)
          max_abs=max(max_abs,err)
          max_rel=max(max_rel,err/max(abs(old_value),1.d-300))
        else
          max_abs=huge(1.d0)
          max_rel=huge(1.d0)
        endif
      endif
    enddo
    enddo
    enddo
    write(*,'(A,3I6,A,I3,A,3F5.1,A,I18,A,ES24.16,A,ES24.16)') &
         'VERIFY_V7_SYMMETRY_NOZERO ex=',extc,' ord=',ord,' SoA=',SoA, &
         ' mismatch_count=',mismatch_count,' max_abs=',max_abs,' max_rel=',max_rel
#endif
    return
  endif
#endif
  ! Preserve the original behavior for other grids, orders and extents.
  call symmetry_bd(ord,extc,func,funcc,SoA)
end subroutine symmetry_bd_deriv_fast
