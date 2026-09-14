!
#include "macrodef.fh"
#include "v7_rhs_profile.h"

#ifndef VERIFY_V7_SHIFT_FUSION
#define VERIFY_V7_SHIFT_FUSION 0
#endif

#ifndef VERIFY_V7_CHI_RICCI_FUSION
#define VERIFY_V7_CHI_RICCI_FUSION 0
#endif

  function compute_rhs_bssn(ex, T,X, Y, Z,                                     &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! calculate constraint violation when co=0               
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

#if ENABLE_V7_RHS_PROFILING
  real*8 :: v7_rhs_tmark,v7_rhs_ttotal,v14_subtime
  real*8 v7_rhs_wtime
  external v7_rhs_wtime
#endif

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz,div_beta,S
  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxa,Gamya,Gamza,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz
  real*8 :: v7_chi0,v7_inv_chi,v7_inv_2chi,v7_chi_f
  real*8 :: v7_cx,v7_cy,v7_cz,v7_cxx,v7_cxy,v7_cxz,v7_cyy,v7_cyz,v7_czz
  integer :: v7_chi_i,v7_chi_j,v7_chi_k
#if VERIFY_V7_CHI_RICCI_FUSION
  real*8 :: v7_chi_old,v7_chi_new,v7_chi_abs,v7_chi_rel,v7_chi_den,v7_chi_old_f
  real*8 :: v7_chi_rxx0,v7_chi_rxy0,v7_chi_rxz0,v7_chi_ryy0,v7_chi_ryz0,v7_chi_rzz0
  integer :: v7_chi_p,v7_chi_sample
  integer, dimension(4) :: v7_chi_si,v7_chi_sj,v7_chi_sk
  integer, save :: v7_chi_verify_calls=0
#endif

#if VERIFY_V7_SHIFT_FUSION
  real*8, dimension(ex(1),ex(2),ex(3)) :: v7_bxxx,v7_bxxy,v7_bxxz,v7_bxyy,v7_bxyz,v7_bxzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: v7_byxx,v7_byxy,v7_byxz,v7_byyy,v7_byyz,v7_byzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: v7_bzxx,v7_bzxy,v7_bzxz,v7_bzyy,v7_bzyz,v7_bzzz
  real*8 :: v7_old_gradx,v7_old_grady,v7_old_gradz
  real*8 :: v7_old_lapx,v7_old_lapy,v7_old_lapz,v7_den
  integer :: v7_p,v7_si,v7_sj,v7_sk
  integer, dimension(5) :: v7_i,v7_j,v7_k
  integer, save :: v7_verify_calls=0
#endif

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO = 0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0

#if ENABLE_V7_RHS_PROFILING
  v7_rhs_tmark=v7_rhs_wtime()
  v7_rhs_ttotal=v7_rhs_tmark
#endif

#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta
#endif

#if (GAUGE == 6 || GAUGE == 7)
  integer :: BHN,i,j,k
  real*8, dimension(9) :: Porg
  real*8, dimension(3) :: Mass
  real*8 :: r1,r2,M,A,w1,w2,C1,C2
  real*8, dimension(ex(1),ex(2),ex(3)) :: reta

  call getpbh(BHN,Porg,Mass)
#endif

!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
     gont = 1
#if ENABLE_V7_RHS_PROFILING
     call v7_rhs_profile_stop(1,v7_rhs_tmark,v7_rhs_ttotal)
#endif
     return
  endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(1,v7_rhs_tmark)
#endif
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  call fderivs_shift3(ex,betax,betay,betaz, &
       betaxx,betaxy,betaxz,betayx,betayy,betayz, &
       betazx,betazy,betazz,X,Y,Z,Symmetry,Lev)
  
#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(15,v14_subtime)
#endif
  div_beta = betaxx + betayy + betazz
 
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(16,v14_subtime)
#endif
  chi_rhs = F2o3 *chin1*( alpn1 * trK - div_beta ) !rhs for chi

#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(17,v14_subtime)
#endif
  gxx_rhs = - TWO * alpn1 * Axx    -  F2o3 * gxx * div_beta          + &
              TWO *(  gxx * betaxx +   gxy * betayx +   gxz * betazx)

  gyy_rhs = - TWO * alpn1 * Ayy    -  F2o3 * gyy * div_beta          + &
              TWO *(  gxy * betaxy +   gyy * betayy +   gyz * betazy)

  gzz_rhs = - TWO * alpn1 * Azz    -  F2o3 * gzz * div_beta          + &
              TWO *(  gxz * betaxz +   gyz * betayz +   gzz * betazz)

  gxy_rhs = - TWO * alpn1 * Axy    +  F1o3 * gxy    * div_beta       + &
                      gxx * betaxy                  +   gxz * betazy + &
                                       gyy * betayx +   gyz * betazx   &
                                                    -   gxy * betazz

  gyz_rhs = - TWO * alpn1 * Ayz    +  F1o3 * gyz    * div_beta       + &
                      gxy * betaxz +   gyy * betayz                  + &
                      gxz * betaxy                  +   gzz * betazy   &
                                                    -   gyz * betaxx
 
  gxz_rhs = - TWO * alpn1 * Axz    +  F1o3 * gxz    * div_beta       + &
                      gxx * betaxz +   gxy * betayz                  + &
                                       gyz * betayx +   gzz * betazx   &
                                                    -   gxz * betayy     !rhs for gij

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(2,v7_rhs_tmark)
#endif
! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(3,v7_rhs_tmark)
#endif
  if(co == 0)then
! Gam^i_Res = Gam^i + gup^ij_,j
  Gmx_Res = Gamx - (gupxx*(gupxx*gxxx+gupxy*gxyx+gupxz*gxzx)&
                   +gupxy*(gupxx*gxyx+gupxy*gyyx+gupxz*gyzx)&
                   +gupxz*(gupxx*gxzx+gupxy*gyzx+gupxz*gzzx)&
                   +gupxx*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupxy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupxz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxx*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupxy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupxz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmy_Res = Gamy - (gupxx*(gupxy*gxxx+gupyy*gxyx+gupyz*gxzx)&
                   +gupxy*(gupxy*gxyx+gupyy*gyyx+gupyz*gyzx)&
                   +gupxz*(gupxy*gxzx+gupyy*gyzx+gupyz*gzzx)&
                   +gupxy*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupyy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupyz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxy*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupyz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmz_Res = Gamz - (gupxx*(gupxz*gxxx+gupyz*gxyx+gupzz*gxzx)&
                   +gupxy*(gupxz*gxyx+gupyz*gyyx+gupzz*gyzx)&
                   +gupxz*(gupxz*gxzx+gupyz*gyzx+gupzz*gzzx)&
                   +gupxy*(gupxz*gxxy+gupyz*gxyy+gupzz*gxzy)&
                   +gupyy*(gupxz*gxyy+gupyz*gyyy+gupzz*gyzy)&
                   +gupyz*(gupxz*gxzy+gupyz*gyzy+gupzz*gzzy)&
                   +gupxz*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyz*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupzz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  endif

! second kind of connection
  Gamxxx =HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy ) + gupxz*(TWO*gxzx - gxxz ))
  Gamyxx =HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy ) + gupyz*(TWO*gxzx - gxxz ))
  Gamzxx =HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy ) + gupzz*(TWO*gxzx - gxxz ))
 
  Gamxyy =HALF*( gupxx*(TWO*gxyy - gyyx ) + gupxy*gyyy + gupxz*(TWO*gyzy - gyyz ))
  Gamyyy =HALF*( gupxy*(TWO*gxyy - gyyx ) + gupyy*gyyy + gupyz*(TWO*gyzy - gyyz ))
  Gamzyy =HALF*( gupxz*(TWO*gxyy - gyyx ) + gupyz*gyyy + gupzz*(TWO*gyzy - gyyz ))

  Gamxzz =HALF*( gupxx*(TWO*gxzz - gzzx ) + gupxy*(TWO*gyzz - gzzy ) + gupxz*gzzz)
  Gamyzz =HALF*( gupxy*(TWO*gxzz - gzzx ) + gupyy*(TWO*gyzz - gzzy ) + gupyz*gzzz)
  Gamzzz =HALF*( gupxz*(TWO*gxzz - gzzx ) + gupyz*(TWO*gyzz - gzzy ) + gupzz*gzzz)

  Gamxxy =HALF*( gupxx*gxxy + gupxy*gyyx + gupxz*( gxzy + gyzx - gxyz ) )
  Gamyxy =HALF*( gupxy*gxxy + gupyy*gyyx + gupyz*( gxzy + gyzx - gxyz ) )
  Gamzxy =HALF*( gupxz*gxxy + gupyz*gyyx + gupzz*( gxzy + gyzx - gxyz ) )

  Gamxxz =HALF*( gupxx*gxxz + gupxy*( gxyz + gyzx - gxzy ) + gupxz*gzzx )
  Gamyxz =HALF*( gupxy*gxxz + gupyy*( gxyz + gyzx - gxzy ) + gupyz*gzzx )
  Gamzxz =HALF*( gupxz*gxxz + gupyz*( gxyz + gyzx - gxzy ) + gupzz*gzzx )

  Gamxyz =HALF*( gupxx*( gxyz + gxzy - gyzx ) + gupxy*gyyz + gupxz*gzzy )
  Gamyyz =HALF*( gupxy*( gxyz + gxzy - gyzx ) + gupyy*gyyz + gupyz*gzzy )
  Gamzyz =HALF*( gupxz*( gxyz + gxzy - gyzx ) + gupyz*gyyz + gupzz*gzzy )
! Raise indices of \tilde A_{ij} and store in R_ij

  Rxx =    gupxx * gupxx * Axx + gupxy * gupxy * Ayy + gupxz * gupxz * Azz + &
      TWO*(gupxx * gupxy * Axy + gupxx * gupxz * Axz + gupxy * gupxz * Ayz)

  Ryy =    gupxy * gupxy * Axx + gupyy * gupyy * Ayy + gupyz * gupyz * Azz + &
      TWO*(gupxy * gupyy * Axy + gupxy * gupyz * Axz + gupyy * gupyz * Ayz)

  Rzz =    gupxz * gupxz * Axx + gupyz * gupyz * Ayy + gupzz * gupzz * Azz + &
      TWO*(gupxz * gupyz * Axy + gupxz * gupzz * Axz + gupyz * gupzz * Ayz)

  Rxy =    gupxx * gupxy * Axx + gupxy * gupyy * Ayy + gupxz * gupyz * Azz + &
          (gupxx * gupyy       + gupxy * gupxy)* Axy                       + &
          (gupxx * gupyz       + gupxz * gupxy)* Axz                       + &
          (gupxy * gupyz       + gupxz * gupyy)* Ayz

  Rxz =    gupxx * gupxz * Axx + gupxy * gupyz * Ayy + gupxz * gupzz * Azz + &
          (gupxx * gupyz       + gupxy * gupxz)* Axy                       + &
          (gupxx * gupzz       + gupxz * gupxz)* Axz                       + &
          (gupxy * gupzz       + gupxz * gupyz)* Ayz

  Ryz =    gupxy * gupxz * Axx + gupyy * gupyz * Ayy + gupyz * gupzz * Azz + &
          (gupxy * gupyz       + gupyy * gupxz)* Axy                       + &
          (gupxy * gupzz       + gupyz * gupxz)* Axz                       + &
          (gupyy * gupzz       + gupyz * gupyz)* Ayz

! Right hand side for Gam^i without shift terms...
  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

   Gamx_rhs = - TWO * (   Lapx * Rxx +   Lapy * Rxy +   Lapz * Rxz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxx +   chiy * Rxy +   chiz * Rxz ) - &
              gupxx * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupxy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupxz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamxxx * Rxx + Gamxyy * Ryy + Gamxzz * Rzz   + &
                TWO * ( Gamxxy * Rxy + Gamxxz * Rxz + Gamxyz * Ryz ) )

   Gamy_rhs = - TWO * (   Lapx * Rxy +   Lapy * Ryy +   Lapz * Ryz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxy +  chiy * Ryy +    chiz * Ryz ) - &
              gupxy * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupyz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamyxx * Rxx + Gamyyy * Ryy + Gamyzz * Rzz   + &
                TWO * ( Gamyxy * Rxy + Gamyxz * Rxz + Gamyyz * Ryz ) )

   Gamz_rhs = - TWO * (   Lapx * Rxz +   Lapy * Ryz +   Lapz * Rzz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxz +  chiy * Ryz +    chiz * Rzz ) - &
              gupxz * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyz * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupzz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamzxx * Rxx + Gamzyy * Ryy + Gamzzz * Rzz   + &
                TWO * ( Gamzxy * Rxy + Gamzxz * Rxz + Gamzyz * Ryz ) )

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(4,v7_rhs_tmark)
#endif
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  call fdderivs_shift_fusion(ex,betax,betay,betaz, &
       gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, &
       fxx,fxy,fxz,gxxx,gxxy,gxxz,X,Y,Z,Symmetry,Lev)

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(18,v14_subtime)
#endif
#if VERIFY_V7_SHIFT_FUSION
  if(v7_verify_calls == 0)then
    call fdderivs(ex,betax,v7_bxxx,v7_bxxy,v7_bxxz,v7_bxyy,v7_bxyz,v7_bxzz, &
                  X,Y,Z,ANTI,SYM,SYM,Symmetry,Lev)
    call fdderivs(ex,betay,v7_byxx,v7_byxy,v7_byxz,v7_byyy,v7_byyz,v7_byzz, &
                  X,Y,Z,SYM,ANTI,SYM,Symmetry,Lev)
    call fdderivs(ex,betaz,v7_bzxx,v7_bzxy,v7_bzxz,v7_bzyy,v7_bzyz,v7_bzzz, &
                  X,Y,Z,SYM,SYM,ANTI,Symmetry,Lev)

    v7_i=(/1,min(2,ex(1)),min(3,ex(1)),max(1,ex(1)/2),max(1,ex(1)-2)/)
    v7_j=(/1,min(2,ex(2)),min(3,ex(2)),max(1,ex(2)/2),max(1,ex(2)-2)/)
    v7_k=(/1,min(2,ex(3)),min(3,ex(3)),max(1,ex(3)/2),max(1,ex(3)-2)/)

    do v7_p=1,5
      v7_si=v7_i(v7_p)
      v7_sj=v7_j(v7_p)
      v7_sk=v7_k(v7_p)
      v7_old_gradx=v7_bxxx(v7_si,v7_sj,v7_sk)+v7_byxy(v7_si,v7_sj,v7_sk)+v7_bzxz(v7_si,v7_sj,v7_sk)
      v7_old_grady=v7_bxxy(v7_si,v7_sj,v7_sk)+v7_byyy(v7_si,v7_sj,v7_sk)+v7_bzyz(v7_si,v7_sj,v7_sk)
      v7_old_gradz=v7_bxxz(v7_si,v7_sj,v7_sk)+v7_byyz(v7_si,v7_sj,v7_sk)+v7_bzzz(v7_si,v7_sj,v7_sk)
      v7_old_lapx=gupxx(v7_si,v7_sj,v7_sk)*v7_bxxx(v7_si,v7_sj,v7_sk)+ &
           gupyy(v7_si,v7_sj,v7_sk)*v7_bxyy(v7_si,v7_sj,v7_sk)+gupzz(v7_si,v7_sj,v7_sk)*v7_bxzz(v7_si,v7_sj,v7_sk)+ &
           TWO*(gupxy(v7_si,v7_sj,v7_sk)*v7_bxxy(v7_si,v7_sj,v7_sk)+gupxz(v7_si,v7_sj,v7_sk)*v7_bxxz(v7_si,v7_sj,v7_sk)+ &
           gupyz(v7_si,v7_sj,v7_sk)*v7_bxyz(v7_si,v7_sj,v7_sk))
      v7_old_lapy=gupxx(v7_si,v7_sj,v7_sk)*v7_byxx(v7_si,v7_sj,v7_sk)+ &
           gupyy(v7_si,v7_sj,v7_sk)*v7_byyy(v7_si,v7_sj,v7_sk)+gupzz(v7_si,v7_sj,v7_sk)*v7_byzz(v7_si,v7_sj,v7_sk)+ &
           TWO*(gupxy(v7_si,v7_sj,v7_sk)*v7_byxy(v7_si,v7_sj,v7_sk)+gupxz(v7_si,v7_sj,v7_sk)*v7_byxz(v7_si,v7_sj,v7_sk)+ &
           gupyz(v7_si,v7_sj,v7_sk)*v7_byyz(v7_si,v7_sj,v7_sk))
      v7_old_lapz=gupxx(v7_si,v7_sj,v7_sk)*v7_bzxx(v7_si,v7_sj,v7_sk)+ &
           gupyy(v7_si,v7_sj,v7_sk)*v7_bzyy(v7_si,v7_sj,v7_sk)+gupzz(v7_si,v7_sj,v7_sk)*v7_bzzz(v7_si,v7_sj,v7_sk)+ &
           TWO*(gupxy(v7_si,v7_sj,v7_sk)*v7_bzxy(v7_si,v7_sj,v7_sk)+gupxz(v7_si,v7_sj,v7_sk)*v7_bzxz(v7_si,v7_sj,v7_sk)+ &
           gupyz(v7_si,v7_sj,v7_sk)*v7_bzyz(v7_si,v7_sj,v7_sk))

      v7_den=max(dabs(v7_old_gradx),1.d-300)
      write(*,'(A,3I6,A,4ES25.16)') 'V7_SHIFT gradx ',v7_si,v7_sj,v7_sk,' old/new/abs/rel ', &
           v7_old_gradx,fxx(v7_si,v7_sj,v7_sk),dabs(v7_old_gradx-fxx(v7_si,v7_sj,v7_sk)), &
           dabs(v7_old_gradx-fxx(v7_si,v7_sj,v7_sk))/v7_den
      v7_den=max(dabs(v7_old_grady),1.d-300)
      write(*,'(A,3I6,A,4ES25.16)') 'V7_SHIFT grady ',v7_si,v7_sj,v7_sk,' old/new/abs/rel ', &
           v7_old_grady,fxy(v7_si,v7_sj,v7_sk),dabs(v7_old_grady-fxy(v7_si,v7_sj,v7_sk)), &
           dabs(v7_old_grady-fxy(v7_si,v7_sj,v7_sk))/v7_den
      v7_den=max(dabs(v7_old_gradz),1.d-300)
      write(*,'(A,3I6,A,4ES25.16)') 'V7_SHIFT gradz ',v7_si,v7_sj,v7_sk,' old/new/abs/rel ', &
           v7_old_gradz,fxz(v7_si,v7_sj,v7_sk),dabs(v7_old_gradz-fxz(v7_si,v7_sj,v7_sk)), &
           dabs(v7_old_gradz-fxz(v7_si,v7_sj,v7_sk))/v7_den
      v7_den=max(dabs(v7_old_lapx),1.d-300)
      write(*,'(A,3I6,A,4ES25.16)') 'V7_SHIFT lapx  ',v7_si,v7_sj,v7_sk,' old/new/abs/rel ', &
           v7_old_lapx,gxxx(v7_si,v7_sj,v7_sk),dabs(v7_old_lapx-gxxx(v7_si,v7_sj,v7_sk)), &
           dabs(v7_old_lapx-gxxx(v7_si,v7_sj,v7_sk))/v7_den
      v7_den=max(dabs(v7_old_lapy),1.d-300)
      write(*,'(A,3I6,A,4ES25.16)') 'V7_SHIFT lapy  ',v7_si,v7_sj,v7_sk,' old/new/abs/rel ', &
           v7_old_lapy,gxxy(v7_si,v7_sj,v7_sk),dabs(v7_old_lapy-gxxy(v7_si,v7_sj,v7_sk)), &
           dabs(v7_old_lapy-gxxy(v7_si,v7_sj,v7_sk))/v7_den
      v7_den=max(dabs(v7_old_lapz),1.d-300)
      write(*,'(A,3I6,A,4ES25.16)') 'V7_SHIFT lapz  ',v7_si,v7_sj,v7_sk,' old/new/abs/rel ', &
           v7_old_lapz,gxxz(v7_si,v7_sj,v7_sk),dabs(v7_old_lapz-gxxz(v7_si,v7_sj,v7_sk)), &
           dabs(v7_old_lapz-gxxz(v7_si,v7_sj,v7_sk))/v7_den
    enddo
    v7_verify_calls=v7_verify_calls+1
  endif
#endif

#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  Gamxa =       gupxx * Gamxxx + gupyy * Gamxyy + gupzz * Gamxzz + &
          TWO*( gupxy * Gamxxy + gupxz * Gamxxz + gupyz * Gamxyz )
  Gamya =       gupxx * Gamyxx + gupyy * Gamyyy + gupzz * Gamyzz + &
          TWO*( gupxy * Gamyxy + gupxz * Gamyxz + gupyz * Gamyyz )
  Gamza =       gupxx * Gamzxx + gupyy * Gamzyy + gupzz * Gamzzz + &
          TWO*( gupxy * Gamzxy + gupxz * Gamzxz + gupyz * Gamzyz )

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(19,v14_subtime)
#endif
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(20,v14_subtime)
#endif
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  Gamx_rhs =               Gamx_rhs +  F2o3 *  Gamxa * div_beta        - &
                     Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz  + &
             F1o3 * (gupxx * fxx    + gupxy * fxy    + gupxz * fxz    ) + &
                     gxxx

  Gamy_rhs =               Gamy_rhs +  F2o3 *  Gamya * div_beta        - &
                     Gamxa * betayx - Gamya * betayy - Gamza * betayz  + &
             F1o3 * (gupxy * fxx    + gupyy * fxy    + gupyz * fxz    ) + &
                     gxxy

  Gamz_rhs =               Gamz_rhs +  F2o3 *  Gamza * div_beta        - &
                     Gamxa * betazx - Gamya * betazy - Gamza * betazz  + &
             F1o3 * (gupxz * fxx    + gupyz * fxy    + gupzz * fxz    ) + &
                     gxxz                                                   !rhs for Gam^i

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(21,v14_subtime)
#endif
#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(5,v7_rhs_tmark)
#endif
!first kind of connection stored in gij,k
  gxxx = gxx * Gamxxx + gxy * Gamyxx + gxz * Gamzxx
  gxyx = gxx * Gamxxy + gxy * Gamyxy + gxz * Gamzxy
  gxzx = gxx * Gamxxz + gxy * Gamyxz + gxz * Gamzxz
  gyyx = gxx * Gamxyy + gxy * Gamyyy + gxz * Gamzyy
  gyzx = gxx * Gamxyz + gxy * Gamyyz + gxz * Gamzyz
  gzzx = gxx * Gamxzz + gxy * Gamyzz + gxz * Gamzzz

  gxxy = gxy * Gamxxx + gyy * Gamyxx + gyz * Gamzxx
  gxyy = gxy * Gamxxy + gyy * Gamyxy + gyz * Gamzxy
  gxzy = gxy * Gamxxz + gyy * Gamyxz + gyz * Gamzxz
  gyyy = gxy * Gamxyy + gyy * Gamyyy + gyz * Gamzyy
  gyzy = gxy * Gamxyz + gyy * Gamyyz + gyz * Gamzyz
  gzzy = gxy * Gamxzz + gyy * Gamyzz + gyz * Gamzzz

  gxxz = gxz * Gamxxx + gyz * Gamyxx + gzz * Gamzxx
  gxyz = gxz * Gamxxy + gyz * Gamyxy + gzz * Gamzxy
  gxzz = gxz * Gamxxz + gyz * Gamyxz + gzz * Gamzxz
  gyyz = gxz * Gamxyy + gyz * Gamyyy + gzz * Gamzyy
  gyzz = gxz * Gamxyz + gyz * Gamyyz + gzz * Gamzyz
  gzzz = gxz * Gamxzz + gyz * Gamyzz + gzz * Gamzzz

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(6,v7_rhs_tmark)
#endif
!compute Ricci tensor for tilted metric
   call fdderivs(ex,dxx,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rxx =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,dyy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Ryy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,dzz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rzz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gxy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI, ANTI,SYM ,symmetry,Lev)
   Rxy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gxz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI ,SYM ,ANTI,symmetry,Lev)
   Rxz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gyz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,ANTI ,ANTI,symmetry,Lev)
   Ryz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO
  Rxx =     - HALF * Rxx                                   + &
               gxx * Gamxx+ gxy * Gamyx   +    gxz * Gamzx + &
             Gamxa * gxxx +  Gamya * gxyx +  Gamza * gxzx  + &
   gupxx *(                                                  &
       TWO*(Gamxxx * gxxx + Gamyxx * gxyx + Gamzxx * gxzx) + &
            Gamxxx * gxxx + Gamyxx * gxxy + Gamzxx * gxxz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxx * gxyx + Gamyxx * gyyx + Gamzxx * gyzx  + &
            Gamxxy * gxxx + Gamyxy * gxyx + Gamzxy * gxzx) + &
            Gamxxy * gxxx + Gamyxy * gxxy + Gamzxy * gxxz  + &
            Gamxxx * gxyx + Gamyxx * gxyy + Gamzxx * gxyz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxx * gxzx + Gamyxx * gyzx + Gamzxx * gzzx  + &
            Gamxxz * gxxx + Gamyxz * gxyx + Gamzxz * gxzx) + &
            Gamxxz * gxxx + Gamyxz * gxxy + Gamzxz * gxxz  + &
            Gamxxx * gxzx + Gamyxx * gxzy + Gamzxx * gxzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxxy * gxyx + Gamyxy * gyyx + Gamzxy * gyzx) + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz )+ &
   gupyz *(                                                  &
       TWO*(Gamxxy * gxzx + Gamyxy * gyzx + Gamzxy * gzzx  + &
            Gamxxz * gxyx + Gamyxz * gyyx + Gamzxz * gyzx) + &
            Gamxxz * gxyx + Gamyxz * gxyy + Gamzxz * gxyz  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxxz * gxzx + Gamyxz * gyzx + Gamzxz * gzzx) + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz )

  Ryy =     - HALF * Ryy                                   + &
               gxy * Gamxy+  gyy * Gamyy  +  gyz * Gamzy   + &
             Gamxa * gxyy +  Gamya * gyyy +  Gamza * gyzy  + &
   gupxx *(                                                  &
       TWO*(Gamxxy * gxxy + Gamyxy * gxyy + Gamzxy * gxzy) + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxy * gxyy + Gamyxy * gyyy + Gamzxy * gyzy  + &
            Gamxyy * gxxy + Gamyyy * gxyy + Gamzyy * gxzy) + &
            Gamxyy * gxyx + Gamyyy * gxyy + Gamzyy * gxyz  + &
            Gamxxy * gyyx + Gamyxy * gyyy + Gamzxy * gyyz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxy * gxzy + Gamyxy * gyzy + Gamzxy * gzzy  + &
            Gamxyz * gxxy + Gamyyz * gxyy + Gamzyz * gxzy) + &
            Gamxyz * gxyx + Gamyyz * gxyy + Gamzyz * gxyz  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxyy * gxyy + Gamyyy * gyyy + Gamzyy * gyzy) + &
            Gamxyy * gyyx + Gamyyy * gyyy + Gamzyy * gyyz )+ &
   gupyz *(                                                  &
       TWO*(Gamxyy * gxzy + Gamyyy * gyzy + Gamzyy * gzzy  + &
            Gamxyz * gxyy + Gamyyz * gyyy + Gamzyz * gyzy) + &
            Gamxyz * gyyx + Gamyyz * gyyy + Gamzyz * gyyz  + &
            Gamxyy * gyzx + Gamyyy * gyzy + Gamzyy * gyzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxyz * gxzy + Gamyyz * gyzy + Gamzyz * gzzy) + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz )

  Rzz =     - HALF * Rzz                                   + &
               gxz * Gamxz+ gyz * Gamyz  +    gzz * Gamzz  + &
             Gamxa * gxzz +  Gamya * gyzz +  Gamza * gzzz  + &
   gupxx *(                                                  &
       TWO*(Gamxxz * gxxz + Gamyxz * gxyz + Gamzxz * gxzz) + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxz * gxyz + Gamyxz * gyyz + Gamzxz * gyzz  + &
            Gamxyz * gxxz + Gamyyz * gxyz + Gamzyz * gxzz) + &
            Gamxyz * gxzx + Gamyyz * gxzy + Gamzyz * gxzz  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxz * gxzz + Gamyxz * gyzz + Gamzxz * gzzz  + &
            Gamxzz * gxxz + Gamyzz * gxyz + Gamzzz * gxzz) + &
            Gamxzz * gxzx + Gamyzz * gxzy + Gamzzz * gxzz  + &
            Gamxxz * gzzx + Gamyxz * gzzy + Gamzxz * gzzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxyz * gxyz + Gamyyz * gyyz + Gamzyz * gyzz) + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz )+ &
   gupyz *(                                                  &
       TWO*(Gamxyz * gxzz + Gamyyz * gyzz + Gamzyz * gzzz  + &
            Gamxzz * gxyz + Gamyzz * gyyz + Gamzzz * gyzz) + &
            Gamxzz * gyzx + Gamyzz * gyzy + Gamzzz * gyzz  + &
            Gamxyz * gzzx + Gamyyz * gzzy + Gamzyz * gzzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxzz * gxzz + Gamyzz * gyzz + Gamzzz * gzzz) + &
            Gamxzz * gzzx + Gamyzz * gzzy + Gamzzz * gzzz )

  Rxy = HALF*(     - Rxy                                   + &
               gxx * Gamxy +    gxy * Gamyy + gxz * Gamzy  + &
               gxy * Gamxx +    gyy * Gamyx + gyz * Gamzx  + &
             Gamxa * gxyx +  Gamya * gyyx +  Gamza * gyzx  + &
             Gamxa * gxxy +  Gamya * gxyy +  Gamza * gxzy )+ &
   gupxx *(                                                  &
            Gamxxx * gxxy + Gamyxx * gxyy + Gamzxx * gxzy  + &
            Gamxxy * gxxx + Gamyxy * gxyx + Gamzxy * gxzx  + &
            Gamxxx * gxyx + Gamyxx * gxyy + Gamzxx * gxyz )+ &
   gupxy *(                                                  &
            Gamxxx * gxyy + Gamyxx * gyyy + Gamzxx * gyzy  + &
            Gamxxy * gxyx + Gamyxy * gyyx + Gamzxy * gyzx  + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz  + &
            Gamxxy * gxxy + Gamyxy * gxyy + Gamzxy * gxzy  + &
            Gamxyy * gxxx + Gamyyy * gxyx + Gamzyy * gxzx  + &
            Gamxxx * gyyx + Gamyxx * gyyy + Gamzxx * gyyz )+ &
   gupxz *(                                                  &
            Gamxxx * gxzy + Gamyxx * gyzy + Gamzxx * gzzy  + &
            Gamxxy * gxzx + Gamyxy * gyzx + Gamzxy * gzzx  + &
            Gamxxz * gxyx + Gamyxz * gxyy + Gamzxz * gxyz  + &
            Gamxxz * gxxy + Gamyxz * gxyy + Gamzxz * gxzy  + &
            Gamxyz * gxxx + Gamyyz * gxyx + Gamzyz * gxzx  + &
            Gamxxx * gyzx + Gamyxx * gyzy + Gamzxx * gyzz )+ &
   gupyy *(                                                  &
            Gamxxy * gxyy + Gamyxy * gyyy + Gamzxy * gyzy  + &
            Gamxyy * gxyx + Gamyyy * gyyx + Gamzyy * gyzx  + &
            Gamxxy * gyyx + Gamyxy * gyyy + Gamzxy * gyyz )+ &
   gupyz *(                                                  &
            Gamxxy * gxzy + Gamyxy * gyzy + Gamzxy * gzzy  + &
            Gamxyy * gxzx + Gamyyy * gyzx + Gamzyy * gzzx  + &
            Gamxxz * gyyx + Gamyxz * gyyy + Gamzxz * gyyz  + &
            Gamxxz * gxyy + Gamyxz * gyyy + Gamzxz * gyzy  + &
            Gamxyz * gxyx + Gamyyz * gyyx + Gamzyz * gyzx  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupzz *(                                                  &
            Gamxxz * gxzy + Gamyxz * gyzy + Gamzxz * gzzy  + &
            Gamxyz * gxzx + Gamyyz * gyzx + Gamzyz * gzzx  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz )

  Rxz = HALF*(     - Rxz                                   + &
               gxx * Gamxz +  gxy * Gamyz + gxz * Gamzz    + &
               gxz * Gamxx +  gyz * Gamyx + gzz * Gamzx    + &
             Gamxa * gxzx +  Gamya * gyzx +  Gamza * gzzx  + &
             Gamxa * gxxz +  Gamya * gxyz +  Gamza * gxzz )+ &
   gupxx *(                                                  &
            Gamxxx * gxxz + Gamyxx * gxyz + Gamzxx * gxzz  + &
            Gamxxz * gxxx + Gamyxz * gxyx + Gamzxz * gxzx  + &
            Gamxxx * gxzx + Gamyxx * gxzy + Gamzxx * gxzz )+ &
   gupxy *(                                                  &
            Gamxxx * gxyz + Gamyxx * gyyz + Gamzxx * gyzz  + &
            Gamxxz * gxyx + Gamyxz * gyyx + Gamzxz * gyzx  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz  + &
            Gamxxy * gxxz + Gamyxy * gxyz + Gamzxy * gxzz  + &
            Gamxyz * gxxx + Gamyyz * gxyx + Gamzyz * gxzx  + &
            Gamxxx * gyzx + Gamyxx * gyzy + Gamzxx * gyzz )+ &
   gupxz *(                                                  &
            Gamxxx * gxzz + Gamyxx * gyzz + Gamzxx * gzzz  + &
            Gamxxz * gxzx + Gamyxz * gyzx + Gamzxz * gzzx  + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz  + &
            Gamxxz * gxxz + Gamyxz * gxyz + Gamzxz * gxzz  + &
            Gamxzz * gxxx + Gamyzz * gxyx + Gamzzz * gxzx  + &
            Gamxxx * gzzx + Gamyxx * gzzy + Gamzxx * gzzz )+ &
   gupyy *(                                                  &
            Gamxxy * gxyz + Gamyxy * gyyz + Gamzxy * gyzz  + &
            Gamxyz * gxyx + Gamyyz * gyyx + Gamzyz * gyzx  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupyz *(                                                  &
            Gamxxy * gxzz + Gamyxy * gyzz + Gamzxy * gzzz  + &
            Gamxyz * gxzx + Gamyyz * gyzx + Gamzyz * gzzx  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz  + &
            Gamxxz * gxyz + Gamyxz * gyyz + Gamzxz * gyzz  + &
            Gamxzz * gxyx + Gamyzz * gyyx + Gamzzz * gyzx  + &
            Gamxxy * gzzx + Gamyxy * gzzy + Gamzxy * gzzz )+ &
   gupzz *(                                                  &
            Gamxxz * gxzz + Gamyxz * gyzz + Gamzxz * gzzz  + &
            Gamxzz * gxzx + Gamyzz * gyzx + Gamzzz * gzzx  + &
            Gamxxz * gzzx + Gamyxz * gzzy + Gamzxz * gzzz )

  Ryz = HALF*(     - Ryz                                   + &
               gxy * Gamxz + gyy * Gamyz + gyz * Gamzz     + &
               gxz * Gamxy + gyz * Gamyy + gzz * Gamzy     + &
             Gamxa * gxzy +  Gamya * gyzy +  Gamza * gzzy  + &
             Gamxa * gxyz +  Gamya * gyyz +  Gamza * gyzz )+ &
   gupxx *(                                                  &
            Gamxxy * gxxz + Gamyxy * gxyz + Gamzxy * gxzz  + &
            Gamxxz * gxxy + Gamyxz * gxyy + Gamzxz * gxzy  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz )+ &
   gupxy *(                                                  &
            Gamxxy * gxyz + Gamyxy * gyyz + Gamzxy * gyzz  + &
            Gamxxz * gxyy + Gamyxz * gyyy + Gamzxz * gyzy  + &
            Gamxyy * gxzx + Gamyyy * gxzy + Gamzyy * gxzz  + &
            Gamxyy * gxxz + Gamyyy * gxyz + Gamzyy * gxzz  + &
            Gamxyz * gxxy + Gamyyz * gxyy + Gamzyz * gxzy  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupxz *(                                                  &
            Gamxxy * gxzz + Gamyxy * gyzz + Gamzxy * gzzz  + &
            Gamxxz * gxzy + Gamyxz * gyzy + Gamzxz * gzzy  + &
            Gamxyz * gxzx + Gamyyz * gxzy + Gamzyz * gxzz  + &
            Gamxyz * gxxz + Gamyyz * gxyz + Gamzyz * gxzz  + &
            Gamxzz * gxxy + Gamyzz * gxyy + Gamzzz * gxzy  + &
            Gamxxy * gzzx + Gamyxy * gzzy + Gamzxy * gzzz )+ &
   gupyy *(                                                  &
            Gamxyy * gxyz + Gamyyy * gyyz + Gamzyy * gyzz  + &
            Gamxyz * gxyy + Gamyyz * gyyy + Gamzyz * gyzy  + &
            Gamxyy * gyzx + Gamyyy * gyzy + Gamzyy * gyzz )+ &
   gupyz *(                                                  &
            Gamxyy * gxzz + Gamyyy * gyzz + Gamzyy * gzzz  + &
            Gamxyz * gxzy + Gamyyz * gyzy + Gamzyz * gzzy  + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz  + &
            Gamxyz * gxyz + Gamyyz * gyyz + Gamzyz * gyzz  + &
            Gamxzz * gxyy + Gamyzz * gyyy + Gamzzz * gyzy  + &
            Gamxyy * gzzx + Gamyyy * gzzy + Gamzyy * gzzz )+ &
   gupzz *(                                                  &
            Gamxyz * gxzz + Gamyyz * gyzz + Gamzyz * gzzz  + &
            Gamxzz * gxzy + Gamyzz * gyzy + Gamzzz * gzzy  + &
            Gamxyz * gzzx + Gamyyz * gzzy + Gamzyz * gzzz )

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(7,v7_rhs_tmark)
#endif
!covariant second derivative of chi respect to tilted metric
  call fdderivs(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

  fxx = fxx - Gamxxx * chix - Gamyxx * chiy - Gamzxx * chiz
  fxy = fxy - Gamxxy * chix - Gamyxy * chiy - Gamzxy * chiz
  fxz = fxz - Gamxxz * chix - Gamyxz * chiy - Gamzxz * chiz
  fyy = fyy - Gamxyy * chix - Gamyyy * chiy - Gamzyy * chiz
  fyz = fyz - Gamxyz * chix - Gamyyz * chiy - Gamzyz * chiz
  fzz = fzz - Gamxzz * chix - Gamyzz * chiy - Gamzzz * chiz
! Fuse the scalar chi contraction and all six Ricci corrections.
#if VERIFY_V7_CHI_RICCI_FUSION
  if(v7_chi_verify_calls == 0)then
    v7_chi_si=(/1,min(2,ex(1)),max(1,ex(1)/2),max(1,ex(1)-1)/)
    v7_chi_sj=(/1,min(2,ex(2)),max(1,ex(2)/2),max(1,ex(2)-1)/)
    v7_chi_sk=(/1,min(2,ex(3)),max(1,ex(3)/2),max(1,ex(3)-1)/)
  endif
#endif
  do v7_chi_k=1,ex(3)
  do v7_chi_j=1,ex(2)
  do v7_chi_i=1,ex(1)
    v7_chi0=chin1(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_inv_chi=ONE/v7_chi0
    v7_inv_2chi=HALF*v7_inv_chi
    v7_cx=chix(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_cy=chiy(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_cz=chiz(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_cxx=fxx(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_cxy=fxy(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_cxz=fxz(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_cyy=fyy(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_cyz=fyz(v7_chi_i,v7_chi_j,v7_chi_k)
    v7_czz=fzz(v7_chi_i,v7_chi_j,v7_chi_k)

    v7_chi_f= &
         gupxx(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cxx-F3o2*v7_inv_chi*v7_cx*v7_cx)+ &
         gupyy(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cyy-F3o2*v7_inv_chi*v7_cy*v7_cy)+ &
         gupzz(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_czz-F3o2*v7_inv_chi*v7_cz*v7_cz)+ &
         TWO*gupxy(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cxy-F3o2*v7_inv_chi*v7_cx*v7_cy)+ &
         TWO*gupxz(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cxz-F3o2*v7_inv_chi*v7_cx*v7_cz)+ &
         TWO*gupyz(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cyz-F3o2*v7_inv_chi*v7_cy*v7_cz)

#if VERIFY_V7_CHI_RICCI_FUSION
    v7_chi_sample=0
    if(v7_chi_verify_calls == 0)then
      do v7_chi_p=1,4
        if(v7_chi_i == v7_chi_si(v7_chi_p) .and. v7_chi_j == v7_chi_sj(v7_chi_p) .and. &
           v7_chi_k == v7_chi_sk(v7_chi_p)) v7_chi_sample=1
      enddo
    endif
    if(v7_chi_sample == 1)then
      v7_chi_rxx0=Rxx(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_rxy0=Rxy(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_rxz0=Rxz(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_ryy0=Ryy(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_ryz0=Ryz(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_rzz0=Rzz(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_old_f= &
           gupxx(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cxx-F3o2/v7_chi0*v7_cx*v7_cx)+ &
           gupyy(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cyy-F3o2/v7_chi0*v7_cy*v7_cy)+ &
           gupzz(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_czz-F3o2/v7_chi0*v7_cz*v7_cz)+ &
           TWO*gupxy(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cxy-F3o2/v7_chi0*v7_cx*v7_cy)+ &
           TWO*gupxz(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cxz-F3o2/v7_chi0*v7_cx*v7_cz)+ &
           TWO*gupyz(v7_chi_i,v7_chi_j,v7_chi_k)*(v7_cyz-F3o2/v7_chi0*v7_cy*v7_cz)
    endif
#endif

    Rxx(v7_chi_i,v7_chi_j,v7_chi_k)=Rxx(v7_chi_i,v7_chi_j,v7_chi_k)+ &
         (v7_cxx-v7_cx*v7_cx*v7_inv_2chi+gxx(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_f)*v7_inv_2chi
    Ryy(v7_chi_i,v7_chi_j,v7_chi_k)=Ryy(v7_chi_i,v7_chi_j,v7_chi_k)+ &
         (v7_cyy-v7_cy*v7_cy*v7_inv_2chi+gyy(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_f)*v7_inv_2chi
    Rzz(v7_chi_i,v7_chi_j,v7_chi_k)=Rzz(v7_chi_i,v7_chi_j,v7_chi_k)+ &
         (v7_czz-v7_cz*v7_cz*v7_inv_2chi+gzz(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_f)*v7_inv_2chi
    Rxy(v7_chi_i,v7_chi_j,v7_chi_k)=Rxy(v7_chi_i,v7_chi_j,v7_chi_k)+ &
         (v7_cxy-v7_cx*v7_cy*v7_inv_2chi+gxy(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_f)*v7_inv_2chi
    Rxz(v7_chi_i,v7_chi_j,v7_chi_k)=Rxz(v7_chi_i,v7_chi_j,v7_chi_k)+ &
         (v7_cxz-v7_cx*v7_cz*v7_inv_2chi+gxz(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_f)*v7_inv_2chi
    Ryz(v7_chi_i,v7_chi_j,v7_chi_k)=Ryz(v7_chi_i,v7_chi_j,v7_chi_k)+ &
         (v7_cyz-v7_cy*v7_cz*v7_inv_2chi+gyz(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_f)*v7_inv_2chi

#if VERIFY_V7_CHI_RICCI_FUSION
    if(v7_chi_sample == 1)then
      v7_chi_old=v7_chi_rxx0+(v7_cxx-v7_cx*v7_cx/v7_chi0/TWO+ &
           gxx(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_old_f)/v7_chi0/TWO
      v7_chi_new=Rxx(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_abs=dabs(v7_chi_old-v7_chi_new); v7_chi_den=max(dabs(v7_chi_old),1.d-300); v7_chi_rel=v7_chi_abs/v7_chi_den
      write(*,'(A,3I6,A,4ES25.16)') 'V7_CHI Rxx ',v7_chi_i,v7_chi_j,v7_chi_k, &
           ' old/new/abs/rel ',v7_chi_old,v7_chi_new,v7_chi_abs,v7_chi_rel
      v7_chi_old=v7_chi_rxy0+(v7_cxy-v7_cx*v7_cy/v7_chi0/TWO+ &
           gxy(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_old_f)/v7_chi0/TWO
      v7_chi_new=Rxy(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_abs=dabs(v7_chi_old-v7_chi_new); v7_chi_den=max(dabs(v7_chi_old),1.d-300); v7_chi_rel=v7_chi_abs/v7_chi_den
      write(*,'(A,3I6,A,4ES25.16)') 'V7_CHI Rxy ',v7_chi_i,v7_chi_j,v7_chi_k, &
           ' old/new/abs/rel ',v7_chi_old,v7_chi_new,v7_chi_abs,v7_chi_rel
      v7_chi_old=v7_chi_rxz0+(v7_cxz-v7_cx*v7_cz/v7_chi0/TWO+ &
           gxz(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_old_f)/v7_chi0/TWO
      v7_chi_new=Rxz(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_abs=dabs(v7_chi_old-v7_chi_new); v7_chi_den=max(dabs(v7_chi_old),1.d-300); v7_chi_rel=v7_chi_abs/v7_chi_den
      write(*,'(A,3I6,A,4ES25.16)') 'V7_CHI Rxz ',v7_chi_i,v7_chi_j,v7_chi_k, &
           ' old/new/abs/rel ',v7_chi_old,v7_chi_new,v7_chi_abs,v7_chi_rel
      v7_chi_old=v7_chi_ryy0+(v7_cyy-v7_cy*v7_cy/v7_chi0/TWO+ &
           gyy(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_old_f)/v7_chi0/TWO
      v7_chi_new=Ryy(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_abs=dabs(v7_chi_old-v7_chi_new); v7_chi_den=max(dabs(v7_chi_old),1.d-300); v7_chi_rel=v7_chi_abs/v7_chi_den
      write(*,'(A,3I6,A,4ES25.16)') 'V7_CHI Ryy ',v7_chi_i,v7_chi_j,v7_chi_k, &
           ' old/new/abs/rel ',v7_chi_old,v7_chi_new,v7_chi_abs,v7_chi_rel
      v7_chi_old=v7_chi_ryz0+(v7_cyz-v7_cy*v7_cz/v7_chi0/TWO+ &
           gyz(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_old_f)/v7_chi0/TWO
      v7_chi_new=Ryz(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_abs=dabs(v7_chi_old-v7_chi_new); v7_chi_den=max(dabs(v7_chi_old),1.d-300); v7_chi_rel=v7_chi_abs/v7_chi_den
      write(*,'(A,3I6,A,4ES25.16)') 'V7_CHI Ryz ',v7_chi_i,v7_chi_j,v7_chi_k, &
           ' old/new/abs/rel ',v7_chi_old,v7_chi_new,v7_chi_abs,v7_chi_rel
      v7_chi_old=v7_chi_rzz0+(v7_czz-v7_cz*v7_cz/v7_chi0/TWO+ &
           gzz(v7_chi_i,v7_chi_j,v7_chi_k)*v7_chi_old_f)/v7_chi0/TWO
      v7_chi_new=Rzz(v7_chi_i,v7_chi_j,v7_chi_k)
      v7_chi_abs=dabs(v7_chi_old-v7_chi_new); v7_chi_den=max(dabs(v7_chi_old),1.d-300); v7_chi_rel=v7_chi_abs/v7_chi_den
      write(*,'(A,3I6,A,4ES25.16)') 'V7_CHI Rzz ',v7_chi_i,v7_chi_j,v7_chi_k, &
           ' old/new/abs/rel ',v7_chi_old,v7_chi_new,v7_chi_abs,v7_chi_rel
    endif
#endif
  enddo
  enddo
  enddo
#if VERIFY_V7_CHI_RICCI_FUSION
  if(v7_chi_verify_calls == 0) v7_chi_verify_calls=v7_chi_verify_calls+1
#endif
#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(8,v7_rhs_tmark)
#endif
! covariant second derivatives of the lapse respect to physical metric
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  call fdderivs(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                SYM,SYM,SYM,symmetry,Lev)

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(22,v14_subtime)
#endif
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  gxxx = (gupxx * chix + gupxy * chiy + gupxz * chiz)/chin1
  gxxy = (gupxy * chix + gupyy * chiy + gupyz * chiz)/chin1
  gxxz = (gupxz * chix + gupyz * chiy + gupzz * chiz)/chin1
! now get physical second kind of connection
  Gamxxx = Gamxxx - ( (chix + chix)/chin1 - gxx * gxxx )*HALF
  Gamyxx = Gamyxx - (                     - gxx * gxxy )*HALF
  Gamzxx = Gamzxx - (                     - gxx * gxxz )*HALF
  Gamxyy = Gamxyy - (                     - gyy * gxxx )*HALF
  Gamyyy = Gamyyy - ( (chiy + chiy)/chin1 - gyy * gxxy )*HALF
  Gamzyy = Gamzyy - (                     - gyy * gxxz )*HALF
  Gamxzz = Gamxzz - (                     - gzz * gxxx )*HALF
  Gamyzz = Gamyzz - (                     - gzz * gxxy )*HALF
  Gamzzz = Gamzzz - ( (chiz + chiz)/chin1 - gzz * gxxz )*HALF
  Gamxxy = Gamxxy - (  chiy        /chin1 - gxy * gxxx )*HALF
  Gamyxy = Gamyxy - (         chix /chin1 - gxy * gxxy )*HALF
  Gamzxy = Gamzxy - (                     - gxy * gxxz )*HALF
  Gamxxz = Gamxxz - (  chiz        /chin1 - gxz * gxxx )*HALF
  Gamyxz = Gamyxz - (                     - gxz * gxxy )*HALF
  Gamzxz = Gamzxz - (         chix /chin1 - gxz * gxxz )*HALF
  Gamxyz = Gamxyz - (                     - gyz * gxxx )*HALF
  Gamyyz = Gamyyz - (  chiz        /chin1 - gyz * gxxy )*HALF
  Gamzyz = Gamzyz - (         chiy /chin1 - gyz * gxxz )*HALF

  fxx = fxx - Gamxxx*Lapx - Gamyxx*Lapy - Gamzxx*Lapz
  fyy = fyy - Gamxyy*Lapx - Gamyyy*Lapy - Gamzyy*Lapz
  fzz = fzz - Gamxzz*Lapx - Gamyzz*Lapy - Gamzzz*Lapz
  fxy = fxy - Gamxxy*Lapx - Gamyxy*Lapy - Gamzxy*Lapz
  fxz = fxz - Gamxxz*Lapx - Gamyxz*Lapy - Gamzxz*Lapz
  fyz = fyz - Gamxyz*Lapx - Gamyyz*Lapy - Gamzyz*Lapz

! store D^i D_i Lap in trK_rhs upto chi
  trK_rhs =    gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
#if 1        
!! follow bam code
  S =  chin1 * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )
  f = F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) )) -1.6d1*PI*rho + EIGHT * PI * S
  f = - F1o3 *(  gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + alpn1/chin1*f)
  
  fxx = alpn1 * (Rxx - EIGHT * PI * Sxx) - fxx
  fxy = alpn1 * (Rxy - EIGHT * PI * Sxy) - fxy
  fxz = alpn1 * (Rxz - EIGHT * PI * Sxz) - fxz
  fyy = alpn1 * (Ryy - EIGHT * PI * Syy) - fyy
  fyz = alpn1 * (Ryz - EIGHT * PI * Syz) - fyz
  fzz = alpn1 * (Rzz - EIGHT * PI * Szz) - fzz
#else        
! Add lapse and S_ij parts to Ricci tensor:

  fxx = alpn1 * (Rxx - EIGHT * PI * Sxx) - fxx
  fxy = alpn1 * (Rxy - EIGHT * PI * Sxy) - fxy
  fxz = alpn1 * (Rxz - EIGHT * PI * Sxz) - fxz
  fyy = alpn1 * (Ryy - EIGHT * PI * Syy) - fyy
  fyz = alpn1 * (Ryz - EIGHT * PI * Syz) - fyz
  fzz = alpn1 * (Rzz - EIGHT * PI * Szz) - fzz

! Compute trace-free part (note: chi^-1 and chi cancel!):

  f = F1o3 *(  gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) )
#endif

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(23,v14_subtime)
#endif
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
  Axx_rhs = fxx - gxx * f
  Ayy_rhs = fyy - gyy * f
  Azz_rhs = fzz - gzz * f
  Axy_rhs = fxy - gxy * f
  Axz_rhs = fxz - gxz * f
  Ayz_rhs = fyz - gyz * f

! Now: store A_il A^l_j into fij:

  fxx =       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz)
  fyy =       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz)
  fzz =       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz)
  fxy =       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
              gupxy *(Axx * Ayy + Axy * Axy)                            + &
              gupxz *(Axx * Ayz + Axz * Axy)                            + &
              gupyz *(Axy * Ayz + Axz * Ayy)
  fxz =       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
              gupxy *(Axx * Ayz + Axy * Axz)                            + &
              gupxz *(Axx * Azz + Axz * Axz)                            + &
              gupyz *(Axy * Azz + Axz * Ayz)
  fyz =       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
              gupxy *(Axy * Ayz + Ayy * Axz)                            + &
              gupxz *(Axy * Azz + Ayz * Axz)                            + &
              gupyz *(Ayy * Azz + Ayz * Ayz)

  f = chin1
! store D^i D_i Lap in trK_rhs
  trK_rhs = f*trK_rhs

  Axx_rhs =           f * Axx_rhs+ alpn1 * (trK * Axx - TWO * fxx)  + &
           TWO * (  Axx * betaxx +   Axy * betayx +   Axz * betazx )- &
             F2o3 * Axx * div_beta

  Ayy_rhs =           f * Ayy_rhs+ alpn1 * (trK * Ayy - TWO * fyy)  + &
           TWO * (  Axy * betaxy +   Ayy * betayy +   Ayz * betazy )- &
             F2o3 * Ayy * div_beta

  Azz_rhs =           f * Azz_rhs+ alpn1 * (trK * Azz - TWO * fzz)  + &
           TWO * (  Axz * betaxz +   Ayz * betayz +   Azz * betazz )- &
             F2o3 * Azz * div_beta

  Axy_rhs =           f * Axy_rhs+ alpn1 *( trK * Axy  - TWO * fxy )+ &
                    Axx * betaxy                  +   Axz * betazy  + &
                                     Ayy * betayx +   Ayz * betazx  + &
             F1o3 * Axy * div_beta                -   Axy * betazz

  Ayz_rhs =           f * Ayz_rhs+ alpn1 *( trK * Ayz  - TWO * fyz )+ &
                    Axy * betaxz +   Ayy * betayz                   + &
                    Axz * betaxy                  +   Azz * betazy  + &
             F1o3 * Ayz * div_beta                -   Ayz * betaxx
 
  Axz_rhs =           f * Axz_rhs+ alpn1 *( trK * Axz  - TWO * fxz )+ &
                    Axx * betaxz +   Axy * betayz                   + &
                                     Ayz * betayx +   Azz * betazx  + &
             F1o3 * Axz * div_beta                -   Axz * betayy      !rhs for Aij

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(24,v14_subtime)
#endif
#if ENABLE_V7_RHS_PROFILING
  v14_subtime=v7_rhs_wtime()
#endif
! Compute trace of S_ij

  S =  f * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )

  trK_rhs = - trK_rhs + alpn1 *( F1o3 * trK * trK         + &
                gupxx * fxx + gupyy * fyy + gupzz * fzz   + &
        TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + &
       FOUR * PI * ( rho + S ))                                !rhs for trK
  
#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(25,v14_subtime)
#endif
!!!! gauge variable part

  Lap_rhs = -TWO*alpn1*trK
#if (GAUGE == 0)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - eta*dtSfx
  dtSfy_rhs = Gamy_rhs - eta*dtSfy
  dtSfz_rhs = Gamz_rhs - eta*dtSfz
#elif (GAUGE == 1)
  betax_rhs = Gamx - eta*betax
  betay_rhs = Gamy - eta*betay
  betaz_rhs = Gamz - eta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 2)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 3)
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 4)
  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-dsqrt(chin1))**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 5)
  call fderivs(ex,chi,dtSfx_rhs,dtSfy_rhs,dtSfz_rhs,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  reta = gupxx * dtSfx_rhs * dtSfx_rhs + gupyy * dtSfy_rhs * dtSfy_rhs + gupzz * dtSfz_rhs * dtSfz_rhs + &
       TWO * (gupxy * dtSfx_rhs * dtSfy_rhs + gupxz * dtSfx_rhs * dtSfz_rhs + gupyz * dtSfy_rhs * dtSfz_rhs)
  reta = 1.31d0/2*dsqrt(reta/chin1)/(1-chin1)**2
  betax_rhs = FF*Gamx - reta*betax
  betay_rhs = FF*Gamy - reta*betay
  betaz_rhs = FF*Gamz - reta*betaz

  dtSfx_rhs = ZEO
  dtSfy_rhs = ZEO
  dtSfz_rhs = ZEO
#elif (GAUGE == 6)
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1/(ONE+w1*r1) + C2/(ONE+w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 1",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#elif (GAUGE == 7)
  if(BHN==2)then
   M = Mass(1)+Mass(2)
   A = 2.d0/M
   w1 = 1.2d1
   w2 = w1
   C1 = 1.d0/Mass(1) - A
   C2 = 1.d0/Mass(2) - A

   do k=1,ex(3)
   do j=1,ex(2)
   do i=1,ex(1)
     r1 = ((Porg(1)-X(i))**2+(Porg(2)-Y(j))**2+(Porg(3)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     r2 = ((Porg(4)-X(i))**2+(Porg(5)-Y(j))**2+(Porg(6)-Z(k))**2)/ &
          ((Porg(1)-Porg(4))**2+(Porg(2)-Porg(5))**2+(Porg(3)-Porg(6))**2)
     reta(i,j,k) = A + C1*dexp(-w1*r1) + C2*dexp(-w2*r2)
    enddo
    enddo
    enddo
  else
    write(*,*) "not support BH_num in Jason's form 2",BHN
  endif
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - reta*dtSfx
  dtSfy_rhs = Gamy_rhs - reta*dtSfy
  dtSfz_rhs = Gamz_rhs - reta*dtSfz
#endif  

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI

!!!!!!!!!advection term part

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(9,v7_rhs_tmark)
#endif
  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)
!!
  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)

#if (GAUGE == 0 || GAUGE == 1 || GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)
  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(10,v7_rhs_tmark)
#endif
  if(eps>0)then
! usual Kreiss-Oliger dissipation
  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
#if 0
#define i 42
#define j 40
#define k 40
if(Lev == 1)then
write(*,*) X(i),Y(j),Z(k)
write(*,*) "before",Axx_rhs(i,j,k)
endif
#undef i
#undef j
#undef k
!!stop
#endif
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
#if 0
#define i 42
#define j 40
#define k 40
if(Lev == 1)then
write(*,*) X(i),Y(j),Z(k)
write(*,*) "after",Axx_rhs(i,j,k)
endif
#undef i
#undef j
#undef k
!!stop
#endif
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)

#if 1 
!! bam does not apply dissipation on gauge variables
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
#endif
#endif

  endif

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(11,v7_rhs_tmark)
#endif
  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = chin1*ham_Res + F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) ))- F16 * PI * rho

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  gxxx = gxxx - (  Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz &
                 + Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz) - chix*Axx/chin1
  gxyx = gxyx - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxx * Axy + Gamyxx * Ayy + Gamzxx * Ayz) - chix*Axy/chin1
  gxzx = gxzx - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxx * Axz + Gamyxx * Ayz + Gamzxx * Azz) - chix*Axz/chin1
  gyyx = gyyx - (  Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chix*Ayy/chin1
  gyzx = gyzx - (  Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chix*Ayz/chin1
  gzzx = gzzx - (  Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chix*Azz/chin1
  gxxy = gxxy - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz) - chiy*Axx/chin1
  gxyy = gxyy - (  Gamxyy * Axx + Gamyyy * Axy + Gamzyy * Axz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chiy*Axy/chin1
  gxzy = gxzy - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chiy*Axz/chin1
  gyyy = gyyy - (  Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz &
                 + Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz) - chiy*Ayy/chin1
  gyzy = gyzy - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyy * Axz + Gamyyy * Ayz + Gamzyy * Azz) - chiy*Ayz/chin1
  gzzy = gzzy - (  Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiy*Azz/chin1
  gxxz = gxxz - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz) - chiz*Axx/chin1
  gxyz = gxyz - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz) - chiz*Axy/chin1
  gxzz = gxzz - (  Gamxzz * Axx + Gamyzz * Axy + Gamzzz * Axz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chiz*Axz/chin1
  gyyz = gyyz - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz) - chiz*Ayy/chin1
  gyzz = gyzz - (  Gamxzz * Axy + Gamyzz * Ayy + Gamzzz * Ayz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiz*Ayz/chin1
  gzzz = gzzz - (  Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz &
                 + Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz) - chiz*Azz/chin1
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

movx_Res = movx_Res - F2o3*Kx - F8*PI*sx
movy_Res = movy_Res - F2o3*Ky - F8*PI*sy
movz_Res = movz_Res - F2o3*Kz - F8*PI*sz
  endif

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_mark(12,v7_rhs_tmark)
#endif
#if (ABV == 1)
  call ricci_gamma(ex, X, Y, Z,                                      &
               chi,                                                  &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,&
               Gamx   ,  Gamy    ,  Gamz    , &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,&
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,&
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,&
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,&
               Symmetry)
  call constraint_bssn(ex, X, Y, Z,&
               chi,trK, &
               dxx,gxy,gxz,dyy,gyz,dzz, &
               Axx,Axy,Axz,Ayy,Ayz,Azz, &
               Gamx,Gamy,Gamz,&
               Lap,betax,betay,betaz,rho,Sx,Sy,Sz,&
               Gamxxx, Gamxxy, Gamxxz,Gamxyy, Gamxyz, Gamxzz, &
               Gamyxx, Gamyxy, Gamyxz,Gamyyy, Gamyyz, Gamyzz, &
               Gamzxx, Gamzxy, Gamzxz,Gamzyy, Gamzyz, Gamzzz, &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz, &
               ham_Res,movx_Res,movy_Res,movz_Res,Gmx_Res,Gmy_Res,Gmz_Res, &
               Symmetry)
#endif 
#if 0
#define i 2
if(Lev == 1)then
write(*,*) X(i),Y(i),Z(i)
write(*,*) Axx(i,i,i),Axy(i,i,i),Axz(i,i,i),Ayy(i,i,i),Ayz(i,i,i),Azz(i,i,i)
write(*,*) 1+Lap(i,i,i),dtSfx(i,i,i),dtSfy(i,i,i),dtSfz(i,i,i)
write(*,*) betax(i,i,i),betay(i,i,i),betaz(i,i,i)
write(*,*) 1+chi(i,i,i),Gamx(i,i,i),Gamy(i,i,i),Gamz(i,i,i)
write(*,*) gxx(i,i,i),gxy(i,i,i),gxz(i,i,i),gyy(i,i,i),gyz(i,i,i),gzz(i,i,i)
write(*,*) trK(i,i,i)
write(*,*) "====="
write(*,*) Axx_rhs(i,i,i),Axy_rhs(i,i,i),Axz_rhs(i,i,i),Ayy_rhs(i,i,i),Ayz_rhs(i,i,i),Azz_rhs(i,i,i)
write(*,*) Lap_rhs(i,i,i),dtSfx_rhs(i,i,i),dtSfy_rhs(i,i,i),dtSfz_rhs(i,i,i)
write(*,*) betax_rhs(i,i,i),betay_rhs(i,i,i),betaz_rhs(i,i,i)
write(*,*) chi_rhs(i,i,i),Gamx_rhs(i,i,i),Gamy_rhs(i,i,i),Gamz_rhs(i,i,i)
write(*,*) gxx_rhs(i,i,i),gxy_rhs(i,i,i),gxz_rhs(i,i,i),gyy_rhs(i,i,i),gyz_rhs(i,i,i),gzz_rhs(i,i,i)
write(*,*) trK_rhs(i,i,i)
endif
#undef i
!!stop
#endif

  gont = 0

#if ENABLE_V7_RHS_PROFILING
  call v7_rhs_profile_stop(13,v7_rhs_tmark,v7_rhs_ttotal)
#endif

  return
  end function compute_rhs_bssn

#if ENABLE_V7_RHS_PROFILING
  subroutine v7_rhs_profile_mark(phase,tmark)
  implicit none
  integer, intent(in) :: phase
  real*8, intent(inout) :: tmark
  real*8 :: now,v7_rhs_times(27)
  integer*8 :: v7_rhs_counts(27)
  real*8 v7_rhs_wtime
  external v7_rhs_wtime
  common /v7_rhs_profile_data/ v7_rhs_times,v7_rhs_counts
  now=v7_rhs_wtime()
  v7_rhs_times(phase)=v7_rhs_times(phase)+now-tmark
  v7_rhs_counts(phase)=v7_rhs_counts(phase)+1
  tmark=now
  end subroutine v7_rhs_profile_mark

  subroutine v7_rhs_profile_stop(last_phase,tmark,ttotal)
  implicit none
  integer, intent(in) :: last_phase
  real*8, intent(in) :: tmark,ttotal
  real*8 :: now,v7_rhs_times(27)
  integer*8 :: v7_rhs_counts(27)
  real*8 v7_rhs_wtime
  external v7_rhs_wtime
  common /v7_rhs_profile_data/ v7_rhs_times,v7_rhs_counts
  now=v7_rhs_wtime()
  v7_rhs_times(last_phase)=v7_rhs_times(last_phase)+now-tmark
  v7_rhs_counts(last_phase)=v7_rhs_counts(last_phase)+1
  v7_rhs_times(14)=v7_rhs_times(14)+now-ttotal
  v7_rhs_counts(14)=v7_rhs_counts(14)+1
  end subroutine v7_rhs_profile_stop

  block data v7_rhs_profile_init
  implicit none
  real*8 :: v7_rhs_times(27)
  integer*8 :: v7_rhs_counts(27)
  common /v7_rhs_profile_data/ v7_rhs_times,v7_rhs_counts
  data v7_rhs_times /27*0.d0/
  data v7_rhs_counts /27*0/
  end block data v7_rhs_profile_init
#endif
