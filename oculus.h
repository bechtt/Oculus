#ifdef HAVE_CONFIGF
#include "config.f"
#else
#define HAVE_NAG
#endif

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

module oculus
  
  implicit none
  
! constants; 25 Mar 15;
    
  REAL, parameter :: zero          =          0.0
  REAL, parameter :: one           =          1.0
  REAL, parameter :: two           =          2.0
  REAL, parameter :: three         =          3.0
  REAL, parameter :: four          =          4.0
  REAL, parameter :: five          =          5.0
  REAL, parameter :: six           =          6.0
  REAL, parameter :: seven         =          7.0
  REAL, parameter :: eight         =          8.0
  REAL, parameter :: nine          =          9.0
  REAL, parameter :: ten           =         10.0
  REAL, parameter :: twelve        =         12.0
  REAL, parameter :: hundred       =        100.0
  REAL, parameter :: thousand      =       1000.0

  REAL, parameter :: pi            =          3.141592653589793238462643383279502884197
  REAL, parameter :: pi2           =          pi * two

  REAL, parameter :: goldenmean    =          1.618033988749895 ! golden mean = ( one + sqrt(five) ) / two ;

  REAL, parameter :: small         =       1.0e-13 ! this should really be machine precision; 25 Mar 15;

  REAL, parameter :: half          =       one / two
  REAL, parameter :: third         =       one / three
  REAL, parameter :: quart         =       one / four
  REAL, parameter :: fifth         =       one / five
  REAL, parameter :: sixth         =       one / six

! structures; 25 Mar 15;
  
  type biotsavart
     LOGICAL              :: LB=.false., LA=.false., LL=.false.
     INTEGER              :: N, ixyz, ifail
     REAL                 :: tol, x, y, z, Bx, By, Bz, Ax, Ay, Az, length
  end type biotsavart
  
  type(biotsavart)        :: bsfield
  !$OMP threadprivate(bsfield)
  
  type magneticaxis
     INTEGER              :: Nfp, Ntor, maxits, its, Lallocated, nbfield(0:1)
     REAL                 :: R, Z, odetol, tol, error, residue, iota, rzf(0:2,0:31)
     REAL                 :: tangent(1:2,1:2), wr(1:2), wi(1:2), vr(1:2,1:2), vi(1:2,1:2)
     REAL   , allocatable :: Ri(:), Zi(:), Rnc(:), Rns(:), Znc(:), Zns(:)
  end type magneticaxis

  type(magneticaxis)      :: axis
  !$OMP threadprivate(axis)

  type homoclinictangle 
     INTEGER              :: Nfp, xits, maxits, hits, ilobe(1:2), Lallocated, maxilobe, nbfield(0:1)
     REAL                 :: odetol, xtol, htol, ltol, R, Z, dU, dS, xerror, herror, lerror(1:2), residue
     REAL                 :: tangent(1:2,1:2), wr(1:2), wi(1:2), vr(1:2,1:2), vi(1:2,1:2)
     REAL, allocatable    :: hpoints(:,:,:)
  end type homoclinictangle

  type(homoclinictangle)  :: tangle
  !$OMP threadprivate(tangle)

  type extremizingcurve 
     INTEGER              :: Nfp, Ntor, emethod, Lrz, maxits, its, nbfield(0:1), Lallocated, itau, ibfield
     REAL                 :: etol, ftol, err, sN, odetol, tauend, dtau, FF, epsilon
     REAL, allocatable    :: Rnc(:), Rns(:), Znc(:), Zns(:)
     REAL, allocatable    :: iRnc(:,:), iRns(:,:), iZnc(:,:), iZns(:,:) ! to be deleted; 28 Aug 15;
     REAL, allocatable    :: tt(:), ct(:,:), st(:,:), trig(:), rr(:), rd(:), zz(:), zd(:)
     CHARACTER            :: init
  end type extremizingcurve

  type(extremizingcurve)  :: curve
  !$OMP threadprivate(curve)

#ifdef FH00AA
  type backgroundcoords
     INTEGER              :: Nfp, pplobe, nsep(1:2), Nx, Ny, Nh, Lallocated, nbfield(0:1)
     INTEGER, allocatable :: ii(:), jj(:)
     REAL                 :: odetol, svdtol, sigma
     REAL                 :: wr(1:2), vr(1:2,1:2)
     REAL, allocatable    :: sep(:,:,:), hi(:)
  end type backgroundcoords  
  type(backgroundcoords)  :: bcoords
  !$OMP threadprivate(bcoords)
#endif

  type coordinates 
     INTEGER              :: Lrad, Nfp, Mpol, Ntor, mn, mm
     INTEGER, allocatable :: im(:), in(:)
     REAL   , allocatable :: ss(:), Rbc(:,:), Rbs(:,:), Zbc(:,:), Zbs(:,:)
     REAL   , allocatable :: Xbc(:,:,:), Xbs(:,:,:), Ybc(:,:,:), Ybs(:,:,:)
  end type coordinates

  type(coordinates)       :: rzmn
  !$OMP threadprivate(rzmn)

  type vectorpotential
     INTEGER              :: Nfp, Lrad, Mpol, Ntor, mn, Nt, Nz, Ntz
     INTEGER, allocatable :: im(:), in(:)
     REAL   , allocatable :: ss(:), gBtc(:,:,:), gBts(:,:,:), gBzc(:,:,:), gBzs(:,:,:), gBstz(:,:,:)
  end type vectorpotential

  type(vectorpotential)   :: Atzmn
  !$OMP threadprivate(Atzmn)

  type irrationalsurface
     INTEGER              :: Nfp, Mpol, Ntor, Mits, Nxdof, Npts, Nc, Nk
     INTEGER, allocatable :: im(:), in(:)
     REAL                 :: odetol, irrtol
     REAL   , allocatable :: rcm(:), rsm(:), tcm(:), tsm(:), sap(:,:), tap(:,:), xdof(:), err(:), derr(:,:)
  end type irrationalsurface

  type(irrationalsurface) :: isurface
  !$OMP threadprivate(isurface)

  type qfminsurface
     INTEGER              :: Nfp, pp, qq, nn, mm, Ntor, Mpol, Np, qNp, Nd, qNd, Lrestart, Id, mn, ok
     INTEGER, allocatable :: im(:), in(:)
     REAL                 :: rr, nu, odetol, offset, pqtol
     REAL   , allocatable :: t(:,:), r(:,:), n(:), X(:,:), Y(:,:), Rbc(:), Rbs(:), Zbc(:), Zbs(:)
  end type qfminsurface

  type(qfminsurface)      :: qfms
  !$OMP threadprivate(qfms)

  type pressurerelax
     LOGICAL, allocatable :: F(:,:,:)
     INTEGER              :: Nfp, NR, Np, NZ, ifb, ii, Ntime, Ldiff, Lode, Lcontrol
     REAL                 :: Dp, dtime, kpara, kperp, Rmin, Rmax, Zmin, Zmax, odetol, lB(-1:1)
     INTEGER, allocatable :: I(:,:,:,:), J(:,:,:,:)
     REAL   , allocatable :: s(:,:,:), p(:,:,:), B(:,:,:,:), x(:,:,:,:), y(:,:,:,:)
  end type pressurerelax

  type(pressurerelax)     :: pressure
  !$OMP threadprivate(pressure)

  type transformdata
     INTEGER              :: Nfp, Ppts, Lallocated, nbfield(0:1)
     REAL                 :: Ra, Za, R, Z, iota, odetol
     REAL   , allocatable :: RZ(:,:)
  end type transformdata

  type(transformdata)     :: transform
  !$OMP threadprivate(transform)

  type bnormal
     INTEGER              :: cmn, Nfp, Mpol, Ntor, mn, Nt, Nz, Ntz
     INTEGER, allocatable :: cim(:), cin(:), im(:), in(:)
     REAL                 :: tol, Itor, Gpol
     REAL   , allocatable :: Rcc(:), Rcs(:), Zcc(:), Zcs(:), gBi(:), gBc(:), gBs(:)
  end type bnormal

  type(bnormal)           :: bn
  !$OMP threadprivate(bn)

  type rzcoordsdata
     INTEGER              :: Nfp, Ppts, Lallocated
     REAL                 :: R, Z, odetol
     REAL   , allocatable :: RZ(:,:)
  end type rzcoordsdata

  type(rzcoordsdata)      :: rzdata
  !$OMP threadprivate(rzdata)

  type poincaredata
     INTEGER              :: Nfp = 1, Ppts = 100, idirection = 1, iLyapunov = 0, flparameter = 0
     INTEGER              :: Lallocated, ipts, nbfield(0:1)
     REAL                 :: phi = zero, R, Z, odetol = 1.0e-08, Ltol = 1.0e-08, Lyapunov = zero, phistart
     REAL   , allocatable :: RZ(:,:), Ly(:)
  end type poincaredata

  type(poincaredata)      :: poincare
  !$OMP threadprivate(poincare)

! miscellaneous internal varibles; 25 Mar 15;

  ! internal copy of "ifail"'s
  INTEGER                 :: ibs00aa
  INTEGER                 :: iga00aa, iec00aa, iho00aa, itr00aa, ipp00aa, ibc00aa, iqf00aa, iad00aa, iaa00aa, iaa00ba
  INTEGER                 :: ibn00aa
#ifdef FH00AA
  INTEGER                 :: ifh00aa
  !$OMP threadprivate(ifh00aa)
#endif
  INTEGER                 :: iir00aa, irz00aa
  
  INTEGER                 :: itangent  ! input required for  user-supplied subroutine bfield                       ; 19 Jun 15;
  LOGICAL                 :: Lbfieldok ! logical flag: check user-supplied subroutine bfield was correctly executed; 19 Jun 15;
  INTEGER                 :: nbfield(0:1) ! counts how many times the subroutine bfield is called; 19 Jun 15;

  REAL                    :: actiongradient
  INTEGER                 :: izeta, iteta, fNtor, ipip, niter(0:1)
  REAL                    :: lxx(1:2), lff(1:2) ! user-termination of C05PBF; 19 Jun 15;

  !$OMP threadprivate(ibs00aa, iga00aa, iec00aa, iho00aa, itr00aa, ipp00aa, ibc00aa, iqf00aa, iad00aa, iaa00aa, iaa00ba)
  !$OMP threadprivate(ibn00aa, iir00aa, irz00aa)
  !$OMP threadprivate(itangent, Lbfieldok, nbfield, actiongradient, izeta, iteta, fNtor, ipip, niter, lxx, lff)


! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

contains

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
#ifdef HAVE_NAG
  subroutine bf00aa( zeta, RZ, BRZ ) ! format constrained by the NAG ode integration routines.  
    implicit none
    INTEGER, parameter  :: Node = 6
#else    
  subroutine bf00aa( Node, zeta, RZ, BRZ ) ! format constrained by the lsode integration routines. 
    implicit none
    INTEGER, intent(in) :: Node
#endif
        
    REAL  , intent(in)  :: zeta, RZ(1:Node)
    REAL  , intent(out) :: BRZ(1:Node)
    
    INTEGER             :: ifail
    REAL                :: RpZ(1:3), dBRpZ(1:3,0:3), TM(1:2,1:2)
    
    BRZ(1:Node) = zero ! default intent out; 25 Mar 15;
    
    if( .not.Lbfieldok ) return
    
    RpZ(1:3) = (/ RZ(1), zeta, RZ(2) /)

    nbfield(itangent) = nbfield(itangent) + 1 ; ifail = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), ifail )
    
    if( ifail.ne.0 ) then ; Lbfieldok = .false. ; return
    endif
    
    if( abs(dBRpZ(2,0)).lt.small ) then
     write(0,'("bf00aa : ! bfield ! : ifail.eq.0, but |B^z| < small ; divide-by-zero error ;")')
     write(0,'("bf00aa : ! bfield ! : B("es23.15","es23.15","es23.15")=("es23.15","es23.15","es23.15")")') &
           RpZ(1:3), dBRpZ(1:3,0)
     Lbfieldok = .false. ; return
    endif

    BRZ(1:2) = (/ dBRpZ(1,0), dBRpZ(3,0) /) /  dBRpZ(2,0) ! normalize to toroidal field; 25 Mar 15;
    
    if( itangent.eq.0 ) return

    TM(1,1) = ( dBRpZ(1,1) - BRZ(1) * dBRpZ(2,1) ) / dBRpZ(2,0)
    TM(1,2) = ( dBRpZ(1,3) - BRZ(1) * dBRpZ(2,3) ) / dBRpZ(2,0)
    TM(2,1) = ( dBRpZ(3,1) - BRZ(2) * dBRpZ(2,1) ) / dBRpZ(2,0)
    TM(2,2) = ( dBRpZ(3,3) - BRZ(2) * dBRpZ(2,3) ) / dBRpZ(2,0)
     
    BRZ(3) = TM(1,1) * RZ(3) + TM(1,2) * RZ(5) ! tangent map obtained by matrix multiplication;
    BRZ(4) = TM(1,1) * RZ(4) + TM(1,2) * RZ(6)
    BRZ(5) = TM(2,1) * RZ(3) + TM(2,2) * RZ(5)
    BRZ(6) = TM(2,1) * RZ(4) + TM(2,2) * RZ(6)
     
9999 continue
    
    return
    
  end subroutine bf00aa

#ifdef HAVE_NAG
  subroutine bf00ab( zeta, RZ, BRZ ) ! format constrained by the NAG ode integration routines.  
    implicit none
    INTEGER, parameter  :: Node = 5
#else    
  subroutine bf00ab( Node, zeta, RZ, BRZ ) ! format constrained by the lsode integration routines. 
    implicit none
    INTEGER, intent(in) :: Node
#endif
    
    REAL  , intent(in)  :: zeta, RZ(1:Node)
    REAL  , intent(out) :: BRZ(1:Node)
    
    INTEGER             :: ifail, jfail
    REAL                :: RpZa(1:3), RpZb(1:3), dBRpZa(1:3,0:3), dBRpZb(1:3,0:3), dR, dZ, length
    
    BRZ(1:Node) = zero ! default intent out; 25 Mar 15;
    
    if( .not.Lbfieldok ) return

    RpZa(1:3) = (/ RZ(1), zeta, RZ(2) /)
    RpZb(1:3) = (/ RZ(3), zeta, RZ(4) /)
    
    nbfield(itangent) = nbfield(itangent) + 1 ; ifail = -9 ; call bfield( RpZa(1:3), itangent, dBRpZa(1:3,0:3), ifail )
    nbfield(itangent) = nbfield(itangent) + 1 ; jfail = -9 ; call bfield( RpZb(1:3), itangent, dBRpZb(1:3,0:3), jfail )

    if( ifail.ne.0 .or. jfail.ne.0 ) then ; Lbfieldok = .false. ; return
    endif
    
    if( abs(dBRpZa(2,0)).lt.small ) then
     write(0,'("bf00ab : ! bfield ! : ifail.eq.0, but |B^z| < small ; divide-by-zero error ;")')
     write(0,'("bf00ab : ! bfield ! : B("es23.15","es23.15","es23.15")=("es23.15","es23.15","es23.15")")') &
           RpZa(1:3), dBRpZa(1:3,0)
     Lbfieldok = .false. ; return
    endif
    
    if( abs(dBRpZb(2,0)).lt.small ) then
     write(0,'("bf00ab : ! bfield ! : ifail.eq.0, but |B^z| < small ; divide-by-zero error ;")')
     write(0,'("bf00ab : ! bfield ! : B("es23.15","es23.15","es23.15")=("es23.15","es23.15","es23.15")")') &
           RpZb(1:3), dBRpZb(1:3,0)
     Lbfieldok = .false. ; return
    endif
    
    BRZ(1:2) = (/ dBRpZa(1,0), dBRpZa(3,0) /) /  dBRpZa(2,0) ! normalize to toroidal field; 25 Mar 15;
    BRZ(3:4) = (/ dBRpZb(1,0), dBRpZb(3,0) /) /  dBRpZb(2,0) ! normalize to toroidal field; 25 Mar 15;
    
    dR = RpZb(1) - RpZa(1)
    dZ = RpZb(3) - RpZa(3) ; length = dR*dR + dZ*dZ

    FATAL(bf00ab, length.lt.small, divide by zero in bf00ab -- request soft fail)

    BRZ(5) = ( dR * ( BRZ(4)-BRZ(2) ) - dZ * ( BRZ(3)-BRZ(1) ) ) / length

9999 continue

    return
    
  end subroutine bf00ab

#ifdef HAVE_NAG
  subroutine bf00ac( phi, RZ, BRZ ) ! format constrained by the NAG ode integration routines. 
    implicit none
    INTEGER, parameter  :: Node = 2
#else    
  subroutine bf00ac( Node, phi, RZ, BRZ ) ! format constrained by the lsode integration routines. 
    implicit none
    INTEGER, intent(in) :: Node
#endif
    
    REAL  , intent(in)  :: phi, RZ(1:Node)
    REAL  , intent(out) :: BRZ(1:Node)
    
    INTEGER             :: ibfield
    REAL                :: RpZ(1:3), dBRpZ(1:3,0:3), Bphi, BB
    
    BRZ(1:Node) = zero ! default intent out; 25 Mar 15;
    
    if( .not.Lbfieldok ) return
    
    nbfield(itangent) = nbfield(itangent) + 1
    
    RpZ(1:3) = (/ RZ(1), phi, RZ(2) /) ; ibfield = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), ibfield )
    
    Bphi = dBRpZ(2,0) ; BB = sqrt( dBRpZ(1,0)**2 + (RpZ(1)*dBRpZ(2,0))**2 + dBRpZ(3,0)**2 )

    if( ibfield.ne.0 .or. abs(Bphi).lt.small .or. BB.lt.small ) then ; Lbfieldok = .false. ; return
    endif

    BRZ(1:2) = (/ dBRpZ(1,0), dBRpZ(3,0) /) / Bphi
      
    return
    
  end subroutine bf00ac

  subroutine bf00ba( time, RZp, BRZp ) ! this subroutine is currently unused
    
    implicit none

    INTEGER, parameter  :: Node = 7
    
    REAL  , intent(in)  :: time, RZp(1:Node)
    REAL  , intent(out) :: BRZp(1:Node)
    
    INTEGER             :: ifail
    REAL                :: RpZ(1:3), dBRpZ(1:3,0:3), TM(1:2,1:2), denominator
    
    BRZp(1:Node) = zero ! default intent out; 25 Mar 15;
    
    if( .not.Lbfieldok ) return
    
    RpZ(1:3) = (/ RZp(1), RZp(3), RZp(2) /)

    nbfield(itangent) = nbfield(itangent) + 1 ; ifail = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), ifail )
    
    if( ifail.ne.0 ) then ; Lbfieldok = .false. ; return
    endif
    
    if( poincare%flparameter.eq.0 ) then ; denominator =                               dBRpZ(2,0)
    else                                 ; denominator = sqrt( dBRpZ(1,0)**2 + (RpZ(1)*dBRpZ(2,0))**2 + dBRpZ(3,0)**2 )
    endif
     
    if( abs(denominator).lt.small ) then
     write(0,'("bf00aa : ! bfield ! : ifail.eq.0, but |B^z| or |B| < small ; divide-by-zero error ;")')
     write(0,'("bf00aa : ! bfield ! : B("es23.15","es23.15","es23.15")=("es23.15","es23.15","es23.15")")') &
           RpZ(1:3), dBRpZ(1:3,0)
     Lbfieldok = .false. ; return
    endif

    BRZp(1:3) = (/ dBRpZ(1,0), dBRpZ(3,0), dBRpZ(2,0) /) / denominator
    
    if( itangent.eq.0 ) return

    TM(1,1) = ( dBRpZ(1,1) - BRZp(1) * dBRpZ(2,1) ) / dBRpZ(2,0)
    TM(1,2) = ( dBRpZ(1,3) - BRZp(1) * dBRpZ(2,3) ) / dBRpZ(2,0)
    TM(2,1) = ( dBRpZ(3,1) - BRZp(2) * dBRpZ(2,1) ) / dBRpZ(2,0)
    TM(2,2) = ( dBRpZ(3,3) - BRZp(2) * dBRpZ(2,3) ) / dBRpZ(2,0)
     
    BRZp(4) = TM(1,1) * RZp(4) + TM(1,2) * RZp(6) ! tangent map obtained by matrix multiplication;
    BRZp(5) = TM(1,1) * RZp(5) + TM(1,2) * RZp(7) ! only correct if \phi is integration parameter;
    BRZp(6) = TM(2,1) * RZp(4) + TM(2,2) * RZp(6)
    BRZp(7) = TM(2,1) * RZp(5) + TM(2,2) * RZp(7)
     
9999 continue
    
    return
    
  end subroutine bf00ba
    
#ifdef HAVE_NAG
  subroutine pf00aa( zeta, st, Bst ) ! format constrained by the NAG ode integration routines. 
    implicit none
    INTEGER, parameter  :: Node = 6
#else    
  subroutine pf00aa( Node, zeta, st, Bst ) ! format constrained by the lsode integration routines. 
    implicit none
    INTEGER, intent(in) :: Node
#endif
    
    REAL  , intent(in)  :: zeta, st(1:Node)
    REAL  , intent(out) :: Bst(1:Node)
    
    INTEGER             :: ifail, Lderiv, id
    REAL                :: stz(1:3), dRpZ(1:3,0:3,0:3), sqrtg, dBRpZ(1:3,0:3), dBstz(1:3,0:3)
    REAL                :: det, VM(1:3,1:3), TM(1:2,1:2), ddet, DM(1:3,1:3), DB(1:3,1:3)
        
    Bst(1:Node) = zero ; stz(1:3) = zero ; dRpZ(1:3,0:3,0:3) = zero ; dBRpZ(1:3,0:3) = zero ; dBstz(1:3,0:3) = zero
    
    if( .not.Lbfieldok ) goto 9999
    
    Lderiv = itangent + 1 ; stz(1:3) = (/ st(1), st(2), zeta /)
    call bc00ab( rzmn, stz(1:3), Lderiv, dRpZ(1:3,0:3,0:3), sqrtg, ifail )
    
    if( ifail.ne.0 ) then ; Lbfieldok = .false. ; goto 9999
    endif

    nbfield(itangent) = nbfield(itangent) + 1 ; ifail = -9
    call bfield( dRpZ(1:3,0,0), itangent, dBRpZ(1:3,0:3), ifail )
    
    if( ifail.ne.0 ) then ; Lbfieldok = .false. ; goto 9999
    endif
    
    if( abs(dBRpZ(2,0)).lt.small ) then
     write(0,'("pf00aa : ! bfield ! : ifail.eq.0, but |B^z| < small ; divide-by-zero error ;")')
     write(0,'("pf00aa : ! bfield ! : B("es23.15","es23.15","es23.15")=("es23.15","es23.15","es23.15")")') &
           dRpZ(1:3,0,0), dBRpZ(1:3,0)
     Lbfieldok = .false. ; goto 9999
    endif

    VM(1,1:3) = (/ - dRpZ(3,2, 0), dRpZ(1,3, 0) * dRpZ(3,2, 0) - dRpZ(1,2, 0) * dRpZ(3,3, 0), + dRpZ(1,2, 0) /)
    VM(2,1:3) = (/ + dRpZ(3,1, 0), dRpZ(1,1, 0) * dRpZ(3,3, 0) - dRpZ(1,3, 0) * dRpZ(3,1, 0), - dRpZ(1,1, 0) /)
    VM(3,1:3) = (/   zero        , dRpZ(1,2, 0) * dRpZ(3,1, 0) - dRpZ(1,1, 0) * dRpZ(3,2, 0),   zero         /)
    
    dBstz(1:3,0) = matmul( VM(1:3,1:3), dBRpZ(1:3,0) )
    
    Bst(1:2) = (/ dBstz(1,0), dBstz(2,0) /) / dBstz(3,0)
    
    Bst(1) = Bst(1) - actiongradient
    
    if( itangent.eq.0 ) goto 9999
        
    DB(1,1:3) = (/ dBRpZ(1,1), dBRpZ(1,2), dBRpZ(1,3) /)
    DB(2,1:3) = (/ dBRpZ(2,1), dBRpZ(2,2), dBRpZ(2,3) /)
    DB(3,1:3) = (/ dBRpZ(3,1), dBRpZ(3,2), dBRpZ(3,3) /)
    
    do id = 1, 3
    DM(1,1:3) = (/ - dRpZ(3,2,id), dRpZ(1,3,id) * dRpZ(3,2, 0) - dRpZ(1,2,id) * dRpZ(3,3, 0)                    &
                                 + dRpZ(1,3, 0) * dRpZ(3,2,id) - dRpZ(1,2, 0) * dRpZ(3,3,id), + dRpZ(1,2,id) /) 
    DM(2,1:3) = (/ + dRpZ(3,1,id), dRpZ(1,1,id) * dRpZ(3,3, 0) - dRpZ(1,3,id) * dRpZ(3,1, 0)                    &
                                 + dRpZ(1,1, 0) * dRpZ(3,3,id) - dRpZ(1,3, 0) * dRpZ(3,1,id), - dRpZ(1,1,id) /) 
    DM(3,1:3) = (/   zero        , dRpZ(1,2,id) * dRpZ(3,1, 0) - dRpZ(1,1,id) * dRpZ(3,2, 0)                    &
                                 + dRpZ(1,2, 0) * dRpZ(3,1,id) - dRpZ(1,1, 0) * dRpZ(3,2,id),   zero         /)

    dBstz(1:3,id) = matmul( DM(1:3,1:3), dBRpZ(1:3,0) ) + matmul( VM(1:3,1:3), matmul( DB(1:3,1:3), dRpZ(1:3,id,0) ) )

    enddo
    
    TM(1,1) = ( dBstz(1,1) - dBstz(1,0) * dBstz(3,1) / dBstz(3,0) ) / dBstz(3,0)
    TM(1,2) = ( dBstz(1,2) - dBstz(1,0) * dBstz(3,2) / dBstz(3,0) ) / dBstz(3,0)
    TM(2,1) = ( dBstz(2,1) - dBstz(2,0) * dBstz(3,1) / dBstz(3,0) ) / dBstz(3,0)
    TM(2,2) = ( dBstz(2,2) - dBstz(2,0) * dBstz(3,2) / dBstz(3,0) ) / dBstz(3,0)
    
    Bst(3) = TM(1,1) * st(3) + TM(1,2) * st(5) ! tangent map obtained by matrix multiplication;
    Bst(4) = TM(1,1) * st(4) + TM(1,2) * st(6) - one
    Bst(5) = TM(2,1) * st(3) + TM(2,2) * st(5)
    Bst(6) = TM(2,1) * st(4) + TM(2,2) * st(6)
    
9999 continue
    
!   write(0,'("pf00aa :          "2x" : stz="3es13.5" ; RpZ="3es13.5" ; BRpZ="3es13.5" ; Bstz="3es13.5" ;")') &
!         stz(1:3), dRpZ(1:3,0,0), dBRpZ(1:3,0), Bstz(1:3)

    return
    
  end subroutine pf00aa

!latex \newpage \section{Biot-Savart subroutines}

!latex In this section are described subroutines for computing the magnetic field produced by a current distribution.

!latex \subroutine{bs00aa}{compute the magnetic field produced by a filamentary current loop of arbitrary shape;}
!latex \bi
!latex \item[1.] Given a closed (i.e. periodic), one-dimensional loop embedded in three-dimensional space, 
!latex with position described by $\bar {\bf x}(t) \equiv \bar x(t) {\bf i} + \bar y(t) {\bf j} + \bar z(t) {\bf k}$, 
!latex with the arbitrary curve parameter $t\in[0,2\pi]$, and $\bar {\bf x}(t+2\pi) = \bar {\bf x}(t)$,
!latex assumed to carry unit current, i.e. $I=1$,  
!latex the magnetic field at ${\bf x} \equiv x {\bf i} + y {\bf j} + z {\bf k}$ is given by the Biot-Savart integral,
!latex \be {\bf B} \equiv \int_{\cal C} \frac{d{\bf l}\times {\bf r}}{r^3}, \label{eq:BiotSavart}
!latex \ee
!latex where ${\bf r}\equiv {\bf x}-\bar{\bf x}$.
!latex \item[2.] In component form, \Eqn{BiotSavart} is
!latex \be B^x & \equiv & \int_0^{2\pi} \frac{\dot {\bar y} (z-\bar z) - \dot {\bar z} (y-\bar y) }{r^3} dt, \\
!latex     B^y & \equiv & \int_0^{2\pi} \frac{\dot {\bar z} (x-\bar x) - \dot {\bar x} (z-\bar z) }{r^3} dt, \\
!latex     B^z & \equiv & \int_0^{2\pi} \frac{\dot {\bar x} (y-\bar y) - \dot {\bar y} (x-\bar x) }{r^3} dt,
!latex \ee
!latex where $\dot {\bar x} \equiv d\bar x/dt$, etc.
!latex \item[3.] The magnetic vector potential is
!latex \be {\bf A} \equiv \int_{\cal C} \frac{d{\bf l}}{r}. \label{eq:BiotSavartA}
!latex \ee
!latex \item[4.] The total length of the curve is
!latex \be      L  \equiv \int_{\cal C}       d     l     , \label{eq:BiotSavartL}
!latex \ee
!latex          where $dl \equiv (\dot {\bar x}^2 + \dot {\bar y}^2 + \dot {\bar z}^2 )^{1/2}$.
!latex \item[5.] The user must supply a subroutine, \verb+iccoil+, that returns $\bar x$, $\bar y$ \& $\bar z$, 
!latex           and $d \bar x /dt$, $d \bar y /dt$ \& $d \bar z /dt$, given $t$:
!latex           \verb+subroutine iccoil( t, x(0:1), y(0:1), z(0:1), ifail )+ \\
!latex           where \verb+t+, \verb+x(0:1)+, \verb+y(0:1)+ and \verb+z(0:1)+ are real and ifail is an integer,
!latex           \verb+x(0)+$\equiv \bar x(t)$ and \verb+x(1)+$\equiv \dot {\bar x}(t)$, and similarly for $y$ and $z$.
!latex \item[6.] The integration is performed using \nag{www.nag.co.uk/numeric/FL/manual19/pdf/D01/d01ajf_fl19.pdf}{D01AJF}.
!latex           (This routine is based upon the {\footnotesize QUADPACK} routine QAGS, which is freely available.)
!latex \item[7.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : biotsavart, bs00aa+ \\ \\
!latex           \verb+type(biotsavart) :: bsfield+ \\ \\ in their source that calls \verb+bs00aa+,
!latex           where \verb+biotsavart+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+bsfield+, is arbitrary.
!latex \item[8.] \underline{\bf Required inputs}
!latex \item[  ] \verb+bsfield%x                   : real    ;+
!latex \item[  ] \verb+bsfield%y                   : real    ;+
!latex \item[  ] \verb+bsfield%z                   : real    ;+
!latex \bi
!latex \item[i.] position ${\bf x} \equiv x \, {\bf i} + y \, {\bf j} + z \, {\bf k}$ at which magnetic field is required;
!latex \ei
!latex \item[  ] \verb+bsfield%tol                 : real    ;+
!latex \item[  ] \verb+bsfield%N                   : integer ;+
!latex \bi
!latex \item[i.] integration accuracy parameters provided to \verb+D01AJF+; (\verb+EPSABS=tol+, \verb+EPSREL=zero+ and \verb+LW=4*N+);
!latex \ei
!latex \item[  ] \verb+bsfield%LB                  : logical ;+
!latex \bi
!latex \item[i.] set \verb+LB = .true.+ to compute the magnetic field;
!latex \ei
!latex \item[  ] \verb+bsfield%LA                  : logical ;+
!latex \bi
!latex \item[i.] set \verb+LA = .true.+ to compute the magnetic vector potential;
!latex \ei
!latex \item[  ] \verb+bsfield%LL                  : logical ;+
!latex \bi
!latex \item[i.] set \verb+LL = .true.+ to compute the length of the curve;
!latex \ei
!latex \item[9.] \underline{\bf Execution}
!latex \item[  ] \verb+call bs00aa( bsfield, ifail )+
!latex \item[10.] \underline{\bf Outputs}
!latex \item[  ] \verb+bsfield%Bx                  : real    ;+
!latex \item[  ] \verb+bsfield%By                  : real    ;+
!latex \item[  ] \verb+bsfield%Bz                  : real    ;+
!latex \bi
!latex \item[i.] only if \verb+LB = .true.+;
!latex \ei
!latex \item[  ] \verb+bsfield%Ax                  : real    ;+
!latex \item[  ] \verb+bsfield%Ay                  : real    ;+
!latex \item[  ] \verb+bsfield%Az                  : real    ;+
!latex \bi
!latex \item[i.] only if \verb+LA = .true.+;
!latex \ei
!latex \item[  ] \verb+bsfield%length              : real    ;+
!latex \bi
!latex \item[i.] only if \verb+LL = .true.+;
!latex \ei
!latex \item[  ] \verb+ifail                       : integer ;+
!latex \bi
!l!tex \item[i.] on output: 
!latex     \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex     \item[] \verb+ifail=2+ : \verb+D01AJF+ encountered a divide-by-zero;
!latex                              this is only possible if $\exists t \in [0,2\pi]$ such that ${\bf x}=\bar{\bf x}(t)$.
!latex \ei
!latex \ei

  subroutine bs00aa( lbsfield, ifail )
    
    implicit none
    
    type(biotsavart)     :: lbsfield
    INTEGER              :: ifail
    
    INTEGER              :: ixyz, id01ajf, LRwork, LIwork, astat
    INTEGER, allocatable :: Iwork(:)
    REAL                 :: lowlimit, upplimit, epsabs, epsrel
    REAL                 :: abserrbx, abserrby, abserrbz, abserrax, abserray, abserraz, abserr
    REAL   , allocatable :: Rwork(:)

#ifndef HAVE_NAG
    INTEGER              :: neval, last
#endif

    
    ibs00aa = ifail

    CHECKINPUT( bs00aa, lbsfield%tol.le.small, 9999 )
    CHECKINPUT( bs00aa, lbsfield%N  .le.    0, 9999 )
    
     bsfield%x      = lbsfield%x
     bsfield%y      = lbsfield%y
     bsfield%z      = lbsfield%z

    lbsfield%Bx     = zero
    lbsfield%By     = zero
    lbsfield%Bz     = zero
    lbsfield%Ax     = zero
    lbsfield%Ay     = zero
    lbsfield%Az     = zero
    lbsfield%length = zero

    abserrbx = - one ; abserrby = - one ; abserrbz = - one
    abserrax = - one ; abserray = - one ; abserraz = - one
    abserr   = - one

    LRwork = lbsfield%N * 4
    SALLOCATE(Rwork,(1:LRwork),zero)

    LIwork = LRwork
    SALLOCATE(Iwork,(1:LIwork),zero)
    
    do ixyz = 1, 3
     
     bsfield%ixyz = ixyz
     
     if( lbsfield%LB ) then ! compute magnetic field; 09 Nov 15;
      id01ajf = 1 ; lowlimit = zero ; upplimit = pi2 ; epsabs = lbsfield%tol ; epsrel = zero
      select case( ixyz) 

#ifdef HAVE_NAG
      case( 1 )
        call D01AJF( bs00bx, lowlimit, upplimit, epsabs, epsrel, lbsfield%Bx, abserrbx, Rwork(1:LRwork), LRwork, &
                     Iwork(1:LIwork), LIwork, id01ajf )
      case( 2 )
        call D01AJF( bs00bx, lowlimit, upplimit, epsabs, epsrel, lbsfield%By, abserrby, Rwork(1:LRwork), LRwork, &
                     Iwork(1:LIwork), LIwork, id01ajf )
      case( 3 )
        call D01AJF( bs00bx, lowlimit, upplimit, epsabs, epsrel, lbsfield%Bz, abserrbz, Rwork(1:LRwork), LRwork, &
                     Iwork(1:LIwork), LIwork, id01ajf )
#else
      case( 1 )
        call DQAGS( bs00bx, lowlimit, upplimit, epsabs, epsrel, lbsfield%Bx, abserrbx, neval, id01ajf, LIwork, LRwork, &
                    last, Iwork(1:LIwork), Rwork(1:LRwork) )
      case( 2 )
        call DQAGS( bs00bx, lowlimit, upplimit, epsabs, epsrel, lbsfield%By, abserrby, neval, id01ajf, LIwork, LRwork, &
                    last, Iwork(1:LIwork), Rwork(1:LRwork) )
      case( 3 )
        call DQAGS( bs00bx, lowlimit, upplimit, epsabs, epsrel, lbsfield%Bz, abserrbz, neval, id01ajf, LIwork, LRwork, &
                    last, Iwork(1:LIwork), Rwork(1:LRwork) )
#endif

      end select
      if( id01ajf.ne.0 ) ibs00aa = 2
     endif
          
     if( lbsfield%LA ) then ! compute magnetic vector potential; 09 Nov 15;
      id01ajf = 1 ; lowlimit = zero ; upplimit = pi2 ; epsabs = lbsfield%tol ; epsrel = zero
      select case( ixyz) 

#ifdef HAVE_NAG
      case( 1 )
        call D01AJF( bs00ax, lowlimit, upplimit, epsabs, epsrel, lbsfield%Ax, abserrax, Rwork(1:LRwork), LRwork, &
                     Iwork(1:LIwork), LIwork, id01ajf )
      case( 2 )
        call D01AJF( bs00ax, lowlimit, upplimit, epsabs, epsrel, lbsfield%Ay, abserray, Rwork(1:LRwork), LRwork, &
                     Iwork(1:LIwork), LIwork, id01ajf )
      case( 3 )
        call D01AJF( bs00ax, lowlimit, upplimit, epsabs, epsrel, lbsfield%Az, abserraz, Rwork(1:LRwork), LRwork, &
                     Iwork(1:LIwork), LIwork, id01ajf )
#else
      case( 1 )
        call DQAGS( bs00ax, lowlimit, upplimit, epsabs, epsrel, lbsfield%Ax, abserrbx, neval, id01ajf, LIwork, LRwork, &
                    last, Iwork(1:LIwork), Rwork(1:LRwork) )
      case( 2 )
        call DQAGS( bs00ax, lowlimit, upplimit, epsabs, epsrel, lbsfield%Ay, abserrby, neval, id01ajf, LIwork, LRwork, &
                    last, Iwork(1:LIwork), Rwork(1:LRwork) )
      case( 3 )
        call DQAGS( bs00ax, lowlimit, upplimit, epsabs, epsrel, lbsfield%Az, abserrbz, neval, id01ajf, LIwork, LRwork, &
                    last, Iwork(1:LIwork), Rwork(1:LRwork) )
#endif
      
      end select
      if( id01ajf.ne.0 ) ibs00aa = 2
     endif
     
    enddo ! end of do ixyz; 13 Oct 15;
               
     if( lbsfield%LL ) then ! compute length; 09 Nov 15;
      lowlimit = zero ; upplimit = pi2 ; epsabs = lbsfield%tol ; epsrel = zero

#ifdef HAVE_NAG
      id01ajf = 1
      call D01AJF( bs00lx, lowlimit, upplimit, epsabs, epsrel, lbsfield%length, abserr, Rwork(1:LRwork), LRwork, &
                   Iwork(1:LIwork), LIwork, id01ajf )
#else
      call DQAGS( bs00lx, lowlimit, upplimit, epsabs, epsrel, lbsfield%length, abserr, neval, id01ajf, LIwork, LRwork, &
                  last, Iwork(1:LIwork), Rwork(1:LRwork) )
#endif
      
      if( id01ajf.ne.0 ) ibs00aa = 2
     endif

    DALLOCATE(Iwork)
    DALLOCATE(Rwork)

!   lbsfield%its = bsfield%its

   !ibs00aa = 0 ! success; 02 Jun 15;
    
9999 continue
    
    if( ibs00aa.ne. 0 ) write(0,9000) ibs00aa, lbsfield%tol, lbsfield%x, lbsfield%y, lbsfield%z, &
                                      abserrbx, abserrby, abserrbz
    
9000 format("bs00aa : d01ajf "i3" : tol="es9.1" ; (x,y,z) = ("f23.15" ,"f23.15" ,"f23.15" ) ; " &
            :"dB(x,y,z) = ("es23.15" ,"es23.15" ,"es23.15" ) ; "i9" ; ")
    
    ifail = ibs00aa
    
    return
    
  end subroutine bs00aa
  
  REAL function bs00bx( tt )
    
    implicit none
    
    REAL                 :: tt
    
    REAL                 :: xx(0:1), yy(0:1), zz(0:1), dx, dy, dz, rr, r3
    INTEGER              :: ifail
    
    ifail = 0
    call iccoil( tt, xx(0:1), yy(0:1), zz(0:1), ifail )
    
    dx = xx(0) - bsfield%x
    dy = yy(0) - bsfield%y
    dz = zz(0) - bsfield%z
    
    rr = sqrt( dx**2 + dy**2 + dz**2 ) ; r3 = rr**3
    
    if( rr.lt.small ) then ; ibs00aa = 2 ; bs00bx = ten**6 ; return
    endif

    select case( bsfield%ixyz)
    case( 1 ) ; bs00bx = ( dy * zz(1) - dz * yy(1) ) / r3
    case( 2 ) ; bs00bx = ( dz * xx(1) - dx * zz(1) ) / r3
    case( 3 ) ; bs00bx = ( dx * yy(1) - dy * xx(1) ) / r3
    end select
    
!   bsfield%its = bsfield%its + 1

    return
    
  end function bs00bx
  
  REAL function bs00ax( tt )
    
    implicit none
    
    REAL                 :: tt
    
    REAL                 :: xx(0:1), yy(0:1), zz(0:1), dx, dy, dz, rr
    INTEGER              :: ifail
    
    ifail = 0
    call iccoil( tt, xx(0:1), yy(0:1), zz(0:1), ifail )
    
    dx = xx(0) - bsfield%x
    dy = yy(0) - bsfield%y
    dz = zz(0) - bsfield%z
    
    rr = sqrt( dx**2 + dy**2 + dz**2 )
    
    if( rr.lt.small ) then ; ibs00aa = 2 ; bs00ax = ten**6 ; return
    endif

    select case( bsfield%ixyz)
    case( 1 ) ; bs00ax = xx(1) / rr
    case( 2 ) ; bs00ax = yy(1) / rr
    case( 3 ) ; bs00ax = zz(1) / rr
    end select
    
!   bsfield%its = bsfield%its + 1

    return
    
  end function bs00ax
  
  REAL function bs00lx( tt )
    
    implicit none
    
    REAL                 :: tt
    
    REAL                 :: xx(0:1), yy(0:1), zz(0:1)
    INTEGER              :: ifail
    
    ifail = 0
    call iccoil( tt, xx(0:1), yy(0:1), zz(0:1), ifail )
    
    bs00lx = sqrt( xx(1)**2 + yy(1)**2 + zz(1)**2 )
    
    return
    
  end function bs00lx
  
!latex \newpage \section{``cylindrical'' subroutines}

!latex In this section are described subroutines that do not depend explicitly on a pre-defined, `background', toroidal coordinate framework.

!latex \subroutine{ga00aa}{find the magnetic axis;}
!latex \bi
!latex \item[1.] Iterative fieldline tracing methods are used to find the magnetic axis,
!latex           defined as the magnetic fieldline that closes on itself after a toroidal distance of $\Delta \p = 2\pi$/\verb+Nfp+,
!latex           i.e. ${\bf x}(\Delta\p)={\bf x}(0)$,
!latex           where \verb+Nfp+ is the field periodicity.
!latex \item[ *] The fieldline mapping is defined by integrating along the magnetic field,
!latex           and is constructed numerically in cylindrical coordinates by integrating the o.d.e.'s
!latex           \be \frac{dR(\p)}{d\p} & = & \frac{B^R(R,\p,Z)}{B^\p(R,\p,Z)} \equiv \dot R(R,\p,Z), \label{eq:BR} \\
!latex               \frac{dZ(\p)}{d\p} & = & \frac{B^Z(R,\p,Z)}{B^\p(R,\p,Z)} \equiv \dot Z(R,\p,Z), \label{eq:BZ}
!latex           \ee
!latex           from an initial, user-supplied starting point, $(R_0,0,Z_0)$.
!latex           The toroidal angle, $\p$, is used as the integration parameter, and so $B^\p$ cannot be zero.
!latex           Upon request, this routine will be modified in order to follow field lines in regions where $B^\p=0$.
!latex \item[ *] A Newton-iterative method is used to find the zero of
!latex           \be {\bf f}\left(\begin{array}{c}R_0\\Z_0\end{array}\right) \equiv \left(\begin{array}{c}R_1-R_0\\Z_1-Z_0\end{array}\right)
!latex           \ee
!latex           where $R_1 \equiv R(\Delta\p)$ and $Z_1 \equiv Z(\Delta\p)$.
!latex \item[ *] Given an initial guess, ${\bf x}\equiv(R_0, Z_0)^T$,
!latex           a better guess for the location of the axis, $(R_0,Z_0)^T+(\delta R,\delta Z)^T$, is given by the linear approximation
!latex           \be {\bf f}\left(\begin{array}{c}R_0+\delta R_0\\Z_0+\delta Z_0\end{array}\right) = {\bf f}\left(\begin{array}{c}R_0\\Z_0\end{array}\right)
!latex               + \underbrace{
!latex               \left( \begin{array}{lcl}\partial_{R_0}R_1 -1&, & \partial_{Z_0}R_1 \\ 
!latex               \partial_{R_0}Z_1&, & \partial_{Z_0}Z_1-1\end{array} \right)}_{\nabla {\bf f}}  \cdot
!latex               \left(\begin{array}{c}\delta R_0\\ \delta Z_0\end{array}\right) + {\cal O}(\delta^2) = 0,
!latex           \ee
!latex           and the correction is given by $\delta{\bf x} = - (\nabla {\bf f})^{-1} \cdot {\bf f}({\bf x})$.
!latex \item[ *] The derivatives, $\partial_{R_0}R_1$, $\partial_{Z_0}R_1$, etc. are determined by fieldline integration,
!latex           \be \frac{d}{d\p} \left( \begin{array}{cc}\partial_{R_0}R(\p), & \partial_{Z_0}R(\p) \\ 
!latex            \partial_{R_0}Z(\p), & \partial_{Z_0}Z(\p) \end{array}\right) = 
!latex                             \left( \begin{array}{cc}\partial_{R  }\dot R, & \partial_{Z  }\dot R\\ 
!latex           \partial_{R  }\dot Z, & \partial_{Z  }\dot Z \end{array}\right) \cdot
!latex                             \left( \begin{array}{cc}\partial_{R_0}R(\p), & \partial_{Z_0}R(\p) \\
!latex           \partial_{R_0}Z(\p), & \partial_{Z_0}Z(\p) \end{array}\right),
!latex           \ee
!latex           from an initial starting point being the identity matrix,
!latex           \be               \left( \begin{array}{cc}\partial_{R_0}R( 0), & \partial_{Z_0}R( 0) \\ \partial_{R_0}Z( 0), & \partial_{Z_0}Z( 0) 
!latex                                      \end{array}\right) = 
!latex                             \left( \begin{array}{cc} 1, & 0 \\ 0,& 1                                                              
!latex                                      \end{array}\right).
!latex           \ee
!l!tex           i.e. $\partial_{R_0}R_1(0)=1$, $\partial_{Z_0}R_1(0)=0$, $\partial_{R_0}Z_1(0)=0$ and $\partial_{Z_0}Z_1(0)=1$.
!latex \item[ *] The iterative search is enabled by \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pbf_fl19.pdf}{C05PBF}.
!latex \item[ *] If \verb+ifail=-4+, then instead the axis search is provided by
!latex           \nag{http://www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05nbf_fl19.pdf}{C05NBF}, which does not require derivatives.
!latex           Note that for this option, the tangent map, the transform on axis, and the residue will not be calculated.
!latex \item[ *] The above definition of the magnetic axis does not have a unique solution:
!latex           an $\iotabar=1/1$ fieldline also satisfies this definition, 
!latex           as does the $\iotabar=2/1$, $3/1$, etc., as also does the ``X'' point at the separatrix.
!latex           Furthermore, during a sawteeth cycle, the $\iotabar=1/1$ fieldline and the original magnetic axis swap places.
!latex           If there is a continuous family of ``magnetic axes", e.g. there exists an intact $q=1$ surface,
!latex           then $\nabla {\bf f}$ will not be invertible (unless singular value decomposition methods are used).
!latex           Thus, this routine should be used with care.
!latex           Which closed field line that \verb+ga00aa+ locates is determined by the initial guess provided.
!l!tex           Only the ``true'' magnetic axis should be located with \verb+ga00aa+.
!l!tex           Frequently, the coordinate harmonics of the magnetic axis returned by \verb+ga00aa+
!l!tex           will be used as the coordinate axis of toroidal ``chaotic'' coordinates.
!l!tex           Other routines will be provided for locating the $X$ point etc.
!l!tex \item[ *] The definition of the ``true'' magnetic axis may be made precise by examining the Fourier harmonic content of the closed fieldline,
!l!tex           and defining the true magnetic axis may be that closed fieldline with minimal Fourier content.
!latex \item[ *] The returned information includes: 
!latex \bi
!latex \item[i.] the Fourier representation of $R(\p)$ and $Z(\p)$;
!latex \item[ii.] the tangent-mapping near the axis, which allows the rotational-transform on axis to be determined;
!latex \item[iii.] Greene's residue [Greene, \link{dx.doi.org/10.1063/1.524170}{J. Math. Phys. 20, 1183 (1979)}]
!latex             calculated at the magnetic axis (determines stability).
!latex \ei
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : magneticaxis, ga00aa+ \\ \\
!latex           \verb+type(magneticaxis) :: axis+ \\ \\ in their source that calls \verb+ga00aa+,
!latex           where \verb+axis+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+axis+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+axis%Nfp              : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, e.g. \verb+Nfp=1+;
!latex \ei
!latex \item[  ] \verb+axis%Ntor             : integer ;+
!latex \bi
!latex \item[i.] the desired Fourier resolution of the magnetic axis, 
!latex \item[ii.] if it is not required to have a Fourier decomposition of the magnetic axis, or the magnetic field is axisymmetric, choose \verb+Ntor=0+;
!latex \ei
!latex \item[  ] \verb+axis%R               : real    ;+
!latex \bi
!latex \item[i.] guess for the $R$ location of the magnetic axis on the $\phi=0$ plane;
!latex \ei
!latex \item[  ] \verb+axis%Z               : real    ;+
!latex \bi
!latex \item[i.] guess for the $Z$ location of the magnetic axis on the $\phi=0$ plane;
!latex \ei
!latex \item[  ] \verb+axis%maxits           : integer ;+
!latex \bi
!latex \item[i.] max. iterations allowed in search;
!latex \item[ii.] e.g. \verb+maxits=16+;
!latex \ei
!latex \item[  ] \verb+axis%tol              : real    ;+
!latex \bi
!latex \item[i.] required accuracy to which the position of the magnetic axis on the $\p=0$ plane is required
!latex \item[ii.] e.g. \verb+tol=1.0e-06+;
!latex \ei
!latex \item[  ] \verb+axis%odetol           : real    ;+
!latex \bi
!latex \item[i.] o.d.e. integration tolerance; 
!latex \item[ii.] e.g. \verb+odetol=1.0e-08+;
!latex \ei
!latex \item[  ] \verb+ifail                      : integer ;+
!latex \bi
!latex \item[i.]   if \verb+ifail.ge. 1+ : there is no screen output, {\em except} if there is an input error;
!latex \item[ii.]  if \verb+ifail.le. 0+ : a one-line summary is provided, giving the $(R,Z)$ location of the axis, $\iotabar_{axis}$, etc.
!latex \item[iii.] if \verb+ifail.le.-1+ : the Fourier harmonics of the magnetic axis are displayed;
!latex \item[iv.]  if \verb+ifail.le.-2+ : the eigenvalues and eigenvectors of the tangent map at the axis are displayed;
!latex \item[iv.]  if \verb+ifail.le.-3+ : information detailing the progress of the iterative search is provided;
!latex \item[iv.]  if \verb+ifail.le.-4+ : C05NBF (which does not require derivatives) is used instead of C05PBF;
!latex \ei
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call ga00aa( axis, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+axis%R                : real    ;+
!latex \bi
!latex \item[i.] updated;
!latex \ei
!latex \item[  ] \verb+axis%Z                : real    ;+
!latex \bi
!latex \item[i.] updated;
!latex \ei
!latex \item[  ] \verb+axis%tangent(1:2,1:2) : real    ;+ 
!latex \bi
!latex \item[i.]  the tangent mapping at axis;
!latex \item[ii.] if the eigenvalues of the tangent map are imaginary, e.g. $\lambda\equiv\alpha+\beta i$, 
!latex            then the rotational-transform on axis satisfies $\tan(||\iotabar||)=\beta/\alpha$, where $||\iotabar||\equiv\iotabar \mod 2\pi$.
!latex \item[iii.] if the eigenvalues of the tangent map are real, then the eigenvalues give the direction of the stable and unstable manifolds.
!latex \ei
!latex \item[  ] \verb+axis%wr(1:2)          : real    ;+ 
!latex \item[  ] \verb+axis%wi(1:2)          : real    ;+ 
!latex \item[  ] \verb+axis%vr(1:2,1:2)      : real    ;+ 
!latex \item[  ] \verb+axis%vi(1:2,1:2)      : real    ;+ 
!latex \bi
!latex \item[i.]  the eigenvalues and eigenvectors of the tangent mapping at axis;
!latex \ei
!latex \item[  ] \verb+axis%iota             : real    ;+
!latex \bi
!latex \item[i.] rotational-transform on axis
!latex \item[ii.] will only be meaningful if axis is stable, which is indicated by both (i) the sign of residue, 
!latex            and (ii) whether the eigenvalues and eigenvectors are real or imaginary.
!latex \ei
!latex \item[  ] \verb+axis%Lallocated       : integer ;+
!latex \bi
!latex \item[i.] if \verb+Lallocated=1+, the \verb+Ri+, \verb+Zi+, \verb+Rnc+, \verb+Rns+, \verb+Zns+, \verb+Znc+ have           been allocated;
!latex \item[ii.] if \verb+Lallocated=0+, the \verb+Ri+, \verb+Zi+, \verb+Rnc+, \verb+Rns+, \verb+Zns+, \verb+Znc+ have {\bf not} been allocated;
!latex \ei
!latex \item[  ] \verb+axis%Ri(0:4*Ntor)     : real    ;+
!latex \bi
!latex \item[i.] the magnetic axis, $R(i\Delta\varphi)$, for $i=0,4*$\verb+Ntor+, where $\Delta\varphi=\Delta\p/(4*$\verb+Ntor+$)$;
!latex \item[ii.] \verb+Ri+ is allocated internally; if on input \verb+Ri+ is already allocated it will first be deallocated; similarly for \verb+Zi+
!latex \ei
!latex \item[  ] \verb+axis%Zi(0:4*Ntor)     : real    ;+
!latex \bi
!latex \item[i.] the magnetic axis, $Z(i\Delta\varphi)$, for $i=0,4*$\verb+Ntor+, where $\Delta\varphi=\Delta\p/(4*$\verb+Ntor+$)$;
!latex \ei
!latex \item[  ] \verb+axis%Rnc(0:Ntor)      : real    ;+
!latex \bi
!latex \item[i.] the Fourier harmonics, $R(\phi)=\sum_n [R_{n,c} \cos(-n\phi)+R_{n,s} \sin(-n\phi)]$;
!latex \item[ii.] \verb+Rnc+ is allocated internally; if on input \verb+Rnc+ is already allocated it will first be deallocated;
!latex            similarly for \verb+Zns+, \verb+Rns+ and \verb+Znc+.
!latex \ei
!latex \item[  ] \verb+axis%Zns(0:Ntor)      : real    ;+
!latex \bi
!latex \item[i.] the Fourier harmonics, $Z(\phi)=\sum_n [Z_{n,c} \cos(-n\phi)+Z_{n,s} \sin(-n\phi)]$;
!latex \ei
!latex \item[  ] \verb+axis%Rns(0:Ntor)      : real    ;+
!latex \bi
!latex \item[i.] the Fourier harmonics, $R(\phi)=\sum_n [R_{n,c} \cos(-n\phi)+R_{n,s} \sin(-n\phi)]$;
!latex \ei
!latex \item[  ] \verb+axis%Znc(0:Ntor)      : real    ;+
!latex \bi
!latex \item[i.] the Fourier harmonics, $Z(\phi)=\sum_n [Z_{n,c} \cos(-n\phi)+Z_{n,s} \sin(-n\phi)]$;
!latex \ei
!latex \item[  ] \verb+axis%error            : real    ;+
!latex \bi
!latex \item[i.] the error, $\sqrt{\Delta R^2 + \Delta Z^2}$, where $\Delta R \equiv R(\Delta\p)-R_0$ and $\Delta Z \equiv Z(\Delta\p)-Z_0$.
!latex \ei
!latex \item[  ] \verb+axis%its              : integer ;+
!latex \bi
!latex \item[i.] the number of iterations required;
!latex \ei
!latex \item[  ] \verb+axis%residue          : real    ;+
!latex \bi
!latex \item[i.] Greene's residue of the magnetic axis; [Greene, \link{dx.doi.org/10.1063/1.524170}{J. Math. Phys. 20, 1183 (1979)}];
!latex \ei
!latex \item[  ] \verb+axis%rzf(0:2,0:31)    : real    ;+
!latex \bi
!latex \item[i.] Information regarding the progress of the iterations;
!latex \item[ii.] $R_i = $ \verb+axis%rzf(0,0:axis%its)+ are the $R$ values used in the iterations;
!latex \item[ii.] $Z_i = $ \verb+axis%rzf(1,0:axis%its)+ are the $Z$ values used in the iterations;
!latex \item[ii.] $F_i = $ \verb+axis%rzf(2,0:axis%its)+ are the $|f|$ at each iteration;
!latex \ei
!latex \item[  ] \verb+ifail                 : integer ;+
!latex \bi
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex     \item[] \verb+ifail=2+ : the routine \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pbf_fl19.pdf}{C05PBF}
!latex                              failed to locate the zero of the function,
!latex                              perhaps because of a failure in integrating along the fieldlines;
!latex     \item[] \verb+ifail=3+ : the NAG routine \verb+C06EAF+ failed to construct the Fourier harmonics of the axis, for either $R$ or $Z$;
!latex     \item[] \verb+ifail=4+ : the NAG routine \verb+F02EBF+ failed to construct the eigenvalues/vectors of the tangent mapping;
!latex \ei
!latex \ei
!latex \item[6.] Comments: 
!latex \item[* ] The NAG routine \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pbf_fl19.pdf}{C05PBF}
!latex           is used for the nonlinear root find, and \verb+tol+ is given directly to
!latex           \nag{www.nag.co.uk/numeric/FL/manual19/pdf/C05/c05pbf_fl19.pdf}{C05PBF}.
!latex \item[* ] The NAG routine \verb+D02BJF+ is used for the o.d.e. integration, and \verb+odetol+ is supplied directly to \verb+D02BJF+.
!latex \item[* ] If a good initial guess is given, this number should be small, as Newton methods should converge rapidly; 
!latex           however, if there are multiple magnetic axes (as during a sawtooth event) then the Newton method may encounter problems;
!latex           also, numerical errors in the magnetic field (perhaps $\nabla\cdot{\bf B}$ is not exactly zero)
!latex           can cause the fieldline integration to be inaccurate, and so it may be difficult to find the solution to the desired accuracy. 
!latex \item[* ] Please also consider using \verb+ec00aa+;
!latex \ei

  subroutine ga00aa( laxis, ifail )

    implicit none
    
    type(magneticaxis)   :: laxis
    INTEGER              :: ifail
   
    INTEGER, parameter   :: NN = 2, Ldf = NN, Lwk = 3 * NN + 13, NM = 2, Lda = NM, Ldvr = NM, Ldvi = NM, Lwka = 4 * NM
    INTEGER              :: ic05xbf, ic06eaf, if02ebf, iev, astat, iflag, Ntor
    REAL                 :: tol, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), wk(1:Lwk), wka(1:Lwka) ! eigenvalues/vectors of tangent;
    CHARACTER            :: job

#ifndef HAVE_NAG
    INTEGER              :: Lwkb
    REAL, allocatable    :: eigv(:,:), wkb(:), wkc(:)
#endif

    iga00aa = ifail ; nbfield(0:1) = 0

    laxis%Lallocated = 0 ; laxis%error = one ; laxis%residue = zero ; laxis%iota = zero
    
    laxis%tangent(1,1) = one ; laxis%tangent(1,2) = zero ; laxis%tangent(2,1) = zero ; laxis%tangent(2,2) = one
    
    laxis%wr(1:2) = zero ; laxis%wi(1:2) = zero ; laxis%vr(1:2,1:2) = zero ; laxis%vi(1:2,1:2) = zero
    
    CHECKINPUT(ga00aa, laxis%Nfp   .le.0    , 9999 )
    CHECKINPUT(ga00aa, laxis%Ntor  .lt.0    , 9999 )
    CHECKINPUT(ga00aa, laxis%R     .le.zero , 9999 )
    CHECKINPUT(ga00aa, laxis%tol   .le.zero , 9999 )
    CHECKINPUT(ga00aa, laxis%odetol.le.zero , 9999 )
    CHECKINPUT(ga00aa, laxis%maxits.le.0    , 9999 )

    axis%Nfp    = laxis%Nfp    ! required in ga00ab to define integration endpoint ; 19 Jun 15;
    axis%Ntor   = laxis%Ntor   ! required in ga00ba to define intermediate output  ; 19 Jun 15;
    Ntor        = laxis%Ntor   ! shorthand                                         ; 19 Jun 15;
    axis%R      = laxis%R      ! required in ??                                    ; 19 Jun 15;
    axis%Z      = laxis%Z      ! required in ??                                    ; 19 Jun 15;
    axis%tol    = laxis%tol    ! required in ga00ab for user-termination           ; 19 Jun 15;
    axis%odetol = laxis%odetol ! required in ga00ab for o.d.e. integration         ; 19 Jun 15;
    axis%maxits = laxis%maxits ! required in ga00ab for user-termination           ; 19 Jun 15;
    
    fNtor = 4 * max(Ntor,1) ! shorthand; 02 Mar 15;
    
    SALLOCATE(axis%Ri,(0:fNtor),zero)
    SALLOCATE(axis%Zi,(0:fNtor),zero)
    
    if( allocated(laxis%Ri) ) deallocate(laxis%Ri)
    if( allocated(laxis%Zi) ) deallocate(laxis%Zi)
    
    SALLOCATE(laxis%Ri,(0:fNtor),zero)
    SALLOCATE(laxis%Zi,(0:fNtor),zero)
    
    if( allocated(laxis%Rnc) ) deallocate(laxis%Rnc)
    if( allocated(laxis%Rns) ) deallocate(laxis%Rns)
    if( allocated(laxis%Znc) ) deallocate(laxis%Znc)
    if( allocated(laxis%Zns) ) deallocate(laxis%Zns)
    
    SALLOCATE(laxis%Rnc,(0:Ntor),zero)
    SALLOCATE(laxis%Rns,(0:Ntor),zero)
    SALLOCATE(laxis%Znc,(0:Ntor),zero)
    SALLOCATE(laxis%Zns,(0:Ntor),zero)
    
    laxis%Lallocated = 1 ! laxis%Ri, . . , laxis%RRnc, . . have been allocated; 02 Jun 15;
    
    axis%its = -1 ; xx(1:2) = (/ axis%R, axis%Z /) ; ff(1:NN) = one ; axis%rzf(0:2,0:31) = zero
    
    select case( ifail )
    case( -3: )

#ifdef HAVE_NAG
     ic05xbf = 1
     if( iga00aa.le.-3 ) write(0,'("ga00aa :         "2x" : calling C05PBF, which calls ga00ab ;")')
     call C05PBF( ga00ab, NN, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), Ldf, laxis%tol, wk(1:Lwk), Lwk, ic05xbf ) ! NAG; 13 Oct 15;
#else 
     if( iga00aa.le.-3 ) write(0,'("ga00aa :         "2x" : calling HYBRJ1, which calls ga00ab ;")')
     call HYBRJ1( ga00ab, NN, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), Ldf, laxis%tol, ic05xbf, wk(1:Lwk), Lwk )
     if (ic05xbf==1) ic05xbf=0 ! because the MINPACK error code is stupid
#endif

    case( -4  )

#ifdef HAVE_NAG    
     ic05xbf = 1
     if( iga00aa.le.-3 ) write(0,'("ga00aa :         "2x" : calling C05NBF, which calls ga00ac ;")')
     call C05NBF( ga00ac, NN, xx(1:NN), ff(1:NN), laxis%tol, wk(1:Lwk), Lwk, ic05xbf ) ! NAG; 13 Oct 15;
#else 
     if( iga00aa.le.-3 ) write(0,'("ga00aa :         "2x" : calling HYBRJ1, which calls ga00ac ;")')
     call HYBRJ1( ga00ac, NN, xx(1:NN), ff(1:NN), laxis%tol, ic05xbf, wk(1:Lwk), Lwk )
     if (ic05xbf==1) ic05xbf=0 ! because the MINPACK error code is stupid
#endif

    end select

    if( ic05xbf.le.-2 ) then ; xx(1:2) = lxx(1:2) ; ff(1:2) = lff(1:2)
    endif

    laxis%error = sqrt( sum(ff(1:NN)**2) ) ; laxis%its = axis%its ; laxis%rzf(0:2,0:31) = axis%rzf(0:2,0:31)
    
    if( iga00aa.le.-3) then
     select case( ic05xbf ) ! 02 Jun 15;                                012345678901234567890
     case( -3 )   ; write(0,1000) ic05xbf, axis%its, xx(1:2), ff(1:2), "exceeded maxits ;  "
     case( -2 )   ; write(0,1000) ic05xbf, axis%its, xx(1:2), ff(1:2), "user accepts ;     "
     case( -1 )   ; write(0,1000) ic05xbf, axis%its, xx(1:2), ff(1:2), "integration error ;"
     case(  0 )   ; write(0,1000) ic05xbf, axis%its, xx(1:2), ff(1:2), "success ;          "
     case(  1 )   ; write(0,1000) ic05xbf, axis%its, xx(1:2), ff(1:2), "input error ;      "
     case(  2 )   ; write(0,1000) ic05xbf, axis%its, xx(1:2), ff(1:2), "consider restart ; "
     case(  3 )   ; write(0,1000) ic05xbf, axis%its, xx(1:2), ff(1:2), "xtol too small ;   "
     case(  4 )   ; write(0,1000) ic05xbf, axis%its, xx(1:2), ff(1:2), "bad progress ;     "
     end select
    endif

1000 format("ga00aa : ic05xbf="i2" ; its="i3" ; (R,Z)="2es15.7" ; F(R,Z)="2es12.4" ; "a19)

    if( ic05xbf.eq.-3 .and. laxis%error.lt.laxis%tol*ten ) ic05xbf = 0
    if( ic05xbf.eq.-2                                    ) ic05xbf = 0
    if( ic05xbf.ge. 2 .and. laxis%error.lt.laxis%tol*ten ) ic05xbf = 0
    if( ic05xbf.ne. 0                                    ) then ; iga00aa = 2 ; goto 9999
    endif
    
    select case( ifail )
    case( -3: )
     ! make sure axis and tangent map are constructed at solution;
     iflag = 2 ; call ga00ab( NN, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), Ldf, iflag )
    case( -4  )
     ! make sure axis and tangent map are constructed at solution;
     iflag = 2 ; call ga00ac( NN, xx(1:NN), ff(1:NN),                      iflag )
    end select
    
    laxis%Ri(0:fNtor) = axis%Ri(0:fNtor)
    laxis%Zi(0:fNtor) = axis%Zi(0:fNtor)
    
    laxis%error = sqrt( sum(ff(1:NN)**2) )
    
    laxis%R = xx(1) ; laxis%Z = xx(2) ; laxis%tangent(1:2,1:2) = axis%tangent(1:2,1:2) ; laxis%residue = axis%residue
    
    if( laxis%residue.gt.zero ) laxis%iota = axis%iota

#ifdef HAVE_NAG
    ic06eaf = 1
    call C06EAF( axis%Ri(0:fNtor-1), fNtor, ic06eaf ) ! this is destructive; 02 Mar 15; ! NAG; 13 Oct 15;
    if( ic06eaf.ne.0 ) then ; iga00aa = 3 ; goto 9999
    endif
     
    ic06eaf = 1
    call C06EAF( axis%Zi(0:fNtor-1), fNtor, ic06eaf ) ! this is destructive; 02 Mar 15; ! NAG; 13 Oct 15;
#else
    Lwkb = fNtor + INT(LOG(real(fNtor))/LOG(2.)) + 4
    SALLOCATE(wkb,(1:Lwkb),zero)
    call rfft1i( fNtor, wkb, Lwkb, ic06eaf )
    if( ic06eaf.ne.0 ) then ; iga00aa = 3 ; goto 9999
    endif

    SALLOCATE(wkc,(1:fNtor),zero)
    call rfft1f( fNtor, 1, axis%Ri(0:fNtor-1), fNtor, wkb, Lwkb, wkc, fNtor, ic06eaf ) ! this is destructive
    if( ic06eaf.ne.0 ) then ; iga00aa = 3 ; goto 9999
    endif

    call rfft1f( fNtor, 1, axis%Zi(0:fNtor-1), fNtor, wkb, Lwkb, wkc, fNtor, ic06eaf ) ! this is destructive
    DALLOCATE(wkb,wkc)
#endif

    if( ic06eaf.ne.0 ) then ; iga00aa = 3 ; goto 9999
    endif

    axis%Ri(0:fNtor-1) = axis%Ri(0:fNtor-1) / sqrt(one*fNtor) ; axis%Ri(1:fNtor-1) = axis%Ri(1:fNtor-1) * two
    axis%Zi(0:fNtor-1) = axis%Zi(0:fNtor-1) / sqrt(one*fNtor) ; axis%Zi(1:fNtor-1) = axis%Zi(1:fNtor-1) * two

    laxis%Rnc(0:Ntor) =   axis%Ri(      0:      Ntor   )
    laxis%Rns(1:Ntor) = - axis%Ri(fNtor-1:fNtor-Ntor:-1)
    laxis%Znc(0:Ntor) =   axis%Zi(      0:      Ntor   )
    laxis%Zns(1:Ntor) = - axis%Zi(fNtor-1:fNtor-Ntor:-1)
     
    DALLOCATE(axis%Ri)
    DALLOCATE(axis%Zi)
    
    job = 'V' ! construct eigenvalues etc of tangent map at magnetic axis; 02 Mar 15;

#ifdef HAVE_NAG
    if02ebf = 1
    call F02EBF( job, NM, axis%tangent(1:Lda,1:NM), Lda, & ! NAG; 13 Oct 15;
                 laxis%wr(1:NM), laxis%wi(1:NM), laxis%vr(1:Ldvr,1:NM), Ldvr, laxis%vi(1:Ldvi,1:NM), Ldvi, &
                 wka(1:Lwka), Lwka, if02ebf )
#else
    SALLOCATE(eigv,(1:Lda,1:NM),zero)
    call DGEEV( 'N', job, NM, axis%tangent(1:Lda,1:NM), Lda, laxis%wr(1:NM), laxis%wi(1:NM), &
                eigv, Ldvi, eigv, Ldvi, &
                wka(1:Lwka), Lwka, if02ebf )
    laxis%vr(1,:)=eigv(:,1)
    laxis%vr(2,:)=eigv(:,2)
    laxis%vi=0. ! Assume real eigenvectors
    DALLOCATE(eigv)
#endif
    
    if( if02ebf.ne.0 ) then ; iga00aa = 4 ; goto 9999
    endif

    if( laxis%residue.gt.zero ) laxis%iota = atan2( laxis%wi(1), laxis%wr(1) ) / (pi2/laxis%Nfp)
    
    iga00aa = 0 ! success; 02 Jun 15;
    
9999 continue

    laxis%nbfield(0:1) = nbfield(0:1)

    if( ifail.lt. 0 ) write(0,9000) iga00aa, laxis%R, laxis%Z, laxis%its, laxis%error, laxis%residue, laxis%iota, &
                                    nbfield(0)+nbfield(1)

    if( ifail.le.-1 ) write(0,9001)          "Rnc", laxis%Rnc(0:Ntor)
    if( ifail.le.-1 ) write(0,9002)          "Rns", laxis%Rns(1:Ntor)
    if( ifail.le.-1 ) write(0,9001)          "Znc", laxis%Znc(0:Ntor)
    if( ifail.le.-1 ) write(0,9002)          "Zns", laxis%Zns(1:Ntor)
    
    if( ifail.le.-2) then
     do iev = 1, 2 ; write(0,9003)            laxis%wr(iev), laxis%wi(iev), laxis%vr(1,iev), laxis%vi(1,iev), &
                                              laxis%vr(2,iev), laxis%vi(2,iev)
     enddo
    endif

9000 format("ga00aa : ifail ="i3" : (R ,Z )=("f14.10" ,"f15.10" ) ; its="i3" ; err="es9.2" ; residue="f18.13 &
            " ; iota="f18.13" ; nbfield="i11" ;")
9001 format("ga00aa :        "3x" : "a3"="    99es14.06)
9002 format("ga00aa :        "3x" : "a3"="14x,99es14.06)
9003 format("ga00aa :        "3x" : evalue =",es13.05," +",es13.05" i ; evector="es13.05" +"es13.05" i,"es13.05 &
            " +"es13.05" i ;")

    ifail = iga00aa
    
    return
    
  end subroutine ga00aa
  
    
  subroutine ga00ab( NN, xx, ff, df, Ldf, iflag )
    
    implicit none
    
    INTEGER, intent(in)    :: NN, Ldf
    REAL                   :: xx(1:NN), ff(1:NN), df(1:Ldf,1:NN)
    
    INTEGER, intent(inout) :: iflag
    
    INTEGER, parameter     :: Node = 6, Lwk = 20 * Node
    
    INTEGER                :: id02bjf
    REAL                   :: RZ(1:Node), phistart, phiend
    CHARACTER              :: relabs

#ifdef HAVE_NAG
    REAL                   :: wk(1:Lwk)
    external               :: D02BJW
#else
    INTEGER, parameter     :: liw=20, lrw=20+16*Node
    INTEGER                :: iwork(liw)
    REAL                   :: rwork(lrw)
    INTEGER                :: iopt,istate,itask,itol,mf,ii
    REAL                   :: atol,rtol,phic,phie
#endif
        
    relabs = 'D'
    
    phistart = zero ; phiend = phistart + pi2/axis%Nfp ! integration endpoints ; 05 Mar 14;
    
    ! initial guess, intialize tangent map integration;
    RZ(1:Node) = (/ xx(1), xx(2),  one, zero, zero, one /) ; lxx(1:2) = xx(1:2)
    
    izeta = 0 ! this counter is incremented in ga00ba; 31 Jul 13;
    
    select case( iflag )
    case( 1 )
      itangent = 0 ! derivatives (i.e. tangent map) is not required; 21 Mar 15; this is passed through to user-supplied bfield;
    case( 2 )
      itangent = 1 ! derivatives (i.e. tangent map) is     required; 21 Mar 15; this is passed through to user-supplied bfield;
    end select
    
    Lbfieldok = .true.

#ifdef HAVE_NAG
    id02bjf = 1  ! NAG ode integration;
    call D02BJF( phistart, phiend, Node, RZ(1:Node), bf00aa, axis%odetol, relabs, ga00ba, D02BJW, wk(1:Lwk), id02bjf )
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=transform%odetol
!    atol :  absolute tolerance
     atol=transform%odetol
!    initializations for loop
     phic=phistart
     phie=phic
     call ga00ba(phie,RZ)
     do ii=1,fNtor
       call lsode(bf00aa,(/Node/),RZ,phic,phie, &
                  itol,(/rtol/),(/atol/),itask,istate,iopt, &
                  rwork,lrw,iwork,liw,du00aa,mf)
       phic=phie
       call ga00ba(phie,RZ)
     enddo
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif
    
    if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
    
    if( id02bjf.ne.0 ) then

     if( iga00aa.le.-3 ) write(0,'("ga00ab : id02bjf="i2" ; its="i3" ; (R,Z)="2es15.7" -> (R,Z)="2es12.4" !")') &
                               id02bjf, axis%its, xx(1:2), RZ(1:2)

     iflag = -1
     
    else
     
     select case( iflag ) 
     case( 1 ) ; ff(1  ) = RZ(1) - xx(1)                         ! must return function;  5 Jun 13;
      ;        ; ff(2  ) = RZ(2) - xx(2)
     case( 2 ) ; df(1,1) = RZ(3) - one   ; df(1,2) = RZ(4)       ! must return Jacobian;  5 Jun 13;
      ;        ; df(2,1) = RZ(5)         ; df(2,2) = RZ(6) - one
     end select
     
     axis%tangent(1,1) = RZ(3) ; axis%tangent(1,2) = RZ(4)
     axis%tangent(2,1) = RZ(5) ; axis%tangent(2,2) = RZ(6) 
     
     axis%residue = ( two - ( RZ(3) + RZ(6) ) ) / four
    !determinant = RZ(3)*RZ(6) - RZ(4)*RZ(5) ! this should be equal to unity at fixed points; 21 Mar 15;
     
     if( iflag.eq.1 ) then

      axis%its = axis%its + 1

      if( axis%its.lt.32 ) axis%rzf(0:2,axis%its) = (/ xx(1), xx(2), sqrt(sum(ff(1:2)**2)) /) ! record iterations; 19 Jun 15;

      lff(1:2) = ff(1:2)

      if( iga00aa.le.-3 ) write(0,'("ga00ab : id02bjf="i2" ; its="i3" ; (R,Z)="2es15.7" ; F(R,Z)="2es12.4" ;")') &
                                id02bjf, axis%its, xx(1:2), ff(1:2)

      if( sqrt( sum(ff(1:2)**2)).lt.axis%tol    ) iflag = -2 ! user termination; accuracy satisfied;      

      if( axis%its              .ge.axis%maxits ) iflag = -3

     endif ! end of if( iflag.eq.1) ; 19 Jun 15;
     
    endif
    
    return
    
  end subroutine ga00ab
  
    
  subroutine ga00ac( NN, xx, ff,          iflag )
    
    implicit none
    
    INTEGER, intent(in)    :: NN
    REAL                   :: xx(1:NN), ff(1:NN)
    
    INTEGER, intent(inout) :: iflag
    
    INTEGER, parameter     :: Node = 6, Lwk = 20 * Node
    
    INTEGER                :: id02bjf
    REAL                   :: RZ(1:Node), phistart, phiend
    CHARACTER              :: relabs
    
#ifdef HAVE_NAG
    REAL                   :: wk(1:Lwk)
    external               :: D02BJW
#else
    INTEGER, parameter     :: liw=20, lrw=20+16*Node
    INTEGER                :: iwork(liw)
    REAL                   :: rwork(lrw)
    INTEGER                :: iopt,istate,itask,itol,mf,ii
    REAL                   :: atol,rtol,phic,phie
#endif

    relabs = 'D'
    
    phistart = zero ; phiend = phistart + pi2/axis%Nfp ! integration endpoints ; 05 Mar 14;
    
    ! initial guess, intialize tangent map integration;
    RZ(1:Node) = (/ xx(1), xx(2), one, zero, zero, one /) ; lxx(1:2) = xx(1:2) 
    
    izeta = 0 ! this counter is incremented in ga00ba; 31 Jul 13;
    
   !select case( iflag )
   !case( 1 )
   !  itangent = 0 ! derivatives (i.e. tangent map) not required; 21 Mar 15; this is passed through to user-supplied bfield;
   !case( 2 )
   !  itangent = 1 ! derivatives (i.e. tangent map)     required; 21 Mar 15; this is passed through to user-supplied bfield;
   !end select
    
    itangent = 0 ! derivatives (i.e. tangent map) not required; 21 Mar 15; this is passed through to user-supplied bfield;

    Lbfieldok = .true.

#ifdef HAVE_NAG
    id02bjf = 1  ! NAG ode integration;
    call D02BJF( phistart, phiend, Node, RZ(1:Node), bf00aa, axis%odetol, relabs, ga00ba, D02BJW, wk(1:Lwk), id02bjf )
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=axis%odetol
!    atol :  absolute tolerance
     atol=axis%odetol
!    initializations for loop
     phic=phistart
     phie=phic
     call ga00ba(phie,RZ)
     do ii=1,fNtor
       call lsode(bf00aa,(/Node/),RZ,phic,phie, &
                  itol,(/rtol/),(/atol/),itask,istate,iopt, &
                  rwork,lrw,iwork,liw,du00aa,mf)
       phic=phie
       call ga00ba(phie,RZ)
     enddo
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif      
    
    if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
    
    if( id02bjf.ne.0 ) then

     if( iga00aa.le.-3 ) write(0,'("ga00ac : id02bjf="i2" ; its="i3" ; (R,Z)="2es15.7" -> (R,Z)="2es12.4" !")') &
                               id02bjf, axis%its, xx(1:2), RZ(1:2)

     iflag = -1
     
    else
     
    !select case( iflag ) 
      ;        ; ff(1  ) = RZ(1) - xx(1)                         ! must return function;  5 Jun 13;
      ;        ; ff(2  ) = RZ(2) - xx(2)
    !case( 2 ) ; df(1,1) = RZ(3) - one   ; df(1,2) = RZ(4)       ! must return Jacobian;  5 Jun 13;
    ! ;        ; df(2,1) = RZ(5)         ; df(2,2) = RZ(6) - one
    !end select
     
     axis%tangent(1,1) =  one  ; axis%tangent(1,2) = zero
     axis%tangent(2,1) = zero  ; axis%tangent(2,2) =  one  
     
     axis%residue = zero 
    !determinant = RZ(3)*RZ(6) - RZ(4)*RZ(5) ! this should be equal to unity at fixed points; 21 Mar 15;
     
    !if( iflag.eq.1 ) then

      axis%its = axis%its + 1

      if( axis%its.lt.32 ) axis%rzf(0:2,axis%its) = (/ xx(1), xx(2), sqrt(sum(ff(1:2)**2)) /) ! record iterations; 19 Jun 15;

      lff(1:2) = ff(1:2)

      if( iga00aa.le.-3 ) write(0,'("ga00ac : id02bjf="i2" ; its="i3" ; (R,Z)="2es15.7" ; F(R,Z)="2es12.4" ;")') &
                                id02bjf, axis%its, xx(1:2), ff(1:2)

      if( sqrt( sum(ff(1:2)**2)).lt.axis%tol    ) iflag = -2 ! user termination; accuracy satisfied;      

      if( axis%its              .ge.axis%maxits ) iflag = -3

    !endif ! end of if( iflag.eq.1) ; 19 Jun 15;
     
    endif
    
    return
    
  end subroutine ga00ac

  
  subroutine ga00ba( zeta, RZ )
    
    implicit none
    
    INTEGER, parameter :: Node = 6 ! 01 Jun 15;
    
    REAL               :: zeta, RZ(1:Node)
    
    axis%Ri(izeta) = RZ(1) ! save information along magnetic axis in preparation for Fourier decomposition; 05 Mar 14;
    axis%Zi(izeta) = RZ(2) ! save information along magnetic axis in preparation for Fourier decomposition; 05 Mar 14;
    
    izeta = izeta + 1
    
    zeta = izeta * ( pi2/axis%Nfp ) / fNtor
    
    return
    
  end subroutine ga00ba

!latex \subroutine{ho00aa}{find the homoclinic points (of the stable/unstable manifold);}
!latex \bi
!latex \item[1.] This subroutine performs two iterative searches.
!latex \item[* ] The first search is for the unstable fixed point, $\bar {\bf x}$. Iterative fieldline tracing methods are used,
!latex           and the numerical method is identical to that used in \verb+ga00aa+, see \Sec{ga00aa}.
!latex \item[* ] The second search is to find the homoclinic points.
!latex           Homoclinic points are defined as those points that approach the unstable fixed point forwards in time {\em and} backwards in time,
!latex           i.e. $T^n({\bf x})\rightarrow \bar {\bf x}$ as $n\rightarrow \pm \infty$, 
!latex           where $T({\bf x})$ is the image of ${\bf x}$ under the \Poincare map, and $T^{-1}({\bf x})$ is the pre-image, and where
!latex           time is analogous to the toroidal angle.
!latex \item[* ] The tangent mapping at the fixed point is constructed. 
!latex           The unstable manifold is identified by the unstable eigenvector, ${\bf v}_u$,
!latex           which is that eigenvector with a real eigenvalue, $\lambda_u$, with magnitude greater than unity.
!latex           The   stable manifold is identified by the   stable eigenvector, ${\bf v}_s$, 
!latex           which is that eigenvector with a real eigenvalue, $\lambda_s$, with magnitude less    than unity. 
!latex           Note that $\lambda_u \lambda_s = 1$.
!latex \item[* ] The numerical task is then to find $(d_u,d_s)$ such that 
!latex           \be {\bf f}(d_u,d_s) \equiv T^{+i}(\bar {\bf x} + d_u {\bf v}_u) - T^{-j}(\bar {\bf x} + d_s {\bf v}_s) = 0, \label{eq:homoclinic}
!latex           \ee
!latex           where the $(d_u,d_s)$ must be sufficiently small so that the linear approximation is valid
!latex           (required as the eigenvectors of the tangent map are used to identify the stable and unstable manifolds).
!latex \item[* ] There are a countable infinity of homoclinic points
!latex           (and they come in families, and they are increasingly close together near the unstable fixed point).
!latex           The integers $i$ and $j$ are used to identify which homoclinic points are located
!latex           and are determined as part of the calculation. The solution for $(d_u,d_s)$ depends on the $i$ and $j$.
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : homoclinictangle, ho00aa+ \\ \\
!latex           \verb+type(homoclinictangle) :: tangle+ \\ \\ in their source that calls \verb+ho00aa+,
!latex           where \verb+tangle+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+tangle+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+tangle%Nfp              : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, 
!latex \item[ii.] e.g. \verb+Nfp=1+;
!latex \ei
!latex \item[  ] \verb+tangle%R               : real    ;+
!latex \bi
!latex \item[i.] guess for the $R$ location of the unstable fixed point on the $\phi=0$ plane;
!latex \ei
!latex \item[  ] \verb+tangle%Z               : real    ;+
!latex \bi
!latex \item[i.] guess for the $Z$ location of the unstable on the $\phi=0$ plane;
!latex \ei
!latex \item[  ] \verb+tangle%maxits           : integer ;+
!latex \bi
!latex \item[i.] max. iterations allowed in search;
!latex \item[ii.] e.g. \verb+maxits=16+;
!latex \ei
!latex \item[  ] \verb+tangle%xtol              : real    ;+
!latex \bi
!latex \item[i.] required accuracy to which the position of the unstable fixed point on the $\p=0$ plane is required
!latex \item[ii.] e.g. \verb+xtol=1.0e-06+;
!latex \ei
!latex \item[  ] \verb+tangle%odetol           : real    ;+
!latex \bi
!latex \item[i.] o.d.e. integration tolerance; 
!latex \item[ii.] e.g. \verb+odetol=1.0e-08+;
!latex \ei
!latex \item[  ] \verb+tangle%dU              : real    ;+
!latex \bi
!latex \item[i.] guess for the displacement along the unstable manifold of a homoclinic point;
!latex \item[ii.] e.g. \verb+dU+ = $10^{-2}$;
!latex \ei
!latex \item[  ] \verb+tangle%dS              : real    ;+
!latex \bi
!latex \item[i.] guess for the displacement along the   stable manifold of a homoclinic point;
!latex \item[ii.] e.g. \verb+dS+ = $10^{-2}$;
!latex \ei
!latex \item[  ] \verb+tangle%htol              : real    ;+
!latex \bi
!latex \item[i.] required accuracy to which the homoclinic point is required;
!latex \item[ii.] e.g. \verb+htol=1.0e-05+;
!latex \ei
!latex \item[  ] \verb+tangle%ltol              : real    ;+
!latex \bi
!latex \item[i.] required accuracy to which the linear approximation is required;
!latex \item[ii.] e.g. \verb+ltol=1.0e-03+;
!latex \ei
!latex \item[  ] \verb+ifail                      : integer ;+
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call ho00aa( tangle, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+tangle%R                : real    ;+
!latex \bi
!latex \item[i.] updated;
!latex \ei
!latex \item[  ] \verb+tangle%Z                : real    ;+
!latex \bi
!latex \item[i.] updated;
!latex \ei
!latex \item[  ] \verb+tangle%tangent(1:2,1:2) : real    ;+ 
!latex \bi
!latex \item[i.]  the tangent mapping at unstable fixed point;
!latex \item[ii.] if the eigenvalues of the tangent map are imaginary, e.g. $\lambda\equiv\alpha+\beta i$, 
!latex            then the rotational-transform on tangle satisfies $\tan(||\iotabar||)=\beta/\alpha$, where $||\iotabar||\equiv\iotabar \mod 2\pi$.
!latex \item[iii.] if the eigenvalues of the tangent map are real, then the eigenvalues give the direction of the stable and unstable manifolds.
!latex \ei
!latex \item[  ] \verb+tangle%wr(1:2)          : real    ;+ 
!latex \item[  ] \verb+tangle%wi(1:2)          : real    ;+ 
!latex \item[  ] \verb+tangle%vr(1:2,1:2)      : real    ;+ 
!latex \item[  ] \verb+tangle%vi(1:2,1:2)      : real    ;+ 
!latex \bi
!latex \item[i.]  the eigenvalues and eigenvectors of the tangent mapping at the fixed point;
!latex \ei
!latex \item[  ] \verb+tangle%xerror            : real    ;+
!latex \bi
!latex \item[i.] the error, $\sqrt{\Delta R^2 + \Delta Z^2}$, where $\Delta R \equiv R(\Delta\p)-R_X$ and $\Delta Z \equiv Z(\Delta\p)-Z_X$.
!latex \ei
!latex \item[  ] \verb+tangle%xits              : integer ;+
!latex \bi
!latex \item[i.] the number of iterations required;
!latex \ei
!latex \item[  ] \verb+tangle%residue          : real    ;+
!latex \bi
!latex \item[i.] Greene's residue of the fixed point; [Greene, \link{dx.doi.org/10.1063/1.524170}{J. Math. Phys. 20, 1183 (1979)}];
!latex \ei
!latex \item[  ] \verb+tangle%hits              : integer ;+
!latex \bi
!latex \item[i.] the number of iterations required to locate the homoclinic points;
!latex \ei
!latex \item[  ] \verb+tangle%herror            : real    ;+
!latex \bi
!latex \item[i.] the error in locating the homoclinic points;
!latex \ei
!latex \item[  ] \verb+tangle%ilobe(1:2)       : integer ;+
!latex \bi
!latex \item[i.] the $i$ and $j$ as defined above, \Eqn{homoclinic};
!latex \ei
!latex \item[  ] \verb+tangle%maxilobe         : integer ;+
!latex \bi
!latex \item[i.] just for convenience; \verb+maxilobe=max(ilobe(1:2))+;
!latex \ei
!latex \item[  ] \verb+tangle%Lallocated       : integer ;+
!latex \bi
!latex \item[i.] if \verb+Lallocated=1+, the \verb+hpoints+ array has           been allocated;
!latex \item[ii.] if \verb+Lallocated=0+, the \verb+hpoints+ array has {\bf not} been allocated;
!latex \ei
!latex \item[  ] \verb+tangle%hpoints(1:2,0:ltangle%maxilobe,1:2) : real    ;+
!latex \bi
!latex \item[i.] the locations of the homoclinic points;
!latex \item[ii.] the homoclinic points on the unstable branch are $R_{i,u}\equiv$ \verb+hpoints(1,i,1)+ and $Z_{i,u}\equiv$ \verb+hpoints(2,i,1)+,
!latex            for $i=0,$ \verb+ilobe(1)+;
!latex \item[iii.] the homoclinic points on the   stable branch are $R_{i,s}\equiv$ \verb+hpoints(1,i,1)+ and $Z_{i,s}\equiv$ \verb+hpoints(2,i,2)+,
!latex            for $i=0,$ \verb+ilobe(2)+;
!latex \item[iv.] \verb+hpoints+ need not be allocated on input; if it is, it is immediately deallocated;
!latex \ei
!latex \item[  ] \verb+tangle%lerror(1:2)       : real    ;+
!latex \bi
!latex \item[i.] the error in the linear approximation;
!latex \ei
!latex \item[  ] \verb+ifail                 : integer ;+
!latex \bi
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex     \item[] \verb+ifail=2+ : the NAG routine \verb+C05PBF+ failed to locate the zero of the `fixed-point' function,
!latex                              perhaps because of a failure in integrating along the fieldlines;
!latex     \item[] \verb+ifail=3+ : the fixed point is not `unstable': the residue must be negative;
!latex     \item[] \verb+ifail=4+ : the NAG routine \verb+F02EBF+ failed to construct the eigenvalues/vectors of the tangent mapping;
!latex     \item[] \verb+ifail=5+ : the fixed point is not `unstable': the eigenvalues must be real;
!latex     \item[] \verb+ifail=6+ : at least one eigenvector has illegal magnitude;                  
!latex     \item[] \verb+ifail=7+ : failed to fieldline integration to determine \verb+ilobe+;
!latex     \item[] \verb+ifail=8+ : \verb+ilobe(1).le.0+ or \verb+ilobe(2).eq.0+ ;
!latex     \item[] \verb+ifail=9+ : failed to locate the homoclinic points;
!latex \ei
!latex \ei
!latex \item[6.] Comments: 
!latex \item[* ] The NAG routine \verb+C05PBF+ is used for the nonlinear root find, and \verb+tol+ is given directly to \verb+C05PBF+.
!latex \item[* ] The NAG routine \verb+D02BJF+ is used for the o.d.e. integration, and \verb+odetol+ is supplied directly to \verb+D02BJF+.
!latex \item[* ] If a good initial guess is given, this number should be small, as Newton methods should converge rapidly; 
!latex           however, if there are multiple magnetic axes (as during a sawtooth event) then the Newton method may encounter problems;
!latex           also, numerical errors in the magnetic field (perhaps $\nabla\cdot{\bf B}$ is not exactly zero)
!latex           can cause the fieldline integration to be inaccurate, 
!latex           and so it may be difficult to find the solution to the desired accuracy.
!latex \item[* ] Please also consider using \verb+ec00aa+ to find the unstable fixed point;
!latex \ei

  subroutine ho00aa( ltangle, ifail )
    
    implicit none
    
    type(homoclinictangle) :: ltangle
    INTEGER                :: ifail
    
    INTEGER, parameter     :: NN = 2, Ldf = NN, Lwk = NN + 13, NM = 2, Lda = NM, Ldvr = NM, Ldvi = NM, Lwka = 4 * NM
    INTEGER, parameter     :: Node = 6, Lwkb = 20 * Node
    
    INTEGER                :: ic05pbf, if02ebf, iev, astat, iflag, ic05nbf, id02bjf, ius
    REAL                   :: xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), dd(1:NN), wk(1:Lwk), wka(1:Lwka), RZ(1:Node,0:1,1:2)
    REAL                   :: phistart, phiend, dUS
    CHARACTER              :: job, relabs
        
#ifdef HAVE_NAG
    REAL                   :: wkb(1:Lwkb)
    external               :: D02BJX, D02BJW
#else
    REAL, allocatable      :: eigv(:,:)
    INTEGER, parameter     :: liw=20, lrw=20+16*Node
    INTEGER                :: iwork(liw)
    REAL                   :: rwork(lrw)
    INTEGER                :: iopt,istate,itask,itol,mf
    REAL                   :: atol,rtol
#endif

    iho00aa = ifail ; nbfield(0:1) = 0

    ltangle%xits = 0 ; ltangle%xerror = one

    ltangle%tangent(1,1) = one ; ltangle%tangent(1,2) = zero ; ltangle%tangent(2,1) = zero ; ltangle%tangent(2,2) = one

    ltangle%residue = zero

    tangle%wr(1:NM) = zero ; tangle%wi(1:NM) = zero ; tangle%vr(1:Ldvr,1:NM) = zero ; tangle%vi(1:Ldvi,1:NM) = zero
     
    ltangle%herror = one

    if( allocated(ltangle%hpoints) ) deallocate(ltangle%hpoints)
    
    ltangle%Lallocated = 0

    ltangle%ilobe(1:2) = 0 ; ltangle%maxilobe = 0

    ltangle%lerror(1:2) = one
    
    if( ltangle%Nfp   .le.0     ) then ; write(0,'("ho00aa : input error ;    Nfp.le.0     ;")') ; iho00aa = 1 ; goto 9999
    endif
    if( ltangle%R     .le.zero  ) then ; write(0,'("ho00aa : input error ;      R.le.zero  ;")') ; iho00aa = 1 ; goto 9999
    endif
    if( ltangle%odetol.le.zero  ) then ; write(0,'("ho00aa : input error ; odetol.le.zero  ;")') ; iho00aa = 1 ; goto 9999
    endif
    if( ltangle%xtol  .le.zero  ) then ; write(0,'("ho00aa : input error ;   xtol.le.zero  ;")') ; iho00aa = 1 ; goto 9999
    endif
    if( ltangle%htol  .le.zero  ) then ; write(0,'("ho00aa : input error ;   htol.le.zero  ;")') ; iho00aa = 1 ; goto 9999
    endif
    if( ltangle%ltol  .le.zero  ) then ; write(0,'("ho00aa : input error ;   ltol.le.zero  ;")') ; iho00aa = 1 ; goto 9999
    endif
    if( ltangle%dU    .le.small ) then ; write(0,'("ho00aa : input error ;     dU.le.small ;")') ; iho00aa = 1 ; goto 9999
    endif
    if( ltangle%dS    .le.small ) then ; write(0,'("ho00aa : input error ;     dS.le.small ;")') ; iho00aa = 1 ; goto 9999
    endif
    if( ltangle%maxits.le.0     ) then ; write(0,'("ho00aa : input error ; maxits.le.0     ;")') ; iho00aa = 1 ; goto 9999
    endif
    
    tangle%Nfp = ltangle%Nfp ; tangle%odetol = ltangle%odetol ; tangle%htol = ltangle%htol ; tangle%maxits = ltangle%maxits
    
    tangle%xits = 0 ; xx(1:2) = (/ ltangle%R, ltangle%Z /)

#ifdef HAVE_NAG
    ic05pbf = 1
    if( iho00aa.le.-3 ) write(0,'("ho00aa :         "2x" : calling C05PBF, which calls ho00ab ;")')
    call C05PBF( ho00ab, NN, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), Ldf, ltangle%xtol, wk(1:Lwk), Lwk, ic05pbf ) ! NAG; 13 Oct 15;
#else 
    if( iho00aa.le.-3 ) write(0,'("ho00aa :         "2x" : calling HYBRJ1, which calls ho00ab ;")')
    call HYBRJ1( ho00ab, NN, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), Ldf, ltangle%xtol, ic05pbf, wk(1:Lwk), Lwk )
    if (ic05pbf==1) ic05pbf=0 ! because the MINPACK error code is stupid
#endif

    if( ic05pbf.le.-2 ) then ; xx(1:2) = lxx(1:2) ; ff(1:2) = lff(1:2)
    endif
   
    ltangle%xerror = sqrt( sum(ff(1:NN)**2) ) ; ltangle%xits = tangle%xits
    
    if( iho00aa.le.-3) then
     select case( ic05pbf ) ! 02 Jun 15;                                012345678901234567890
     case( -3 )   ; write(0,1000) ic05pbf, tangle%xits, xx(1:2), ff(1:2), "exceeded maxits ;  "
     case( -2 )   ; write(0,1000) ic05pbf, tangle%xits, xx(1:2), ff(1:2), "user accepts ;     "
     case( -1 )   ; write(0,1000) ic05pbf, tangle%xits, xx(1:2), ff(1:2), "integration error ;"
     case(  0 )   ; write(0,1000) ic05pbf, tangle%xits, xx(1:2), ff(1:2), "success ;          "
     case(  1 )   ; write(0,1000) ic05pbf, tangle%xits, xx(1:2), ff(1:2), "input error ;      "
     case(  2 )   ; write(0,1000) ic05pbf, tangle%xits, xx(1:2), ff(1:2), "consider restart ; "
     case(  3 )   ; write(0,1000) ic05pbf, tangle%xits, xx(1:2), ff(1:2), "xtol too small ;   "
     case(  4 )   ; write(0,1000) ic05pbf, tangle%xits, xx(1:2), ff(1:2), "bad progress ;     "
     end select
    endif

1000 format("ho00aa : ic05pbf="i2" ; its="i3" ; (R ,Z ) = ("f19.15" ,"f19.15" ) ; F(R ,Z )="2es13.5" ; "a19)

    if( ic05pbf.eq.-3 .and. ltangle%xerror.lt.ltangle%xtol*ten ) ic05pbf = 0 ! override fail; 02 Jun 15;
    if( ic05pbf.eq.-2                                          ) ic05pbf = 0 ! override fail; 02 Jun 15;
    if( ic05pbf.ge. 2 .and. ltangle%xerror.lt.ltangle%xtol*ten ) ic05pbf = 0 ! override fail; 02 Jun 15;
    if( ic05pbf.ne.0                                           ) then ; iho00aa = 2 ; goto 9999
    endif
    
    iflag = 2
    call ho00ab( NN, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), Ldf, iflag ) ! ensure tangent map is constructed at X point; 01 Jun 15;
    
    ltangle%xerror = sqrt(sum(ff(1:NN)**2))

    ltangle%R = xx(1) ; ltangle%Z = xx(2)
    ltangle%tangent(1:Lda,1:NM) = tangle%tangent(1:Lda,1:NM) ; ltangle%residue = tangle%residue
    ;tangle%R = xx(1) ;  tangle%Z = xx(2)
    
    if( ltangle%residue.ge.zero ) then ; iho00aa = 3 ; goto 9999 ! not unstable; 02 Jun 15;
    endif
    
    job = 'V' ! construct eigenvalues etc of tangent map; 02 Mar 15;
    
#ifdef HAVE_NAG
    if02ebf = 1
    call F02EBF( job, NM, tangle%tangent(1:Lda,1:NM), Lda, & ! NAG; 13 Oct 15;
                 tangle%wr(1:NM), tangle%wi(1:NM), tangle%vr(1:Ldvr,1:NM), Ldvr, tangle%vi(1:Ldvi,1:NM), Ldvi, &
                 wka(1:Lwka), Lwka, if02ebf )
#else
    SALLOCATE(eigv,(1:Lda,1:NM),zero)
     call DGEEV( 'N', job, NM, tangle%tangent(1:Lda,1:NM), Lda, tangle%wr(1:NM), tangle%wi(1:NM), &
                 eigv, Ldvr, eigv, Ldvi, wka(1:Lwka), Lwka, if02ebf )
     ltangle%vr(1,:)=eigv(:,1)
     ltangle%vr(2,:)=eigv(:,2)
     ltangle%vi=0. ! Assume real eigenvectors
    DALLOCATE(eigv)
#endif
    
    if( if02ebf.ne.0 ) then ; iho00aa = 4 ; goto 9999
    endif
    
    if( tangle%wr(1) .lt. tangle%wr(2) ) then
      ltangle%wr(       1) = tangle%wr(       2) ; ltangle%wr(       2) = tangle%wr(       1)
      ltangle%wi(       1) = tangle%wi(       2) ; ltangle%wi(       2) = tangle%wi(       1)
      ltangle%vr(1:Ldvr,1) = tangle%vr(1:Ldvr,2) ; ltangle%vr(1:Ldvr,2) = tangle%vr(1:Ldvr,1)
      ltangle%vi(1:Ldvi,1) = tangle%vi(1:Ldvi,2) ; ltangle%vi(1:Ldvi,2) = tangle%vi(1:Ldvi,1)
    else
      ltangle%wr(       1) = tangle%wr(       1) ; ltangle%wr(       2) = tangle%wr(       2)
      ltangle%wi(       1) = tangle%wi(       1) ; ltangle%wi(       2) = tangle%wi(       2)
      ltangle%vr(1:Ldvr,1) = tangle%vr(1:Ldvr,1) ; ltangle%vr(1:Ldvr,2) = tangle%vr(1:Ldvr,2)
      ltangle%vi(1:Ldvi,1) = tangle%vi(1:Ldvi,1) ; ltangle%vi(1:Ldvi,2) = tangle%vi(1:Ldvi,2)     
    endif
    
    if( ltangle%Z * ltangle%vr(2,1) .gt. zero ) ltangle%vr(1:2,1) = - ltangle%vr(1:2,1)
    if( ltangle%Z * ltangle%vr(2,2) .gt. zero ) ltangle%vr(1:2,2) = - ltangle%vr(1:2,2)
    
    if( abs(ltangle%wi(1)).gt.small*ten .or. abs(ltangle%wi(2)).gt.small*ten ) then ; iho00aa = 5 ; goto 9999 ! not unstable; 02 Jun 15;
    endif

    tangle%wr(       1:NM) = ltangle%wr(       1:NM) ! it is tangle, and not ltangle, that is used in dependent routines; 02 Jun 15;
    tangle%wi(       1:NM) = ltangle%wi(       1:NM)
    tangle%vr(1:Ldvr,1:NM) = ltangle%vr(1:Ldvr,1:NM)
    tangle%vi(1:Ldvi,1:NM) = ltangle%vi(1:Ldvi,1:NM)

    if( sqrt(sum(ltangle%vr(1:2,1)**2)).lt.half .or. &
        sqrt(sum(ltangle%vr(1:2,2)**2)).lt.half      ) then ; iho00aa = 6 ; goto 9999
    endif
    
    phistart = zero
    
    do ius = 1, 2 ! determine suitable ilobe; loop over unstable and stable branches; 02 Jun 15;
     
     ;                   ; RZ(1:   2, 0, ius) = (/ ltangle%R, ltangle%Z /)
     if( ius.eq.1 ) then ; RZ(1:   2, 1, ius) = (/ ltangle%R, ltangle%Z /) + ltangle%dU * ltangle%vr(1:2,ius) ! unstable direction; 01 Jun 15;
     else                ; RZ(1:   2, 1, ius) = (/ ltangle%R, ltangle%Z /) + ltangle%dS * ltangle%vr(1:2,ius) !   stable direction; 01 Jun 15;
     endif
     ;                   ; RZ(3:Node, 1, ius) = (/ one, zero, zero, one /) ! initialize tangent map; not actually used; 02 Jun 15
     
     tangle%ilobe(ius) = 0 ! initialize counter; 01 Jun 15;
     
     do ; tangle%ilobe(ius) = tangle%ilobe(ius) + 1
      
      if( ius.eq.1 ) then ; phiend = phistart + pi2/tangle%Nfp !  forwards in `time' ; 01 Jun 15;
      else                ; phiend = phistart - pi2/tangle%Nfp ! backwards in `time' ; 01 Jun 15;
      endif
      
      relabs = 'D'; Lbfieldok = .true. ; itangent = 0
      
#ifdef HAVE_NAG
      id02bjf = 1
      call D02BJF( phistart, phiend, Node, RZ(1:Node,1,ius), bf00aa, ltangle%odetol, relabs, D02BJX, D02BJW, &
                  wkb(1:Lwkb), id02bjf ) ! NAG; 13 Oct 15;
#else
!     set the integrator parameters.
      rwork=0.
      iwork=0
!     istate=1 :  indicates the first lsode call
      istate=1
!     itask=4 :  normal integration with limited over-shoot (set by
!                rwork(1) in the i_lsode loop
      itask=1
!     iopt=1 :  optional inputs are used
      iopt=1
!     rwork(6) :  set maximum lsode-internal step size
      iwork(6)=400000
!     rwork(7) :  set minimum lsode-internal step size
!     iwork(6) :  set maximum lsode-internal steps
!     iwork(7) :  set maximum lsode-internal error messages printed
!     mf=10 :  non-stiff Adams method of integration
      mf=10
!     itol=1 :  indicates absolute tolerance is just a scalar
      itol=1
!     rtol :  relative tolerance
      rtol=ltangle%odetol
!     atol :  absolute tolerance
      atol=ltangle%odetol
!     do integration
      call lsode(bf00aa,(/Node/),RZ(1:Node,1,ius),phistart,phiend, &
                 itol,(/rtol/),(/atol/),itask,istate,iopt, &
                 rwork,lrw,iwork,liw,du00aa,mf)
      id02bjf=istate
      if (istate>0) id02bjf=0
#endif
      
      if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;

      if( id02bjf.ne.0 ) then ; iho00aa = 7 ; goto 9999
      endif
      
      if( (RZ(1,1,ius)-ltangle%R)**2+(RZ(2,1,ius)-ltangle%Z)**2 &
           .lt.(RZ(1,0,ius)-ltangle%R)**2+(RZ(2,0,ius)-ltangle%Z)**2 ) exit
      
      RZ(1:2,0,ius) = RZ(1:2,1,ius) ! update; 02 Jun 15;
      
     enddo ! end of do ; increment tangle%lobe loop; 01 Jun 15;
     
    enddo ! end of do ius; 01 Jun 15;
    
    ltangle%ilobe(1:2) = tangle%ilobe(1:2) + (/ +0, -0 /) ! can include offsets here; 28 Aug 15;
   
    if( ltangle%ilobe(1).le.0 .or. ltangle%ilobe(2).le.0 ) then ; iho00aa = 8 ; goto 9999
    endif
    
8000 continue
    
    tangle%ilobe(1:2) = ltangle%ilobe(1:2) !; tangle%dU = ltangle%dU ; tangle%dS = ltangle%dS
    
    ltangle%maxilobe = maxval(tangle%ilobe(1:2))
    
    SALLOCATE( tangle%hpoints,(1:2,0:ltangle%maxilobe,1:2),zero)
    SALLOCATE(ltangle%hpoints,(1:2,0:ltangle%maxilobe,1:2),zero)
    
    ltangle%Lallocated = 1 ! indicates that tangle%hpoints has been allocated; 02 Jun 15;
    
    do ipip = 1, 1 ! loop over primary intersection points; 19 Nov 14;
     
     select case( ipip )
     case( 1 ) ; dd(1:NN) = (/ ltangle%dU, ltangle%dS /)
     case( 2 ) ; dd(1:NN) = (/ ltangle%dU, ltangle%dS /) * sqrt( ltangle%wr(1:2) )
     end select
     
     xx(1:NN) = log( (/ dd(1), dd(2) /) )
     
     tangle%hits = 0

     if( iho00aa.le.-3 ) write(0,'("ho00aa :         "2x" : calling C05PBF, which calls ho00bb ;")')

#ifdef HAVE_NAG
     ic05pbf = 1 ! NAG; 13 Oct 15;
     call C05PBF( ho00bb, NN, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), Ldf, ltangle%htol, wk(1:Lwk), Lwk, ic05pbf )
#else 
     call HYBRJ1( ho00bb, NN, xx(1:NN), ff(1:NN), df(1:Ldf,1:NN), Ldf, ltangle%htol, ic05pbf, wk(1:Lwk), Lwk )
     if (ic05pbf==1) ic05pbf=0 ! because the MINPACK error code is stupid
#endif
     
     if( ic05pbf.le.-2 ) then ; xx(1:2) = lxx(1:2) ; ff(1:2) = lff(1:2)
     endif

     ltangle%hits = tangle%hits ; ltangle%herror = sqrt(sum(ff(1:NN)**2))
     
     dd(1:NN) = exp( xx(1:NN) )
     
     if( iho00aa.le.-3) then
      select case( ic05pbf ) ! 02 Jun 15;                                  012345678901234567890
      case( -3 )   ; write(0,1010) ic05pbf, tangle%hits, dd(1:2), ff(1:2), "exceeded maxits;   "
      case( -2 )   ; write(0,1010) ic05pbf, tangle%hits, dd(1:2), ff(1:2), "user accepts;      "
      case( -1 )   ; write(0,1010) ic05pbf, tangle%hits, dd(1:2), ff(1:2), "integration error ;"
      case(  0 )   ; write(0,1010) ic05pbf, tangle%hits, dd(1:2), ff(1:2), "success ;          "
      case(  1 )   ; write(0,1010) ic05pbf, tangle%hits, dd(1:2), ff(1:2), "input error ;      "
      case(  2 )   ; write(0,1010) ic05pbf, tangle%hits, dd(1:2), ff(1:2), "consider restart ; "
      case(  3 )   ; write(0,1010) ic05pbf, tangle%hits, dd(1:2), ff(1:2), "xtol too small ;   "
      case(  4 )   ; write(0,1010) ic05pbf, tangle%hits, dd(1:2), ff(1:2), "bad progress ;     "
      end select
     endif
     
1010 format("ho00aa : ic05pbf="i2" ; its="i3" ; (dU,dS) = ("f19.15" ,"f19.15" ) ; F(dU,dS)="2es13.5" ; "a19)
     
     if( ic05pbf.eq.-3 .and. ltangle%herror.lt.ltangle%htol*ten ) ic05pbf = 0 ! overrule error flag; 02 Jun 15;     
     if( ic05pbf.eq.-2                                          ) ic05pbf = 0 ! overrule error flag; 02 Jun 15;     
     if( ic05pbf.ge. 2 .and. ltangle%herror.lt.ltangle%htol*ten ) ic05pbf = 0 ! overrule error flag; 02 Jun 15;     
     if( ic05pbf.ne. 0                                          ) then ; iho00aa = 9 ; goto 9998
     endif
     
     if( ipip.eq.1 ) ltangle%hpoints(1:2,0:ltangle%maxilobe,1:2) = tangle%hpoints(1:2,0:ltangle%maxilobe,1:2)
     
     if( ipip.eq.1 ) then ; ltangle%dU = dd(1) ; ltangle%dS = dd(2) ! update input guesses; 25 Mar 15;
     endif
     
    enddo ! end of do ipip; 19 Nov 14;
    
    tangle%dU = ltangle%dU ; tangle%dS = ltangle%dS ; ltangle%lerror(1:2) = one ! 13 Jun 15;
    
    do ius = 1, 2
     
     ;                   ; phistart = zero
     if( ius.eq.1 ) then ; phiend   = phistart + pi2/tangle%Nfp ; dUS = ltangle%dU / ltangle%wr(ius)
     else                ; phiend   = phistart - pi2/tangle%Nfp ; dUS = ltangle%dS * ltangle%wr(ius)
     endif
     
     RZ(1:   2, 0, ius) = (/ ltangle%R, ltangle%Z /) + dUS * ltangle%vr(1:2,ius)
     RZ(3:Node, 0, ius) = (/ one, zero, zero , one /)
     
     relabs = 'D'; Lbfieldok = .true. ; itangent = 0
     
#ifdef HAVE_NAG
     id02bjf = 1
     call D02BJF( phistart, phiend, Node, RZ(1:Node,0,ius), bf00aa, ltangle%odetol, relabs, D02BJX, D02BJW, &
                  wkb(1:Lwkb), id02bjf ) ! NAG; 13 Oct 15;
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=ltangle%odetol
!    atol :  absolute tolerance
     atol=ltangle%odetol
!    do integration
     call lsode(bf00aa,(/Node/),RZ(1:Node,0,ius),phistart,phiend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,du00aa,mf)
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif
     
     FATAL(ho00aa, id02bjf.ne.0 .or. .not.Lbfieldok, failed to check lin approx -- soft fail) ! 13 Jun 15;
     
    !if( ius.eq.1) then
    !  lerror = (/ ltangle%R, ltangle%Z /) + ltangle%dU * ltangle%wr(ius) * ltangle%vr(1:2,ius) - ltangle%hpoints(1:2,1,ius)
    !else               ;
    !endif

     ltangle%lerror(ius) = sqrt( sum( ( RZ(1:2,0,ius)-ltangle%hpoints(1:2,0,ius) )**2 ) )
     
     if( ltangle%lerror(ius).gt.ltangle%ltol ) then
      if( ius.eq.1 ) then ; ltangle%dU = dUS ; ltangle%ilobe(ius) = ltangle%ilobe(ius) + 1
      else                ; ltangle%dS = dUS ; ltangle%ilobe(ius) = ltangle%ilobe(ius) + 1
      endif
     endif
     
    enddo ! end of do ius; 02 Jun 15;

    if( iho00aa.le.-3 ) then
      write(0,'("ho00aa :         "2x" :     "3x" ; (dU,dS) = ("f19.15" ,"f19.15" ) ; lerror=("es10.2" ,"es10.2" ) ;")') &
        tangle%dU, tangle%dS, ltangle%lerror(1:2)
    endif
    
    if( ltangle%lerror(1).gt.ltangle%ltol .or. ltangle%lerror(2).gt.ltangle%ltol ) then
     deallocate( tangle%hpoints)
     deallocate(ltangle%hpoints)
     ltangle%Lallocated = 0
     goto 8000
    endif

    iho00aa = 0
    
9998 continue
    
    DALLOCATE(tangle%hpoints)
    
9999 continue

    ltangle%nbfield(0:1) = nbfield(0:1)
    
    if( ifail.le. 0 ) write(0,9000) iho00aa, ltangle%R , ltangle%Z , ltangle%xits, ltangle%xerror, ltangle%residue, &
                                    nbfield(0)+nbfield(1)
    if( ifail.le.-1 ) then
     do iev = 1, 2 ; write(0,9002)           ltangle%wr(iev), ltangle%wi(iev), ltangle%vr(1,iev), ltangle%vi(1,iev), &
                                             ltangle%vr(2,iev), ltangle%vi(2,iev)
     enddo
    endif
    if( ifail.le. 0 ) write(0,9001)          ltangle%dU, ltangle%dS, ltangle%hits, ltangle%herror, ltangle%ilobe(1:2), &
                                             ltangle%Lallocated, ltangle%lerror    
    if( ltangle%Lallocated.eq. 1 ) then
     if( ifail.le.-2 ) write(0,9003)         ltangle%hpoints(1:2,ltangle%ilobe(1)-1:0:-1,1)
     if( ifail.le.-2 ) write(0,9003)         ltangle%hpoints(1:2,ltangle%ilobe(2)-1:0:-1,2)
    endif
    
9000 format("ho00aa : ifail ="i3" : (R ,Z )=("f14.10" ,"f15.10" ) ; xits="i3" ; xerror="es10.2" ; residue="f18.13 &
            " ; nbfield="i11" ;")
9001 format("ho00aa :        "3x" : (dU,dS)=("f14.10" ,"f15.10" ) ; hits="i3" ; herror="es10.2" ; ilobe="2i3 &
            " ; Lallocated="i2" ; lerror="2es10.2" ;")
9002 format("ho00aa :        "3x" : evalue =",es13.05," +",es13.05" i ; evector="es13.05" +"es13.05 &
            " i,"es13.05" +"es13.05" i ;")
9003 format("ho00aa :        "3x" : homoclinic (R,Z) =("f7.4","f7.4"),"99("("f5.2","f5.2"),":))
    
    ifail = iho00aa
    
    return
    
  end subroutine ho00aa
  
  subroutine ho00ab( NN, xx, ff, df, Ldf, iflag )
    
    implicit none
    
    INTEGER, intent(in)    :: NN, Ldf
    REAL                   :: xx(1:NN), ff(1:NN), df(1:Ldf,1:NN)
    
    INTEGER, intent(inout) :: iflag
    
    INTEGER, parameter     :: Node = 6, Lwk = 20 * Node ! 01 Jun 15;
    
    INTEGER                :: id02bjf
    REAL                   :: RZ(1:Node), phistart, phiend
    CHARACTER              :: relabs
    
#ifdef HAVE_NAG
    REAL                   :: wk(1:Lwk)
    external               :: D02BJX, D02BJW
#else
    INTEGER, parameter     :: liw=20, lrw=20+16*Node
    INTEGER                :: iwork(liw)
    REAL                   :: rwork(lrw)
    INTEGER                :: iopt,istate,itask,itol,mf
    REAL                   :: atol,rtol
#endif
    
    relabs = 'D'
    
    phistart = zero ; phiend = phistart + pi2/tangle%Nfp ! integration endpoints ; 05 Mar 14;
    
    RZ(1:Node) = (/ xx(1), xx(2),  one, zero, zero, one /) ; lxx(1:2) = xx(1:NN)
    
    select case( iflag )
    case( 1 ) ; itangent = 0 ! derivatives (i.e. tangent map) is not required; 21 Mar 15; this is passed through to user-supplied bfield;
    case( 2 ) ; itangent = 1 ! derivatives (i.e. tangent map) is     required; 21 Mar 15; this is passed through to user-supplied bfield;
    end select
    
    Lbfieldok = .true.

#ifdef HAVE_NAG
    id02bjf = 1 ! NAG ode integration;
    call D02BJF( phistart, phiend, Node, RZ(1:Node), bf00aa, tangle%odetol, relabs, D02BJX, D02BJW, wk(1:Lwk), id02bjf )
#else
!   set the integrator parameters.
    rwork=0.
    iwork=0
!   istate=1 :  indicates the first lsode call
    istate=1
!   itask=4 :  normal integration with limited over-shoot (set by
!              rwork(1) in the i_lsode loop
    itask=1
!   iopt=1 :  optional inputs are used
    iopt=1
!   rwork(6) :  set maximum lsode-internal step size
    iwork(6)=400000
!   rwork(7) :  set minimum lsode-internal step size
!   iwork(6) :  set maximum lsode-internal steps
!   iwork(7) :  set maximum lsode-internal error messages printed
!   mf=10 :  non-stiff Adams method of integration
    mf=10
!   itol=1 :  indicates absolute tolerance is just a scalar
    itol=1
!   rtol :  relative tolerance
    rtol=tangle%odetol
!   atol :  absolute tolerance
    atol=tangle%odetol
!   do integration
    call lsode(bf00aa,(/Node/),RZ,phistart,phiend, &
               itol,(/rtol/),(/atol/),itask,istate,iopt, &
               rwork,lrw,iwork,liw,du00aa,mf)
    id02bjf=istate
    if (istate>0) id02bjf=0
#endif
    
    if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
  
    if( id02bjf.ne.0 ) then

    !if( iho00aa.le.-3 ) write(0,'("ho00ab : id02bjf="i2" ; its="i3" ; (R,Z)="2es13.5" -> (R,Z)="2es13.5" !")') &
    !                      id02bjf, tangle%xits, xx(1:2), RZ(1:2)

     iflag = -1 ! tell NAG that an error has occured; 02 Jun 15;
     
    else
     
     select case( iflag ) 
     case( 1 ) ; ff(1  ) = RZ(1) - xx(1)                            ! must return function;  5 Jun 13;
      ;        ; ff(2  ) = RZ(2) - xx(2)
     case( 2 ) ; df(1,1) = RZ(3) - one   ; df(1,2) = RZ(4)       ! must return Jacobian;  5 Jun 13;
      ;        ; df(2,1) = RZ(5)         ; df(2,2) = RZ(6) - one
     end select
     
     tangle%tangent(1,1) = RZ(3) ; tangle%tangent(1,2) = RZ(4)
     tangle%tangent(2,1) = RZ(5) ; tangle%tangent(2,2) = RZ(6) 

     tangle%residue = ( two - ( RZ(3) + RZ(6) ) ) / four
    !determinant = RZ(3)*RZ(6) - RZ(4)*RZ(5) ! this should be equal to unity at fixed points; 21 Mar 15;
     
     if( iflag.eq.1 ) then
      tangle%xits = tangle%xits + 1 ; lff(1:2) = ff(1:2)
      if( iho00aa.le.-3 ) write(0,1000) id02bjf, tangle%xits, xx(1:2), ff(1:2)
      if( sqrt(sum(ff(1:2)**2 )) .lt. axis%tol      ) iflag = -2
      if( tangle%xits            .ge. tangle%maxits ) iflag = -3
     endif

1000 format("ho00ab : id02bjf="i2" ; its="i3" ; (R ,Z ) = ("f19.15" ,"f19.15" ) ; F(R ,Z )="2es13.5" ; ")
     
    endif
    
    return
    
  end subroutine ho00ab
  
  subroutine ho00bb( NN, xx, ff, df, Ldf, iflag )
    
    implicit none

    INTEGER            :: NN, Ldf, iflag
    REAL               :: xx(1:NN), ff(1:NN), df(1:Ldf,1:NN)
    
    INTEGER, parameter :: Node = 6, Lwk = 20 * Node
    INTEGER            :: ius, ii, id02bjf
    REAL               :: dd(1:NN), RZ(1:Node,1:2), odetol, phistart, phiend
    CHARACTER          :: relabs
    
#ifdef HAVE_NAG
    REAL                   :: wk(1:Lwk)
    external               :: D02BJX, D02BJW
#else
    INTEGER, parameter     :: liw=20, lrw=20+16*Node
    INTEGER                :: iwork(liw)
    REAL                   :: rwork(lrw)
    INTEGER                :: iopt,istate,itask,itol,mf
    REAL                   :: atol,rtol
#endif
    
    dd(1:NN) = exp( xx(1:NN) ) ; lxx(1:2) = xx(1:2)

    odetol = tangle%odetol / ten ; relabs = 'D'
    
    if( iflag.eq.1 ) then ; itangent = 0 ! derivatives are not required; 02 Jun 15;
    else                  ; itangent = 1 ! derivatives are     required; 02 Jun 15;
    endif
    
    do ius = 1, 2
     
     ii = 0 ; phistart = zero
     
     RZ(1:2,ius) = (/ tangle%R, tangle%Z /) + dd(ius) * tangle%vr(1:2,ius) ! displace along unstable direction; 19 Nov 14;
     
     if( ipip.eq.1 ) tangle%hpoints(1:2,ii,ius) = RZ(1:2,ius) ! should eliminate ipip; 02 Jun 15;
     
     RZ(3:Node,ius) = (/ one, zero, zero, one /) ! initialize tangent map; 01 Jun 15;
     
     do ii = 1, tangle%ilobe(ius)-1
      
      if( ius.eq.1 ) then ; phiend = phistart + pi2/tangle%Nfp !  forwards; 19 Nov 14;
      else                ; phiend = phistart - pi2/tangle%Nfp ! backwards; 19 Nov 14;
      endif
      
      Lbfieldok = .true.
      
#ifdef HAVE_NAG
      id02bjf = 1 ! NAG; 13 Oct 15;
      call D02BJF( phistart, phiend, Node, RZ(1:Node,ius), bf00aa, odetol, relabs, D02BJX, D02BJW, wk(1:Lwk), id02bjf )
#else
!     set the integrator parameters.
      rwork=0.
      iwork=0
!     istate=1 :  indicates the first lsode call
      istate=1
!     itask=4 :  normal integration with limited over-shoot (set by
!                rwork(1) in the i_lsode loop
      itask=1
!     iopt=1 :  optional inputs are used
      iopt=1
!     rwork(6) :  set maximum lsode-internal step size
      iwork(6)=400000
!     rwork(7) :  set minimum lsode-internal step size
!     iwork(6) :  set maximum lsode-internal steps
!     iwork(7) :  set maximum lsode-internal error messages printed
!     mf=10 :  non-stiff Adams method of integration
      mf=10
!     itol=1 :  indicates absolute tolerance is just a scalar
      itol=1
!     rtol :  relative tolerance
      rtol=odetol
!     atol :  absolute tolerance
      atol=odetol
!     do integration
      call lsode(bf00aa,(/Node/),RZ(1:Node,ius),phistart,phiend, &
                 itol,(/rtol/),(/atol/),itask,istate,iopt, &
                 rwork,lrw,iwork,liw,du00aa,mf)
      id02bjf=istate
      if (istate>0) id02bjf=0
#endif
      
      if( .not.Lbfieldok .or. id02bjf.ne. 0 ) then ; iflag = -1 ; goto 9999
      endif
      
      if( ipip.eq.1 ) tangle%hpoints(1:2,ii,ius) = RZ(1:2,ius)
      
     enddo ! end of do ii; 19 Nov 14;
     
    enddo ! end of do ius; 19 Nov 14;
    
    select case( iflag )
    case( 1 ) ; ff(1  ) = +   RZ(1,1)                  - RZ(1,2)
     ;        ; ff(2  ) = +   RZ(2,1)                  - RZ(2,2)
    case( 2 ) ; df(1,1) = + ( RZ(3,1) * tangle%vr(1,1) + RZ(4,1) * tangle%vr(2,1) ) * dd(1)
     ;        ; df(2,1) = + ( RZ(5,1) * tangle%vr(1,1) + RZ(6,1) * tangle%vr(2,1) ) * dd(1)
     ;        ; df(1,2) = - ( RZ(3,2) * tangle%vr(1,2) + RZ(4,2) * tangle%vr(2,2) ) * dd(2)
     ;        ; df(2,2) = - ( RZ(5,2) * tangle%vr(1,2) + RZ(6,2) * tangle%vr(2,2) ) * dd(2)
    end select
    
    if( iflag.eq.1 ) then ! only in this case is the "function" calculated; 02 Jun 15;
     tangle%hits = tangle%hits + 1 ; lff(1:2) = ff(1:2)
     if( iho00aa.le.-3 ) write(0,1000) tangle%hits, dd(1:NN), ff(1:NN)
     if( sqrt(sum(ff(1:NN)**2)) .lt. tangle%htol   ) iflag = -2 ! user termination; accuracy satisfied;
     if( tangle%hits            .ge. tangle%maxits ) iflag = -3
    endif
    
1000 format("ho00bb :         "2x" ; its="i3" ; (dU,dS) = ("f19.15" ,"f19.15" ) ; F(dU,dS)="2es13.5" ;")
    
9999 continue

    return
    
  end subroutine ho00bb

!latex \subroutine{ec00aa}{find action extremizing curves using global integration;}
!latex \bi
!latex \item[1.] Find curves that extremize the action integral,
!latex           \be S[{\cal C}] \equiv \int_{\cal C} {\bf A}\cdot d{\bf l}.
!latex           \ee
!latex           Note: for a given vector potential, ${\bf B}=\nabla\times{\bf A}$, 
!latex           the action, $S$, is considered to be a function of the `trial-curve', ${\cal C}$.
!latex \item[* ] From variational calculus, the variation in $S$ due to variations, $\delta d{\bf l}$, in the curve is
!latex           \be \delta S = \int_{\cal C} {\bf B} \times d{\bf l} \cdot \delta d{\bf l},
!latex           \ee
!latex           from which we see that curves that extremize $S$ satisfy ${\bf B}\times d{\bf l}=0$, i.e. the extremal curves are parallel to ${\bf B}$.
!latex \item[* ] This method of detemining magnetic fieldlines is called Lagrangian or global integration.
!latex \item[* ] Working in cylindrical coordinates, an arbitrary trial curve is represented by 
!latex           \be {\bf x}(\phi) = R(\phi) {\bf e}_R(\phi) + Z(\phi) {\bf e}_Z,
!latex           \ee
!latex           where ${\bf e}_R(\phi) \equiv \cos(\phi) {\bf i} + \sin(\phi) {\bf j}$ and ${\bf e}_Z \equiv {\bf k}$.
!latex \item[* ] The infinitesimal change in ${\bf x}$ due to an infinitesimal increase in $\phi$ is
!latex           \be d{\bf l} = \left( \dot R \, {\bf e}_R + {\bf e}_\phi + \dot Z \, {\bf e}_Z \right) d\phi,
!latex           \ee
!latex           where the `dot' denotes the derivative with respect to $\phi$.
!latex \item[* ] It is assumed that the magnetic vector potential is in the form 
!latex           \be {\bf A} \equiv A_R \nabla R + A_\phi \nabla \phi + A_Z \nabla Z.
!latex           \ee
!latex           Interestingly, and fortunately, only the curl of the vector potential will be required.
!latex \item[* ] The `Lagrangian', i.e. the integrand ${\bf A}\cdot d{\bf l} /d \phi$, is $A_R \dot R + A_\phi + A_Z \dot Z$.
!latex
!latex \item[* ] A Fourier representation for $R(\phi)$ and $Z(\phi)$ is employed:
!latex           \be R & \equiv & \sum_{n=0}^{N} \left[ R_{n,c} \cos(n \phi) + R_{n,c} \sin(n \phi) \right], \\
!latex               Z & \equiv & \sum_{n=0}^{N} \left[ Z_{n,c} \cos(n \phi) + Z_{n,c} \sin(n \phi) \right].
!latex           \ee
!latex \item[* ] Extremal curves are curves for which the derivative of $S$ with respect to to the $R_{n,c}$, $R_{n,s}$ etc. are zero. 
!latex           We have, for example, 
!latex           \be \frac{\partial S}{\partial R_{n,c}} = \int_{0}^{2\pi} \left(
!latex           \frac{\partial A_R}{\partial R} \frac{\partial R}{\partial R_{n,c}} \dot R + A_R \frac{\partial}{\partial R_{n,c}} \dot R 
!latex           + \frac{\partial A_\phi}{\partial R} \frac{\partial R}{\partial R_{n,c}} 
!latex           + \frac{\partial A_Z}{\partial R} \frac{\partial R}{\partial R_{n,c}} \dot Z \right) d\phi.
!latex           \ee
!latex \item[* ] Consider the term
!latex           \be A_R \frac{\partial}{\partial R_{n,c}} \dot R = A_R \frac{\partial}{\partial R_{n,c}} \frac{d}{d\phi} R =
!latex           A_R \frac{d}{d\phi} \frac{\partial}{\partial R_{n,c}} R = A_R \frac{d}{d\phi} \cos(n \phi)
!latex           \ee
!latex           Rather than writing $d_\phi \cos(n \phi) = - n \sin(n \phi)$, a very neat trick is to
!latex           to instead use integration-by-parts to write $A_R \, d_\phi \cos(n \phi) \equiv - d_\phi A_R \, \cos(n \phi)$.
!latex           Note that $d_\phi$ is the total derivative with respect to $\phi$, 
!latex           so that $d_\phi A_R \equiv \partial_R A_R \; \dot R + \partial_\phi A_R + \partial_Z A_R \; \dot Z$.
!latex           Several terms cancel, and only the derivatives of ${\bf A}$ are now required, rather than ${\bf A}$ itself; in fact,
!latex           it is the components of $\nabla \times {\bf A}$ that appear!
!latex \item[* ] The derivatives of the action with respect to the parameters defining the curve are
!latex           \be \frac{\partial S}{\partial R_{n,c}} & = & \int_{0}^{2\pi} \cos(n\phi) \left( B^Z - B^\phi \dot Z \right) R \, d\phi, \\
!latex               \frac{\partial S}{\partial R_{n,s}} & = & \int_{0}^{2\pi} \sin(n\phi) \left( B^Z - B^\phi \dot Z \right) R \, d\phi, \\
!latex               \frac{\partial S}{\partial Z_{n,c}} & = & \int_{0}^{2\pi} \cos(n\phi) \left( B^\phi \dot R - B^R \right) R \, d\phi, \\
!latex               \frac{\partial S}{\partial Z_{n,s}} & = & \int_{0}^{2\pi} \sin(n\phi) \left( B^\phi \dot R - B^R \right) R \, d\phi.
!latex           \ee
!latex \item[* ] Various numerical methods may be employed to find $F_R \equiv \left( B^Z - B^\phi \dot Z \right) R = 0$ and 
!latex                                                             $F_Z \equiv \left( B^\phi \dot R - B^R \right) R = 0$.
!latex          %For example, evolving the curve according to $\partial_\tau {\bf x} \equiv - \partial_{\bf x} S$ ensures the action decreases,
!latex          %and so the curve will eventually approach the action-minimizing curve.
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : extremizingcurve, ec00aa+ \\ \\
!latex           \verb+type(extremizingcurve) :: curve+ \\ \\ in their source that calls \verb+ec00aa+.
!latex           The variable name, \verb+curve+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+curve%Nfp              : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, 
!latex \item[ii.] e.g. \verb+Nfp=1+;
!latex \ei
!latex \item[  ] \verb+curve%Ntor            : integer ;+
!latex \bi
!latex \item[i.] the required Fourier resolution in the toroidal direction;
!latex \item[ii.] constraint: \verb+Ntor.ge.0+;
!latex \ei
!latex \item[  ] \verb+curve%Rnc(0:Ntor)     : real    ;+
!latex \item[  ] \verb+curve%Rns(0:Ntor)     : real    ;+
!latex \item[  ] \verb+curve%Znc(0:Ntor)     : real    ;+
!latex \item[  ] \verb+curve%Zns(0:Ntor)     : real    ;+
!latex \bi
!latex \item[i.] an initial guess for the Fourier harmonics of the extremizing curve;
!latex \item[ii.] these are allocatable, and must be allocated {\em before} calling \verb+ec00aa+;
!latex \ei
!latex \item[  ] \verb+curve%etol            : real    ;+
!latex \bi
!latex \item[i.] the accuracy to which the extremal curve is required;
!latex \item[ii.] e.g. \verb+etol+$=10^{-6}$;
!latex \ei
!latex \item[  ] \verb+curve%ftol            : real    ;+
!latex \bi
!latex \item[i.] the accuracy to which the extremal curve is required; only if gradient flow is used;
!latex \item[ii.] e.g. \verb+ftol+$=10^{-3}$;
!latex \ei
!latex \item[  ] \verb+curve%maxits          : integer ;+
!latex \bi
!latex \item[i.] maximum number of iterations allowed;
!latex \ei
!latex \item[  ] \verb+curve%emethod         : integer ;+
!latex \bi
!latex \item[i.] determines the numerical method used to locate extremal curves:
!latex \item[ii.] \verb+emethod+ $=0$ : uses Newton method to find $F_R=0$ and $F_Z=0$;
!latex \ei
!latex \item[  ] \verb+curve%odetol          : real    ;+
!latex \bi
!latex \item[i.] the o.d.e. integration tolerance;
!latex \item[ii.] e.g. \verb+odetol=+$10^{-6}$;
!latex \item[iii.] only used if \verb+emethod=1,2,3,4+;
!latex \ei
!latex \item[  ] \verb+curve%tauend          : real    ;+
!latex \bi
!latex \item[i.] the upper limit on the o.d.e. integration;
!latex \item[ii.] e.g. \verb+tauend=1.00+;
!latex \item[iii.] only used if \verb+emethod=1,2,3,4+;
!latex \ei
!latex \item[  ] \verb+curve%dtau            : real    ;+
!latex \bi
!latex \item[i.] the intermediate output on the o.d.e. integration;
!latex \item[ii.] e.g. \verb+dtau=0.05+;
!latex \item[iii.] only used if \verb+emethod=1,2,3,4+;
!latex \ei
!latex \item[  ] \verb+ifail                      : integer ;+
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call ec00aa( curve, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+curve%Rnc(0:Ntor)     : real    ;+
!latex \item[  ] \verb+curve%Rns(0:Ntor)     : real    ;+
!latex \item[  ] \verb+curve%Znc(0:Ntor)     : real    ;+
!latex \item[  ] \verb+curve%Zns(0:Ntor)     : real    ;+
!latex \bi
!latex \item[i.] updated;
!latex \ei
!latex \item[  ] \verb+curve%its             : integer ;+
!latex \bi
!latex \item[i.] required iterations;
!latex \ei
!latex \item[  ] \verb+curve%err             : real    ;+
!latex \bi
!latex \item[i.] accuracy achieved;
!latex \ei
!latex \item[  ] \verb+ifail                 : integer ;+
!latex \bi
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex     \item[] \verb+ifail=2+ : the routine \verb+C05NBF+ failed to find a solution;
!latex \ei
!latex \ei
!latex \item[6.] Comments:
!latex \item[* ] The NAG routine \verb+C05NBF+ is used;
!latex \ei

  subroutine ec00aa( lcurve, ifail )
    
    implicit none
    
    type(extremizingcurve) :: lcurve
    INTEGER                :: ifail
    
    INTEGER                :: Nfp, astat, ii, jj, Ntor, tN, NXmn, LXmn, Nfmn, Lfmn, Lwk, ic05nbf, ic05pbf
    INTEGER                :: NN, iflag, id02bjf, itau
    REAL                   :: tau, odetol
    REAL,    allocatable   :: XXmn(:), FXmn(:), DXmn(:,:), Xfmn(:), Ffmn(:), Dfmn(:,:), wk(:)
    CHARACTER              :: relabs

#ifdef HAVE_NAG
    external               :: D02BJX, D02BJW
#else              
    INTEGER                :: liw, lrw
    INTEGER, allocatable   :: iwork(:)
    REAL,    allocatable   :: rwork(:)
    INTEGER                :: iopt,istate,itask,itol,mf,jt,ng,jroot
    REAL                   :: atol,rtol,tauc,taue
#endif

    iec00aa = ifail ; nbfield(0:1) = 0
    
    if( allocated(lcurve%iRnc) ) deallocate(lcurve%iRnc)
    if( allocated(lcurve%iRns) ) deallocate(lcurve%iRns)
    if( allocated(lcurve%iZnc) ) deallocate(lcurve%iZnc)
    if( allocated(lcurve%iZns) ) deallocate(lcurve%iZns)
    
    lcurve%Lallocated = 0 ; Nfp = lcurve%Nfp
    
    CHECKINPUT(ec00aa,                       Nfp   .le.0    , 9999 )
    CHECKINPUT(ec00aa,                lcurve%Ntor  .lt.0    , 9999 )
    CHECKINPUT(ec00aa, .not.allocated(lcurve%Rnc)           , 9999 )
    CHECKINPUT(ec00aa, .not.allocated(lcurve%Rns)           , 9999 )
    CHECKINPUT(ec00aa, .not.allocated(lcurve%Znc)           , 9999 )
    CHECKINPUT(ec00aa, .not.allocated(lcurve%Zns)           , 9999 )
    CHECKINPUT(ec00aa,                lcurve%Rnc(0).lt.small, 9999 )
    CHECKINPUT(ec00aa,                lcurve%etol  .lt.small, 9999 )
    CHECKINPUT(ec00aa,                lcurve%maxits.le.    0, 9999 )
    
    Ntor = lcurve%Ntor ; tN = 1+Ntor + Ntor ; NXmn = 2 * tN ; LXmn = NXmn ; Nfmn = tN ; Lfmn = Nfmn
    
    curve%Ntor = Ntor ; curve%Nfp = Nfp ; curve%etol = lcurve%etol
    curve%maxits = lcurve%maxits ; curve%emethod = lcurve%emethod
    
    SALLOCATE(curve%Rnc,(0:Ntor),lcurve%Rnc(0:Ntor)) ! lcurve is local input; curve is global/common;  8 Jul 15;
    SALLOCATE(curve%Rns,(0:Ntor),lcurve%Rns(0:Ntor))
    SALLOCATE(curve%Znc,(0:Ntor),lcurve%Znc(0:Ntor))
    SALLOCATE(curve%Zns,(0:Ntor),lcurve%Zns(0:Ntor))

    NN = 4*max(1,Ntor) ; curve%sN = sqrt(one*NN) ; curve%init = 'I'
    
    SALLOCATE(curve%ct,(1:NN,0:Ntor),zero)
    SALLOCATE(curve%st,(1:NN,0:Ntor),zero)

#ifdef HAVE_NAG
    SALLOCATE(curve%trig,(1:2*NN),zero)
#endif

    SALLOCATE(curve%tt,(1:NN),zero)    
    curve%tt(1:NN) = (/ ( ii, ii = 0, NN-1 ) /) * (pi2/Nfp) / NN
    do jj = 0, Ntor ; curve%ct(1:NN,jj) = cos( jj * Nfp * curve%tt(1:NN) )
     ;              ; curve%st(1:NN,jj) = sin( jj * Nfp * curve%tt(1:NN) )
    enddo
    DALLOCATE(curve%tt)

    SALLOCATE(curve%rr,(1:NN),zero)
    SALLOCATE(curve%rd,(1:NN),zero)
    SALLOCATE(curve%zz,(1:NN),zero)
    SALLOCATE(curve%zd,(1:NN),zero)
    
    select case( lcurve%emethod )
     
    case( 0   )
     
     if( ifail.le.-4 ) write(0,9001)          "Rnc", lcurve%Rnc(0:lcurve%Ntor)
     if( ifail.le.-4 ) write(0,9002)          "Rns", lcurve%Rns(1:lcurve%Ntor)
     if( ifail.le.-4 ) write(0,9001)          "Znc", lcurve%Znc(0:lcurve%Ntor)
     if( ifail.le.-4 ) write(0,9002)          "Zns", lcurve%Zns(1:lcurve%Ntor)
     
     SALLOCATE(XXmn,(1:NXmn),zero)
     ;  ii = 0       ; XXmn(   ii+1) = lcurve%Rnc(ii) ; XXmn(tN+   ii+1) = lcurve%Znc(ii) ! pack independent variable;  8 Jul 15;
     do ii = 1, Ntor ; XXmn(   ii+1) = lcurve%Rnc(ii) ; XXmn(tN+   ii+1) = lcurve%Znc(ii)
      ;              ; XXmn(tN-ii+1) = lcurve%Rns(ii) ; XXmn(tN+tN-ii+1) = lcurve%Zns(ii)
     enddo
     
     SALLOCATE(FXmn,(1:NXmn       ),zero) ! "function"   vector;  8 Jul 15;
     SALLOCATE(DXmn,(1:LXmn,1:NXmn),zero) ! "derivative" vector;  8 Jul 15;
        
     Lwk = NXmn * ( 3 * NXmn + 13 ) / 2
     SALLOCATE(wk,(1:Lwk),zero) ! work array required by NAG;  8 Jul 15;
          
     curve%its = -1

#ifdef HAVE_NAG
     ic05pbf = 1 ! NAG; 13 Oct 15;
     if( iec00aa.le.-3 ) write(0,'("ec00aa :         "2x" : calling C05PBF, which calls ec00ab ;")')
     call C05PBF( ec00ab, NXmn, XXmn(1:NXmn), FXmn(1:NXmn), DXmn(1:LXmn,1:NXmn), LXmn, lcurve%etol, wk(1:Lwk), Lwk, ic05pbf )
#else 
     if( iec00aa.le.-3 ) write(0,'("ec00aa :         "2x" : calling HYBRJ1, which calls ec00ab ;")')
     call HYBRJ1( ec00ab, NXmn, XXmn(1:NXmn), FXmn(1:NXmn), DXmn(1:LXmn,1:NXmn), LXmn, lcurve%etol, ic05pbf, wk(1:Lwk), Lwk )
     if (ic05pbf==1) ic05pbf=0 ! because the MINPACK error code is stupid
#endif

     lcurve%its = curve%its ; lcurve%err = sqrt( sum( FXmn(1:NXmn)**2 ) ) ! need to pass solution for when iflag=-2 or iflag=-1 ; 17 Jun 15;
     
     DALLOCATE(wk)
     
     DALLOCATE(FXmn)
     DALLOCATE(DXmn)
         
     if( iec00aa.le.-2 ) then ! input value of iec00aa=ifail controls screen output;  8 Jul 15;
      select case( ic05pbf )
      case( -3 )
        write(0,'("ec00aa : ic05pbf="i2" ; max. iterations  ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
     !case( -2 )
     !  write(0,'("ec00aa : ic05pbf="i2" ; user terminated  ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
      case( -1 )
        write(0,'("ec00aa : ic05pbf="i2" ; ibfield.ne.0     ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
      case(  0 )
        write(0,'("ec00aa : ic05pbf="i2" ; success          ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
      case(  1 )
        write(0,'("ec00aa : ic05pbf="i2" ; input error      ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
      case(  2 )
        write(0,'("ec00aa : ic05pbf="i2" ; consider restart ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
      case(  3 )
        write(0,'("ec00aa : ic05pbf="i2" ; xtol too small   ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
      case(  4 )
        write(0,'("ec00aa : ic05pbf="i2" ; bad progress     ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
      case default
        write(0,'("ec00aa : ic05pbf="i2" ; illegal ifail    ; its="i3" ; err="es12.5" ;")') ic05pbf, lcurve%its, lcurve%err
      end select
     endif

     if( ic05pbf.eq. 0 ) then ! NAG routine successfully located "x" s.t. "f(x)=0";  8 Jul 15;
      ;  ii = 0       ; lcurve%Rnc(ii) = XXmn(   ii+1) ; lcurve%Znc(ii) = XXmn(tN+   ii+1) ! unpack solution; 17 Jun 15;
      do ii = 1, Ntor ; lcurve%Rnc(ii) = XXmn(   ii+1) ; lcurve%Znc(ii) = XXmn(tN+   ii+1)
       ;              ; lcurve%Rns(ii) = XXmn(tN-ii+1) ; lcurve%Zns(ii) = XXmn(tN+tN-ii+1)
      enddo
     endif
     
     if( ic05pbf.eq.-3 ) then ! NAG routine successfully located "x" s.t. "f(x)=0";  8 Jul 15;
      lcurve%Rnc(0:Ntor) = curve%Rnc(0:Ntor) ; lcurve%Znc(0:Ntor) = curve%Znc(0:Ntor)
      lcurve%Rns(0:Ntor) = curve%Rns(0:Ntor) ; lcurve%Zns(0:Ntor) = curve%Zns(0:Ntor)
      ic05pbf = 0 ! over-rule error flag;  8 Jul 15;
     endif
     
     DALLOCATE(XXmn)
     
     if( ic05pbf.ne.0 ) then ; iec00aa = 2 ; goto 9998 ! if ic05pbf is returned by NAG as non-zero then an error has ocurred;  8 Jul 15;
     endif
     
    case( 1:4 ) ! gradient flow is under construction;  8 Jul 15;
     
     CHECKINPUT(ec00aa, lcurve%ftol  .lt.small, 9998 )
     CHECKINPUT(ec00aa, lcurve%tauend.lt.small, 9998 )
     CHECKINPUT(ec00aa, lcurve%dtau  .lt.small, 9998 )
     CHECKINPUT(ec00aa, lcurve%odetol.lt.small, 9998 )
     
     curve%ftol = lcurve%ftol ; curve%dtau = lcurve%dtau

     itau = lcurve%tauend / lcurve%dtau + 1
     
     SALLOCATE( curve%iRnc,(0:Ntor,0:itau),zero)
     SALLOCATE( curve%iRns,(0:Ntor,0:itau),zero)
     SALLOCATE( curve%iZnc,(0:Ntor,0:itau),zero)
     SALLOCATE( curve%iZns,(0:Ntor,0:itau),zero)

     SALLOCATE(Xfmn,(1:Nfmn),zero)    
     if    ( curve%emethod.eq.1 .or. curve%emethod.eq.2 ) then ! independent variable is X;  8 Jul 15;
      ;  ii = 0       ; Xfmn(   ii+1) = lcurve%Znc(ii)
      do ii = 1, Ntor ; Xfmn(   ii+1) = lcurve%Znc(ii)
       ;              ; Xfmn(tN-ii+1) = lcurve%Zns(ii)
      enddo
     elseif( curve%emethod.eq.3 .or. curve%emethod.eq.4 ) then ! independent variable is R;  8 Jul 15;
      ;  ii = 0       ; Xfmn(   ii+1) = lcurve%Rnc(ii)
      do ii = 1, Ntor ; Xfmn(   ii+1) = lcurve%Rnc(ii)
       ;              ; Xfmn(tN-ii+1) = lcurve%Rns(ii)
      enddo
     endif
     
     SALLOCATE(Ffmn,(1:Nfmn       ),zero)
!    SALLOCATE(Dfmn,(1:Lfmn,1:Nfmn),zero)
        
     tau = zero ; relabs = 'D' ; odetol = lcurve%odetol * ten ; curve%itau = 0 ; curve%ibfield = 0
 
#ifdef HAVE_NAG
     id02bjf = 1  ! NAG; 13 Oct 15;
     Lwk = 20 * Nfmn
     SALLOCATE(wk,(1:Lwk),zero)
1000 format("ec00aa :         "2x" : calling D02BJF, which calls ec00bb ; tau="es13.5" ; tauend="es13.5" ; N="i3 &
            " ; relabs="a2" ;")
     call D02BJF( tau, lcurve%tauend, Nfmn, Xfmn(1:Nfmn), ec00bb, odetol, relabs, ec00bc, ec00bd, wk(1:Lwk), id02bjf )
     DALLOCATE(wk)
     call ec00bb( tau, Xfmn(1:Nfmn), Ffmn(1:Nfmn) ) ; lcurve%its = 0 ; lcurve%err = curve%FF
#else
1000 format("ec00aa :         "2x" : calling lsodar, which calls ec00bb ; tau="es13.5" ; tauend="es13.5" ; N="i3 &
            " ; relabs="a2" ;")
!    set the integrator parameters.
     liw=20
     lrw=20+16*Nfmn
     SALLOCATE(iwork,(1:liw),zero)
     SALLOCATE(rwork,(1:lrw),zero)
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    jt=2 :  internally generated full Jacobian
     jt=2
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=odetol
!    atol :  absolute tolerance
     atol=odetol
!    ng : number of roots being found
     ng=1
!    initialize for loop
     tauc=tau
     taue=tauc
     call ec00bc(taue,Xfmn)
     do while( curve%ibfield.ne.1 )
       call lsodar(ec00bb,(/Nfmn/),Xfmn,taue,tauc, &
                  itol,(/rtol/),(/atol/),itask,istate,iopt, &
                  rwork,lrw,iwork,liw,du00aa,jt,ec00bd,ng,jroot)
       tauc=taue
       call ec00bc(taue,Xfmn)
     enddo
     id02bjf=istate
     if (istate>0) id02bjf=0
     DALLOCATE(iwork)
     DALLOCATE(rwork)
     call ec00bb( Nfmn, tau, Xfmn(1:Nfmn), Ffmn(1:Nfmn) ) ; lcurve%its = 0 ; lcurve%err = curve%FF
#endif
     
     DALLOCATE(Ffmn)
!    DALLOCATE(Dfmn)
       
     if( iec00aa.le.-1 ) then
      select case( id02bjf )
      case(  0 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; success ;         ")') id02bjf, tau, curve%FF
      case(  1 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; input error ;     ")') id02bjf, tau, curve%FF
      case(  2 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; no more progress ;")') id02bjf, tau, curve%FF
        id02bjf = 0
      case(  3 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; tol too small ;   ")') id02bjf, tau, curve%FF
      case(  4 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; xsol not reset ;  ")') id02bjf, tau, curve%FF
      case(  5 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; xsol too large ;  ")') id02bjf, tau, curve%FF
      case(  6 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; increase tauend ; ")') id02bjf, tau, curve%FF
      case(  7 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; serious error ;   ")') id02bjf, tau, curve%FF
      case default
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; illegal ifail ;   ")') id02bjf, tau, curve%FF
      end select
     endif
     
    !if( id02bjf.eq.0 ) then
     if    ( curve%emethod.eq.1 .or. curve%emethod.eq.2 ) then
      ;  ii = 0       ; lcurve%Znc(ii) = Xfmn(   ii+1) ; lcurve%Rnc(ii) = curve%Rnc(ii)
      do ii = 1, Ntor ; lcurve%Znc(ii) = Xfmn(   ii+1) ; lcurve%Rnc(ii) = curve%Rnc(ii)
       ;              ; lcurve%Zns(ii) = Xfmn(tN-ii+1) ; lcurve%Rns(ii) = curve%Rns(ii)
      enddo
     elseif( curve%emethod.eq.3 .or. curve%emethod.eq.4 ) then
      ;  ii = 0       ; lcurve%Rnc(ii) = Xfmn(   ii+1) ; lcurve%Znc(ii) = curve%Znc(ii)
      do ii = 1, Ntor ; lcurve%Rnc(ii) = Xfmn(   ii+1) ; lcurve%Znc(ii) = curve%Znc(ii)
       ;              ; lcurve%Rns(ii) = Xfmn(tN-ii+1) ; lcurve%Zns(ii) = curve%Zns(ii)
      enddo
     endif

     if( id02bjf.eq.6 ) curve%itau = curve%itau - 1
     
     FATAL(ec00aa, curve%itau.gt.itau, exceeds allocation) 

     itau = curve%itau
     
     ;  ii = 0       ;  curve%iRnc(ii,  itau) = lcurve%Rnc(ii       ) ;  curve%iZnc(ii,  itau) = lcurve%Znc(ii       )
     do ii = 1, Ntor ;  curve%iRnc(ii,  itau) = lcurve%Rnc(ii       ) ;  curve%iZnc(ii,  itau) = lcurve%Znc(ii       )
      ;              ;  curve%iRns(ii,  itau) = lcurve%Rns(ii       ) ;  curve%iZns(ii,  itau) = lcurve%Zns(ii       )
     enddo
     
     lcurve%itau = itau
     
     SALLOCATE(lcurve%iRnc,(0:Ntor,0:itau),zero)
     SALLOCATE(lcurve%iRns,(0:Ntor,0:itau),zero)
     SALLOCATE(lcurve%iZnc,(0:Ntor,0:itau),zero)
     SALLOCATE(lcurve%iZns,(0:Ntor,0:itau),zero)
     
     lcurve%Lallocated = 1
     
     ;  ii = 0       ; lcurve%iRnc(ii,0:itau) = curve%iRnc(ii,0:itau) ; lcurve%iZnc(ii,0:itau) = curve%iZnc(ii,0:itau)
     do ii = 1, Ntor ; lcurve%iRnc(ii,0:itau) = curve%iRnc(ii,0:itau) ; lcurve%iZnc(ii,0:itau) = curve%iZnc(ii,0:itau)
      ;              ; lcurve%iRns(ii,0:itau) = curve%iRns(ii,0:itau) ; lcurve%iZns(ii,0:itau) = curve%iZns(ii,0:itau)
     enddo
    !endif

     DALLOCATE(curve%iRnc)
     DALLOCATE(curve%iRns)
     DALLOCATE(curve%iZnc)
     DALLOCATE(curve%iZns)
     
     DALLOCATE(Xfmn)
     
     if( id02bjf.ne.0 ) then ; iec00aa = 2 ; goto 9998
     endif

    case( 5 ) ! gradient flow is under construction;  8 Jul 15;
     
     CHECKINPUT(ec00aa, lcurve%tauend .lt.small, 9998 )
     CHECKINPUT(ec00aa, lcurve%odetol .lt.small, 9998 )
     CHECKINPUT(ec00aa, lcurve%epsilon.lt.small, 9998 )
     
     curve%epsilon = lcurve%epsilon
     
     SALLOCATE(XXmn,(1:NXmn),zero) 
     ! pack independent variable;  8 Jul 15;
     ;  ii = 0       ; XXmn(   ii+1) = lcurve%Rnc(ii) ; XXmn(tN+   ii+1) = lcurve%Znc(ii)
     do ii = 1, Ntor ; XXmn(   ii+1) = lcurve%Rnc(ii) ; XXmn(tN+   ii+1) = lcurve%Znc(ii)
      ;              ; XXmn(tN-ii+1) = lcurve%Rns(ii) ; XXmn(tN+tN-ii+1) = lcurve%Zns(ii)
     enddo
     
     SALLOCATE(FXmn,(1:NXmn       ),zero)
!    SALLOCATE(DXmn,(1:LXmn,1:NXmn),zero)
             
     tau = zero ; relabs = 'D' ; odetol = lcurve%odetol * ten ; curve%ibfield = 0
 
#ifdef HAVE_NAG
     id02bjf = 1 ! NAG; 13 Oct 15;
     Lwk = 20 * NXmn
     SALLOCATE(wk,(1:Lwk),zero)
     call D02BJF( tau, lcurve%tauend, NXmn, XXmn(1:NXmn), ec00da, odetol, relabs, D02BJX, D02BJW, wk(1:Lwk), id02bjf )
     DALLOCATE(wk)
#else
!    set the integrator parameters.
     liw=20
     lrw=20+16*NXmn
     SALLOCATE(iwork,(1:liw),zero)
     SALLOCATE(rwork,(1:lrw),zero)
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=odetol
!    atol :  absolute tolerance
     atol=odetol
!    do integration
     call lsode(ec00da,(/NXmn/),XXmn,tau,lcurve%tauend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,du00aa,mf)
     id02bjf=istate
     if (istate>0) id02bjf=0
     DALLOCATE(iwork)
     DALLOCATE(rwork)
#endif

     DALLOCATE(FXmn)
!    DALLOCATE(DXmn)
       
     if( iec00aa.le.-1 ) then
      select case( id02bjf )
      case(  0 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; success ;         ")') id02bjf, tau, curve%FF
      case(  1 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; input error ;     ")') id02bjf, tau, curve%FF
      case(  2 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; no more progress ;")') id02bjf, tau, curve%FF
        id02bjf = 0
      case(  3 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; tol too small ;   ")') id02bjf, tau, curve%FF
      case(  4 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; xsol not reset ;  ")') id02bjf, tau, curve%FF
      case(  5 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; xsol too large ;  ")') id02bjf, tau, curve%FF
      case(  6 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; increase tauend ; ")') id02bjf, tau, curve%FF
      case(  7 )
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; serious error ;   ")') id02bjf, tau, curve%FF
      case default
        write(0,'("ec00aa : id02bjf="i2" : tau="f09.06" ; FF="es12.5" ; illegal ifail ;   ")') id02bjf, tau, curve%FF
      end select
     endif
     
    !if( id02bjf.eq.0 ) then
      ;  ii = 0       ; lcurve%Rnc(ii) = XXmn(   ii+1) ; lcurve%Znc(ii) = XXmn(tN+   ii+1) ! unpack solution; 17 Jun 15;
      do ii = 1, Ntor ; lcurve%Rnc(ii) = XXmn(   ii+1) ; lcurve%Znc(ii) = XXmn(tN+   ii+1)
       ;              ; lcurve%Rns(ii) = XXmn(tN-ii+1) ; lcurve%Zns(ii) = XXmn(tN+tN-ii+1)
      enddo
    !endif

     DALLOCATE(XXmn)
      
    case default
      
     write(0,'("ec00aa : input error ; emethod.eq."i2" is not supported ; ")') lcurve%emethod ; iec00aa = 1 ; goto 9998
     
    end select
    
    iec00aa = 0
    
9998 continue

#ifdef HAVE_NAG    
    DALLOCATE(curve%trig)
#endif    

    DALLOCATE(curve%rr)
    DALLOCATE(curve%rd)
    DALLOCATE(curve%zz)
    DALLOCATE(curve%zd)
    
    DALLOCATE(curve%ct)
    DALLOCATE(curve%st)
    
    DALLOCATE(curve%Rnc)
    DALLOCATE(curve%Rns)
    DALLOCATE(curve%Znc)
    DALLOCATE(curve%Zns)
    
9999 continue
    
    lcurve%nbfield(0:1) = nbfield(0:1)
    
    if( ifail.le. 0 ) write(0,9000) iec00aa, lcurve%Ntor, lcurve%its, lcurve%err, lcurve%Lallocated, nbfield(0)+nbfield(1)
    
    if( ifail.le.-1 ) write(0,9001)          "Rnc", lcurve%Rnc(0:lcurve%Ntor)
    if( ifail.le.-1 ) write(0,9002)          "Rns", lcurve%Rns(1:lcurve%Ntor)
    if( ifail.le.-1 ) write(0,9001)          "Znc", lcurve%Znc(0:lcurve%Ntor)
    if( ifail.le.-1 ) write(0,9002)          "Zns", lcurve%Zns(1:lcurve%Ntor)
    
9000 format("ec00aa : ifail ="i3" : Ntor="i3" ; its="i3" ; err="es12.5" ; Lallocated="i2" ; nbfield="i11" ;")
9001 format("ec00aa :        "3x" : "a3"="    99es14.06)
9002 format("ec00aa :        "3x" : "a3"="14x,99es14.06)
    
    ifail = iec00aa
    
    return
    
  end subroutine ec00aa
  
  subroutine ec00ab( NXmn, XXmn, FXmn, DXmn, LXmn, iflag )
    
    implicit none
    
    INTEGER                :: NXmn, LXmn, iflag
    REAL                   :: XXmn(1:NXmn), FXmn(1:NXmn), DXmn(1:LXmn,1:NXmn)
    
    INTEGER                :: astat, ii, jj, irz, Ntor, MM, NN, ic06fpf, ic06fqf, ifail, tN, tA, Nfp
    REAL                   :: phi, RpZ(1:3)
    REAL, allocatable      :: xx(:,:), ff(:,:), bb(:,:,:)
    
    REAL, allocatable      :: wk(:), gg(:,:)

#ifndef HAVE_NAG
    INTEGER                :: Lwk
    REAL, allocatable      :: wkb(:)
#endif

    Ntor = curve%Ntor ; tN = 1+Ntor + Ntor ; MM = 4 ; NN = 4*max(1,Ntor) ; Nfp = curve%Nfp

    curve%its = curve%its + 1

    if( iflag.eq.1 .and. iec00aa.le.-3 ) then
     write(0,'("ec00ab :       "4x" :     "  13x " ; Rnc="    99f11.6" ;")') (  XXmn(      ii+1), ii = 0, Ntor  )
     write(0,'("ec00ab :       "4x" :     "  13x " ; Rns="11x,99f11.6" ;")') (  XXmn(   tN-ii+1), ii = 1, Ntor  )
     write(0,'("ec00ab :       "4x" :     "  13x " ; Znc="    99f11.6" ;")') (  XXmn(tN+   ii+1), ii = 0, Ntor  )
     write(0,'("ec00ab :       "4x" :     "  13x " ; Zns="11x,99f11.6" ;")') (  XXmn(tN+tN-ii+1), ii = 1, Ntor  )
    endif
    
    SALLOCATE(xx,(1:MM,1:NN),zero)
    
    ! unpack independent variable;  8 Jul 15;
    ;   ii = 0       ; xx(1,   ii+1) =        XXmn(   ii+1)        ; xx(3,   ii+1) =        XXmn(tN+   ii+1)
    ;                ; xx(2,   ii+1) =        zero                 ; xx(4,   ii+1) =        zero
    if( Ntor.ge.1 ) then
     do ii = 1, Ntor ; xx(1,   ii+1) =        XXmn(   ii+1) * half ; xx(3,   ii+1) =        XXmn(tN+   ii+1) * half
      ;              ; xx(1,NN-ii+1) =        XXmn(tN-ii+1) * half ; xx(3,NN-ii+1) =        XXmn(tN+tN-ii+1) * half
      ;              ; xx(2,   ii+1) =   ii * XXmn(tN-ii+1) * half ; xx(4,   ii+1) =   ii * XXmn(tN+tN-ii+1) * half
      ;              ; xx(2,NN-ii+1) = - ii * XXmn(   ii+1) * half ; xx(4,NN-ii+1) = - ii * XXmn(tN+   ii+1) * half
     enddo
    endif

    xx(1:MM,1:NN) = xx(1:MM,1:NN) * curve%sN ! FFT scaling factor;  8 Jul 15;
        
#ifdef HAVE_NAG
    SALLOCATE(wk,(1:MM*NN),zero)
    ic06fqf = 0
    call C06FQF( MM, NN, xx(1:MM,1:NN), curve%init, curve%trig(1:2*NN), wk(1:MM*NN), ic06fqf ) ; curve%init = 'S' ! NAG; 13 Oct 15;
#else
    Lwk = NN + INT(LOG(real(NN))/LOG(2.)) + 4
    SALLOCATE(wk,(1:Lwk),zero)
    call rfft1i( NN, wk, Lwk, ic06fqf )
    ! apply transform to each row
    SALLOCATE(wkb,(1:NN),zero)
    do ii = 1, MM
      call rfft1b( NN, 1, xx(ii,1:NN), NN, wk, Lwk, wkb, NN, ic06fqf ) ! this is destructive
    enddo
#endif

    curve%rr(1:NN) = xx(1,1:NN)       !  R(\phi)      ;  8 Jul 15;
    curve%rd(1:NN) = xx(2,1:NN) * Nfp ! dR(\phi)/d\phi;  8 Jul 15;
    curve%zz(1:NN) = xx(3,1:NN)       !  Z(\phi)      ;  8 Jul 15;
    curve%zd(1:NN) = xx(4,1:NN) * Nfp ! dZ(\phi)/d\phi;  8 Jul 15;
    
    DALLOCATE(xx)
    
    if( iflag.eq.0 ) then ; deallocate(wk) ; goto 9999
    endif

    SALLOCATE(bb,(1:3,0:3,1:NN),zero)
    
    if( iflag.eq.1 ) then ; itangent = 0
    else                  ; itangent = 1
    endif
    
    do ii = 1, NN

     phi = (ii-1) * (pi2/Nfp) / NN ; RpZ(1:3) = (/ curve%rr(ii), phi, curve%zz(ii) /)
     
     nbfield(itangent) = nbfield(itangent) + 1 ; ifail = -9 ; call bfield( RpZ(1:3), itangent, bb(1:3,0:3,ii), ifail )
     
     if( ifail.ne.0 ) then
      if( iec00aa.le.-5 ) then
        bb(1:3,0,ii) = ten
      else
        write(0,'("ec00ab : ifail ="i3" : (R,p,Z)=("3f15.10") ;")') ifail, RpZ(1:3)
        deallocate(wk,bb) ; iflag = -1 ; goto 9999
      endif
     endif
     
    enddo ! end of do ii = 1, NN ; 17 Jun 15;
    
    SALLOCATE(ff,(1:MM,1:NN),zero)
    
    select case( iflag )
     
    case( 1 )
     
     ff(1,1:NN) = ( + bb(3,0,1:NN) - bb(2,0,1:NN) * curve%zd(1:NN) ) * curve%rr(1:NN)  !  FR  ;
     ff(3,1:NN) = ( - bb(1,0,1:NN) + bb(2,0,1:NN) * curve%rd(1:NN) ) * curve%rr(1:NN)  !  FZ  ;
     
    case( 2 )
     
     ff(1,1:NN) = ( + bb(3,1,1:NN) - bb(2,1,1:NN) * curve%zd(1:NN) ) * curve%rr(1:NN) &
                  + bb(3,0,1:NN) - bb(2,0,1:NN) * curve%zd(1:NN)                       ! dFRdR;
     ff(2,1:NN) = ( + bb(3,3,1:NN) - bb(2,3,1:NN) * curve%zd(1:NN) ) * curve%rr(1:NN)  ! dFRdZ;
     ff(3,1:NN) = ( - bb(1,1,1:NN) + bb(2,1,1:NN) * curve%rd(1:NN) ) * curve%rr(1:NN) &
                  - bb(1,0,1:NN) + bb(2,0,1:NN) * curve%rd(1:NN)                       ! dFZdR;
     ff(4,1:NN) = ( - bb(1,3,1:NN) + bb(2,3,1:NN) * curve%rd(1:NN) ) * curve%rr(1:NN)  ! dFZdZ;
     
     SALLOCATE(gg,(1:MM,1:NN),zero)
     
    end select
    
    select case( iflag )
     
    case( 1 )
     
#ifdef HAVE_NAG
       ic06fpf = 0
       call C06FPF( MM, NN, ff(1:MM,1:NN), curve%init, curve%trig(1:2*NN), wk(1:MM*NN), ic06fpf ) ! NAG; 13 Oct 15;
#else
       ! apply transform to each row
       do ii = 1, MM
         call rfft1f( NN, 1, ff(ii,1:NN), NN, wk, Lwk, wkb, NN, ic06fpf ) ! this is destructive
       enddo
#endif
     
       ff(1:MM,1:NN) = ff(1:MM,1:NN) / curve%sN ; ff(1:MM,2:NN) = ff(1:MM,2:NN) * two
     
       ;   ii = 0       ; FXmn(      ii+1           ) =   ff(1,   ii+1) ; FXmn(tN+   ii+1           ) =   ff(3,   ii+1)
       if( Ntor.ge.1 ) then
        do ii = 1, Ntor ; FXmn(      ii+1           ) =   ff(1,   ii+1) ; FXmn(tN+   ii+1           ) =   ff(3,   ii+1)
         ;              ; FXmn(   tN-ii+1           ) = - ff(1,NN-ii+1) ; FXmn(tN+tN-ii+1           ) = - ff(3,NN-ii+1)
        enddo
       endif
     
    case( 2 )
     
     do jj = 0, Ntor ! labels Fourier harmonic that the derivative is with respect to; 17 Jun 15;
      
      do irz = 0, 1 ; tA = tN * irz
       
       if( irz.eq.0 ) then 
         gg(1,1:NN) = ff(1,1:NN) * curve%ct(1:NN,jj)                                                                ! dFRdRc ;
         gg(2,1:NN) = ff(1,1:NN) * curve%st(1:NN,jj)                                                                ! dFRdRs ;
         gg(3,1:NN) = ff(3,1:NN) * curve%ct(1:NN,jj) - bb(2,0,1:NN) * curve%rr(1:NN) * curve%st(1:NN,jj) * jj * Nfp ! dFZdRc ;
         gg(4,1:NN) = ff(3,1:NN) * curve%st(1:NN,jj) + bb(2,0,1:NN) * curve%rr(1:NN) * curve%ct(1:NN,jj) * jj * Nfp ! dFZdRs ;
       else
         gg(1,1:NN) = ff(2,1:NN) * curve%ct(1:NN,jj) + bb(2,0,1:NN) * curve%rr(1:NN) * curve%st(1:NN,jj) * jj * Nfp ! dFRdZc ;
         gg(2,1:NN) = ff(2,1:NN) * curve%st(1:NN,jj) - bb(2,0,1:NN) * curve%rr(1:NN) * curve%ct(1:NN,jj) * jj * Nfp ! dFRdZs ;
         gg(3,1:NN) = ff(4,1:NN) * curve%ct(1:NN,jj)                                                                ! dFZdZc ;
         gg(4,1:NN) = ff(4,1:NN) * curve%st(1:NN,jj)                                                                ! dFZdZs ;
       endif
       
#ifdef HAVE_NAG
       ic06fpf = 0
       call C06FPF( MM, NN, gg(1:MM,1:NN), curve%init, curve%trig(1:2*NN), wk(1:MM*NN), ic06fpf ) ! NAG; 13 Oct 15;
#else
       ! apply transform to each row
       do ii = 1, MM
         call rfft1f( NN, 1, gg(ii,1:NN), NN, wk, Lwk, wkb, NN, ic06fpf ) ! this is destructive
       enddo
#endif
       
       gg(1:MM,1:NN) = gg(1:MM,1:NN) / curve%sN ; gg(1:MM,2:NN) = gg(1:MM,2:NN) * two
       
       ;   ii = 0       ; DXmn(      ii+1,tA+   jj+1) =   gg(1,   ii+1) ; DXmn(tN+   ii+1,tA+   jj+1) =   gg(3,   ii+1) ! dFRcdXc ; dFZcdXc ; 17 Jun 15;
       ;   if( jj.gt.0 ) then
       ;                ; DXmn(      ii+1,tA+tN-jj+1) =   gg(2,   ii+1) ; DXmn(tN+   ii+1,tA+tN-jj+1) =   gg(4,   ii+1) ! dFRcdXs ; dFZcdXs ; 17 Jun 15;
       ;   endif
       if( Ntor.gt.0 ) then
        do ii = 1, Ntor ; DXmn(      ii+1,tA+   jj+1) =   gg(1,   ii+1) ; DXmn(tN+   ii+1,tA+   jj+1) =   gg(3,   ii+1) ! dFRcdXc ; dFZcdXc ; 17 Jun 15;
         ; if( jj.gt.0 ) then
         ;              ; DXmn(      ii+1,tA+tN-jj+1) =   gg(2,   ii+1) ; DXmn(tN+   ii+1,tA+tN-jj+1) =   gg(4,   ii+1) ! dFRcdXs ; dFZcdXs ; 17 Jun 15;
         ; endif
         ;              ; DXmn(   tN-ii+1,tA+   jj+1) = - gg(1,NN-ii+1) ; DXmn(tN+tN-ii+1,tA+   jj+1) = - gg(3,NN-ii+1) ! dFRsdXc ; dFZsdXc ; 17 Jun 15;
         ; if( jj.gt.0 ) then
         ;              ; DXmn(   tN-ii+1,tA+tN-jj+1) = - gg(2,NN-ii+1) ; DXmn(tN+tN-ii+1,tA+tN-jj+1) = - gg(4,NN-ii+1) ! dFRsdXs ; dFZsdXs ; 17 Jun 15;
         ; endif
        enddo ! end of do ii; 17 Jun 15;
       endif
       
      enddo ! end of do irz; 17 Jun 15;
     enddo ! end of do jj; 17 Jun 15;
     
     DALLOCATE(gg)
     
    end select
    
    DALLOCATE(bb)
    DALLOCATE(wk)
    
    DALLOCATE(ff)
    
9998 continue
    
    if( iflag.eq.1 ) then
     curve%err = sqrt( sum(FXmn(1:NXmn)**2) )
     if( iec00aa.le.-3 ) then
      write(0,'("ec00ab : its = "i4" : |F|="es13.5" ; ")') curve%its, curve%err
      ;  ii = 0       ; curve%Rnc(ii) = XXmn(   ii+1) ; curve%Znc(ii) = XXmn(tN+   ii+1)
      do ii = 1, Ntor ; curve%Rnc(ii) = XXmn(   ii+1) ; curve%Znc(ii) = XXmn(tN+   ii+1)
       ;              ; curve%Rns(ii) = XXmn(tN-ii+1) ; curve%Zns(ii) = XXmn(tN+tN-ii+1)
      enddo
      if( iec00aa.le.-5 ) pause
     endif
!    if( curve%err.le.curve%etol   ) iflag = -2 ! can save solution in curve%Rnc, etc. ; 19 Jun 15;
     if( curve%its.ge.curve%maxits ) iflag = -3
    else
     if( iec00aa.le.-3 ) then
      write(0,'("ec00ab : its = "i4" :  D           ;")') curve%its
     endif
    endif
    
9999 continue
    
    return
    
  end subroutine ec00ab
  
#ifdef HAVE_NAG
  subroutine ec00bb( tau, Xfmn, Ffmn ) ! format constrained by NAG; 19 Jun 15; ! returns negative gradient for flow minimization; 29 Jul 15;
    implicit none
#else    
  subroutine ec00bb( Node, tau, Xfmn, Ffmn ) ! format constrained by lsode
    implicit none
    INTEGER, intent(in) :: Node
#endif
        
    INTEGER             :: Nfmn, Lfmn
    REAL                :: tau, Xfmn(*), Ffmn(*)
    
    INTEGER             :: astat, Ntor, tN, NXmn, LXmn, Ncmn, Lcmn, Lwk, iflag, ic05pbf, ii
    REAL, allocatable   :: Dfmn(:,:), XXmn(:), FXmn(:), DXmn(:,:), Xcmn(:), Fcmn(:), Dcmn(:,:), wk(:)
    
    Ntor = curve%Ntor ; tN = 1+Ntor + Ntor ; Nfmn = tN ; Lfmn = Nfmn
    NXmn = 2 * Nfmn ; LXmn = NXmn ; Ncmn = Nfmn ; Lcmn = Nfmn

    if( curve%ibfield.eq.1 ) then ; curve%FF = one ; Ffmn(1:Nfmn) = zero ; goto 9999
    endif

    Lwk = Ncmn * ( Ncmn + 13 ) / 2

    SALLOCATE(XXmn,(1:NXmn),zero)
    if    ( curve%emethod.eq.1 .or. curve%emethod.eq.2 ) then
     ! Zmn are independent freedom; need to give initial guess for Xc = R ; 19 Jun 15;
     ;  ii = 0       ; XXmn(   ii+1) = curve%Rnc(ii) ; XXmn(tN+   ii+1) = Xfmn(   ii+1)
     do ii = 1, Ntor ; XXmn(   ii+1) = curve%Rnc(ii) ; XXmn(tN+   ii+1) = Xfmn(   ii+1)
      ;              ; XXmn(tN-ii+1) = curve%Rns(ii) ; XXmn(tN+tN-ii+1) = Xfmn(tN-ii+1)
     enddo
    elseif( curve%emethod.eq.3 .or. curve%emethod.eq.4 ) then
     ! Rmn are independent freedom; need to give initial guess for Xc = R ; 19 Jun 15;
     ;  ii = 0       ; XXmn(   ii+1) = Xfmn(   ii+1) ; XXmn(tN+   ii+1) = curve%Znc(ii)
     do ii = 1, Ntor ; XXmn(   ii+1) = Xfmn(   ii+1) ; XXmn(tN+   ii+1) = curve%Znc(ii)
      ;              ; XXmn(tN-ii+1) = Xfmn(tN-ii+1) ; XXmn(tN+tN-ii+1) = curve%Zns(ii)
     enddo
    endif
        
    SALLOCATE(FXmn,(1:NXmn       ),zero)
    SALLOCATE(DXmn,(1:LXmn,1:NXmn),zero)
    
    iflag = 0 ; curve%its = 0
    call ec00ab( NXmn, XXmn(1:NXmn), FXmn(1:NXmn), DXmn(1:LXmn,1:NXmn), LXmn, iflag ) ! to compute curve%rr, curve%rd, etc. ;
    
    DALLOCATE(FXmn)
    DALLOCATE(DXmn)

    DALLOCATE(XXmn)

    SALLOCATE(Xcmn,(1:Ncmn),zero)
    
    if    ( curve%emethod.eq.1 .or. curve%emethod.eq.2 ) then
     ! Zmn are independent freedom; need to give initial guess for Xc = R ; 19 Jun 15;
     ;  ii = 0       ; Xcmn(   ii+1) = curve%Rnc(ii)
     do ii = 1, Ntor ; Xcmn(   ii+1) = curve%Rnc(ii)
      ;              ; Xcmn(tN-ii+1) = curve%Rns(ii)
     enddo
    elseif( curve%emethod.eq.3 .or. curve%emethod.eq.4 ) then
     ! Rmn are independent freedom; need to give initial guess for Xc = Z ; 19 Jun 15;
     ;  ii = 0       ; Xcmn(   ii+1) = curve%Znc(ii)
     do ii = 1, Ntor ; Xcmn(   ii+1) = curve%Znc(ii)
      ;              ; Xcmn(tN-ii+1) = curve%Zns(ii)
     enddo
    endif
    
    SALLOCATE(Fcmn,(1:Ncmn       ),zero)
    SALLOCATE(Dcmn,(1:Lcmn,1:Lcmn),zero)

    SALLOCATE(wk,(1:Lwk),zero)
    
    curve%Lrz = + curve%emethod

    curve%its = 0

#ifdef HAVE_NAG
    ic05pbf = 1 ! NAG; 13 Oct 15;
    call C05PBF( ec00ca, Ncmn, Xcmn(1:Ncmn), Fcmn(1:Ncmn), Dcmn(1:Lcmn,1:Ncmn), Lcmn, curve%etol, wk(1:Lwk), Lwk, ic05pbf )
#else 
    call HYBRJ1( ec00ca, Ncmn, Xcmn(1:Ncmn), Fcmn(1:Ncmn), Dcmn(1:Lcmn,1:Ncmn), Lcmn, curve%etol, ic05pbf, wk(1:Lwk), Lwk )
    if (ic05pbf==1) ic05pbf=0 ! because the MINPACK error code is stupid
#endif

    DALLOCATE(wk)
    
    curve%err = sqrt( sum( Fcmn(1:Ncmn)**2 ) ) ! need to pass solution for when iflag=-2 or iflag=-1 ; 17 Jun 15;
    
    if( iec00aa.le.-4 ) then
     select case( ic05pbf )
    !case( -3 )
    !  write(0,'("ec00bb : ic05pbf="i2" ; max. iterations  ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
    !case( -2 )
    !  write(0,'("ec00bb : ic05pbf="i2" ; user terminated  ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
     case( -1 )
       write(0,'("ec00bb : ic05pbf="i2" ; ibfield.ne.0     ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
     case(  0 )
       write(0,'("ec00bb : ic05pbf="i2" ; success          ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
     case(  1 )
       write(0,'("ec00bb : ic05pbf="i2" ; input error      ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
     case(  2 )
       write(0,'("ec00bb : ic05pbf="i2" ; consider restart ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
     case(  3 )
       write(0,'("ec00bb : ic05pbf="i2" ; xtol too small   ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
     case(  4 )
       write(0,'("ec00bb : ic05pbf="i2" ; bad progress     ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
     case default
       write(0,'("ec00bb : ic05pbf="i2" ; illegal ifail    ; its="i3" ; err="es12.5" ;")') ic05pbf, curve%its, curve%err
     end select
    endif

    if    ( curve%emethod.eq.1 .or. curve%emethod.eq.2 ) then ! Zmn are independent freedom; 19 Jun 15;
     ;  ii = 0       ; curve%Rnc(ii) = Xcmn(   ii+1)
     do ii = 1, Ntor ; curve%Rnc(ii) = Xcmn(   ii+1)
      ;              ; curve%Rns(ii) = Xcmn(tN-ii+1)
     enddo
    elseif( curve%emethod.eq.3 .or. curve%emethod.eq.4 ) then ! Rmn are independent freedom; 19 Jun 15;
     ;  ii = 0       ; curve%Znc(ii) = Xcmn(   ii+1)
     do ii = 1, Ntor ; curve%Znc(ii) = Xcmn(   ii+1)
      ;              ; curve%Zns(ii) = Xcmn(tN-ii+1)
     enddo
    endif
    
    curve%Lrz = + curve%emethod

    iflag = 0 ; curve%its = 0 ; call ec00ca( Ncmn, Xcmn(1:Ncmn), Fcmn(1:Ncmn), Dcmn(1:Lcmn,1:Ncmn), Lcmn, iflag )
    
    DALLOCATE(Xcmn)
    DALLOCATE(Fcmn)
    DALLOCATE(Dcmn)
   
    if( ic05pbf.ne.0 ) then ; curve%FF = one ; Ffmn(1:Nfmn) = zero ; goto 9999
    endif
    
    SALLOCATE(Dfmn,(1:Lfmn,1:Nfmn),zero)

    curve%Lrz = - curve%emethod

    iflag = 1 ; curve%its = 0 ; call ec00ca( Nfmn, Xfmn(1:Nfmn), Ffmn(1:Nfmn), Dfmn(1:Lfmn,1:Nfmn), Lfmn, iflag )

    curve%FF = sqrt( sum( Ffmn(1:Nfmn)**2 ) )
    
   !Ffmn(1:Nfmn) = - Ffmn(1:Nfmn) / abs( Ffmn(1) ) ! this normalization may not be robust; 29 Jul 15;
    Ffmn(1:Nfmn) = - Ffmn(1:Nfmn) / curve%FF

    DALLOCATE(Dfmn)

!    curve%its = 0    
!    iflag = 2
!     call ec00cb( NZmn, XZmn(1:NZmn), FZmn(1:NZmn), DZmn(1:LZmn,1:NZmn), LZmn, iflag ) ! to re-compute derivatives of FZ wrt Z ; DZmn;
!     LDA = NZmn
!     SALLOCATE(ipivot,(1:NZmn),0)
!     call F07ADF( LZmn, NZmn, DZmn(1:LDA,1:NZmn), LDA, ipivot(1:NZmn), info )
!     FATAL( ec00aa, info.ne.0, error calling F07AFF to solve A = PLU )
!     call F07AJF( NZmn, DZmn(1:LDA,1:NZmn), LDA, ipivot(1:NZmn), wk(1:Lwk), Lwk, info )
!     FATAL( ec00aa, info.ne.0, error calling F07AJF to solve 1 / A )
!     DALLOCATE(ipivot)    
!     if( iec00aa.le.-2 ) write(0,'("ec00aa :         "2x" : constructed inverse of DzFz ; derivatives of Z wrt R are available ;")') 
    
9999 continue
    
    if( iec00aa.le.-3 ) then
      write(0,'("ec00bb :         "2x" ; tau="f09.06" ; FF="es12.5" ; Xfmn="999f12.08)') tau, curve%FF, Xfmn(1:Nfmn)
    endif
   !if( iec00aa.le. 0 ) then
   !  write(0,'("ec00bb :         "2x" ;     " 10x  " ; FF="es12.5" ; Ffmn="999f12.08)')      curve%FF, Ffmn(1:Nfmn)

    if( iec00aa.le.-5 ) pause
    
    return
    
  end subroutine ec00bb
  
#ifdef HAVE_NAG
  subroutine ec00da( tau, XXmn, FXmn ) ! format constrained by NAG; 19 Jun 15; ! returns negative gradient for flow minimization; 29 Jul 15;
    implicit none
#else    
  subroutine ec00da( Node, tau, XXmn, FXmn ) ! format constrained by lsode
    implicit none
    INTEGER, intent(in) :: Node
#endif
        
    REAL              :: tau, XXmn(*), FXmn(*)
    
    INTEGER           :: NXmn, LXmn, iflag, Ntor, tN, astat
    REAL, allocatable :: DXmn(:,:)
    
    Ntor = curve%Ntor ; tN = 1+Ntor + Ntor ; NXmn = 2 * tN ; LXmn = NXmn
    
    SALLOCATE(DXmn,(1:LXmn,1:NXmn),zero)
    
    iflag = 1
    call ec00ab( NXmn, XXmn(1:NXmn), FXmn(1:NXmn), DXmn(1:LXmn,1:NXmn), LXmn, iflag )

    DALLOCATE(DXmn)
    
    FXmn( 1:   tN) = - FXmn( 1:   tN)
    FXmn(tN:tN+tN) = + FXmn(tN:tN+tN) / curve%epsilon

9999 continue
    
!   if( iec00aa.le.-3 ) write(0,'("ec00da :         "2x" ; tau="f09.06" ; XXmn="999f12.08)') tau, XXmn(1:NXmn)
!   if( iec00aa.le. 0 ) write(0,'("ec00da :         "2x" ;     " 09x  " ; FXmn="999f12.08)')      FXmn
    
!   if( iec00aa.le.-5 ) pause
    
    return
    
  end subroutine ec00da
  
  subroutine ec00bc( tau, Xfmn ) ! int. output of action-minimizing gradient flow; 29 Jul 15;

    implicit none

    REAL                   :: tau, Xfmn(*)
    INTEGER                :: astat, Ntor, tN, Nfmn, ii, itau
    REAL, allocatable      :: Ffmn(:)

    Ntor = curve%Ntor ; tN = 1+Ntor + Ntor ; Nfmn = tN

    SALLOCATE(Ffmn,(1:Nfmn),zero)

#ifdef HAVE_NAG
    call ec00bb( tau, Xfmn(1:Nfmn), Ffmn(1:Nfmn) )
#else
    call ec00bb( Nfmn, tau, Xfmn(1:Nfmn), Ffmn(1:Nfmn) )
#endif

    DALLOCATE(Ffmn)

    if( iec00aa.le.-2 ) then 
      write(0,'("ec00bc :         "2x" ; tau="f09.06" ; FF="es12.5" ; Xfmn="999f12.08)') tau, curve%FF, Xfmn(1:Nfmn)
    endif

    itau = curve%itau ! shorthand; 29 Jul 15;

    if    ( curve%emethod.eq.1 .or. curve%emethod.eq.2 ) then
     ;  ii = 0       ; curve%iZnc(ii,itau) = Xfmn(   ii+1) ; curve%iRnc(ii,itau) = curve%Rnc(ii)
     do ii = 1, Ntor ; curve%iZnc(ii,itau) = Xfmn(   ii+1) ; curve%iRnc(ii,itau) = curve%Rnc(ii)
      ;              ; curve%iZns(ii,itau) = Xfmn(tN-ii+1) ; curve%iRns(ii,itau) = curve%Rns(ii)
     enddo
    elseif( curve%emethod.eq.3 .or. curve%emethod.eq.4 ) then
     ;  ii = 0       ; curve%iRnc(ii,itau) = Xfmn(   ii+1) ; curve%iZnc(ii,itau) = curve%Znc(ii)
     do ii = 1, Ntor ; curve%iRnc(ii,itau) = Xfmn(   ii+1) ; curve%iZnc(ii,itau) = curve%Znc(ii)
      ;              ; curve%iRns(ii,itau) = Xfmn(tN-ii+1) ; curve%iZns(ii,itau) = curve%Zns(ii)
     enddo
    endif
        
   !;               curve%iRnc(0:Ntor) = Xfmn(       1:Ntor+1: 1)
   !if( Ntor.gt.0 ) curve%iRns(1:Ntor) = Xfmn(2*Ntor+1:Ntor+2:-1)
    
    curve%itau = curve%itau + 1 ; tau = curve%itau * curve%dtau ! increment output; 29 Jul 15;

9999 continue

    return

  end subroutine ec00bc
  
#ifdef HAVE_NAG
  REAL function ec00bd( tau, Xfmn ) ! termination of action-minimizing gradient flow; 29 Jul 15;

    implicit none

    REAL                   :: tau, Xfmn(*)
    INTEGER                :: astat, Ntor, tN, Nfmn
    REAL, allocatable      :: Ffmn(:)

    Ntor = curve%Ntor ; tN = 1+Ntor + Ntor ; Nfmn = tN
    SALLOCATE(Ffmn,(1:Nfmn),zero)
    call ec00bb( tau, Xfmn(1:Nfmn), Ffmn(1:Nfmn) )
    DALLOCATE(Ffmn)
    ec00bd = curve%FF - curve%ftol ! curve%FF is calculated in ec00bb; 29 Jul 15;

9999 continue

    return
  end function ec00bd
#else    
  subroutine ec00bd( Node, tau, Xfmn, Nofmn, Ofmn ) ! format constrained by lsodar

    implicit none

    INTEGER, intent(in)    :: Node, Nofmn
    REAL, intent(out)      :: Ofmn(1:Nofmn)
    REAL                   :: tau, Xfmn(1:Node)
    INTEGER                :: astat, Ntor, tN, Nfmn
    REAL, allocatable      :: Ffmn(:)

    Ntor = curve%Ntor ; tN = 1+Ntor + Ntor ; Nfmn = tN
    SALLOCATE(Ffmn,(1:Nfmn),zero)
    call ec00bb( Node, tau, Xfmn(1:Nfmn), Ffmn(1:Nfmn) )
    Ofmn = curve%FF - curve%ftol ! curve%FF is calculated in ec00bb; 29 Jul 15;
    DALLOCATE(Ffmn)

9999 continue

    return
  end subroutine ec00bd
#endif

  subroutine ec00ca( Ncmn, Xcmn, Fcmn, Dcmn, Lcmn, iflag ) ! returns FR and dFRdR, or FZ and dFZdR; needs curve%rr, curve%zz & curve%zd ; 19 Jun 15;
    
    implicit none
    
    INTEGER                :: Ncmn, Lcmn, iflag
    REAL                   :: Xcmn(1:Ncmn), Fcmn(1:Ncmn), Dcmn(1:Lcmn,1:Ncmn)
    
    INTEGER                :: astat, ii, jj, kk, Ntor, MM, NN, ic06fpf, ic06fqf, ifail, tN, Nfp
    REAL                   :: phi, RpZ(1:3)
    REAL, allocatable      :: xx(:,:), ff(:,:), bb(:,:,:)
    
    REAL, allocatable      :: wk(:), gg(:,:)

#ifndef HAVE_NAG
    INTEGER                :: Lwk
    REAL, allocatable      :: wkb(:)
#endif
    
    Ntor = curve%Ntor ; tN = 1+Ntor + Ntor ; MM = 4 ; NN = 4*max(1,Ntor) ; Nfp = curve%Nfp

    if( iec00aa.le.-6 ) then
      write(0,'("ec00ca :         "2x" : curve%its="i3" ; iflag="i2" ; curve%Lrz="i3" ;")') curve%its, iflag, curve%Lrz
    endif

    if( curve%Lrz.lt.0 ) then
     SALLOCATE(wk,(1:MM*NN),zero)
     goto 8000 ! this assumes that the curve%rr and curve%zz are correct; 19 Jun 15;
    endif
    
    SALLOCATE(xx,(1:MM,1:NN),zero) ! pack the Fourier harmonics of the dependent variable into FFT format; 19 Jun 15;
    
    ;   ii = 0       ; xx(1,   ii+1) =        Xcmn(   ii+1)       
    ;                ; xx(2,   ii+1) =        zero                
    if( Ntor.ge.1 ) then                                          
     do ii = 1, Ntor ; xx(1,   ii+1) =        Xcmn(   ii+1) * half
      ;              ; xx(1,NN-ii+1) =        Xcmn(tN-ii+1) * half
      ;              ; xx(2,   ii+1) =   ii * Xcmn(tN-ii+1) * half
      ;              ; xx(2,NN-ii+1) = - ii * Xcmn(   ii+1) * half
     enddo
    endif

    xx(1:MM,1:NN) = xx(1:MM,1:NN) * curve%sN
        
#ifdef HAVE_NAG
    SALLOCATE(wk,(1:MM*NN),zero)
    ic06fqf = 1
    call C06FQF( MM, NN, xx(1:MM,1:NN), curve%init, curve%trig(1:2*NN), wk(1:MM*NN), ic06fqf ) ; curve%init = 'S' ! NAG; 13 Oct 15;
#else
    Lwk = NN + INT(LOG(real(NN))/LOG(2.)) + 4
    SALLOCATE(wk,(1:Lwk),zero)
    call rfft1i( NN, wk, Lwk, ic06fqf )
    ! apply transform to each row
    SALLOCATE(wkb,(1:NN),zero)
    do ii = 1, MM
      call rfft1b( NN, 1, xx(ii,1:NN), NN, wk, Lwk, wkb, NN, ic06fqf ) ! this is destructive
    enddo
#endif

    FATAL(ec00ca, ic06fqf.ne.0, error constructing FFT - please ask srh to implement soft fail)
    
    if    ( curve%Lrz.eq. 1 .or. curve%Lrz.eq. 2 ) then ; curve%rr(1:NN) = xx(1,1:NN) ; curve%rd(1:NN) = xx(2,1:NN) * Nfp
    elseif( curve%Lrz.eq. 3 .or. curve%Lrz.eq. 4 ) then ; curve%zz(1:NN) = xx(1,1:NN) ; curve%zd(1:NN) = xx(2,1:NN) * Nfp
    endif

    DALLOCATE(xx)

    if( iflag.eq.0 ) then ; deallocate(wk) ; goto 9999 ! only the real space transformation of R is required; 19 Jun 15;
    endif

8000 continue

    SALLOCATE(bb,(1:3,0:3,1:NN),zero) ! magnetic field along trial curve; 19 Jun 15;
    
    if( iflag.eq.1 ) then ; itangent = 0
    else                  ; itangent = 1
    endif

    if( iec00aa.le.-6 ) write(0,'("ec00ca :         "2x" : rr = "99f8.4)') curve%rr
    if( iec00aa.le.-6 ) write(0,'("ec00ca :         "2x" : zz = "99f8.4)') curve%zz
   !if( iec00aa.le.-6 ) pause
    
    do ii = 1, NN

     phi = (ii-1) * (pi2/Nfp) / NN ; RpZ(1:3) = (/ curve%rr(ii), phi, curve%zz(ii) /) ! note that either curve%rr or curve%zz assumed;

     nbfield(itangent) = nbfield(itangent) + 1 ; ifail = -9 ; call bfield( RpZ(1:3), itangent, bb(1:3,0:3,ii), ifail )

     if( ifail.ne.0 ) then
      write(0,'("ec00ca : ifail ="i3" : B(R,p,Z)=("es22.15","es23.15","es23.15" ) = ( ?, ?, ? ) ;")') ifail, RpZ(1:3)
      deallocate(wk,bb)
      curve%ibfield = 1 ; iflag = -1 ; goto 9999
     endif
     
    enddo ! end of do ii = 1, NN ; 17 Jun 15;
    
    SALLOCATE(ff,(1:MM,1:NN),zero)
    
    select case( iflag )
     
    case( 1 )
     
     ff(1,1:NN) = ( + bb(3,0,1:NN) - bb(2,0,1:NN) * curve%zd(1:NN) ) * curve%rr(1:NN)   ! FR; 19 Jun 15;
     ff(3,1:NN) = ( - bb(1,0,1:NN) + bb(2,0,1:NN) * curve%rd(1:NN) ) * curve%rr(1:NN)   ! FZ; 19 Jun 15;
     
    case( 2 )
     
     ff(1,1:NN) = ( + bb(3,1,1:NN) - bb(2,1,1:NN) * curve%zd(1:NN) ) * curve%rr(1:NN) &
                  + bb(3,0,1:NN) - bb(2,0,1:NN) * curve%zd(1:NN)                        ! dFRdR; 19 Jun 15;
     ff(2,1:NN) = ( + bb(3,3,1:NN) - bb(2,3,1:NN) * curve%zd(1:NN) ) * curve%rr(1:NN)                                                ! dFRdZ; 19 Jun 15;
     ff(3,1:NN) = ( - bb(1,1,1:NN) + bb(2,1,1:NN) * curve%rd(1:NN) ) * curve%rr(1:NN) &
                  - bb(1,0,1:NN) + bb(2,0,1:NN) * curve%rd(1:NN)                        ! dFZdR; 19 Jun 15;
     ff(4,1:NN) = ( - bb(1,3,1:NN) + bb(2,3,1:NN) * curve%rd(1:NN) ) * curve%rr(1:NN)                                                ! dFZdZ; 19 Jun 15;
     
     SALLOCATE(gg,(1:MM,1:NN),zero)
     
    end select ! end select case( iflag ) ; 19 Jun 15;
    
    select case( iflag )
     
    case( 1 )
     
#ifdef HAVE_NAG
      ic06fpf = 0
      call C06FPF( MM, NN, ff(1:MM,1:NN), curve%init, curve%trig(1:2*NN), wk(1:MM*NN), ic06fpf ) ! NAG; 13 Oct 15;
#else
      ! apply transform to each row
      do ii = 1, MM
        call rfft1f( NN, 1, ff(ii,1:NN), NN, wk, Lwk, wkb, NN, ic06fpf ) ! this is destructive
      enddo
#endif

     ff(1:MM,1:NN) = ff(1:MM,1:NN) / curve%sN ; ff(1:MM,2:NN) = ff(1:MM,2:NN) * two
     
     if    ( curve%Lrz.eq. 1 .or. curve%Lrz.eq. 3 .or. curve%Lrz.eq.-2 .or. curve%Lrz.eq. -4) then
       kk = 1 ! need to return FR; 19 Jun 15;
     elseif( curve%Lrz.eq. 2 .or. curve%Lrz.eq. 4 .or. curve%Lrz.eq.-1 .or. curve%Lrz.eq. -3) then
       kk = 3 ! need to return FZ; 19 Jun 15;
     endif
     
     ;    ii = 0       ; Fcmn(   ii+1        ) =   ff(kk,   ii+1)
     ;if( Ntor.ge.1 ) then                                        
     ; do ii = 1, Ntor ; Fcmn(   ii+1        ) =   ff(kk,   ii+1)
     ;  ;              ; Fcmn(tN-ii+1        ) = - ff(kk,NN-ii+1)
     ; enddo
     ;endif
     
    case( 2 )
     
     FATAL(ec00ca, curve%Lrz.le.0, need to implement derivatives)
     
     do jj = 0, Ntor ! labels Fourier harmonic; 17 Jun 15;
      
      if    ( curve%Lrz.eq.1 .or. curve%Lrz.eq.2 ) then
       gg(1,1:NN) = ff(1,1:NN) * curve%ct(1:NN,jj)                                                                ! dFRdRc ; 29 Jul 15;
       gg(2,1:NN) = ff(1,1:NN) * curve%st(1:NN,jj)                                                                ! dFRdRs ; 29 Jul 15;
       gg(3,1:NN) = ff(3,1:NN) * curve%ct(1:NN,jj) - bb(2,0,1:NN) * curve%rr(1:NN) * curve%st(1:NN,jj) * jj * Nfp ! dFZdRc ; 29 Jul 15;
       gg(4,1:NN) = ff(3,1:NN) * curve%st(1:NN,jj) + bb(2,0,1:NN) * curve%rr(1:NN) * curve%ct(1:NN,jj) * jj * Nfp ! dFZdRs ; 29 Jul 15;
      elseif( curve%Lrz.eq.3 .or. curve%Lrz.eq.4 ) then
       gg(1,1:NN) = ff(2,1:NN) * curve%ct(1:NN,jj) + bb(2,0,1:NN) * curve%rr(1:NN) * curve%st(1:NN,jj) * jj * Nfp ! dFRdZc ; 29 Jul 15;
       gg(2,1:NN) = ff(2,1:NN) * curve%st(1:NN,jj) - bb(2,0,1:NN) * curve%rr(1:NN) * curve%ct(1:NN,jj) * jj * Nfp ! dFRdZs ; 29 Jul 15;
       gg(3,1:NN) = ff(4,1:NN) * curve%ct(1:NN,jj)                                                                ! dFZdZc ; 29 Jul 15;
       gg(4,1:NN) = ff(4,1:NN) * curve%st(1:NN,jj)                                                                ! dFZdZs ; 29 Jul 15;
      endif
      
#ifdef HAVE_NAG
      ic06fpf = 0
      call C06FPF( MM, NN, gg(1:MM,1:NN), curve%init, curve%trig(1:2*NN), wk(1:MM*NN), ic06fpf ) ! NAG; 13 Oct 15;
#else
      ! apply transform to each row
      do ii = 1, MM
        call rfft1f( NN, 1, gg(ii,1:NN), NN, wk, Lwk, wkb, NN, ic06fpf ) ! this is destructive
      enddo
#endif

      gg(1:MM,1:NN) = gg(1:MM,1:NN) / curve%sN ; gg(1:MM,2:NN) = gg(1:MM,2:NN) * two
      
      if    ( curve%Lrz.eq. 1 .or. curve%Lrz.eq. 3 .or. curve%Lrz.eq.-2 .or. curve%Lrz.eq. -4) then
        kk = 0 ! need to return FR; 19 Jun 15;
      elseif( curve%Lrz.eq. 2 .or. curve%Lrz.eq. 4 .or. curve%Lrz.eq.-1 .or. curve%Lrz.eq. -3) then 
        kk = 2 ! need to return FZ; 19 Jun 15;
      endif
     
      ;   ii = 0       ; Dcmn(   ii+1,   jj+1) =   gg(1+kk,   ii+1)
      ;   if( jj.gt.0 ) then
      ;                ; Dcmn(   ii+1,tN-jj+1) =   gg(2+kk,   ii+1)
      ;   endif
      if( Ntor.gt.0 ) then
       do ii = 1, Ntor ; Dcmn(   ii+1,   jj+1) =   gg(1+kk,   ii+1)
        ; if( jj.gt.0 ) then
        ;              ; Dcmn(   ii+1,tN-jj+1) =   gg(2+kk,   ii+1)
        ; endif
        ;              ; Dcmn(tN-ii+1,   jj+1) = - gg(1+kk,NN-ii+1)
        ; if( jj.gt.0 ) then
        ;              ; Dcmn(tN-ii+1,tN-jj+1) = - gg(2+kk,NN-ii+1)
        ; endif
       enddo ! end of do ii; 17 Jun 15;
      endif
      
     enddo ! end of do jj; 17 Jun 15;
     
     DALLOCATE(gg)
     
    end select
    
    DALLOCATE(bb)
    DALLOCATE(wk)
    
    DALLOCATE(ff)
    
    if( iflag.eq.1 ) then
     curve%its = curve%its + 1 ; curve%err = sqrt( sum( Fcmn(1:Ncmn)**2 ) )
     if( iec00aa.le.-4 ) write(0,'("ec00ca : its = "i4" : |F|="es13.5" ;")') curve%its, curve%err
!    if( curve%err.le.curve%etol   ) iflag = -2
!    if( curve%its.ge.curve%maxits ) iflag = -3
    else
     if( iec00aa.le.-4 ) write(0,'("ec00ca : its = "i4" :  D           ;")') curve%its
    endif
    
9999 continue
    
    return
    
  end subroutine ec00ca

!latex \subroutine{tr00aa}{measure rotational-transform;}
!latex \bi
!latex \item[1.] Fieldline tracing methods are used to determine the relative rotational-transform of one fieldline about a given ``reference'' fieldline,
!latex           which will usually be a magnetic axis.
!latex \item[ *] The equations governing the fieldlines are the same as that given in \Eqn{BR} and \Eqn{BZ}.
!latex \item[ *] The user must supply {\em two} starting points, $(R_a,Z_a)$ and $(R,Z)$.
!latex           The location of the magnetic axis, $(R_a,Z_a)$, can be obtained from a previous call to \verb+ga00aa+, see \Sec{ga00aa}.
!latex           however, the values of $(R_a,Z_a)$ and $(R,Z)$ are completely arbitrary.
!latex           What is really measured by this routine is average ``linking'' of one fieldline about the other.
!latex \item[ *] A poloidal angle, $\t(\phi)$, is introduced as 
!latex           \be \tan\t(\phi)=\frac{\delta Z(\phi)}{\delta R(\phi)}, \ee
!latex           where $\delta R(\phi) = R(\phi) - R_a(\phi)$ and $\delta Z(\phi) = Z(\phi) - Z_a(\phi)$.
!latex \item[ *] This angle varies with $\phi$ according to 
!latex           \be \frac{d\t}{d\phi} = \frac{\delta R\,(Z^\prime-Z_a^\prime) - \delta Z\,(R^\prime-R_a^\prime)}{\delta R^2 + \delta Z^2}, 
!latex           \label{eq:diotadphi}
!latex           \ee
!latex           where $\prime$ denotes total derivative with respect to $\phi$.
!latex \item[ *] The o.d.e. integration defined in \Eqn{diotadphi} may be initialized with $\t(0)=0$,
!latex           and after a sufficiently large distance, $\Delta\phi$,
!latex           this angle satisfies $\Delta \t \approx \iotabar \Delta\phi$, where $\iotabar$ is the rotational-transform.
!latex \item[ *] A more accurate calculation of $\iotabar$ is enabled by fitting a straight line to $\t(\phi)$, rather than just subtracting the endpoints.
!latex           This will be implemented in time, upon request, . . . 
!latex \item[ *] Formally, the rotational-transform is defined as the limit
!latex           \be \iotabar \equiv \lim_{\Delta \phi \rightarrow \infty} \frac{\Delta \t}{\Delta \phi}.
!latex           \ee
!latex           This limit {\em only} converges on regular fieldlines. 
!latex           For irregular, or chaotic, fieldlines, this limit does not converge and the rotational-transform is not defined!
!latex \item[ *] Note that \Eqn{diotadphi} requires knowledge of $R_a(\phi)$ and $Z_a(\phi)$, and $R(\phi)$ and $Z(\phi)$,
!latex           and these are obtained by integrating
!latex           \Eqn{BR} and \Eqn{BZ}. So, in total there are $5$ coupled o.d.e.s.
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : transformdata, tr00aa+ \\ \\
!latex           \verb+type(transformdata)    :: transform+ \\ \\ in their source that calls \verb+tr00aa+,
!latex           where \verb+transform+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+transform+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+transform%Nfp              : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, e.g. \verb+Nfp=1+ for tokamaks, \verb+Nfp=5+ for LHD, . . ;
!latex \ei
!latex \item[  ] \verb+transform%Ppts             : integer ;+
!latex \bi
!latex \item[i.] the number of toroidal transits that the fieldlines will be followed, e.g. \verb+Ppts>100+;
!latex \ei
!latex \item[  ] \verb+transform%odetol           : real    ;+
!latex \bi
!latex \item[i.] o.d.e. integration tolerance; e.g. \verb+odetol=1.0e-08+;
!latex \ei
!latex \item[  ] \verb+transform%Ra               : real    ;+
!latex \item[  ] \verb+transform%Za               : real    ;+
!latex \item[  ] \verb+transform%R                : real    ;+
!latex \item[  ] \verb+transform%Z                : real    ;+
!latex \bi
!latex \item[i.] starting points for fieldline integrations;
!latex \item[ii.] usually, $(R_a,Z_a)$ will be the location of the magnetic axis on the $\phi=0$ plane, and $(R,Z)$ is arbitrary.
!latex \item[iii.] the starting points must be distinct: in particular, \verb+length+ $\equiv(R-R_a)^2+(Z-Z_a)^2$ must exceed $10^{-12}$.
!latex \ei
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call tr00aa( transform, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+transform%iota             : real    ;+
!latex \bi
!latex \item[i.] the ``rotational-transform'' of the fieldline starting at $(R,Z)$ relative to the fieldline starting at $(R_a,Z_a)$.
!latex \ei
!latex \item[  ] \verb+transform%Lallocated       : integer ;+
!latex \bi
!latex \item[i.] if \verb+Lallocated+=1, the \verb+transform%RZ+ array has           been allocated;
!latex \item[ii.] if \verb+Lallocated+=0, the \verb+transform%RZ+ array has {\bf not} been allocated;
!latex \ei
!latex \item[  ] \verb+transform%RZ(1:2,0:Ppts)   : real ;+
!latex \bi
!latex \item [i.] the \Poincare plot data: $R_i\equiv$\verb+transform%RZ(1,i)+, $Z_i\equiv$\verb+transform%RZ(2,i)+
!latex \item [ii.] the \verb+transform%RZ+ need not be allocated on input; if it is, it is first deallocated;
!latex \ei
!latex \item[  ] \verb+ifail                       : integer ;+
!latex \bi
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error; 
!latex     \item[] \verb+ifail=2+ : the NAG routine \verb+D02BJF+ failed to integrate along the fieldline.
!latex                              perhaps because the fieldline left the computational domain, . . .
!latex \ei
!latex \ei
!latex \ei
  
  subroutine tr00aa( ltransform, ifail )
    
    implicit none
    
    type(transformdata) :: ltransform
    INTEGER             :: ifail
    
    INTEGER, parameter  :: Node = 5, Lwk = 20 * Node
    INTEGER             :: id02bjf, ipoint, astat
    REAL                :: RZRZt(1:Node), phistart, phiend, length, leastfit(1:5), dzeta
    CHARACTER           :: relabs
    
#ifdef HAVE_NAG
    REAL                :: wk(1:Lwk)
    external            :: D02BJX, D02BJW
#else
    INTEGER, parameter  :: liw=20, lrw=20+16*Node
    INTEGER             :: iwork(liw)
    REAL                :: rwork(lrw)
    INTEGER             :: iopt,istate,itask,itol,mf
    REAL                :: atol,rtol
#endif

    itr00aa = ifail ; nbfield(0:1) = 0

    if( allocated(ltransform%RZ) ) deallocate(ltransform%RZ)
    
    ltransform%Lallocated = 0 ; ltransform%iota = zero
    
    CHECKINPUT(tr00aa, ltransform%Nfp   .le.0       , 9999 )
    CHECKINPUT(tr00aa, ltransform%Ppts  .le.0       , 9999 )
    CHECKINPUT(tr00aa, ltransform%odetol.le.zero    , 9999 )
    CHECKINPUT(tr00aa, ltransform%Ra    .le.zero    , 9999 )
    CHECKINPUT(tr00aa, ltransform%R     .le.zero    , 9999 )

    length = (ltransform%R-ltransform%Ra)**2 + (ltransform%Z-ltransform%Za)**2

    CHECKINPUT(tr00aa, length           .lt.small   , 9999 )

    relabs = 'D' ; itangent = 0 ; dzeta = pi2/ltransform%Nfp
    
    SALLOCATE(ltransform%RZ,(1:2,0:ltransform%Ppts),zero) ! for block writing to file; 25 Mar 15;
    
    ltransform%Lallocated = 1

    ipoint = 0
    
    RZRZt(1:5) = (/ ltransform%Ra, ltransform%Za, ltransform%R, ltransform%Z, zero /) ! is this the correct initial angle; 20 Jan 16;
    
    ltransform%RZ(1:2,ipoint) = RZRZt(3:4)

    leastfit(1:5) = (/ zero, zero, zero, zero, one /)

    do ipoint = 1, ltransform%Ppts
     
     phistart = zero ; phiend = dzeta
     
     Lbfieldok = .true.
     
#ifdef HAVE_NAG
     id02bjf = 1
     call D02BJF( phistart, phiend, Node, RZRZt(1:Node), bf00ab, ltransform%odetol, relabs, D02BJX, D02BJW, &
                  wk(1:Lwk), id02bjf ) ! NAG ode integration;
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=ltransform%odetol
!    atol :  absolute tolerance
     atol=ltransform%odetol
!    do integration
     call lsode(bf00ab,(/Node/),RZRZt,phistart,phiend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,du00aa,mf)
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif
     
     if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
     
     if( id02bjf.ne.0 ) then ; itr00aa = 2 ; goto 9999
     endif
     
     ltransform%RZ(1:2,ipoint) = RZRZt(3:4)
     
     leastfit(1:5) = leastfit(1:5) + (/ (ipoint*dzeta)**2, (ipoint*dzeta), (ipoint*dzeta)*RZRZt(5), RZRZt(5), one /)

    enddo ! end of do ipoint; 20 Jan 16;
    
   !ltransform%iota = RZRZt(5) / ( ltransform%Ppts * dzeta )
    ltransform%iota = ( leastfit(5)*leastfit(3)-leastfit(2)*leastfit(4) ) &
                    / ( leastfit(5)*leastfit(1)-leastfit(2)*leastfit(2) )

    itr00aa = 0

9999 continue

    ltransform%nbfield(0:1) = nbfield(0:1)
    
    if( ifail.le.-1 ) then
      write(0,9000) itr00aa, ltransform%Ra, ltransform%Za, ltransform%R, ltransform%Z, ltransform%Ppts, ltransform%iota
    endif
    
9000 format("tr00aa : ifail ="i3" ; (Ra,Za)=("es23.15" ,"es23.15" ) ; (R ,Z )=("es23.15" ,"es23.15" ) ; Ppts="i6 &
            " ; iota="es23.15" ;")
    
    ifail = itr00aa
    
    return
    
  end subroutine tr00aa

!latex \subroutine{pp00aa}{fieldline tracing for \Poincare plot, calculate Lyapunov exponent;}
!latex \bi
!latex \item[1.] This subroutine follows a magnetic fieldline from a given starting point, for a given number of toroidal periods, 
!latex           either forwards or backwards in the toroidal angle, $\phi$.
!latex           The intersection points of the fieldline with the \Poincare section $\phi=0$ are returned.
!latex \item[* ] The (maximum) Lyapunov exponent, $\lambda$, may also be calculated 
!latex           [Benettin, Galgani \& Strelcyn, \link{dx.doi.org/10.1103/PhysRevA.14.2338}{Phys. Rev. A, 14:2338 (1976)}].
!latex           This measures the average exponential rate of separation of a nearby trajectory, 
!latex           \be |\delta {\bf x}(\phi)| = e^{\lambda \phi} |\delta {\bf x}(0)|,
!latex           \ee
!latex           as $\phi\rightarrow\infty$ and as $|\delta {\bf x}(0)|\rightarrow 0$.
!latex \item[* ] The limit $|{\delta \bf x}(0)|\rightarrow 0$ is best treated by linearizing the fieldline equations about the given reference fieldline, 
!latex           as is herafter assumed.
!latex \item[* ] To calculate $\lambda$, assuming that $\delta {\bf x}$ lies in the tangent space and $|{\delta\bf x}(0)|=1$, as
!latex           \be \lambda(\phi) = \frac{1}{\phi} \log(|\delta {\bf x}(\phi)|). \label{eq:Lyapunov}
!latex           \ee
!latex \item[* ] In \Eqn{Lyapunov}, $\lambda$ has been expressed as a function of $\phi$ 
!latex           so that the user may examine whether the limit $\lim_{\phi\rightarrow \infty} \lambda(\phi)$ has converged.
!latex \item[* ] If $\lambda>0$, the long-time trajectory of the fieldline is infinitesimally sensitive to the initial position,
!latex           and this is a defining characteristic of chaos; however, not all fieldlines are chaotic.
!latex           Consider a nearby fieldline that separates linearly, i.e. $|\delta {\bf x}(\phi)| = 1 + c \phi$ for some constant, $c$; for example,
!latex           a fieldline that lies on a nearby flux surface with slightly-different rotational-transform.
!latex           The value of $\lambda$ given by \Eqn{Lyapunov} gives
!latex           \be \lambda \approx \frac{\log c}{\phi} + \frac{\log \phi}{\phi}. \label{eq:linearLyapunov}
!latex           \ee
!latex \item[* ] To distinguish a weakly-exponentially-separating fieldline from a linearly-separating fieldline may require a very long integration in $\phi$.
!latex           It is recommended to plot $\log(|\lambda(\phi)|)$ against $\log(\phi)$ 
!latex           to determine if $\lambda(\phi)$ is approaching a non-zero value as $\phi\rightarrow\infty$; 
!latex           and it may be useful to compare this to what would be expected, \Eqn{linearLyapunov}, for linearly separating fieldlines.
!latex \item[* ] Also, there maybe nearby fieldlines that neither exponentially separate nor linearly separate,
!latex           but instead oscillate about the reference fieldline.
!latex           Such behavior is displayed by fieldlines in the vicinity of a stable fixed point.
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : poincaredata, pp00aa+ \\ \\
!latex           \verb+type(poincaredata) :: poincare+ \\ \\ in their source that calls \verb+pp00aa+,
!latex           where \verb+poincaredata+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+poincare+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+poincare%Nfp            : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, 
!latex \item[ii.] e.g. \verb+Nfp=1+;
!latex \ei
!latex \item[  ] \verb+poincare%R              : real    ;+
!latex \bi
!latex \item[i.] starting point;
!latex \ei
!latex \item[  ] \verb+poincare%Z              : real    ;+
!latex \bi
!latex \item[i.] starting point;
!latex \ei
!latex \item[  ] \verb+poincare%Ppts           : integer ;+
!latex \bi
!latex \item[i.] toroidal periods / iterations;
!latex \ei
!latex \item[  ] \verb+poincare%idirection     : integer ;+
!latex \bi
!latex \item[i.] \verb+idirection= 1+ : follow fieldline in increasing toroidal direction;
!latex \item[i.] \verb+idirection=-1+ : follow fieldline in decreasing toroidal direction;
!latex \ei
!latex \item[  ] \verb+poincare%flparameter    : integer ;+
!latex \bi
!latex \item[i.] \verb+flparameter=0+ : use toroidal angle to parameterize distance along fieldline; i.e. divide equations by $B^\zeta$;
!latex \item[i.] \verb+flparameter=1+ : use length to parameterize distance along fieldline; i.e. divide equations by $|B|$; under construction;
!latex \ei
!latex \item[  ] \verb+poincare%iLyapunov      : integer ;+
!latex \bi
!latex \item[i.]  \verb+iLyapunov = 0+ : the Lyapunov exponent will not be calculated;
!latex \item[ii.] \verb+iLyapunov = 1+ : the Lyapunov exponent      not be calculated;
!latex                                  note that this will slow the fieldline integration, 
!latex                                  as additional o.d.e.s that define the tangent mapping need to be integrated;
!latex \ei
!latex \item[  ] \verb+poincare%Ltol           : real    ;+
!latex \bi
!latex \item[i.] tolerance required for calculation of Lyapunov exponent; 
!latex \item[ii.] e.g. \verb+Ltol=1.0e-03+;
!latex \item[iii.] under construction: I plan to automatically adjust the integration length (i.e. \verb+Ppts+) 
!latex             to ensure that $\lambda$ is converged to within \verb+Ltol+.
!latex \ei
!latex \item[  ] \verb+poincare%odetol         : real    ;+
!latex \bi
!latex \item[i.] o.d.e. integration tolerance; 
!latex \item[ii.] e.g. \verb+odetol=1.0e-08+;
!latex \ei
!latex \item[  ] \verb+ifail                   : integer ;+
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call pp00aa( poincare, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+poincare%RZ(1:2,0:Ppts) : real    ;+
!latex \bi
!latex \item[i.] \Poincare data;
!latex \item[ii.] $R_i = $ \verb+RZ(1,0:Ppts)+; $Z_i = $ \verb+RZ(2,0:Ppts)+;
!latex \item[i1i.] \verb+RZ+ need not be allocated on input; if it is, it is immediately deallocated;
!latex \ei
!latex \item[  ] \verb+poincare%Ly(1:Ppts)      : real    ;+
!latex \bi
!latex \item[i.] ``evolution" of Lyapunov exponent; i.e. the estimate of the Lyapunov exponent as a function of iteration;
!latex \item[iv.] \verb+Ly+ need not be allocated on input; if it is, it is immediately deallocated;
!latex \ei
!latex \item[  ] \verb+poincare%Lallocated     : integer ;+
!latex \bi
!latex \item[i.] if \verb+Lallocated=1+, the \verb+RZ+ and \verb+Ly+ arrays have     been allocated;
!latex \item[ii.] if \verb+Lallocated=0+, the \verb+RZ+ and \verb+Ly+ arrays have not been allocated;
!latex \ei
!latex \item[  ] \verb+poincare%Lyapunov       : real    ;+
!latex \bi
!latex \item[i.] the Lyapunov exponent;
!latex \ei
!latex \item[  ] \verb+poincare%ipts           : integer ;+
!latex \bi
!latex \item[i.] toroidal periods actually followed: if an error is encountered (e.g. the fieldline leaves the computational domain,
!latex           the fieldline integration is terminated;
!latex \ei
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex \ei
!latex \item[6.] Comments:
!latex \item[* ] The NAG routine \nag{www.nag.co.uk/numeric/FL/manual19/pdf/D02/d02bjf_fl19.pdf}{D02BJF} is used to perform the o.d.e. integration.
!latex \ei
  
  subroutine pp00aa( lpoincare, ifail )
    
    implicit none
    
    type(poincaredata) :: lpoincare
    INTEGER            :: ifail
    
    INTEGER, parameter :: Node = 7, Lwk = 20 * Node
    INTEGER            :: id02bjf, ii, astat
    REAL               :: RZp(1:Node), timestart, timeend, tv(1:2), Lyapunovsum ! tv = tangent vector; 10 Apr 16;
    CHARACTER          :: relabs
    
#ifdef HAVE_NAG
    REAL               :: wk(1:Lwk)
    external           :: D02BJX
#else
    INTEGER, parameter :: liw=20, lrw=20+16*Node
    INTEGER            :: iwork(liw)
    REAL               :: rwork(lrw)
    INTEGER            :: iopt,istate,itask,itol,jt,ng,jroot
    REAL               :: atol,rtol
#endif    
    ipp00aa = ifail ; nbfield(0:1) = 0
    
    if( allocated(lpoincare%RZ) ) deallocate(lpoincare%RZ)
    if( allocated(lpoincare%Ly) ) deallocate(lpoincare%Ly)
    
    lpoincare%Lallocated = 0 ; lpoincare%Lyapunov = zero
    
    CHECKINPUT(pp00aa, lpoincare%Nfp        .le. 0                                       , 9999 )
    CHECKINPUT(pp00aa, lpoincare%Ppts       .le. 0                                       , 9999 )
    CHECKINPUT(pp00aa, lpoincare%odetol     .le.small                                    , 9999 )
    CHECKINPUT(pp00aa, lpoincare%idirection .ne.-1   .and. lpoincare%idirection .ne. +1  , 9999 )
    CHECKINPUT(pp00aa, lpoincare%iLyapunov  .ne. 0   .and. lpoincare%iLyapunov  .ne. +1  , 9999 )
    CHECKINPUT(pp00aa, lpoincare%iLyapunov  .eq. 1   .and. lpoincare%Ltol       .le.small, 9999 )
    CHECKINPUT(pp00aa, lpoincare%flparameter.ne. 0   .and. lpoincare%flparameter.ne. 1   , 9999 )

    relabs = 'D' ; itangent = lpoincare%iLyapunov

    poincare%phi = lpoincare%phi ; poincare%idirection = lpoincare%idirection
    poincare%Nfp = lpoincare%Nfp ; poincare%flparameter = lpoincare%flparameter
    
    SALLOCATE(lpoincare%RZ,(1:2,0:lpoincare%Ppts),zero)
    SALLOCATE(lpoincare%Ly,(    1:lpoincare%Ppts),zero)
    
    lpoincare%Lallocated = 1
    
    ii = 0 ; lpoincare%ipts = ii

    RZp(1:3) = (/ lpoincare%R, lpoincare%Z, lpoincare%phi /)
    tv(1:2) = (/ sqrt(half), sqrt(half) /) ; lyapunovsum = zero
    
    lpoincare%RZ(1:2,ii) = RZp(1:2)
    
    do ii = 1, lpoincare%Ppts ; lpoincare%ipts = ii ; Lbfieldok = .true.
     
     if( poincare%flparameter.eq.0 ) then
      timestart = zero ; timeend = poincare%idirection * (pi2/poincare%Nfp)
     else
      timestart = zero ; timeend = poincare%idirection * (pi2/poincare%Nfp) * lpoincare%R * ten
     endif

     Lbfieldok = .true. ; poincare%phistart = RZp(3)
     
     RZp(4:Node) = (/ one, zero, zero, one  /) ! re-set tangent; 10 Apr 16;
     
#ifdef HAVE_NAG
     id02bjf = 1
     call D02BJF( timestart, timeend, Node, RZp(1:Node), bf00aa, lpoincare%odetol, relabs, D02BJX, pp00ab, &
                  wk(1:Lwk), id02bjf ) ! NAG ode integration;
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    jt=2 :  internally generated full Jacobian
     jt=2
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=lpoincare%odetol
!    atol :  absolute tolerance
     atol=lpoincare%odetol
!    ng : number of roots being found
     ng=1
!    do integration
     call lsodar(bf00aa,(/Node/),RZp,timestart,timeend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,du00aa,jt,pp00ab,ng,jroot)
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif

     if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; an error has occured in user-supplied bfield;
     
     if( id02bjf.ne.0 ) then
      lpoincare%RZ(1,ii:lpoincare%Ppts) = lpoincare%RZ(1,ii-1)
      lpoincare%RZ(2,ii:lpoincare%Ppts) = lpoincare%RZ(2,ii-1)
      if( ii.gt.1 ) then
      lpoincare%Ly(  ii:lpoincare%Ppts) = lpoincare%Ly(  ii-1) 
      endif
      ipp00aa = 1 ; goto 9999 ! complete array with dummy data; 02 Jun 15;
     endif
     
     lpoincare%RZ(1:2,ii) = RZp(1:2)
     
     if( lpoincare%flparameter.eq.0 .and. lpoincare%iLyapunov.eq.1 ) then
      tv(1:2) = tv(1:2) / sqrt( tv(1)**2 + tv(2)**2 ) ! re-normalize; 13 Jun 15;
      tv(1:2) = (/ RZp(4)*tv(1) + RZp(5)*tv(2), RZp(6)*tv(1) + RZp(7)*tv(2) /) ! tangent mapping matrix vector multiplication; 13 Jun 15;
      lyapunovsum = lyapunovsum + log( sqrt( tv(1)**2 + tv(2)**2 ) ) 
      lpoincare%Lyapunov = Lyapunovsum / ( ii*(pi2/lpoincare%Nfp) )
      lpoincare%Ly(ii) = lpoincare%Lyapunov
     endif
     
    enddo ! end of do ii ; 02 Jun 15;
    
    ipp00aa = 0
    
9999 continue

    lpoincare%nbfield(0:1) = nbfield(0:1)
    
    if( ifail.le.-1 ) then
      write(0,'("pp00aa : ifail ="i3" : p="f15.10" ; (R,Z)=("f15.10" ,"f15.10" ) ; Ppts="i7" ; lyapunov="es12.5" ;")') &
        ipp00aa, lpoincare%phi, lpoincare%R, lpoincare%Z, lpoincare%Ppts, lpoincare%Lyapunov 
    endif

    ifail = ipp00aa

    return
    
  end subroutine pp00aa

#ifdef HAVE_NAG
  REAL function pp00ab( time, RZp )
    
    implicit none
    
    INTEGER, parameter  :: Node = 7
    REAL                :: time, RZp(1:Node)

   !pp00ab = RZp(3) - ( poincare%phistart + poincare%idirection * (pi2/poincare%Nfp) )
    pp00ab = time -   (                     poincare%idirection * (pi2/poincare%Nfp) ) ! 10 Apr 16;
    
    return
    
  end function pp00ab
#else    
  subroutine pp00ab( Node, time, RZp, NB, BRZp ) ! format constrained by lsodar

    implicit none

    INTEGER, intent(in) :: Node, NB
    REAL, intent(in)    :: time, RZp(1:Node)
    REAL, intent(out)   :: BRZp

   !BRZp = RZp(3) - ( poincare%phistart + poincare%idirection * (pi2/poincare%Nfp) )
    BRZp = time -   (                     poincare%idirection * (pi2/poincare%Nfp) ) ! 10 Apr 16;
    
    return
    
  end subroutine pp00ab
#endif

!latex \subroutine{gc00aa}{follow guiding center;}
!latex \bi
!latex \item[1.] This is under construction. Contact shudson@pppl.gov.
!latex \ei
  
!latex \subroutine{rz00aa}{construct cylindrical Fourier harmonics of flux surface using fieldline tracing;}
!latex \bi
!latex \item[1.] This subroutine follows a magnetic fieldline from a given starting point, for a given number of toroidal periods, 
!latex           to construct the Fourier coefficients that satisfy
!latex           \be R = \sum_i R_i \cos(m_i\t-n_i\z) \\
!latex               Z = \sum_i Z_i \sin(m_i\t-n_i\z)
!latex           \ee
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : poincaredata, rz00aa+ \\ \\
!latex           \verb+type(poincaredata) :: poincare+ \\ \\ in their source that calls \verb+rz00aa+,
!latex           where \verb+poincaredata+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+poincare+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+poincare%Nfp            : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, 
!latex \item[ii.] e.g. \verb+Nfp=1+;
!latex \ei
!latex \item[  ] \verb+poincare%R              : real    ;+
!latex \bi
!latex \item[i.] starting point;
!latex \ei
!latex \item[  ] \verb+poincare%Z              : real    ;+
!latex \bi
!latex \item[i.] starting point;
!latex \ei
!latex \item[  ] \verb+poincare%Ppts           : integer ;+
!latex \bi
!latex \item[i.] toroidal periods / iterations;
!latex \ei
!latex \item[  ] \verb+poincare%odetol         : real    ;+
!latex \bi
!latex \item[i.] o.d.e. integration tolerance; 
!latex \item[ii.] e.g. \verb+odetol=1.0e-08+;
!latex \ei
!latex \item[  ] \verb+ifail                   : integer ;+
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call rz00aa( poincare, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex \ei
!latex \item[6.] Comments:
!latex \item[* ] 
!latex \ei
  
  subroutine rz00aa( lrzdata, ifail )
    
    implicit none
    
    type(rzcoordsdata) :: lrzdata
    INTEGER            :: ifail
    
    INTEGER, parameter :: Node = 6, Lwk = 20 * Node
    INTEGER            :: id02bjf, ii, astat, Nfp
    REAL               :: RZ(1:Node), phistart, phiend
    CHARACTER          :: relabs
    
#ifdef HAVE_NAG
    REAL               :: wk(1:Lwk)
    external           :: D02BJX, D02BJW
#else
    INTEGER, parameter :: liw=20, lrw=20+16*Node
    INTEGER            :: iwork(liw)
    REAL               :: rwork(lrw)
    INTEGER            :: iopt,istate,itask,itol,mf
    REAL               :: atol,rtol
#endif

!    irz00aa = ifail ; nbfield(0:1) = 0 ; Nfp = lpoincare%Nfp
!    
!    if( allocated(lpoincare%RZ) ) deallocate(lpoincare%RZ)
!    if( allocated(lpoincare%Ly) ) deallocate(lpoincare%Ly)
!    
!    lpoincare%Lallocated = 0 ; lpoincare%Lyapunov = zero
!    
!    if( lpoincare%Ppts      .le. 0     ) then 
!      write(0,'("rz00aa :             ;   Ppts.le.0        ;")') irz00aa = 0 ; goto 9999
!    endif
!    if( lpoincare%Nfp       .le. 0     ) then
!      write(0,'("rz00aa : input error ;    Nfp.le.0        ;")') ; irz00aa = 1 ; goto 9999
!    endif
!    if( lpoincare%odetol    .le. small ) then
!      write(0,'("rz00aa :             ; odetol.le.small    ;")') ; irz00aa = 1 ; goto 9999
!    endif
!    if( lpoincare%idirection.ne. -1 .and. &
!        lpoincare%idirection.ne. +1    ) then
!      write(0,'("rz00aa :             ; idirection.ne.\pm 1;")') ; irz00aa = 1 ; goto 9999
!    endif
!    if( lpoincare%iLyapunov .ne.  0 .and. &
!        lpoincare%iLyapunov .ne. +1    ) then
!      write(0,'("rz00aa :             ; iLyapunov .ne. 0, 1;")') ; irz00aa = 1 ; goto 9999
!    endif
!    if( lpoincare%iLyapunov .eq.  1 .and. &
!        lpoincare%Ltol      .le. small ) then
!      write(0,'("rz00aa :             ;   Ltol.le.small    ;")') ; irz00aa = 1 ; goto 9999
!    endif
!    
!    relabs = 'D' ; itangent = lpoincare%iLyapunov
!    
!    SALLOCATE(lpoincare%RZ,(1:2,0:lpoincare%Ppts),zero)
!    SALLOCATE(lpoincare%Ly,(    1:lpoincare%Ppts),zero)
!    
!    lpoincare%Lallocated = 1
!    
!    ii = 0 ; lpoincare%ipts = ii
!    
!    RZ(1:2) = (/ lpoincare%R, lpoincare%Z /) ; tv(1:2) = (/ sqrt(half), sqrt(half) /) ; lyapunovsum = zero
!    
!    lpoincare%RZ(1:2,ii) = RZ(1:2)
!    
!    do ii = 1, lpoincare%Ppts ; lpoincare%ipts = ii
!     
!     phistart = zero ; phiend = lpoincare%idirection * (pi2/Nfp) ; Lbfieldok = .true.
!     
!     Lbfieldok = .true. ; RZ(3:Node) = (/ one, zero, zero, one /)
!     
!#ifdef HAVE_NAG
!     id02bjf = 1
!     call D02BJF( phistart, phiend, Node, RZ(1:Node), bf00aa, lpoincare%odetol, relabs, D02BJX, D02BJW, &
!                  wk(1:Lwk), id02bjf ) ! NAG ode integration;
!#else
!!    set the integrator parameters.
!     rwork=0.
!     iwork=0
!!    istate=1 :  indicates the first lsode call
!     istate=1
!!    itask=4 :  normal integration with limited over-shoot (set by
!!               rwork(1) in the i_lsode loop
!     itask=1
!!    iopt=1 :  optional inputs are used
!     iopt=1
!!    rwork(6) :  set maximum lsode-internal step size
!     iwork(6)=400000
!!    rwork(7) :  set minimum lsode-internal step size
!!    iwork(6) :  set maximum lsode-internal steps
!!    iwork(7) :  set maximum lsode-internal error messages printed
!!    mf=10 :  non-stiff Adams method of integration
!     mf=10
!!    itol=1 :  indicates absolute tolerance is just a scalar
!     itol=1
!!    rtol :  relative tolerance
!     rtol=lpoincare%odetol
!!    atol :  absolute tolerance
!     atol=lpoincare%odetol
!!    do integration
!     call lsode(bf00aa,(/Node/),RZ,phistart,phiend, &
!                itol,(/rtol/),(/atol/),itask,istate,iopt, &
!                rwork,lrw,iwork,liw,du00aa,mf)
!     id02bjf=istate
!     if (istate>0) id02bjf=0
!#endif
!
!     
!     if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; an error has occured in user-supplied bfield;
!     
!     if( id02bjf.ne.0 ) then
!      lpoincare%RZ(1,ii:lpoincare%Ppts) = lpoincare%RZ(1,ii-1)
!      lpoincare%RZ(2,ii:lpoincare%Ppts) = lpoincare%RZ(2,ii-1)
!      if( ii.gt.1 ) then
!      lpoincare%Ly(  ii:lpoincare%Ppts) = lpoincare%Ly(  ii-1) 
!      endif
!      irz00aa = 1 ; goto 9999 ! complete array with dummy data; 02 Jun 15;
!     endif
!     
!     lpoincare%RZ(1:2,ii) = RZ(1:2)
!     
!     if( lpoincare%iLyapunov.eq.1 ) then
!      tv(1:2) = tv(1:2) / sqrt( tv(1)**2 + tv(2)**2 ) ! re-normalize; 13 Jun 15;
!      tv(1:2) = (/ RZ(3)*tv(1) + RZ(4)*tv(2), RZ(5)*tv(1) + RZ(6)*tv(2) /) ! tangent mapping matrix vector multiplication; 13 Jun 15;
!      lyapunovsum = lyapunovsum + log( sqrt( tv(1)**2 + tv(2)**2 ) ) 
!      lpoincare%Lyapunov = Lyapunovsum / ( ii*(pi2/Nfp) )
!      lpoincare%Ly(ii) = lpoincare%Lyapunov
!     endif
!     
!    enddo ! end of do ii ; 02 Jun 15;
    
    irz00aa = 0
    
9999 continue
    
    if( ifail.le.-1 ) write(0,'("rz00aa : ifail ="3x" : (R,Z)=("f15.10" ,"f15.10" ) ; Ppts="i6" ;")') &
                                          irz00aa, lrzdata%R, lrzdata%Z,      lrzdata%Ppts

    ifail = irz00aa

    return
    
  end subroutine rz00aa
  
#ifdef FH00AA

!l tex \subroutine{fh00aa}{fit background coordinates to partial separatrix}
!l tex \bi
!l tex \item[1.] This routine is ready for use; contact \verb+shudson@pppl.gov+ for documentation.
!l tex \ei
  
  subroutine fh00aa( laxis, ltangle, lbcoords, ifail )
    
    implicit none
    
    type(magneticaxis)     :: laxis
    type(homoclinictangle) :: ltangle
    type(backgroundcoords) :: lbcoords

    INTEGER                :: ifail

    LOGICAL                :: svd
    INTEGER                :: pplobe, ns(1:2), ius, ii, jj, kk, ll, astat, ipp00aa, Nh, Nc, NN, MM, NRA
    INTEGER                :: ie, je, irank, Lwk, if04jgf
    REAL                   :: dUS, levalue, ya, yb, xx, yy, dl
    REAL, allocatable      :: Mc(:,:), Md(:), mat(:,:), rhs(:), wk(:)
    
    ifh00aa = ifail ; nbfield(0:1) = 0

    if( allocated(lbcoords%sep) ) deallocate(lbcoords%sep)
    if( allocated(lbcoords%ii ) ) deallocate(lbcoords%ii )
    if( allocated(lbcoords%jj ) ) deallocate(lbcoords%jj )
    if( allocated(lbcoords%hi ) ) deallocate(lbcoords%hi )

    lbcoords%Lallocated = 0 ; lbcoords%Nh = 0 ; lbcoords%nsep(1:2) = 0 ; lbcoords%sigma = zero

    if( laxis%R          .le.0    ) then
      write(0,'("fh00aa : input error ;   axis%R.le.0    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%Nfp      .le.0    ) then
      write(0,'("fh00aa : input error ;      Nfp.le.0    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%odetol   .le.zero ) then 
      write(0,'("fh00aa : input error ;   odetol.le.zero ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%R        .le.zero ) then
      write(0,'("fh00aa : input error ; tangle%R.le.0    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%dU       .le.zero ) then
      write(0,'("fh00aa : input error ;       dU.le.0    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%dS       .le.zero ) then
      write(0,'("fh00aa : input error ;       dS.le.0    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%ilobe(1) .le.0    ) then
      write(0,'("fh00aa : input error ; ilobe(1).le.0    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%ilobe(2) .le.0    ) then
      write(0,'("fh00aa : input error ; ilobe(2).le.0    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%wr(1)    .le.zero ) then
      write(0,'("fh00aa : input error ;    wr(1).le.zero ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( ltangle%wr(2)    .le.zero ) then
      write(0,'("fh00aa : input error ;    wr(2).le.zero ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( lbcoords%pplobe  .le.0    ) then
      write(0,'("fh00aa : input error ;   pplobe.le.0    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( lbcoords%Nx      .le.1    ) then
      write(0,'("fh00aa : input error ;       Nx.le.1    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( lbcoords%Ny      .le.1    ) then
      write(0,'("fh00aa : input error ;       Ny.le.1    ;")') ; ifh00aa = 1 ; goto 9999
    endif
    if( lbcoords%svdtol  .le.zero ) then
      write(0,'("fh00aa : input error ;   svdtol.le.zero ;")') ; ifh00aa = 1 ; goto 9999
    endif
    
    poincare%Nfp = ltangle%Nfp ; poincare%odetol = ltangle%odetol ; pplobe = lbcoords%pplobe ; poincare%iLyapunov = 0
    
    lbcoords%nsep(1:2) = pplobe * ltangle%ilobe(1:2) ; ns(1:2) = lbcoords%nsep(1:2)

    SALLOCATE(lbcoords%sep,(1:2,0:maxval(ns(1:2))-1,1:2),zero)
        
    lbcoords%Nh = ( lbcoords%Nx + 1 ) * ( lbcoords%Ny + 1 ) ; Nh = lbcoords%Nh ! polynomial coefficients; 02 Jun 15;
    
    SALLOCATE(lbcoords%ii,(1:Nh),0)
    SALLOCATE(lbcoords%jj,(1:Nh),0)
    SALLOCATE(lbcoords%hi,(1:Nh),zero)
    
    ll = 0
    do kk = 0, lbcoords%Nx + lbcoords%Ny
     do jj = 0, lbcoords%Ny
      do ii = 0, lbcoords%Nx
       if( ii+jj.eq.kk ) then ; ll = ll + 1 ; lbcoords%ii(ll) = ii ; lbcoords%jj(ll) = jj ! order polynomial coefficients as a list; 13 Jun 15;
       endif
      enddo
     enddo
    enddo

    lbcoords%Lallocated = 1 ! this indicates that ii, jj & hi are allocated; 13 Jun 15;

    do ius = 1, 2
     
     select case( ius )
     case( 1 ) ; poincare%idirection = +1 ; dUS = ltangle%dU ; levalue =       ltangle%wr(ius)
     case( 2 ) ; poincare%idirection = -1 ; dUS = ltangle%dS ; levalue = one / ltangle%wr(ius)
     end select
     
     poincare%Ppts = ltangle%ilobe(ius)-1 ; poincare%odetol = ltangle%odetol / ten
     
     do ii = 1, pplobe
      
      poincare%R = ltangle%R + dUS * ltangle%vr(1,ius) * levalue**(ii*one/pplobe-one)
      poincare%Z = ltangle%Z + dUS * ltangle%vr(2,ius) * levalue**(ii*one/pplobe-one)
      
      call pp00aa( poincare, ipp00aa )
      
      lbcoords%sep(1:2,ii-1:ii-1+poincare%Ppts*pplobe:pplobe,ius) = poincare%RZ(1:2,0:poincare%Ppts)
      
      FATAL(fh00aa, ipp00aa.ne.0, integration around separatrix failed need to implement soft fail)
      
     enddo ! end of do ii; 02 Jun 15;

    enddo ! end of do ius; 02 Jun 15;
    
    if( ifh00aa.le.-4 ) then
     do ius = 1, 2
      write(0,'("fh00aa : separatrix(1,*,1) ="999f6.2)') lbcoords%sep(1,0:ns(ius)-1,ius)
      write(0,'("fh00aa : separatrix(2,*,1) ="999f6.2)') lbcoords%sep(2,0:ns(ius)-1,ius)
     enddo
    endif

    ya = ltangle%vr(2,1) / ltangle%vr(1,1)
    yb = ltangle%vr(2,2) / ltangle%vr(1,2)
    
    Nc = 8 ; kk = 0 ! number of constraints; 02 Jun 15;
    
    SALLOCATE(Mc,(1:Nc,1:Nh),zero)
    SALLOCATE(Md,(1:Nc     ),zero)
    
! O point constraints; 02 Jun 15;
    
!   xx =   laxis%R - laxis%R ; yy =   laxis%Z - laxis%Z ; note: xx & yy are zero, exponents will fail; 02 Jun 15;
    
    kk = kk + 1 ; Md(kk) = zero ; Mc(kk,1) = one
    kk = kk + 1 ; Md(kk) = zero ; Mc(kk,2) = one
    kk = kk + 1 ; Md(kk) = zero ; Mc(kk,3) = one
    
! X point constraints; 02 Jun 15;
    
    xx = ltangle%R - laxis%R ; yy = ltangle%Z - laxis%Z ! note: xx & yy are zero, exponents will fail; 02 Jun 15;
    
    kk = kk + 1 ; Md(kk) =  one
    do jj = 1, Nh ; ie = lbcoords%ii(jj) ; je = lbcoords%jj(jj) ; Mc(kk,jj) =      xx**(ie  ) * yy**(je  )
    enddo
    
    kk = kk + 1 ; Md(kk) = zero
    do jj = 1, Nh ; ie = lbcoords%ii(jj) ; je = lbcoords%jj(jj) ; Mc(kk,jj) = ie * xx**(ie-1) * yy**(je  )
    enddo
    
    kk = kk + 1 ; Md(kk) = zero
    do jj = 1, Nh ; ie = lbcoords%ii(jj) ; je = lbcoords%jj(jj) ; Mc(kk,jj) = je * xx**(ie  ) * yy**(je-1)
    enddo
    
    kk = kk + 1 ; Md(kk) = zero
    do jj = 1, Nh ; ie = lbcoords%ii(jj) ; je = lbcoords%jj(jj) ; Mc(kk,jj) = ie * (ie-1) * xx**(ie-2) * yy**(je  ) * one    &
                                                                            + ie *  je    * xx**(ie-1) * yy**(je-1) *  ya    &
                                                                            + je * (je-1) * xx**(ie  ) * yy**(je-2) *  ya**2
    enddo
    kk = kk + 1 ; Md(kk) = zero
    do jj = 1, Nh ; ie = lbcoords%ii(jj) ; je = lbcoords%jj(jj) ; Mc(kk,jj) = ie * (ie-1) * xx**(ie-2) * yy**(je  ) * one    &
                                                                            + ie *  je    * xx**(ie-1) * yy**(je-1) *  yb    &
                                                                            + je * (je-1) * xx**(ie  ) * yy**(je-2) *  yb**2
    enddo

    FATAL(fh00aa, kk.ne.Nc, counting error)

    if( ifh00aa.le.-3 ) then
      write(0,'("fh00aa : axis ("2f15.10") ; tangle ("2f15.10") ;")') laxis%R, laxis%Z, ltangle%R, ltangle%Z
    endif

    NN = Nh + Nc ; MM = NN ; NRA = NN ! NN = #degrees-of-freedom = hamiltonian coefficients + multipliers; 02 Jun 15;
    
    SALLOCATE(mat,(1:NRA,1:NN),zero)
    SALLOCATE(rhs,(1:MM      ),zero)
    
    do ii = 1, Nh
     
     do jj = 1, Nh
      
      do ius = 1, 2
      !do kk =           pplobe, ns(ius)-1
       do kk =                1, ns(ius)-1
        xx = lbcoords%sep(1,kk,ius) - laxis%R
        yy = lbcoords%sep(2,kk,ius) - laxis%Z
        dl = sqrt( (lbcoords%sep(1,kk,ius)-lbcoords%sep(1,kk-1,ius))**2 &
                 + (lbcoords%sep(2,kk,ius)-lbcoords%sep(2,kk-1,ius))**2 )
        mat(ii,jj) = mat(ii,jj) + xx**(lbcoords%ii(ii)+lbcoords%ii(jj)) * yy**(lbcoords%jj(ii)+lbcoords%jj(jj)) * dl
       enddo ! end of do kk ; 02 Jun 15;
      enddo ! end of do ius; 02 Jun 15;
      
     enddo ! end of do jj; 02 Jun 15;

     do kk = 1, Nc ; mat(ii,Nh+kk) = Mc(kk,ii)
     enddo
          
     do jj = 1, 1
      
      do ius = 1, 2
      !do kk =           pplobe, ns(ius)-1
       do kk =                1, ns(ius)-1
        xx = lbcoords%sep(1,kk,ius) - laxis%R
        yy = lbcoords%sep(2,kk,ius) - laxis%Z
        dl = sqrt( (lbcoords%sep(1,kk,ius)-lbcoords%sep(1,kk-1,ius))**2 &
                 + (lbcoords%sep(2,kk,ius)-lbcoords%sep(2,kk-1,ius))**2 )
        rhs(ii   ) = rhs(ii   ) + xx**lbcoords%ii(ii) * yy**lbcoords%jj(ii) * dl
       enddo ! end of do kk; 02 Jun 15;
      enddo ! end of do ius; 02 Jun 15;
      
     enddo ! end of do jj; 02 Jun 15;
     
    enddo ! end of do ii; 02 Jun 15;
    
    do kk = 1, Nc
     ;             ; rhs(Nh+kk   ) = Md(kk   )
     do jj = 1, Nh ; mat(Nh+kk,jj) = Mc(kk,jj)
     enddo
    enddo
    
    if( ifh00aa.le.-9 ) then
     do ii = 1, NN ; write(0,'("fh00aa : matrix="99f9.4)') mat(ii,1:NN)
     enddo
    endif
    
    DALLOCATE(Mc)
    DALLOCATE(Md)
    
#ifdef HAVE_NAG
    Lwk = 4 * NN
    SALLOCATE(wk,(1:Lwk),zero)

    if04jgf = 1
    call F04JGF( MM, NN, mat(1:NRA,1:NN), NRA, rhs(1:MM), lbcoords%svdtol, svd, lbcoords%sigma, irank, &
                 wk(1:Lwk), Lwk, if04jgf ) ! NAG; 13 Oct 15;

    if( ifh00aa.le.-1 ) then
     select case( if04jgf ) ! 02 Jun 15;                                01234567890123
     case( 0 ) ; write(0,1000) if04jgf, svd, lbcoords%sigma, irank, NN, "success ;    "
     case( 1 ) ; write(0,1000) if04jgf, svd, lbcoords%sigma, irank, NN, "input error ;"
     case( 2 ) ; write(0,1000) if04jgf, svd, lbcoords%sigma, irank, NN, "QR failed ;  "
     end select
    endif

1000 format("fh00aa : if04jgf="i2" ; svd="L2" ; sigma="es13.5" ; irank="i6" ; N="i6" ; "a13)
#else
    call DGELS( 'N', MM, NN, 1, mat(1:NRA,1:NN), NRA, rhs(1:MM), MM, lbcoords%sigma, -1, if04jgf ) ! query optimal workspace size
    Lwk = INT(lbcoords%sigma)
    lbcoords%sigma = -1.
    SALLOCATE(wk,(1:Lwk),zero)

    call DGELS( 'N', MM, NN, 1, mat(1:NRA,1:NN), NRA, rhs(1:MM), MM, wk(1:Lwk), Lwk, if04jgf )
    
    if( ifh00aa.le.-1 ) then
      if( if04jgf.le.-1 ) then
        write(0,1000) if04jgf, NN, "input error ;"
      elseif( if04jgf.ge.1 ) then
        write(0,1000) if04jgf, NN, "failed ;     "
      else
        write(0,1000) if04jgf, NN, "success ;    "
      endif
    endif

1000 format("fh00aa : if04jgf="i2" ; N="i6" ; "a13)
#endif

    DALLOCATE(wk)

    lbcoords%hi(1:Nh) = rhs(1:Nh)
    
    DALLOCATE(mat)
    DALLOCATE(rhs)

    if( ifh00aa.le.-3 ) then
      write(0,'("fh00aa : h["i2","i2"]="es13.5" ;")') ( lbcoords%ii(ii), lbcoords%jj(ii), lbcoords%hi(ii), ii = 1, Nh )
    endif

    if( if04jgf.ne.0 ) then ; ifh00aa = if04jgf ; goto 9999
    endif

    ifh00aa = 0
    
9999 continue

    lbcoords%nbfield(0:1) = nbfield(0:1)
    
   !if( ifail.le. 0 ) then 
   !  write(0,'("fh00aa : ifail ="i3" ; Nx="i3" ; Ny="i3" ; Nh="i6" ; Nconstraints="i6" ; svdtol="es13.5" ; sigma="es13.5" ;")') &
   !    ifh00aa, lbcoords%Nx, lbcoords%Ny, lbcoords%Nh, Nc, lbcoords%svdtol, lbcoords%sigma
   !endif
    
    ifail = ifh00aa
    
    return
    
  end subroutine fh00aa

#endif

!latex \subroutine{ad00aa}{anisotropic diffusion using locally-field-aligned coordinates;}
!latex \bi
!latex \item[1.] This subroutine integrates in time the anisotropic-diffusion equation,
!latex           \be \frac{\partial p}{\partial t} = \nabla \cdot \left( \kappa_\parallel \nabla_{\parallel} p + \kappa_\perp \nabla_{\perp} p \right) + S,
!latex           \label{eq:anisotropicdiffusion}
!latex           \ee
!latex           where $p$ is a scalar function of position, e.g. the pressure; $S$ is a scalar function of position, e.g. a source;
!latex            $t$ is an arbitrary integration parameter, e.g. time;
!latex           the directional derivative along the magnetic field is $\nabla_\parallel p  \equiv {\bf b} \, {\bf b}\cdot\nabla p$;
!latex           and $\nabla_{\perp} p \equiv \nabla p - \nabla_{\parallel} p$.
!latex \item[* ] This equation may be re-written as
!latex           $\partial_t p = \nabla \cdot [ (\kappa_\parallel-\kappa_\perp) \nabla_\parallel p + \kappa_\perp \nabla p) + S$,
!latex           and hereafter will use $\bar \kappa_\parallel \equiv \kappa_\parallel - \kappa_\perp$.
!latex \item[* ] \red{[Presently, $\bar \kappa = \kappa_\parallel$, i.e. the $\kappa_\perp \nabla_\parallel p$ term is ignored.]}
!latex \item[* ] It is more transparent to write the anisotropic-diffusion equation as
!latex           \be \frac{\partial p}{\partial t} = 
!latex           \;\; \bar \kappa_\parallel \;\; {\bf B} \cdot \nabla \;\; \frac{1}{B^2} \;\; {\bf B} \cdot \nabla \;\; p \;\; + \;\;
!latex               \kappa_\perp \;\; \nabla\cdot\nabla p \;\; + \;\; S.
!latex           \ee
!latex \item[* ] The operator ${\bf B}\cdot \nabla$ reduces to $B^\phi \partial_\phi$ if coordinates, $(\alpha,\beta,\phi)$, can be constructed
!latex           so that ${\bf B} = \nabla \alpha \times \nabla \beta$ and $\nabla \alpha \times \nabla \beta \cdot \nabla \phi = B^\phi$.
!latex           Such coordinates can always be constructed locally by fieldline integration, {\em provided} $B^\phi \ne 0$. 
!latex           The benefit of this approach is to minimize the ``parallel-pollution'' of the perpendicular diffusion.
!latex \item[* ] A regular, cylindrical grid is employed, $R \in[R_{min},R_{max}]$, $\phi \in[0,2\pi/N_P]$ and $Z\in[Z_{min}, Z_{max}]$,
!latex           where $N_P\equiv$ \verb+Nfp+ is the field periodicity;  with grid resolution $N_R$, $N_\phi$ and $N_Z$.
!latex           The location of each grid point, ${\bf x}_{i,j,k} = R_i {\bf e}_R + Z_j {\bf e}_Z$ is given by 
!latex           $R_i \equiv R_{min} + i \Delta R$, $Z_j \equiv Z_{min} + j \Delta Z$ and $\phi_k = k \Delta \phi$; 
!latex           where $\Delta R = (R_{max}-R_{min})/N_R$, $\Delta Z = (Z_{max}-Z_{min})/N_Z$ and $\Delta \phi = (2\pi/N_P)/N_\phi$.
!latex \item[* ] A second-order discretization of the {\em second} parallel derivative is obtained by 
!latex           differencing the {\em first} parallel derivatives on the ``forward half-$\Delta\phi$" and the ``backward half-$\Delta\phi$" grids as follows:
!latex           \be \nabla^2_{\parallel} p|_{i,j,k} & \equiv & \left. \frac{B^\phi}{\Delta \phi} \right|_{0} \left( 
!latex           \left. \frac{B^\phi}{B^2} \right|_{+1/2} \frac{p_{+1}-p_{ 0}}{\Delta \phi} 
!latex         - \left. \frac{B^\phi}{B^2} \right|_{-1/2} \frac{p_{ 0}-p_{-1}}{\Delta \phi} 
!latex           \right), \label{eq:secondparallelderivative}
!latex            \ee
!latex           where $p_{0} \equiv p({\bf x}_{i,j,k})=p_{i,j,k}$ and $p_{\pm 1} \equiv p( {\cal M}^{\pm 1}( {\bf x}_{i,j,k})) $, 
!latex           where ${\cal M}^{1}({\bf x}_{i,j,k})$ is the image of ${\bf x}_{i,j,k}$ under the ``forward" fieldline mapping 
!latex           from $\phi=k \Delta \phi$ to $\phi=(k+1) \Delta \phi$, 
!latex           and ${\cal M}^{-1}({\bf x}_{i,j,k})$ is the image of ${\bf x}_{i,j,k}$ under the ``backward" fieldline mapping 
!latex           from $\phi=k \Delta \phi$ to $\phi=(k-1) \Delta \phi$; and the factors $B^\phi/B^2$ 
!latex           are evaluated on the forward and backward half-$\Delta\phi$ grids.
!latex \item[* ] \red{[Presently, $B^\phi = 1$ and $B^\phi/B^2=1$, i.e. the ``metric'' information is ignored.]}
!latex \item[* ] The ${\cal M}^{\pm 1}( {\bf x}_{i,j,k})$ do not generally coincide with grid points, 
!latex           so $ p( {\cal M}^{\pm 1}( {\bf x}_{i,j,k}))$ must be determined by interpolation.
!latex           Given that the ${\cal M}^{\pm 1}( {\bf x}_{i,j,k})$ do however lie on the toroidal planes $(k+1)\Delta \phi$ and $(k-1)\Delta \phi$, 
!latex           an interpolation in $\phi$ is not required.
!latex           The simplest interpolation is the second-order, bi-linear interpolation, e.g.
!latex           \be p_{+1} \equiv ( 1 - y ) \; [ ( 1 - x ) \; p_{I,J,k+1} + x \; p_{I+1,J,k+1} ] + y \; [ ( 1 - x ) \; p_{I,J+1,k+1} + x \; p_{I+1,J+1,k+1} ],
!latex           \ee
!latex           where the $I$ and $J$ label the left-lower grid point, i.e. $R_{I} \le \bar R < R_{I+1}$ 
!latex                                                                   and $Z_{J} \le \bar Z < Z_{J+1}$;
!latex           and $x\equiv (\bar R - R_{I})/\Delta R$ and $y \equiv (\bar Z - Z_{J})/\Delta Z$, 
!latex \item[* ] The bi-linear interpolation is used near the computational boundary, 
!latex           and a fourth-order, bi-cubic interpolation is used in the interior domain.
!latex \item[* ] The fieldline tracing information is contained in the 
!latex           $I_{\pm 1,i,j,k}$, $J_{\pm 1,i,j,k}$, $x_{\pm 1,i,j,k}$ and the $y_{\pm 1,i,j,k}$ arrays;
!latex           and the $B^\phi|_{0}$ and $B^\phi/B^2|_{\pm 1/2}$ information is contained in $B_{-1:1,,i,j,k}$.
!latex \item[* ] The discretization of the second, parallel derivative given in \Eqn{secondparallelderivative}
!latex           is only possible if $B^\phi\ne 0$; as the equations defining the 
!latex           fieldline mapping are $d_\phi R \equiv B^R/B^\phi$ and $d_\phi Z \equiv B^Z/B^\phi$, 
!latex           and the Jacobian of the locally-field-aligned coordinates is $1/B^\phi$.
!latex           An additional logical array, $F_{i,j,k}$ indicates whether the fieldline tracing construction of the locally-field-aligned coordinates
!latex           was successful.
!latex           For all points where the fieldline tracing was {\em not} successful the pressure at the ``mapped'' points is set to zero, i.e. $p_{\pm}=0$.
!latex          %This effectively sets the boundary conditions that must be supplied to the anisotropic-diffusion equation, \Eqn{anisotropicdiffusion}.
!latex \item[* ] All of these arrays, i.e. the $I$, $J$, $x$, $y$, $B$ and $F$, are returned by \verb+ad00aa+.
!latex           Note that these arrays depend only on the magnetic field and the computational grid.
!latex           If the user wishes to continue to relax the pressure {\em with the same magnetic field} and {\em with the same computational grid},
!latex           then this information can be passed back to \verb+ad00aa+ on a subsequent call, and this will eliminate the computational cost of
!latex           re-constructing the locally-field-aligned coordinates.
!latex           (Some computational savings can also be made if the user wishes to refine the grid in only the $R$ and $Z$ directions
!latex           for the same magnetic field, but this option is not yet implemented.)
!latex \item[* ] If, however, on a subsequent call to \verb+ad00aa+, the magnetic field has changed, the $I$, $J$, $x$, $y$, $B$ and $F$ arrays will 
!latex           not be consistent with the new magnetic field.
!latex           The input flag \verb+Lcontrol+ determines how this information is to be used.
!latex \item[* ] The ``non-directional diffusion'' term, $\nabla \cdot \nabla p$, 
!latex           is given in cylindrical coordinates as
!latex           \be \nabla \cdot \nabla p = \frac{1}{R} \left[ 
!latex             \frac{\partial}{\partial R} \left(R \, \frac{\partial p}{\partial R}\right) 
!latex           + \frac{\partial}{\partial \phi} \left(\frac{1}{R}\frac{\partial p}{\partial R}\right)
!latex           + \frac{\partial}{\partial Z} \left(R \, \frac{\partial p}{\partial Z}\right) \right]. \label{eq:nondirectionaldiffusiona}
!latex           \ee
!latex \item[* ] \red{This 
!latex           is simplified to 
!latex           \be \nabla \cdot \nabla p =
!latex             \frac{\partial^2 p}{\partial R^2}
!latex           + \frac{\partial^2 p}{\partial Z^2}, \label{eq:nondirectionaldiffusionb}
!latex           \ee
!latex           i.e., metric information is ignored.}
!latex \item[* ] A fourth-order discretization of \Eqn{nondirectionaldiffusionb} is employed.
!latex \item[* ] After constructing the parallel and perpendicular derivatives, the pressure is advanced explicitly via
!latex           \be p^{n+1}_{i,j,k} = p^{n}_{i,j,k} + \Delta t \left( \bar \kappa_\parallel \nabla^2_\parallel p|^n_{i,j,k} 
!latex            + \kappa_\perp \nabla\cdot\nabla p |^n_{i,j,k} + S_{i,j,k} \right) \label{eq:pressurerelax}
!latex           \ee
!latex           for $n=1,N$ time-steps. %, where the Jacobian $\sqrt g_{i,j,k} \equiv R_i$ for cylindrical coordinates.
!latex \item[* ] It is possible, upon request, to implement an implicit time advance.
!latex \item[* ] The boundary condition is that $p=0$ on the computational boundary defined by 
!latex           $R_{min}$, $R_{max}$, and $Z_{min}$, $Z_{max}$, and this boundary condition is enforced internally.
!latex          %However, given that the condition $p_{i,j,k}=0$ is enforced if an error was encountered during the integration along the magnetic field
!latex          %from the starting point ${\bf x}_{i,j,k}$ using the toroidal angle, $\phi$, as the integration parameter,
!latex          %the boundary condition is effectively $p=0$ outside the ``physical'' domain, 
!latex          %which is defined as the set of all points, ${\bf x}_{i,j,k}$, for which the fieldline starting at ${\bf x}_{i,j,k}$ 
!latex          %reaches the adjacent toroidal planes in the forwards and backwards direction along the magnetic field without encountering 
!latex          %either $|B^\phi|<\epsilon$ or $|B|<\epsilon$, and without leaving 
!latex          %the domain for which the user-provided field, \verb+bfield+, does not return an error.
!latex \item[* ] The iterations will terminate if $p$ becomes either too small or too large, 
!latex           which presumably indicate that the CFL condition on the explicit integration timestep has been violated.
!latex           It is assumed that $p$ will be order unity.
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : pressurerelax, ad00aa+ \\ \\
!latex           \verb+type(pressurerelax) :: pressure+ \\ \\ in their source that calls \verb+ad00aa+,
!latex           where \verb+pressurerelax+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+pressure+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+pressure%Nfp            : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, e.g. \verb+Nfp=1+;
!latex \ei
!latex \item[  ] \verb+pressure%NR             : integer ;+
!latex \item[  ] \verb+pressure%Np             : integer ;+
!latex \item[  ] \verb+pressure%NZ             : integer ;+
!latex \bi
!latex \item[i.] grid resolution in $R$;    constraint \verb+NR+ $\ge 4$;
!latex \item[ii.] grid resolution in $\phi$; constraint \verb+Np+ $\ge 1$;
!latex \item[iii.] grid resolution in $Z$;    constraint \verb+NZ+ $\ge 4$;
!latex \ei
!latex \item[  ] \verb+pressure%Rmin           : real    ;+
!latex \item[  ] \verb+pressure%Rmax           : real    ;+
!latex \item[  ] \verb+pressure%Zmin           : real    ;+
!latex \item[  ] \verb+pressure%Zmax           : real    ;+
!latex \bi
!latex \item[i.] calculation domain in $R$; constraint \verb+Rmax+$>$\verb+Rmin+;   \verb+Rmin+$>\epsilon$;
!latex \item[ii.] calculation domain in $Z$; constraint \verb+Zmax+$>$\verb+Zmin+;
!latex \ei
!latex \item[  ] \verb+pressure%Lcontrol       : integer ;+
!latex \bi
!latex \item[i.] if \verb+Lcontrol = 1,+ : only construct the locally-field-aligned coordinates by fieldline tracing; \\
!latex           the information describing the locally-field-aligned coordinates are returned in arrays 
!latex           \verb+x+, \verb+y+, \verb+I+, \verb+J+, \verb+B+ and \verb+F+; \\
!latex           if these arrays are allocated on input, they are first de-allocated; \\
!latex           note that the pressure and source are not referenced, and so they need not be allocated;
!latex \item[ii.] if \verb+Lcontrol = 2,+ : relax pressure by integrating \Eqn{pressurerelax}, 
!latex            assuming that the locally-field-aligned coordinate information is provided; \\
!latex            the \verb+x+, \verb+y+, \verb+I+, \verb+J+, \verb+B+ and \verb+F+ arrays returned by a previous call to \verb+ad00aa+
!latex            \underline{must} be provided;\\
!latex            note that these arrays depend on the magnetic field and the computational grid, so this option should only be used if these
!latex            have not changed since the last call to \verb+ad00aa+; \\
!latex            the pressure and the source must be allocated on input;
!latex \item[iii.] if \verb+Lcontrol = 3,+ : first construct the locally-field-aligned coordinates as for \verb+Lcontrol=1+, 
!latex             then relax the pressure as for \verb+Lcontrol=2+; \\
!latex             only the pressure and source must be allocated on input;
!latex \ei
!latex \item[  ] \verb+pressure%Lode           : integer ;+
!latex \bi
!latex \item[i.]  \verb+Lode = 0+ : a single $4$-th order, Runge-Kutta o.d.e. integration step will be perfomed 
!latex                              to construct the field-aligned coordinates between adjacent toroidal planes; \\
!latex                              this will invariably be faster than using option \verb+Lode=1+ (i.e. require fewer evaluations of the magnetic field) 
!latex                              but this option will only be reliable if $\Delta \phi$ is sufficiently small, i.e. \verb+Np+ is sufficiently large;
!latex \item[ii.] \verb+Lode = 1+ : the NAG routine \verb+D02BFJ+ will be used to more-accurately integrate along the magnetic field 
!latex                              to construct the field-aligned coordinates;
!latex \ei
!latex \item[  ] \verb+pressure%odetol         : real    ;+
!latex \bi
!latex \item[i.] o.d.e. integration tolerance provided to \verb+D02BJF+ for integrating along the magnetic field;
!latex \item[ii.] only used if \verb+Lode = 1+;
!latex \item[iii.] the accuracy of the calculation is limited by $N_R$ and $N_Z$, so \verb+odetol+ need not be too small; suggested \verb+odetol+$=10^{-5}$.
!latex \ei
!latex \item[  ] \verb+pressure%p(0:NR,0:NZ,-1:Np) : real ;+
!latex \bi
!latex \item[i.] initial state for the pressure, where \verb+p(i,j,k)+$\equiv p(R_{min}+i \Delta R, k \Delta \phi, Z_{min}+j \Delta Z )$;
!latex \item[ii.] this must be allocated {\em before} calling \verb+ad00aa+ if \verb+Lcontrol=2,3+;
!latex \item[iii.] a suitable initization is \be p=(1-x^2)(1-y^2),\ee
!latex            where $x\equiv 2(R-R_{mid})/(R_{max}-R_{min})$ where $R_{mid}\equiv (R_{min}+R_{max})/2$; and similarly for $y$;
!latex            or, more simply, $p=1$.
!latex \item[iv.] on exit, this array will be updated;
!latex \ei
!latex \item[  ] \verb+pressure%s(0:NR,0:NZ, 0:Np-1) : real ;+
!latex \bi
!latex \item[i.] source, where \verb+s(i,j,k)+$\equiv s(R_{min}+i \Delta R, k \Delta \phi, Z_{min}+j \Delta Z )$;
!latex \item[ii.] this must be allocated {\em before} calling \verb+ad00aa+ if \verb+Lcontrol=2,3+;
!l tex \item[iii.] it is strongly recommended that the source be zero outside the physical domain, as defined above;
!latex \ei
!latex \item[  ] \verb+pressure%Ntime          : integer ;+
!latex \bi
!latex \item[i.] number of time steps taken; only if \verb+Lcontrol=2,3+ ; e.g. \verb+Ntime=10000+;
!latex \ei
!latex \item[  ] \verb+pressure%dtime          : real    ;+
!latex \bi
!latex \item[i.] \verb+dtime+$\equiv \Delta t$ appearing in \Eqn{pressurerelax}; only if \verb+Lcontrol=2,3+;
!latex \ei
!latex \item[  ] \verb+pressure%kpara          : real    ;+
!latex \item[  ] \verb+pressure%kperp          : real    ;+
!latex \bi
!latex \item[i.] the parallel and perpendicular diffusion coefficients;
!latex \ei
!l tex \item[  ] \verb+pressure%Ldiff          : integer ;+
!l tex \bi
!l tex \item[i.] if \verb+Ldiff=0+ the diffusion operator, \Eqn{nondirectionaldiffusionb}, is approximated by writing 
!l tex           \be \frac{\partial}{\partial R} \left(R p_R\right) \equiv \frac{[Rp_R]_{+1/2}-[Rp_R]_{-1/2}}{\Delta R} = 
!l tex           \frac{1}{\Delta R}\left( \frac{R_{i+1/2}(p_{i+1,j,k}-p_{i,j,k})}{\Delta R} - \frac{R_{i-1/2}(p_{i,j,k}-p_{i-1,j,k})}{\Delta R} \right),
!l tex           \ee
!l tex           and similarly for the other terms, to obtain
!l tex           \be (\nabla\cdot\nabla p)_{i,j,k} 
!l tex           & = & \frac{1}{\Delta R   } \left( R_{i+1/2}     \frac{ p_{i+1,j  ,k  }-p_{i  ,j  ,k  } }{\Delta R   } 
!l tex                                            - R_{i-1/2}     \frac{ p_{i  ,j  ,k  }-p_{i-1,j  ,k  } }{\Delta R   } \right) \frac{1}{R_i} \nonumber \\
!l tex           & + & \frac{1}{\Delta \phi} \left( \frac{1}{R_i} \frac{ p_{i  ,j  ,k+1}-p_{i  ,j  ,k  } }{\Delta \phi} 
!l tex                                            - \frac{1}{R_i} \frac{ p_{i  ,j  ,k  }-p_{i  ,j  ,k-1} }{\Delta \phi} \right) \frac{1}{R_i} \nonumber \\
!l tex           & + & \frac{1}{\Delta Z   } \left( R_i           \frac{ p_{i  ,j+1,k  }-p_{i  ,j  ,k  } }{\Delta Z   } 
!l tex                                            - R_i           \frac{ p_{i  ,j  ,k  }-p_{i  ,j-1,k  } }{\Delta Z   } \right) \frac{1}{R_i} 
!l tex           \ee
!l tex \item[i.] if \verb+Ldiff=1+ a symmetrized, second-order operator is used;
!l tex \ei
!latex \item[  ] \verb+ifail                   : integer ;+
!latex \bi
!latex \item[i.] if \verb+ifail =  1+, quiet mode;
!latex \item[i.] if \verb+ifail =  0+, screen output is terse;
!latex \item[i.] if \verb+ifail = -1+, 
!latex \item[i.] if \verb+ifail = -2+, some information regarding the relaxation will be shown at every 10000 timesteps;
!latex \item[i.] if \verb+ifail = -3+, some information regarding the relaxation will be shown at every 1000  timesteps;
!latex \item[i.] if \verb+ifail = -4+, some information regarding the relaxation will be shown at every 100   timesteps;
!latex \item[i.] if \verb+ifail = -5+, some information regarding the relaxation will be shown at every 10    timesteps;
!latex \item[i.] if \verb+ifail = -6+, some information regarding the relaxation will be shown at every 1     timesteps;
!latex \ei
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call ad00aa( pressure, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+pressure%x( 0:1,1:NR,1:NZ,0:Np-1) : real    ;+
!latex \item[  ] \verb+pressure%y( 0:1,1:NR,1:NZ,0:Np-1) : real    ;+
!latex \item[  ] \verb+pressure%I( 0:1,1:NR,1:NZ,0:Np-1) : integer ;+
!latex \item[  ] \verb+pressure%J( 0:1,1:NR,1:NZ,0:Np-1) : integer ;+
!latex \item[  ] \verb+pressure%B(-1:1,1:NR,1:NZ,0:Np-1) : real    ;+
!latex \item[  ] \verb+pressure%F(     1:NR,1:NZ,0:Np-1) : logical ;+
!latex \bi
!latex \item[i.] information describing the locally-field-aligned coordinates;
!latex \item[ii.] if on a subsequent call to \verb+ad00aa+ the user chooses to set \verb+Lcontrol+$=2$, 
!latex            then these arrays are required on input to the subsequent call of \verb+ad00aa+;
!latex \ei
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex     \item[] \verb+ifail=2+ : the pressure has exploded, probably because the CFL condition on the timestep has been violated;
!latex \ei
!latex \item[6.] Comments:
!latex \item[* ] The steady state solution has ${\cal O} (p/S) \sim {\cal O}(\bar \kappa_\parallel / \kappa_\perp)$, where $S$ is the source; 
!latex           so that if $p\sim {\cal O}(\kappa_\parallel)$ then $S \sim {\cal O}(\kappa_\perp)$.
!l tex \item[* ] For formal convergence of the solution with respect to grid resolution it is required that that source is zero outside the physical domain.
!l tex           This is because, as the numerical resolution increases, the physical boundary is more accurately resolved; 
!l tex           and as the physical boundary determines the boundary condition that supplements \Eqn{anisotropicdiffusion}, 
!l tex           and as convergence of the solution of a differential equation with respect to grid resolution can only be achieved if the boundary condition 
!l tex           is held constant, it is required that the effective source does not change as the grid resolution is increased.
!latex \item[* ] This subroutine is not yet parallelized, but it is easy to do so.
!latex \ei
  
  subroutine ad00aa( pp, ifail )
    
    implicit none
    
    type(pressurerelax)    :: pp
    
    integer, parameter     :: Node = 2, Lwk = 20 * Node
    integer                :: astat, ifail, ii, jj, kk, id02bjf, Nfp, NR, Np, NZ, ifb, lfb, itime, I, J, ibfield, ix, iy
    real                   :: RZ(1:Node), oRZ(1:2), zstart, zend, phi, phih, phif
    real                   :: k1(1:Node), k2(1:Node), k3(1:Node), k4(1:Node)
    real                   :: RpZ(1:3), dBRpZ(1:3,0:3), BB, Bphi
    real                   :: Rmin, Rmax, Zmin, Zmax, DR, HR, Dp, Hp, DZ, DRDR, DpDp, DZDZ, Rh(-4:4), R0
    real                   :: dRprdR, dppRdp, dRpzdZ, g0, g1, interpress(-1:1), ff(0:3,0:3), xx(0:3), yy(0:3)
    real                   :: dpara, dperp, po, pn
    real                   :: dpdRR, dpdZZ, stencil(1:16), dpdx(-2:1,-2:1), dpdy(-2:1,-2:1), sauce
    character              :: relabs
    
#ifdef HAVE_NAG
    REAL                   :: wk(1:Lwk)
    external               :: D02BJW
#else
    INTEGER, parameter     :: liw=20, lrw=20+16*Node
    INTEGER                :: iwork(liw)
    REAL                   :: rwork(lrw)
    INTEGER                :: iopt,istate,itask,itol,mf,ij
    REAL                   :: atol,rtol,ddz,zc,ze
#endif
    
    iad00aa = ifail ; po = zero ; pn = zero
    
    if( pp%Nfp .le.0                                   ) then
      write(0,'("ad00aa : input error ;     Nfp.le.0     ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( pp%NR  .lt.4                                   ) then
      write(0,'("ad00aa : input error ;      NR.lt.4     ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( pp%Np  .lt.1                                   ) then
      write(0,'("ad00aa : input error ;      Np.lt.1     ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( pp%NZ  .lt.4                                   ) then
      write(0,'("ad00aa : input error ;      NZ.lt.4     ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( pp%Rmin.le.small                               ) then
      write(0,'("ad00aa : input error ;    Rmin.le.small ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( pp%Rmax.le.pp%Rmin                             ) then
      write(0,'("ad00aa : input error ;    Rmax.le.Rmin  ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( pp%Zmax.le.pp%Zmin                             ) then
      write(0,'("ad00aa : input error ;    Zmax.le.Zmin  ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( pp%Lode.eq.1       .and. pp%odetol  .le. small ) then
      write(0,'("ad00aa : input error ;  odetol.le.small ;")') ; iad00aa = 1 ; goto 9999
    endif
    
    Nfp = pp%Nfp
    
    Rmin = pp%Rmin ; Rmax = pp%Rmax ; NR = pp%NR ; DR  = ( Rmax - Rmin ) / NR ; DRDR = DR**2 ; HR = DR * half
    ;              ;                ; Np = pp%Np ; Dp  = (  pi2 / Nfp  ) / Np ; DpDp = Dp**2 ; Hp = Dp * half
    Zmin = pp%Zmin ; Zmax = pp%Zmax ; NZ = pp%NZ ; DZ  = ( Zmax - Zmin ) / NZ ; DZDZ = DZ**2
    pressure%Dp = Dp
    
    if( pp%Lcontrol.eq.1 .or. pp%Lcontrol.eq.3 ) then ! construct locally-field-aligned coordinates; 29 Jul 15;
     
     if( allocated(pp%x) ) deallocate(pp%x)
     if( allocated(pp%y) ) deallocate(pp%y)
     if( allocated(pp%I) ) deallocate(pp%I)
     if( allocated(pp%J) ) deallocate(pp%J)
     if( allocated(pp%B) ) deallocate(pp%B)
     if( allocated(pp%F) ) deallocate(pp%F)
     
      ! macro expansion of sallocate = set allocate;

#ifdef DEBUG
   if( allocated( pp%x ) ) then
    write(0,'("oculus : 0123456789 : pp%x already allocated ;")') 
    stop      'oculus : 0123456789 : pp%x already allocated ;'
   endif
#endif

   allocate( pp%x( 0:1,1:NR-1,1:NZ-1,0:Np-1), stat=astat )

#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("oculus : 0123456789 : error allocating pp%x ;")') 
    stop      'oculus : 0123456789 : error allocating pp%x ;'
   endif
#endif

   pp%x( 0:1,1:NR-1,1:NZ-1,0:Np-1) = zero  

 ! bi-linear ; weights; 29 Jul 15;
      ! macro expansion of sallocate = set allocate;

#ifdef DEBUG
   if( allocated( pp%y ) ) then
    write(0,'("oculus : 0123456789 : pp%y already allocated ;")') 
    stop      'oculus : 0123456789 : pp%y already allocated ;'
   endif
#endif

   allocate( pp%y( 0:1,1:NR-1,1:NZ-1,0:Np-1), stat=astat )

#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("oculus : 0123456789 : error allocating pp%y ;")') 
    stop      'oculus : 0123456789 : error allocating pp%y ;'
   endif
#endif

   pp%y( 0:1,1:NR-1,1:NZ-1,0:Np-1) = zero  

 ! bi-linear ; weights; 29 Jul 15;
      ! macro expansion of sallocate = set allocate;

#ifdef DEBUG
   if( allocated( pp%I ) ) then
    write(0,'("oculus : 0123456789 : pp%I already allocated ;")') 
    stop      'oculus : 0123456789 : pp%I already allocated ;'
   endif
#endif

   allocate( pp%I( 0:1,1:NR-1,1:NZ-1,0:Np-1), stat=astat )

#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("oculus : 0123456789 : error allocating pp%I ;")') 
    stop      'oculus : 0123456789 : error allocating pp%I ;'
   endif
#endif

   pp%I( 0:1,1:NR-1,1:NZ-1,0:Np-1) = 0  

 ! bi-linear ; cell   ; 29 Jul 15;
      ! macro expansion of sallocate = set allocate;

#ifdef DEBUG
   if( allocated( pp%J ) ) then
    write(0,'("oculus : 0123456789 : pp%J already allocated ;")') 
    stop      'oculus : 0123456789 : pp%J already allocated ;'
   endif
#endif

   allocate( pp%J( 0:1,1:NR-1,1:NZ-1,0:Np-1), stat=astat )

#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("oculus : 0123456789 : error allocating pp%J ;")') 
    stop      'oculus : 0123456789 : error allocating pp%J ;'
   endif
#endif

   pp%J( 0:1,1:NR-1,1:NZ-1,0:Np-1) = 0  

 ! bi-linear ; cell   ; 29 Jul 15;
      ! macro expansion of sallocate = set allocate;

#ifdef DEBUG
   if( allocated( pp%B ) ) then
    write(0,'("oculus : 0123456789 : pp%B already allocated ;")') 
    stop      'oculus : 0123456789 : pp%B already allocated ;'
   endif
#endif

   allocate( pp%B(-1:1,1:NR-1,1:NZ-1,0:Np-1), stat=astat )

#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("oculus : 0123456789 : error allocating pp%B ;")') 
    stop      'oculus : 0123456789 : error allocating pp%B ;'
   endif
#endif

   pp%B(-1:1,1:NR-1,1:NZ-1,0:Np-1) = zero  

 ! parallel weights   ; 29 Jul 15;
      ! macro expansion of sallocate = set allocate;

#ifdef DEBUG
   if( allocated( pp%F ) ) then
    write(0,'("oculus : 0123456789 : pp%F already allocated ;")') 
    stop      'oculus : 0123456789 : pp%F already allocated ;'
   endif
#endif

   allocate( pp%F(     1:NR-1,1:NZ-1,0:Np-1), stat=astat )

#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("oculus : 0123456789 : error allocating pp%F ;")') 
    stop      'oculus : 0123456789 : error allocating pp%F ;'
   endif
#endif

   pp%F(     1:NR-1,1:NZ-1,0:Np-1) = .false.

 ! logical            ; 29 Jul 15;
     
     itangent = 0
     
     if( ifail.le.-1 ) write(0,1000) NR, Nz, Np, pp%Lode
1000 format("ad00aa :         "2x" : constructing field-alligned coordinates ; NR="i6" ; NZ="i6" ; Np="i6" ; Lode="i2" ;")
     
     do kk = 0, Np-1
      
      do ii = 1, NR-1
       do jj = 1, NZ-1
        
        Lbfieldok = .true. ; pressure%lB(-1:1) = zero
        
        phi = kk * Dp ; oRZ(1:2) = (/ Rmin, Zmin /) + (/ ii, jj /) * (/ DR, DZ /)
        
        RpZ(1:3) = (/ oRZ(1), phi, oRZ(2) /) ; ibfield = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), ibfield )
        
        Bphi = dBRpZ(2,0) ; BB = sqrt( dBRpZ(1,0)**2 + (RpZ(1)*dBRpZ(2,0))**2 + dBRpZ(3,0)**2 )
        
        pressure%lB(0) = abs(Bphi) ! shall record B^phi at every location for plotting; 29 Jul 15;

        if( ibfield.ne.0 .or. abs(Bphi).lt.small .or. BB.lt.small ) Lbfieldok = .false.
        
        pp%F(ii,jj,kk) = .not.Lbfieldok
        
        if( pp%F(ii,jj,kk) ) cycle ! the fieldline integration has failed; 29 Jul 15;        
        
        if( pp%Lode.eq.0 ) k1(1:2) = (/ dBRpZ(1,0), dBRpZ(3,0) /) / Bphi
        
        do ifb = -1, 1, 2 ; pressure%ifb = ifb ; pressure%ii = 0 ; lfb = max(0,ifb) ; nbfield(0:1) = 0
         
         select case( pp%Lode )
          
         case( 0 ) ! take one Runge-Kutta step to adjacent toroidal planes;  8 Jul 15;

#ifdef HAVE_NAG          
!        !phi  = kk * Dp            ; RZ(1:2) = oRZ(1:2)
!        !call bf00ac( phi , RZ(1:2), k1(1:2) )
          phih = kk * Dp + ifb * Hp ; RZ(1:2) = oRZ(1:2) + ifb * Hp * k1(1:2)
          call bf00ac( phih, RZ(1:2), k2(1:2) )
          ;                         ; RZ(1:2) = oRZ(1:2) + ifb * Hp * k2(1:2)
          call bf00ac( phih, RZ(1:2), k3(1:2) )
#else
!        !phi  = kk * Dp            ; RZ(1:2) = oRZ(1:2)
!        !call bf00ac( Node, phi , RZ(1:2), k1(1:2) )
          phih = kk * Dp + ifb * Hp ; RZ(1:2) = oRZ(1:2) + ifb * Hp * k1(1:2)
          call bf00ac( Node, phih, RZ(1:2), k2(1:2) )
          ;                         ; RZ(1:2) = oRZ(1:2) + ifb * Hp * k2(1:2)
          call bf00ac( Node, phih, RZ(1:2), k3(1:2) )
#endif
          
          RpZ(1:3) = (/ RZ(1), phih, RZ(2) /) ; ibfield = -9
          call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), ibfield )
          
          Bphi = dBRpZ(2,0) ; BB = sqrt( dBRpZ(1,0)**2 + (RpZ(1)*dBRpZ(2,0))**2 + dBRpZ(3,0)**2 )
          
          if( ibfield.ne.0 .or. abs(Bphi).lt.small .or. BB.lt.small ) Lbfieldok = .false.
          
          if( Lbfieldok ) pressure%lB(ifb) = abs(Bphi) / BB**2

#ifdef HAVE_NAG          
          phif = kk * Dp + ifb * Dp  ; RZ(1:2) = oRZ(1:2) + ifb * Dp  * k3(1:2)
          call bf00ac( phif, RZ(1:2), k4(1:2) )
#else
          phif = kk * Dp + ifb * Dp  ; RZ(1:2) = oRZ(1:2) + ifb * Dp  * k3(1:2)
          call bf00ac( Node, phif, RZ(1:2), k4(1:2) )
#endif
          
          RZ(1:2) = oRZ(1:2) + ifb * Dp * ( k1(1:2) + two * ( k2(1:2) + k3(1:2) ) + k4(1:2) ) * sixth         
          
         case( 1 ) ! "exact" o.d.e. integration to adjacent toroidal planes;  8 Jul 15;
          
          RZ(1:2) = oRZ(1:2) ; relabs = 'A' ; zstart = phi ; zend = zstart + ifb * Dp
          
#ifdef HAVE_NAG
          id02bjf = 1 ! NAG; 13 Oct 15;
          call D02BJF( zstart, zend, Node, RZ(1:Node), bf00ac, pp%odetol, relabs, ad00ab, D02BJW, wk(1:Lwk), id02bjf )
#else
!         set the integrator parameters.
          rwork=0.
          iwork=0
!         istate=1 :  indicates the first lsode call
          istate=1
!         itask=4 :  normal integration with limited over-shoot (set by
!                    rwork(1) in the i_lsode loop
          itask=1
!         iopt=1 :  optional inputs are used
          iopt=1
!         rwork(6) :  set maximum lsode-internal step size
          iwork(6)=400000
!         rwork(7) :  set minimum lsode-internal step size
!         iwork(6) :  set maximum lsode-internal steps
!         iwork(7) :  set maximum lsode-internal error messages printed
!         mf=10 :  non-stiff Adams method of integration
          mf=10
!         itol=1 :  indicates absolute tolerance is just a scalar
          itol=1
!         rtol :  relative tolerance
          rtol=pp%odetol
!         atol :  absolute tolerance
          atol=pp%odetol
!         initializations for loop
          zc=zstart
          ze=zc
          call ad00ab(ze,RZ)
          do ij=1,2
            call lsode(bf00ac,(/Node/),RZ,zc,ze, &
                       itol,(/rtol/),(/atol/),itask,istate,iopt, &
                       rwork,lrw,iwork,liw,du00aa,mf)
            zc=ze
            call ad00ab(ze,RZ)
          enddo
          id02bjf=istate
          if (istate>0) id02bjf=0
#endif
          
         case default
          
          write(0,'("ad00aa : input error ; Lode="i2" is not supported ;")') pp%Lode ; iad00aa = 1 ; goto 9999
          
         end select
         
         pp%F(ii,jj,kk) = .not.Lbfieldok
         
         if( pp%F(ii,jj,kk) ) cycle ! the fieldline integration has failed; 29 Jul 15;
         
        !FATAL(ad00aa, pressure%lB(ifb)      .lt.small, illegal Bphi/B^2 factors )
        !FATAL(ad00aa, pressure%lB(ifb)*small.gt.one  , illegal Bphi/B^2 factors )
         
         if( RZ(1).lt.Rmin ) RZ(1) = Rmin
         if( RZ(1).gt.Rmax ) RZ(1) = Rmax
         if( RZ(2).lt.Zmin ) RZ(2) = Zmin
         if( RZ(2).gt.Zmax ) RZ(2) = Zmax
         
         pp%I(lfb,ii,jj,kk) = ( RZ(1)-Rmin ) / DR
         pp%J(lfb,ii,jj,kk) = ( RZ(2)-Zmin ) / DZ
         
         if( pp%I(lfb,ii,jj,kk).eq.NR   ) pp%I(lfb,ii,jj,kk) = NR-1
         if( pp%J(lfb,ii,jj,kk).eq.NZ   ) pp%J(lfb,ii,jj,kk) = NZ-1
         
        !FATAL(ad00aa, pp%I(lfb,ii,jj,kk).lt.0 .or. pp%I(lfb,ii,jj,kk).ge.NR, counting error on R )
        !FATAL(ad00aa, pp%J(lfb,ii,jj,kk).lt.0 .or. pp%J(lfb,ii,jj,kk).ge.NZ, counting error on Z )
         
         pp%x(lfb,ii,jj,kk) = ( RZ(1) - Rmin - pp%I(lfb,ii,jj,kk) * DR ) / DR
         pp%y(lfb,ii,jj,kk) = ( RZ(2) - Zmin - pp%J(lfb,ii,jj,kk) * DZ ) / DZ
         
        !FATAL(ad00aa, pp%x(lfb,ii,jj,kk).lt.zero-small .or. pp%x(lfb,ii,jj,kk).gt.one+small, interpolation error on R )
        !FATAL(ad00aa, pp%y(lfb,ii,jj,kk).lt.zero-small .or. pp%y(lfb,ii,jj,kk).gt.one+small, interpolation error on Z )
         
        enddo ! end of do ifb; 29 Jul 15;
        
        if( pp%F(ii,jj,kk) ) cycle ! the fieldline integration has failed; 29 Jul 15;
        
        pp%B(-1:1,ii,jj,kk) = pressure%lB(-1:1)
        
       enddo ! end of do jj ; 29 Jul 15;
      enddo ! end of do ii ; 29 Jul 15;
      if( ifail.le.-2 ) write(0,'("ad00aa :         "2x" : completed toroidal plane kk = "i3" /"i3" ;")') kk, Np
     enddo ! end of do kk ; 29 Jul 15;
     
    endif ! end of if( pp%Lcontrol.eq.1 .or. pp%Lcontrol.eq.3 ) then ; 29 Jul 15;
    
    
    if( pp%Lcontrol.ne.2 .and. pp%Lcontrol.ne.3 ) goto 9998
    

! can now perform bi-linear interpolation and construct the parallel and perpendicular diffusion; 29 Jul 15;

    
    if( .not.allocated(pp%p) ) then
      write(0,'("ad00aa : input error ; p is not allocated ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( .not.allocated(pp%s) ) then
      write(0,'("ad00aa : input error ; s is not allocated ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( .not.allocated(pp%x) ) then
      write(0,'("ad00aa : input error ; x is not allocated ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( .not.allocated(pp%y) ) then
      write(0,'("ad00aa : input error ; y is not allocated ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( .not.allocated(pp%I) ) then
      write(0,'("ad00aa : input error ; I is not allocated ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( .not.allocated(pp%J) ) then
      write(0,'("ad00aa : input error ; J is not allocated ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( .not.allocated(pp%B) ) then

      write(0,'("ad00aa : input error ; B is not allocated ;")') ; iad00aa = 1 ; goto 9999
    endif
    if( .not.allocated(pp%F) ) then
      write(0,'("ad00aa : input error ; F is not allocated ;")') ; iad00aa = 1 ; goto 9999
    endif
    
    if( pp%Ntime.le.0     ) then
      write(0,'("ad00aa :             ; Ntime.le.0     ;")') ;             ; goto 9998
    endif
    if( pp%dtime.le.small ) then 
      write(0,'("ad00aa : input error ; dtime.le.small ;")') ; iad00aa = 1 ; goto 9999
    endif
   !if( pp%kpara.le.small ) then
   !  write(0,'("ad00aa : input error ; kpara.le.small ;")') ; iad00aa = 1 ; goto 9999
   !endif
   !if( pp%kperp.le.small ) then
   !  write(0,'("ad00aa : input error ; kperp.le.small ;")') ; iad00aa = 1 ; goto 9999
   !endif

    
    pp%p(0:NR,0:NZ,-1   ) =       pp%p(0:NR,0:NZ,   Np-1) ! periodic boundary conditions; 29 Jul 15;
    pp%p(0:NR,0:NZ,   Np) =       pp%p(0:NR,0:NZ, 0     ) ! periodic boundary conditions; 29 Jul 15;
    

     ! macro expansion of sallocate = set allocate;

#ifdef DEBUG
   if( allocated( pressure%p ) ) then
    write(0,'("oculus : 0123456789 : pressure%p already allocated ;")') 
    stop      'oculus : 0123456789 : pressure%p already allocated ;'
   endif
#endif

   allocate( pressure%p(0:NR,0:NZ,-1:Np), stat=astat )

#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("oculus : 0123456789 : error allocating pressure%p ;")') 
    stop      'oculus : 0123456789 : error allocating pressure%p ;'
   endif
#endif

   pressure%p(0:NR,0:NZ,-1:Np) = zero


    

    if( ifail.le.-1 ) write(0,1004) pp%kpara, pp%kperp, pp%NR, pp%NZ, pp%Np, pp%dtime, pp%Ntime
1004 format("ad00aa :         "2x" : relaxing ; kpara="es11.4" ; kperp="es11.4 &
            " ; (NR,NZ,Np)=("i4" ,"i4" ,"i4" ) ; dtime="es11.04" ; Ntime="i9" ;")
    

    po = maxval(abs(pp%p(0:NR,0:NZ, 0:Np-1)))
    pn = maxval(abs(pp%p(0:NR,0:NZ, 0:Np-1)))
    

    itime = 0 
    
    
    if( ifail.le.-2 ) write(0,1005) itime, pn
    
    
    do itime = 1, pp%Ntime
     
     
     pressure%p(0:NR,0:NZ,-1:Np) = zero
     
     
     do kk = 0, Np-1
      do ii = 1, NR-1
       do jj = 1, NZ-1
        
        if( pp%F(ii,jj,kk) ) then ! fieldline integration failed;  8 Jul 15;
          
         !pp%lB(-1:1) = pp%B(-1:1,NR/2,NZ/2,kk) ; interpress(-1:1) = zero
         pp%lB(-1:1) = (/ one ,  one ,  one /) ; interpress(-1:1) = zero
          
        else
          
         !pp%lB(-1:1) = pp%B(-1:1,ii  ,jj  ,kk)
         pp%lB(-1:1) = (/ one ,  one ,  one /)
          
         do ifb = -1, 1, 2 ; lfb = max(0,ifb) ! forward and backward mapped points;  8 Jul 15;
          
          I = pp%I(lfb,ii,jj,kk) ; J = pp%J(lfb,ii,jj,kk)
           
          if( I.ge.1 .and. I.le.NR-2 .and. J.ge.1 .and. J.le.NZ-2 ) then ! bi-cubic interpolation; 8 Jul 15;
            
           stencil(1:16) = (/ pp%p(I-1,J-1,kk), pp%p(I+0,J-1,kk), pp%p(I+1,J-1,kk), pp%p(I+2,J-1,kk), &
                              pp%p(I-1,J+0,kk), pp%p(I+0,J+0,kk), pp%p(I+1,J+0,kk), pp%p(I+2,J+0,kk), &
                              pp%p(I-1,J+1,kk), pp%p(I+0,J+1,kk), pp%p(I+1,J+1,kk), pp%p(I+2,J+1,kk), &
                              pp%p(I-1,J+2,kk), pp%p(I+0,J+2,kk), pp%p(I+1,J+2,kk), pp%p(I+2,J+2,kk) /)
            
           ff(0,0) = sum( (/  0,+ 0,+ 0,+ 0,+ 0,+ 1,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0 /) * stencil(1:16) ) /  1
           ff(0,1) = sum( (/  0,- 2,+ 0,+ 0,+ 0,- 3,+ 0,+ 0,+ 0,+ 6,+ 0,+ 0,+ 0,- 1,+ 0,+ 0 /) * stencil(1:16) ) /  6
           ff(0,2) = sum( (/  0,+ 1,+ 0,+ 0,+ 0,- 2,+ 0,+ 0,+ 0,+ 1,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0 /) * stencil(1:16) ) /  2
           ff(0,3) = sum( (/  0,- 1,+ 0,+ 0,+ 0,+ 3,+ 0,+ 0,+ 0,- 3,+ 0,+ 0,+ 0,+ 1,+ 0,+ 0 /) * stencil(1:16) ) /  6
           ff(1,0) = sum( (/  0,+ 0,+ 0,+ 0,- 2,- 3,+ 6,- 1,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0 /) * stencil(1:16) ) /  6
           ff(1,1) = sum( (/  4,+ 6,-12,+ 2,+ 6,+ 9,-18,+ 3,-12,-18,+36,- 6,+ 2,+ 3,- 6,+ 1 /) * stencil(1:16) ) / 36
           ff(1,2) = sum( (/- 2,- 3,+ 6,- 1,+ 4,+ 6,-12,+ 2,- 2,- 3,+ 6,- 1,+ 0,+ 0,+ 0,+ 0 /) * stencil(1:16) ) / 12
           ff(1,3) = sum( (/  2,+ 3,- 6,+ 1,- 6,- 9,+18,- 3,+ 6,+ 9,-18,+ 3,- 2,- 3,+ 6,- 1 /) * stencil(1:16) ) / 36
           ff(2,0) = sum( (/  0,+ 0,+ 0,+ 0,+ 1,- 2,+ 1,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0 /) * stencil(1:16) ) /  2
           ff(2,1) = sum( (/- 2,+ 4,- 2,+ 0,- 3,+ 6,- 3,+ 0,+ 6,-12,+ 6,+ 0,- 1,+ 2,- 1,+ 0 /) * stencil(1:16) ) / 12
           ff(2,2) = sum( (/  1,- 2,+ 1,+ 0,- 2,+ 4,- 2,+ 0,+ 1,- 2,+ 1,+ 0,+ 0,+ 0,+ 0,+ 0 /) * stencil(1:16) ) /  4
           ff(2,3) = sum( (/- 1,+ 2,- 1,+ 0,+ 3,- 6,+ 3,+ 0,- 3,+ 6,- 3,+ 0,+ 1,- 2,+ 1,+ 0 /) * stencil(1:16) ) / 12
           ff(3,0) = sum( (/  0,+ 0,+ 0,+ 0,- 1,+ 3,- 3,+ 1,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0,+ 0 /) * stencil(1:16) ) /  6
           ff(3,1) = sum( (/  2,- 6,+ 6,- 2,+ 3,- 9,+ 9,- 3,- 6,+18,-18,+ 6,+ 1,- 3,+ 3,- 1 /) * stencil(1:16) ) / 36
           ff(3,2) = sum( (/- 1,+ 3,- 3,+ 1,+ 2,- 6,+ 6,- 2,- 1,+ 3,- 3,+ 1,+ 0,+ 0,+ 0,+ 0 /) * stencil(1:16) ) / 12
           ff(3,3) = sum( (/  1,- 3,+ 3,- 1,- 3,+ 9,- 9,+ 3,+ 3,- 9,+ 9,- 3,- 1,+ 3,- 3,+ 1 /) * stencil(1:16) ) / 36
           
           ! if suff. memory then can save this information;  8 Jul 15;
           xx(0) = one ; xx(1) = pp%x(lfb,ii,jj,kk) ; xx(2) = xx(1)*xx(1) ; xx(3) = xx(2)*xx(1)
           yy(0) = one ; yy(1) = pp%y(lfb,ii,jj,kk) ; yy(2) = yy(1)*yy(1) ; yy(3) = yy(2)*yy(1)
           
           interpress(ifb) = zero
           do ix = 0, 3
            do iy = 0, 3 ; interpress(ifb) = interpress(ifb) + ff(ix,iy) * xx(ix) * yy(iy) ! bi-cubic interpolation;  8 Jul 15;
            enddo
           enddo
           
          else ! use bi-linear interpolation near the computational boundary;  8 Jul 15;
           
           g0 = pp%p(I  ,J  ,kk+ifb ) * ( one - pp%x(lfb,ii,jj,kk) ) + pp%p(I+1,J  ,kk+ifb  ) * pp%x(lfb,ii,jj,kk)
           g1 = pp%p(I  ,J+1,kk+ifb ) * ( one - pp%x(lfb,ii,jj,kk) ) + pp%p(I+1,J+1,kk+ifb  ) * pp%x(lfb,ii,jj,kk)
           interpress(ifb) &  ! bi-linear interpolation;
              = g0                    * ( one - pp%y(lfb,ii,jj,kk) ) + g1                     * pp%y(lfb,ii,jj,kk)
           
          endif ! end of if( I.ge.1 .and. I.le.NR-2 .and . . . . ) ; bi-cubic interpolation;  8 Jul 15;
          
         enddo ! end of do ifb; 29 Jul 15;
         
        endif ! end of if( pp%F(ii,jj,kk) ) ;  8 Jul 15;
                
        dpara = ( ( interpress(+1) - pp%p(ii,jj,kk) ) * pp%lB(+1) &
                - ( pp%p(ii,jj,kk) - interpress(-1) ) * pp%lB(-1) ) * pp%lB(0) / DpDp
        
        
       !R0 = Rmin + ii * DR ; Rh(-4:4) = R0 + HR + (/ -4, -3, -2, -1,  0,  1, 2, 3, 4 /) * DR
        R0 = one            ; Rh(-4:4) = one

        if( ii.ge.3 .and. ii.le.NR-3 .and. jj.ge.3 .and. jj.le.NZ-3 ) then ! 4th-order, symmetric differences;  8 Jul 15;
         
         do ix = -2, 1 ; I = ii + ix
          do iy = -2, 1 ; J = jj + iy
           
           stencil(1:16) = (/ pp%p(I-1,J-1,kk), pp%p(I-1,J  ,kk), pp%p(I-1,J+1,kk), pp%p(I-1,J+2,kk), &
                              pp%p(I  ,J-1,kk), pp%p(I  ,J  ,kk), pp%p(I  ,J+1,kk), pp%p(I  ,J+2,kk), &
                              pp%p(I+1,J-1,kk), pp%p(I+1,J  ,kk), pp%p(I+1,J+1,kk), pp%p(I+1,J+2,kk), &
                              pp%p(I+2,J-1,kk), pp%p(I+2,J  ,kk), pp%p(I+2,J+1,kk), pp%p(I+2,J+2,kk) /)
           
           dpdx(ix,iy) = sum( (/ -  1, +  9, +  9, -  1, &
                                 + 27, -243, -243, + 27, &
                                 - 27, +243, +243, - 27, &
                                 +  1, -  9, -  9, +  1  /) * stencil(1:16) ) / 384 / DR
           
           dpdy(ix,iy) = sum( (/ -  1, + 27, - 27, +  1, &
                                 +  9, -243, +243, -  9, &
                                 +  9, -243, +243, -  9, &
                                 -  1, + 27, - 27, +  1  /) * stencil(1:16) ) / 384 / DZ
           
          enddo ! end of do iy ;  8 Jul 15;
         enddo ! end of do ix ;  8 Jul 15;
         
        !sauce = sum( (/ +  1, -  9, -  9, +  1 /) * Rh(-2) * pp%s(ii-2,jj-2:jj+1,kk) &
        !           + (/ -  9, + 81, + 81, -  9 /) * Rh(-1) * pp%s(ii-1,jj-2:jj+1,kk) &
        !           + (/ -  9, + 81, + 81, -  9 /) * Rh( 0) * pp%s(ii  ,jj-2:jj+1,kk) &
        !           + (/ +  1, -  9, -  9, +  1 /) * Rh(+1) * pp%s(ii+1,jj-2:jj+1,kk) ) / 256      / R0
         
         dpdRR = sum( (/ -  1, +  9, +  9, -  1 /) * Rh(-2) * dpdx(  -2,  -2:   1   ) &
                    + (/ + 27, -243, -243, + 27 /) * Rh(-1) * dpdx(  -1,  -2:   1   ) &
                    + (/ - 27, +243, +243, - 27 /) * Rh( 0) * dpdx(   0,  -2:   1   ) &
                    + (/ +  1, -  9, -  9, +  1 /) * Rh( 1) * dpdx(   1,  -2:   1   ) ) / 384 / DR / R0
         
         dpdZZ = sum( (/ -  1, + 27, - 27, +  1 /) * Rh(-2) * dpdy(  -2,  -2:   1   ) &
                    + (/ +  9, -243, +243, -  9 /) * Rh(-1) * dpdy(  -1,  -2:   1   ) &
                    + (/ +  9, -243, +243, -  9 /) * Rh( 0) * dpdy(   0,  -2:   1   ) &
                    + (/ -  1, + 27, - 27, +  1 /) * Rh( 1) * dpdy(   1,  -2:   1   ) ) / 384 / DZ / R0
         
        else ! 4th-order, non-symmetric differences near the computational boundary;  8 Jul 15;
         
         if    ( ii.eq.1                  ) then
           dpdRR    = sum( (/ +11, -20, + 6, + 4, - 1 /) * pp%p(ii-1:ii+3,jj,kk) ) / twelve / DRDR !
          !         + sum( (/ - 3, -10, +18, - 6, + 1 /) * pp%p(ii-1:ii+3,jj,kk) ) / twelve / DR / R0
          ! ; sauce = sum( (/ + 5, +15, - 5, + 1      /) * pp%s(ii-1:ii+2,jj-1,kk) * Rh(-1:2) &
          !              + (/ + 5, +15, - 5, + 1      /) * pp%s(ii-1:ii+2,jj  ,kk) * Rh(-1:2) ) * half / 16 / R0
         elseif( ii.gt.1 .and. ii.lt.NR-1 ) then
           dpdRR    = sum( (/ - 1, +16, -30, +16, - 1 /) * pp%p(ii-2:ii+2,jj,kk) ) / twelve / DRDR !
          !         + sum( (/ + 1, - 8,   0, + 8, - 1 /) * pp%p(ii-2:ii+2,jj,kk) ) / twelve / DR / R0
          ! ; sauce = sum( (/ - 1, + 9, + 9, - 1      /) * pp%s(ii-2:ii+1,jj-1,kk) * Rh(-2:1) &
          !              + (/ - 1, + 9, + 9, - 1      /) * pp%s(ii-2:ii+1,jj  ,kk) * Rh(-2:1) ) * half / 16 / R0
         elseif(               ii.eq.NR-1 ) then
           dpdRR    = sum( (/ - 1, + 4, + 6, -20, +11 /) * pp%p(ii-3:ii+1,jj,kk) ) / twelve / DRDR !
          !         + sum( (/ - 1, + 6, -18, +10, + 3 /) * pp%p(ii-3:ii+1,jj,kk) ) / twelve / DR / R0
          ! ; sauce = sum( (/ + 1, - 5, +15, + 5      /) * pp%s(ii-3:ii  ,jj-1,kk) * Rh(-3:0) &
          !              + (/ + 1, - 5, +15, + 5      /) * pp%s(ii-3:ii  ,jj  ,kk) * Rh(-3:0) ) * half / 16 / R0
         endif
          
         if    ( jj.eq.1                  ) then
            dpdZZ = sum( (/ +11, -20, + 6, + 4, - 1 /) * pp%p(ii,jj-1:jj+3,kk) ) / twelve / DZDZ
         elseif( jj.gt.1 .and. jj.lt.NZ-1 ) then
            dpdZZ = sum( (/ - 1, +16, -30, +16, - 1 /) * pp%p(ii,jj-2:jj+2,kk) ) / twelve / DZDZ
         elseif(               jj.eq.NZ-1 ) then
            dpdZZ = sum( (/ - 1, + 4, + 6, -20, +11 /) * pp%p(ii,jj-3:jj+1,kk) ) / twelve / DZDZ
         endif
         
        endif ! end of if( ii.ge.3 .and. ii.le.NR-3 .and. jj.ge.3 .and. jj.le.NZ-3 ) then ;  8 Jul 15;
         
        dperp = dpdRR + dpdZZ

       !pressure%p(ii,jj,kk) = pp%p(ii,jj,kk) &
       !                     + pp%dtime * ( ( pp%kpara - pp%kperp ) * dpara + pp%kperp * dperp + pp%s(ii,jj,kk) )
       !pressure%p(ii,jj,kk) = pp%p(ii,jj,kk) &
       !                     + pp%dtime * (   pp%kpara              * dpara + pp%kperp * dperp + pp%s(ii,jj,kk) * R0 )
        pressure%p(ii,jj,kk) = pp%p(ii,jj,kk) &
                             + pp%dtime * (   pp%kpara              * dpara + pp%kperp * dperp + pp%s(ii,jj,kk)      )
       !pressure%p(ii,jj,kk) = pp%p(ii,jj,kk) &
       !                     + pp%dtime * (   pp%kpara              * dpara + pp%kperp * dperp + sauce               )

       enddo ! end of do jj; 29 Jul 15;
      enddo ! end of do ii; 29 Jul 15;
      
     enddo ! end of do kk; 29 Jul 15;
     
     pp%p(0:NR,0:NZ,-1:Np) = pressure%p(0:NR,0:NZ,-1:Np  )     
     
     pp%p(0:NR,0:NZ,-1   ) =       pp%p(0:NR,0:NZ,   Np-1) ! periodic boundary conditions; 29 Jul 15;
     pp%p(0:NR,0:NZ,   Np) =       pp%p(0:NR,0:NZ, 0     ) ! periodic boundary conditions; 29 Jul 15;
     
     pn = maxval(abs(pp%p(0:NR,0:NZ, 0:Np-1)))
     
     select case( ifail )
     case( -2 )
      if( (itime/10000)*10000.eq.itime ) write(0,1005) itime, pn, pn-po
     case( -3 )
      if( (itime/1000 )*1000 .eq.itime ) write(0,1005) itime, pn, pn-po
     case( -4 )
      if( (itime/100  )*100  .eq.itime ) write(0,1005) itime, pn, pn-po
     case( -5 )
      if( (itime/10   )*10   .eq.itime ) write(0,1005) itime, pn, pn-po
     case( -6 )
      ;                                                  write(0,1005) itime, pn, pn-po
     end select
     
1005 format("ad00aa : "i9"  : max(p) ="es24.17" ;":" dif(p) ="es13.5" ;")
     
     if( pn.lt.small .or. pn*small.gt.one ) then
       iad00aa = 2 ; deallocate(pressure%p) ; goto 9999 ! avoid over/under-flow if CFL condition is violated;
     endif
     
    enddo ! end of do itime; 29 Jul 15;
    
     ! macro expansion of dallocate = deallocate;

#ifdef DEBUG
   if( .not. allocated( pressure%p ) ) then
    write(0,'("oculus : 0123456789 : pressure%p not allocated ;")') 
    stop      'oculus : 0123456789 : pressure%p not allocated ;'
   endif
#endif

   deallocate( pressure%p, stat=astat )

#ifdef DEBUG
   if( astat.ne.0 ) then    
    write(0,'("oculus : 0123456789 : error de-allocating pressure%p ;")') 
    stop      'oculus : 0123456789 : error de-allocating pressure%p ;'
   endif
#endif


    
9998 continue
    
    iad00aa = 0
    
9999 continue
    
    if( ifail.le. 0 ) write(0,9000) iad00aa, pp%kpara, pp%kperp, pp%NR, pp%Np, pp%NZ, pp%dtime, pp%Ntime, pn, pn-po
    
    ifail = iad00aa
    
9000 format("ad00aa : ifail =",i3," : kpara=",es9.2," ; kperp=",es9.2," ; (NR,Np,NZ)=(",i5,",",i4,",",i5, &
            " ) ; dtime=",es11.04," ; Ntime=",i8," ; |p|=",es12.05," ; d|p|="es13.05" ;")
    
    return
    
  end subroutine ad00aa  
  
  subroutine ad00ab( phi, RZ )
    
    implicit none
    
    integer, parameter :: Node = 2
    
    integer            :: ibfield
    real               :: phi, RZ(1:Node), RpZ(1:3), dBRpZ(1:3,0:3), BB, Bphi
    
    if( pressure%ii.eq.-1 .or. pressure%ii.eq.+1 ) then
     
     RpZ(1:3) = (/ RZ(1), phi, RZ(2) /) ; ibfield = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), ibfield )
     
     Bphi = dBRpZ(2,0) ; BB = sqrt( dBRpZ(1,0)**2 + (RpZ(1)*dBRpZ(2,0))**2 + dBRpZ(3,0)**2 )
     
     if( ibfield.ne.0 .or. abs(Bphi).lt.small .or. BB.lt.small ) then ; Lbfieldok = .false.
     else                                                             ; pressure%lB(pressure%ii) = abs(Bphi) / BB**2
     endif
     
    endif
    
    pressure%ii = pressure%ii + pressure%ifb ; phi = phi + pressure%ifb * half * pressure%Dp
    
    return
    
  end subroutine ad00ab

!latex \subroutine{bn00aa}{compute $\left({\bf B}\cdot{\bf n}\right)_{m,n}$ on a given toroidal surface;}
!latex \bi
!latex \item Given a toroidal, `control' surface, 
!latex \be {\bf x}(\t,\z) \equiv R(\t,\z) \cos\z \, {\bf i} + R(\t,\z) \sin\z \, {\bf j} + Z(\t,\z) \, {\bf k},
!latex \ee
!latex we may compute
!latex \be B^n(\t,\z) \equiv {\bf B} \cdot {\bf e}_\t \times {\bf e}_\z.
!latex \ee
!latex \item With the magnetic field given in cylindrical coordinates, ${\bf B} \equiv B^R {\bf e}_R + B^\phi {\bf e}_\phi + B^Z {\bf e}_Z$,
!latex \be {\bf B} \cdot {\bf e}_\t \times {\bf e}_\z = \left[ B^R ( -Z_\t R ) + B^\phi ( Z_\t R_\z - R_\t Z_\z ) + B^Z ( R_\t R ) \right] R.
!latex \ee
!latex \item The routine also returns the toroidal and poloidal currents, defined as surface integrals through appropriate surfaces,
!latex \be I \equiv \int_{\cal S} {\bf j}\cdot d{\bf s} = \int_{\partial {\cal S}} \!\! {\bf B} \cdot d{\bf l} = \int_0^{2\pi} \!\!\! B_\t d\t, \\
!latex     G \equiv \int_{\cal S} {\bf j}\cdot d{\bf s} = \int_{\partial {\cal S}} \!\! {\bf B} \cdot d{\bf l} = \int_0^{2\pi} \!\!\! B_\z d\z.
!latex \ee
!latex These integrals are calculated using \nag{www.nag.co.uk/numeric/FL/manual19/pdf/D01/d01ahf_fl19.pdf}{D01AHF}.
!latex \item[7.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : bnormal, bn00aa+ \\ \\
!latex           \verb+type(bnormal) :: bn+ \\ \\ in their source that calls \verb+bn00aa+,
!latex           where \verb+bn+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+bn+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+bn%cmn                  : integer ;+
!latex \bi
!latex \item[i.] total number of Fourier harmonics that describe the input control surface;
!latex \ei
!latex \item[  ] \verb+bn%cim(1:cmn)           : integer array ;+
!latex \item[  ] \verb+bn%cin(1:cmn)           : integer array ;+
!latex \bi
!latex \item[i.] Fourier mode identification;
!latex \ei
!latex \item[  ] \verb+bn%Rcc(1:cmn)           : real array;+
!latex \item[  ] \verb+bn%Rcs(1:cmn)           : real array;+
!latex \item[  ] \verb+bn%Zcc(1:cmn)           : real array;+
!latex \item[  ] \verb+bn%Zcs(1:cmn)           : real array;+
!latex \bi
!latex \item[i.] Fourier harmonics of control surface;
!latex \ei
!latex \item[  ] \verb+bn%Mpol                 : integer;+
!latex \item[  ] \verb+bn%Ntor                 : integer;+
!latex \item[  ] \verb+bn%Nfp                  : integer;+
!latex \bi
!latex \item[i.] Poloidal and toroidal Fourier resolution, and field-periodicity;
!latex \ei
!latex \item[  ] \verb+bn%tol                  : real;+
!latex \bi
!latex \item[i.] Relative accuracy required; 
!latex           \verb+tol+ $\equiv$ \verb+EPSR+ provided to \nag{www.nag.co.uk/numeric/FL/manual19/pdf/D01/d01ahf_fl19.pdf}{D01AHF}.
!latex \ei
!latex \item[  ] \verb+ifail                   : integer ;+
!latex \bi
!latex \item[i.] if \verb+ifail =  1+, quiet mode;
!latex \item[i.] if \verb+ifail =  0+, screen output is terse;
!latex \item[i.] if \verb+ifail = -1+, 
!latex \ei
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call bn00aa( bn, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+bn%mn                    : integer      ;+
!latex \item[  ] \verb+bn%im(1:mn)              : integer array;+
!latex \item[  ] \verb+bn%in(1:mn)              : integer array;+
!latex \item[  ] \verb+bn%gBc(1:mn)             : real    array;+
!latex \item[  ] \verb+bn%gBs(1:mn)             : real    array;+
!latex \bi
!latex \item[i.] Fourier harmonics of normal field to control surface;
!latex \ei
!latex \item[  ] \verb+bn%Itor                  : real;+
!latex \item[  ] \verb+bn%Gpol                  : real;+
!latex \bi
!latex \item[i.] enclosed currents;
!latex \ei
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex     \item[] \verb+ifail=2+ : the integration of the toroidal current failed;
!latex     \item[] \verb+ifail=3+ : the integration of the linking  current failed;
!latex \ei
!latex \ei

subroutine bn00aa( lbn, ifail )
  
  implicit none
  
  type(bnormal)         :: lbn
  
  INTEGER               :: ifail, astat, ii, jj, kk, jk, itangent, lfail, ift02aa
  REAL                  :: teta, zeta, arg, carg, sarg, RR, ZZ, dRt, dZt, dRz, dZz, szeta, czeta
  REAL                  :: xx(1:3), xt(1:3), xz(1:3), ds(1:3)
  REAL                  :: RpZ(1:3), dBRpZ(1:3,0:3), Bx, By, Bz ! input/output from user-supplied bfield; 23 Nov 15;
  INTEGER               :: npts, nlimit, id01ahf
  REAL                  :: lowlimit, upplimit, epsr, relerr, d01ahf
  REAL, allocatable     :: imag(:), even(:), oddd(:)
  
  ibn00aa = ifail

  if( allocated(lbn%im ) ) deallocate(lbn%im ) ! these are output quantities; 08 Aug 16;
  if( allocated(lbn%in ) ) deallocate(lbn%in )
  if( allocated(lbn%gBi) ) deallocate(lbn%gBi) ! "normal" field on real-space grid; 08 Aug 16;
  if( allocated(lbn%gBc) ) deallocate(lbn%gBc)
  if( allocated(lbn%gBs) ) deallocate(lbn%gBs)

  CHECKINPUT(bn00aa, lbn%cmn  .le. 0         , 9999 )

  CHECKINPUT(bn00aa, .not.allocated(lbn%cim ), 9999 )
  CHECKINPUT(bn00aa, .not.allocated(lbn%cin ), 9999 )

  CHECKINPUT(bn00aa, .not.allocated(lbn%Rcc) , 9999 )
  CHECKINPUT(bn00aa, .not.allocated(lbn%Rcs) , 9999 )
  CHECKINPUT(bn00aa, .not.allocated(lbn%Zcc) , 9999 )
  CHECKINPUT(bn00aa, .not.allocated(lbn%Zcs) , 9999 )

  CHECKINPUT(bn00aa, lbn%Mpol.lt. 1          , 9999 )
  CHECKINPUT(bn00aa, lbn%Ntor.lt. 0          , 9999 )
  CHECKINPUT(bn00aa, lbn%Nfp .lt. 1          , 9999 )

  CHECKINPUT(bn00aa, lbn%tol .le. zero       , 9999 )

  ibn00aa = 0

  if( ifail.le.-2 ) write(0,8999) ifail, lbn%cmn, lbn%Mpol, lbn%Ntor, lbn%Nfp

  if( ifail.le.-3 ) then 
    write(0,8998) ( lbn%cim(ii), lbn%cin(ii), lbn%Rcc(ii), lbn%Rcs(ii), lbn%Zcc(ii), lbn%Zcs(ii), ii = 1, lbn%cmn )
  endif

  lbn%mn = lbn%Ntor + 1 + lbn%Mpol * ( 2 * lbn%Ntor + 1 )

  lbn%Nt = 4 * lbn%Mpol ; lbn%Nz = 4 * max( lbn%Ntor, 1 ) ; lbn%Ntz = lbn%Nt * lbn%Nz

  SALLOCATE(lbn%im,(1:lbn%mn),0)
  SALLOCATE(lbn%in,(1:lbn%mn),0)

  SALLOCATE(lbn%gBi,(1:lbn%Ntz),zero)

  SALLOCATE(imag   ,(1:lbn%Ntz),zero) ! internal workspace; deallocated below; 23 Nov 15;
  
  call gi00ab( lbn%Mpol, lbn%Ntor, lbn%Nfp, lbn%mn, lbn%im(1:lbn%mn), lbn%in(1:lbn%mn) ) ! set Fourier modes; 23 Nov 15;
  
  do kk = 0, lbn%Nz-1 ;                         ; zeta = kk * (pi2/lbn%Nfp) / lbn%Nz

   czeta = cos(zeta) ; szeta = sin(zeta)
   
   do jj = 0, lbn%Nt-1 ; jk = 1 + jj + kk*lbn%Nt ; teta = jj *  pi2          / lbn%Nt
    
     RR = zero
     ZZ = zero 
    
    dRt = zero
    dZt = zero
    
    dRz = zero
    dZz = zero
    
    do ii = 1, lbn%cmn ; arg = lbn%cim(ii) * teta - lbn%cin(ii) * zeta ; carg = cos(arg) ; sarg = sin(arg)
     
      RR =  RR +     lbn%Rcc(ii) * carg + lbn%Rcs(ii) * sarg
      ZZ =  ZZ +     lbn%Zcc(ii) * carg + lbn%Zcs(ii) * sarg
     
     dRt = dRt + ( - lbn%Rcc(ii) * sarg + lbn%Rcs(ii) * carg ) * lbn%cim(ii)
     dZt = dZt + ( - lbn%Zcc(ii) * sarg + lbn%Zcs(ii) * carg ) * lbn%cim(ii)
     
     dRz = dRz - ( - lbn%Rcc(ii) * sarg + lbn%Rcs(ii) * carg ) * lbn%cin(ii)
     dZz = dZz - ( - lbn%Zcc(ii) * sarg + lbn%Zcs(ii) * carg ) * lbn%cin(ii)
     
    enddo ! end of do ii; 30 Oct 15;
    
    RpZ(1:3) = (/ RR, zeta, ZZ /) ; itangent = 0 ; lfail = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), lfail )

    Bx = dBRpZ(1,0) * czeta - dBRpZ(2,0) * RR * szeta ! Cartesian components of magnetic field; 23 Nov 15;
    By = dBRpZ(1,0) * szeta + dBRpZ(2,0) * RR * czeta
    Bz = dBRpZ(3,0)
    
   !xx(1:3) = (/  RR * czeta,  RR * szeta,  ZZ  /)
    
    xt(1:3) = (/ dRt * czeta, dRt * szeta, dZt  /)
    
    xz(1:3) = (/ dRz * czeta, dRz * szeta, dZz  /) &
            + (/ -RR * szeta, +RR * czeta, zero /)
  
    ds(1:3) = (/ xt(2) * xz(3) - xt(3) * xz(2), &
                 xt(3) * xz(1) - xt(1) * xz(3), &
                 xt(1) * xz(2) - xt(2) * xz(1) /)

    lbn%gBi(jk) = Bx * ds(1) + By * ds(2) + Bz * ds(3)

    if( ifail.le.-4 ) write(0,8997) kk, zeta, jj, jk, teta, lbn%gBi(jk)
    
   enddo ! end of do jj; 23 Nov 15;
  enddo ! end of do kk; 23 Nov 15;

  SALLOCATE(lbn%gBc,(1:lbn%mn),zero)
  SALLOCATE(lbn%gBs,(1:lbn%mn),zero)

  SALLOCATE(even   ,(1:lbn%mn),zero) ! workspace; 08 Aug 16;
  SALLOCATE(oddd   ,(1:lbn%mn),zero)
  

  ift02aa = 0 ; imag(1:lbn%Ntz) = zero

  call ft02aa( lbn%Nfp, lbn%Nt, lbn%Nz, lbn%gBi(1:lbn%Ntz), imag(1:lbn%Ntz), lbn%mn, lbn%im(1:lbn%mn), & ! not-destructive; 13 Oct 15;
               lbn%in(1:lbn%mn), lbn%gBc(1:lbn%mn), lbn%gBs(1:lbn%mn), even(1:lbn%mn), oddd(1:lbn%mn), ift02aa )
  
  DALLOCATE(even)
  DALLOCATE(oddd)
  DALLOCATE(imag)

  ;         bn%cmn           =lbn%cmn
  SALLOCATE(bn%cim,(1:bn%cmn),lbn%cim(1:bn%cmn))
  SALLOCATE(bn%cin,(1:bn%cmn),lbn%cin(1:bn%cmn))
  SALLOCATE(bn%Rcc,(1:bn%cmn),lbn%Rcc(1:bn%cmn))
  SALLOCATE(bn%Rcs,(1:bn%cmn),lbn%Rcs(1:bn%cmn))
  SALLOCATE(bn%Zcc,(1:bn%cmn),lbn%Zcc(1:bn%cmn))
  SALLOCATE(bn%Zcs,(1:bn%cmn),lbn%Zcs(1:bn%cmn))
  
  lowlimit = zero ; upplimit = pi2 ; epsr = lbn%tol ; nlimit = 0

#ifdef HAVE_NAG
  id01ahf = 1
  lbn%Itor = D01AHF( lowlimit, upplimit, epsr, npts, relerr, dBteta, nlimit, id01ahf ) ! NAG
#else
  call DQNG( dBteta, lowlimit, upplimit, lbn%tol, epsr, lbn%Itor, relerr, npts, id01ahf )
#endif

  lbn%Itor = lbn%Itor - pi2 ! see offset in definition of dBteta below; 08 Aug 16;
  
  if( id01ahf.ne.0 .and. ifail.le.0 ) then 
    write(0,'("bn00aa : Itor : npts="i9" ; relerr="es13.05" ; id01ahf="i2" ;")') npts, relerr, id01ahf
  endif
  if( id01ahf.ne.0                  ) ibn00aa = 2
  
#ifdef HAVE_NAG
  id01ahf = 1
  lbn%Gpol = D01AHF( lowlimit, upplimit, epsr, npts, relerr, dBzeta, nlimit, id01ahf ) ! NAG
#else
  call DQNG( dBzeta, lowlimit, upplimit, lbn%tol, epsr, lbn%Itor, relerr, npts, id01ahf )
#endif
  
  if( id01ahf.ne.0 .and. ifail.le.0 ) then
    write(0,'("bn00aa : Gpol : npts="i9" ; relerr="es13.05" ; id01ahf="i2" ;")') npts, relerr, id01ahf
  endif
  if( id01ahf.ne.0                  ) ibn00aa = 3
  
  DALLOCATE(bn%cim)
  DALLOCATE(bn%cin)
  DALLOCATE(bn%Rcc)
  DALLOCATE(bn%Rcs)
  DALLOCATE(bn%Zcc)
  DALLOCATE(bn%Zcs)
  
9999 continue
  
  if( ifail.le. 0 ) write(0,9000) ibn00aa, lbn%Mpol, lbn%Ntor, lbn%Itor, lbn%Gpol
  if( ifail.le.-1 ) write(0,9001) ( lbn%im(ii), lbn%in(ii), lbn%gBc(ii), lbn%gBs(ii), ii = 1, lbn%mn )
  
  ifail = ibn00aa
  
8997 format("bn00aa :        ",3x," ; kk="i4", zeta="f9.6", jj="i4", jk="i5", teta="f9.6", gBi="es23.15" ;")
8998 format("bn00aa :        ",3x," ; ("i3","i4" ) : Rc = ("es13.05","es13.05" ) ; Zc = ("es13.05","es13.05" ) ; ")
8999 format("bn00aa : ifail =",i3," ; cmn=",i3," ; Mpol=",i3," ; Ntor=",i3," ; Nfp=",i3," ;")
9000 format("bn00aa : ifail =",i3," ; Mpol=",i3," ; Ntor=",i3," ; Itor =",es23.15," ; Gpol =",es23.15," ;")
9001 format("bn00aa :        ",3x," ; ("i3","i4" ) : gB = ("es13.05","es13.05" ) ; ")
  
  return
  
end subroutine bn00aa

REAL function dBteta(teta)
  
  implicit none
  
  REAL                 :: teta
  
  INTEGER              :: ii, itangent, lfail
  REAL                 :: arg, zeta, czeta, szeta, carg, sarg, RR, ZZ, dRt, dZt, dRz, dZz
  REAL                 :: RpZ(1:3), dBRpZ(1:3,0:3), Bx, By, Bz, xx(1:3), xt(1:3), xz(1:3)
  
  zeta = zero
  
  czeta = cos(zeta) ; szeta = sin(zeta)

   RR = zero
   ZZ = zero 
  
  dRt = zero
  dZt = zero
  
  dRz = zero
  dZz = zero
  
  do ii = 1, bn%cmn ; arg = bn%cim(ii) * teta - bn%cin(ii) * zeta ; carg = cos(arg) ; sarg = sin(arg)
   
    RR =  RR +     bn%Rcc(ii) * carg + bn%Rcs(ii) * sarg
    ZZ =  ZZ +     bn%Zcc(ii) * carg + bn%Zcs(ii) * sarg
   
   dRt = dRt + ( - bn%Rcc(ii) * sarg + bn%Rcs(ii) * carg ) * bn%cim(ii)
   dZt = dZt + ( - bn%Zcc(ii) * sarg + bn%Zcs(ii) * carg ) * bn%cim(ii)
   
   dRz = dRz - ( - bn%Rcc(ii) * sarg + bn%Rcs(ii) * carg ) * bn%cin(ii)
   dZz = dZz - ( - bn%Zcc(ii) * sarg + bn%Zcs(ii) * carg ) * bn%cin(ii)
   
  enddo ! end of do ii; 30 Oct 15;
  
  RpZ(1:3) = (/ RR, zeta, ZZ /) ; itangent = 0 ; lfail = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), lfail )
  
  Bx = dBRpZ(1,0) * czeta - dBRpZ(2,0) * RR * szeta ! Cartesian components of magnetic field; 23 Nov 15;
  By = dBRpZ(1,0) * szeta + dBRpZ(2,0) * RR* czeta
  Bz = dBRpZ(3,0)
  
! xx(1:3) = (/  RR * czeta,  RR * szeta,  ZZ  /)
  
  xt(1:3) = (/ dRt * czeta, dRt * szeta, dZt  /)
  
! xz(1:3) = (/ dRz * czeta, dRz * szeta, dZz  /) &
!         + (/ -RR * szeta, +RR * czeta, zero /)

  dBteta = Bx * xt(1) + By * xt(2) + Bz * xt(3) + one ! added "offset" of one to normalize "relative error" to unity; 08 Aug 16;

  return
  
end function dBteta

REAL function dBzeta(zeta)
  
  implicit none
  
  REAL                 :: zeta
  
  INTEGER              :: ii, itangent, lfail
  REAL                 :: arg, teta, czeta, szeta, carg, sarg, RR, ZZ, dRt, dZt, dRz, dZz
  REAL                 :: RpZ(1:3), dBRpZ(1:3,0:3), Bx, By, Bz, xx(1:3), xt(1:3), xz(1:3)
  
  teta = zero
  
  czeta = cos(zeta) ; szeta = sin(zeta)

   RR = zero
   ZZ = zero 
  
  dRt = zero
  dZt = zero
  
  dRz = zero
  dZz = zero
  
  do ii = 1, bn%cmn ; arg = bn%cim(ii) * teta - bn%cin(ii) * zeta ; carg = cos(arg) ; sarg = sin(arg)
   
    RR =  RR +     bn%Rcc(ii) * carg + bn%Rcs(ii) * sarg
    ZZ =  ZZ +     bn%Zcc(ii) * carg + bn%Zcs(ii) * sarg
   
   dRt = dRt + ( - bn%Rcc(ii) * sarg + bn%Rcs(ii) * carg ) * bn%cim(ii)
   dZt = dZt + ( - bn%Zcc(ii) * sarg + bn%Zcs(ii) * carg ) * bn%cim(ii)
   
   dRz = dRz - ( - bn%Rcc(ii) * sarg + bn%Rcs(ii) * carg ) * bn%cin(ii)
   dZz = dZz - ( - bn%Zcc(ii) * sarg + bn%Zcs(ii) * carg ) * bn%cin(ii)
   
  enddo ! end of do ii; 30 Oct 15;
  
  RpZ(1:3) = (/ RR, zeta, ZZ /) ; itangent = 0 ; lfail = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), lfail )
  
  Bx = dBRpZ(1,0) * czeta - dBRpZ(2,0) * RR * szeta ! Cartesian components of magnetic field; 23 Nov 15;
  By = dBRpZ(1,0) * szeta + dBRpZ(2,0) * RR * czeta
  Bz = dBRpZ(3,0)
  
! xx(1:3) = (/  RR * czeta,  RR * szeta,  ZZ  /)
  
! xt(1:3) = (/ dRt * czeta, dRt * szeta, dZt  /)
  
  xz(1:3) = (/ dRz * czeta, dRz * szeta, dZz  /) &
          + (/ -RR * szeta, +RR * czeta, zero /)

  dBzeta = Bx * xz(1) + By * xz(2) + Bz * xz(3)

 !write(0,'("dBzeta : dBzeta, dBRpZ(2,0) = ",2es23.15," ;")') dBzeta, dBRpZ(2,0)

  return
  
end function dBzeta

!latex \newpage \section{toroidal-cylindrical coordinate-vector transformation} \label{sec:coordinatetransformation}

!latex \ben

!latex \item Specializing to coordinate transformations of the form
!latex       \be \begin{array}{rcl} R & = & R(\r,\t,\z), \\ \phi & = & \z, \\ Z & = & Z(\r,\t,\z), \end{array} \label{eq:cylindricaltoroidalcoordinates}
!latex       \ee
!latex       where position is given ${\bf x}\equiv R \cos\z \, {\bf i} + R \sin\z \, {\bf j} + Z        \, {\bf k}$,
!latex       the induced vector transformation is
!latex       \be \left( \begin{array}{c} B^R \\ B^\phi \\ B^Z \end{array} \right) = 
!latex       \left( \begin{array}{ccc}    R_\r, &    R_\t, &    R_\z \\ 
!latex                                       0, &       0, &       1 \\ 
!latex                                    Z_\r, &    Z_\t, &    Z_\z \end{array} \right) 
!latex           \left( \begin{array}{c} B^\r \\ B^\t \\ B^\z \end{array} \right). \label{eq:cylindricaltoroidalvector}
!latex       \ee
!latex \item This matrix is invertible if $\Delta \equiv R_\t Z_\r - R_\r Z_\t \ne 0$, and the coordinate Jacobian is
!latex       $\sqrt g \equiv R (R_\t Z_\r - R_\r Z_\t ) = R \Delta$.
!latex \item The coordinate transformation is written:
!latex       \be \begin{array}{ccc}R &=& \ds \sum_{m,n} \left[R_{m,n}(0) + \lambda_{m,n}(\r) \; X_{m,n}(\r)\right] \cos(m\t-n\z),\\
!latex                             Z &=& \ds \sum_{m,n} \left[Z_{m,n}(0) + \lambda_{m,n}(\r) \; Y_{m,n}(\r)\right] \sin(m\t-n\z);\end{array}
!latex       \label{eq:RZinterpolation}
!latex       \ee
!latex       where the $X_{m,n}(\r)$ and $Y_{m,n}(\r)$ are cubic-splines
!latex       (see \verb+Oculus:nr00aa+ and \verb+Oculus:nr00ab+ for details on the cubic-spline interpolation).
!latex \item The regularization factors are given by
!latex       \be \lambda_{m,n}(\r) = \left\{ \begin{array}{lcccc} \displaystyle 1                     &,& \mbox{\rm if } \; m = 0, \\ 
!latex                                                            \displaystyle \bar \r^{\overline{m}}&,& \mbox{\rm if } \; m > 0, \end{array} \right.
!latex       \ee
!latex       where e.g. $\overline{m} \equiv \min(m,2)$; and $\bar \r \equiv \r / V$, where $V$ is a normalization factor, e.g. $V\equiv $ total volume.
!latex       The following constraints must be enforced:
!latex       \bi 
!latex       \item[i.] for $m \ne 0$: $R_{m,n}(0)=0$ and $Z_{m,n}(0)=0$, and $X_{m,n}(0)$ and $Y_{m,n}(0)$ are arbitrary;
!latex       \item[ii.] for $m   = 0$: $X_{m,n}(0)=0$ and $Y_{m,n}(0)=0$;
!latex       \ei
!latex       The summation over $m$ and $n$ includes only the $\{(m,n):m=0;n=0,N\}$ and $\{(m,n):m=1,M; n=-N,N\}$ harmonics, which may be called the ``VMEC convention''.
!latex \item \Eqn{RZinterpolation} is for ``stellarator-symmetric'' equilibria
!latex       [Dewar \& Hudson, \link{dx.doi.org/10.1016/S0167-2789(97)00216-9}{Physica D 112 (1998) 275}].
!latex       For arbitrary geometry, additional sine and cosine harmonics are added.
!latex \item Inverting the vector transformation yields
!latex       \be \left( \begin{array}{c} \sqrt g B^\r \\ \sqrt g B^\t \\ \sqrt g B^\z \end{array} \right) = 
!latex       \left( \begin{array}{ccc} - Z_\t, & R_\z Z_\t - R_\t Z_\z, & + R_\t \\ 
!latex                                 + Z_\r, & R_\r Z_\z - R_\z Z_\r, & - R_\r \\ 
!latex                                      0, & R_\t Z_\r - R_\r Z_\t, &      0 \\ 
!latex              \end{array} \right) 
!latex           \left( \begin{array}{c} R B^R \\ R B^\phi \\ R B^Z \end{array} \right). \label{eq:toroidalcylindricalvector}
!latex       \ee
!latex \item The coordinates are ``right-handed'' if $\sqrt g >0$, and ``left-handed'' if $\sqrt g < 0$.
!latex       This can have important implications, particularly for \verb+qf00aa+ below.
!latex       If the coordinate transformation is $R=R_{0,0}+\r \cos\t$ and $Z=+\r \sin\t$, then the coordinates are left-handed; and 
!latex       if the coordinate transformation is $R=R_{0,0}+\r \cos\t$ and $Z=-\r \sin\t$, then the coordinates are right-handed.
!latex       Usually, the sign of Jacobian is opposite to the sign of $Z_{1,0}$.
!latex \een

!latex \subroutine{bc00aa}{interpolation of toroidal surfaces; construction;}

!latex \bi

!latex \item[1.] Given a discrete set of (closed) toroidal surfaces, which can be described by a set of Fourier harmonics for $R$ and $Z$ in suitable angles,
!latex           a continuous coordinate framework that is consistent with that described in \Sec{coordinatetransformation}
!latex           can be constructed by a suitable interpolation.

!latex \item     Most of the description in \Sec{coordinatetransformation} applies. Some loose ends are:
!latex \bi
!latex \item     for $m \ne 0$, if \verb+Lrad+$=1$, then $X_{m,n}(0)= X_{m,n}(\r_1)                                       $, and similarly for $Y_{m,n}$ etc.; 
!latex \item     for $m \ne 0$, if \verb+Lrad+$>1$, then $X_{m,n}(0)=(X_{m,n}(\r_1) \r_2 - X_{m,n}(\r_2) \r_1)/(\r_2-\r_1)$, and similarly for $Y_{m,n}$ etc.;
!latex \item     for the cubic-spline interpolations the end-point derivatives are required, and these are approximated using
!latex           one-sided, first-order differences, i.e. 
!latex           $X_{m,n}^\prime(0) = ( X_{m,n}(\rho_1) - X_{m,n}(0) ) / ( \r_1 - \r_0)$, etc.
!latex \ei
!latex \item     The radial coordinate is the volume, which is calculated for each toroidal surface by the integral
!latex     \be V = \int_{{\cal V}} \; dv 
!latex           = \frac{1}{3}\int_{{\cal V}} \; \nabla \cdot {\bf x} \; dv 
!latex           = \frac{1}{3}\int_{{\cal S}} \; {\bf x} \cdot d{\bf s} 
!latex           = \frac{1}{3}        \int_{0}^{2\pi} \!\! d\t \int_{0}^{2\pi/N} \!\! d\z \;\;\; \left. {\bf x   } \cdot {\bf x_\t} \times {\bf x_\z} \right|^s
!latex     \ee
!latex       where we have used $\nabla \cdot {\bf x}=3$, and have assumed that the domain is periodic in the angles.
!latex \item Using ${\bf x}\cdot {\bf e}_\theta \times {\bf e}_\zeta  = R ( Z R_\t - R Z_\t )$, 
!latex       \be V & = & \frac{1}{3} \; \int_{0}^{2\pi}\!\!\!d\t \int_{0}^{2\pi/N}\!\!\!\!\! d\z \; R \left( Z R_{\t} - R Z_{\t}  \right)              
!latex             \nonumber \\
!latex             & = & \frac{1}{3}  \; \sum_i \sum_j \sum_k R_{e,i} \left(Z_{e,j} R_{o,k} - R_{e,j} Z_{o,k} \right) (+m_k) 
!latex             \int \!\!\!\! \int \!\! d\t d\z \; \cos\alpha_i \cos\alpha_j \cos\alpha_k \nonumber \\
!latex             & + & \frac{1}{3}  \; \sum_i \sum_j \sum_k R_{e,i} \left(Z_{o,j} R_{e,k} - R_{o,j} Z_{e,k} \right) (-m_k) 
!latex             \int \!\!\!\! \int \!\! d\t d\z \; \cos\alpha_i \sin\alpha_j \sin\alpha_k \nonumber \\
!latex             & + & \frac{1}{3}  \; \sum_i \sum_j \sum_k R_{o,i} \left(Z_{e,j} R_{e,k} - R_{e,j} Z_{e,k} \right) (-m_k) 
!latex             \int \!\!\!\! \int \!\! d\t d\z \; \sin\alpha_i \cos\alpha_j \sin\alpha_k \nonumber \\
!latex             & + & \frac{1}{3}  \; \sum_i \sum_j \sum_k R_{o,i} \left(Z_{o,j} R_{o,k} - R_{o,j} Z_{o,k} \right) (+m_k) 
!latex             \int \!\!\!\! \int \!\! d\t d\z \; \sin\alpha_i \sin\alpha_j \cos\alpha_k
!latex       \ee
!latex       where $\alpha_i \equiv m_i \t - n_i \z$.
!latex \item Triple angle expansions are used to simplify the trigonometric terms as follows:
!latex       \be \begin{array}{ccccccccccccccccccccccccccccccccccccccccccc}
!latex           4 \; \cos\alpha_i \cos\alpha_j \cos\alpha_k & \!\!\! = & \!\!\!
!latex       + \!\!\! & \cos(\alpha_i+\alpha_j+\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i+\alpha_j-\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i-\alpha_j+\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i-\alpha_j-\alpha_k) \\
!latex           4 \; \cos\alpha_i \sin\alpha_j \sin\alpha_k & \!\!\! = & \!\!\!
!latex       - \!\!\! & \cos(\alpha_i+\alpha_j+\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i+\alpha_j-\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i-\alpha_j+\alpha_k) & \!\!\! - \!\!\! & \cos(\alpha_i-\alpha_j-\alpha_k) \\
!latex           4 \; \sin\alpha_i \cos\alpha_j \sin\alpha_k & \!\!\! = & \!\!\!
!latex       - \!\!\! & \cos(\alpha_i+\alpha_j+\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i+\alpha_j-\alpha_k) & \!\!\! - \!\!\! & \cos(\alpha_i-\alpha_j+\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i-\alpha_j-\alpha_k) \\
!latex           4 \; \sin\alpha_i \sin\alpha_j \cos\alpha_k & \!\!\! = & \!\!\!
!latex       - \!\!\! & \cos(\alpha_i+\alpha_j+\alpha_k) & \!\!\! - \!\!\! & \cos(\alpha_i+\alpha_j-\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i-\alpha_j+\alpha_k) & \!\!\! + \!\!\! & \cos(\alpha_i-\alpha_j-\alpha_k)
!latex       \end{array}
!latex       \ee
!latex \item     The cubic-spline interpolations are performed using \verb+Oculus:nr00aa+ and \verb+Oculus:nr00ab+, which will be described elsewhere.

!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : coordinates, bc00aa+   \\ \\
!latex           \verb+type(coordinates)   :: rzmn+   \\ \\ in their source that calls \verb+bc00aa+,
!latex           where \verb+coordinates+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+rzmn+, is arbitrary.

!latex \item[3.] \underline{\bf Required inputs}

!latex \item[  ] \verb+rzmn%mn               : integer ;+
!latex \bi
!latex \item[i.] Number of Fourier harmonics used to describe the toroidal surfaces;
!latex \ei

!latex \item[  ] \verb+rzmn%im(1:mn)         : integer ;+
!latex \bi
!latex \item[i.] Poloidal mode numbers: cannot be negative;
!latex \ei

!latex \item[  ] \verb+rzmn%in(1:mn)         : integer ;+
!latex \bi
!latex \item[i.] Toroidal mode numbers: for \verb+im(j)=0+ must have \verb+in(j).ge.0+ ;
!latex \ei

!latex \item[  ] \verb+rzmn%Lrad             : integer ;+
!latex \bi
!latex \item[i.] Number of surfaces to be interpolated; must have \verb+Lrad.ge.1+ ;
!latex \ei

!latex \item[  ] \verb+rzmn%Rbc(0:Lrad,1:mn) : integer ;+
!latex \item[  ] \verb+rzmn%Rbs(0:Lrad,1:mn) : integer ;+
!latex \item[  ] \verb+rzmn%Zbc(0:Lrad,1:mn) : integer ;+
!latex \item[  ] \verb+rzmn%Zbs(0:Lrad,1:mn) : integer ;+
!latex \bi
!latex \item[ i.] Fourier harmonics of the toroidal surfaces; these arrays must be allocated and assigned before calling \verb+bc00aa+.
!latex \item[ii.] The surfaces must be provided in increasing volume, and they must be ``nested".
!latex \item[iii.] The coordinate transformation is given by \Eqn{cylindricaltoroidalvector}, with the $j$-th coordinate surface being
!latex \be R_j(\t,\z) & \equiv & \sum_{i=1}^{\type{mn}} \type{Rbc[}j,i\type{]} \cos( \type{im[}i\type{]} \t - \type{in[}i\type{]} \z ), \\
!latex     Z_j(\t,\z) & \equiv & \sum_{i=1}^{\type{mn}} \type{Zbs[}j,i\type{]} \sin( \type{im[}i\type{]} \t - \type{in[}i\type{]} \z ),
!latex \ee
!latex and similarly for the non-stellarator-symmetric terms.
!latex \item[iv.] Note that the $0$-th surface is the degenerate surface $\equiv$ the coordinate axis,
!latex            for which \type{Rbc[}$0,i$\type{]}$=0$ and \type{Zbs[}$0,i$\type{]}$=0$ for \type{im[}$i$\type{]}$\ne0$.
!latex \ei

!latex \item[  ] \verb+rzmn%mm               : integer ;+
!latex \bi
!latex \item[i.] Maximum regularization factor; e.g. \verb+mm=2+;
!latex \ei

!latex \item[  ] \verb+ifail                 : integer ;+
!latex \bi
!latex \item[i.] if \verb+ifail =  1+, quiet mode;
!latex \item[i.] if \verb+ifail =  0+, screen output is terse;
!latex \item[i.] if \verb+ifail = -1+, 
!latex \ei

!latex \item[4.] \underline{\bf Execution}

!latex \item[  ] \verb+call bc00aa( rzmn, ifail )+

!latex \item[5.] \underline{\bf Outputs}

!latex \item[i.] on output: 
!latex \bi
!latex \item[  ] \verb+rzmn%ss(0:Lrad)            : real ;+
!latex \bi
!latex \item[i.] Volume of each surface; will be used as radial coordinate; this is under construction: perhaps a factor of $2\pi$ is missing.
!latex \ei

!latex \item[  ] \verb+rzmn%Xbc(0:Lrad,-1:2,1:mn) : real ;+
!latex \item[  ] \verb+rzmn%Xbs(0:Lrad,-1:2,1:mn) : real ;+
!latex \item[  ] \verb+rzmn%Ybc(0:Lrad,-1:2,1:mn) : real ;+
!latex \item[  ] \verb+rzmn%Ybs(0:Lrad,-1:2,1:mn) : real ;+
!latex \bi
!latex \item[i.] interpolating cubic-splines consistent with \Eqn{RZinterpolation};
!latex \item[ii.] if these are allocated on input they will first be de-allocated, and then re-allocated.
!latex \item[iii.] note that the spline {\em extrapolation} past the outermost surface is particularly unreliable if
!latex             (i) there are more than \verb+Lrad+$=1$ surface to be interpolated; or
!latex             (ii) the extrapolation is approaching the separatrix, and recall that a separatrix bounds all confined plasmas.
!latex \ei

!latex \item[  ] \verb+ifail                      : integer ; +
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex     \item[] \verb+ifail=2+ : the $m=0$ harmonics are not zero at the origin;
!latex     \item[] \verb+ifail=3+ : volume is not an increasing function of the surface label;
!latex \ei
!latex \item[6.] Comments:
!latex \bi \item[i.] The subroutine \verb+bc00ab+ will, given the interpolation coefficients,
!latex     calculate the coordinate transformation and derivatives as required.
!latex \ei
!latex \ei
!latex \ei

subroutine bc00aa( rzmn, ifail )
  
  implicit none
  
  type(coordinates)     :: rzmn
  INTEGER               :: ifail
  
  INTEGER               :: ibc00aa, astat, mn, Lrad, pLrad, ll, ii, jj, kk, mi, mj, mk, ni, nj, nk, mm
  REAL                  :: Rei, Roi, Zei, Zoi, Rej, Roj, Zej, Zoj, Rek, Rok, Zek, Zok
  REAL                  :: AA, BB, CC, DD, volume, VV, oV, lrho, rm
  REAL, allocatable     :: swork(:,:), rwork(:)
  
  ibc00aa = ifail ; mn = rzmn%mn ; Lrad = rzmn%Lrad ; pLrad = Lrad+1 ! internal shorthand; 22 Sep 15;
  
  CHECKINPUT( bc00aa, mn.le.1                  , 9999 )
  CHECKINPUT( bc00aa, .not.allocated(rzmn%im ) , 9999 )
  CHECKINPUT( bc00aa, .not.allocated(rzmn%in ) , 9999 )
  
  do ii = 1, mn
   CHECKINPUT(bc00aa,  rzmn%im(ii).lt.0                        , 9999 )
   CHECKINPUT(bc00aa,  rzmn%im(ii).eq.0 .and. rzmn%in(ii).lt.0 , 9999 )
  enddo
  
  CHECKINPUT( bc00aa, Lrad .le.0               , 9999)
  CHECKINPUT( bc00aa, .not.allocated(rzmn%Rbc) , 9999)
  CHECKINPUT( bc00aa, .not.allocated(rzmn%Rbs) , 9999)
  CHECKINPUT( bc00aa, .not.allocated(rzmn%Zbc) , 9999)
  CHECKINPUT( bc00aa, .not.allocated(rzmn%Zbs) , 9999)

  CHECKINPUT( bc00aa, rzmn%mm.lt.1             , 9999 )
  
  do ii = 1, mn ; mi = rzmn%im(ii)
   if( mi.eq.0 ) cycle ! need only check the m.ne.0 harmonics at the origin; 23 Nov 15;
   if( abs(rzmn%Rbc(0,ii)).gt.small ) then
     write(0,'("bc00aa : input error ; Rbc angle dependence at origin ;")') ; ibc00aa = 2 ; goto 9999
   endif
   if( abs(rzmn%Rbs(0,ii)).gt.small ) then
     write(0,'("bc00aa : input error ; Rbs angle dependence at origin ;")') ; ibc00aa = 2 ; goto 9999
   endif
   if( abs(rzmn%Zbc(0,ii)).gt.small ) then
     write(0,'("bc00aa : input error ; Zbc angle dependence at origin ;")') ; ibc00aa = 2 ; goto 9999
   endif
   if( abs(rzmn%Zbs(0,ii)).gt.small ) then
     write(0,'("bc00aa : input error ; Zbs angle dependence at origin ;")') ; ibc00aa = 2 ; goto 9999
   endif
  enddo ! end of do ii; 23 Nov 15;
  
  if( .not.allocated(rzmn%ss) ) allocate(rzmn%ss(0:Lrad))
  
  rzmn%ss(0:Lrad) = zero
  
  do ll = 0, Lrad
   
   volume = zero
   
   do ii = 1, mn ; mi = rzmn%im(ii) ; ni = rzmn%in(ii)
    
    Rei = rzmn%Rbc(ll,ii) ; Roi = rzmn%Rbs(ll,ii) ; Zei = rzmn%Zbc(ll,ii) ; Zoi = rzmn%Zbs(ll,ii)
    
    do jj = 1, mn ; mj = rzmn%im(jj) ; nj = rzmn%in(jj)
     
     Rej = rzmn%Rbc(ll,jj) ; Roj = rzmn%Rbs(ll,jj) ; Zej = rzmn%Zbc(ll,jj) ; Zoj = rzmn%Zbs(ll,jj)
     
     do kk = 1, mn ; mk = rzmn%im(kk) ; nk = rzmn%in(kk) 
      
      Rek = rzmn%Rbc(ll,kk) ; Rok = rzmn%Rbs(ll,kk) ; Zek = rzmn%Zbc(ll,kk) ; Zok = rzmn%Zbs(ll,kk)
      
      AA = Rei * ( Zej * Rok - Rej * Zok ) * (+mk)
      BB = Rei * ( Zoj * Rek - Roj * Zek ) * (-mk)
      CC = Roi * ( Zej * Rek - Rej * Zek ) * (-mk)
      DD = Roi * ( Zoj * Rok - Roj * Zok ) * (+mk)
      
      if( mi+mj+mk.eq.0 .and. ni+nj+nk.eq.0 ) volume = volume + AA - BB - CC - DD
      if( mi+mj-mk.eq.0 .and. ni+nj-nk.eq.0 ) volume = volume + AA + BB + CC - DD
      if( mi-mj+mk.eq.0 .and. ni-nj+nk.eq.0 ) volume = volume + AA + BB - CC + DD
      if( mi-mj-mk.eq.0 .and. ni-nj-nk.eq.0 ) volume = volume + AA - BB + CC + DD
      
     enddo ! end of do kk; 13 Sep 13;
    enddo ! end of do jj; 13 Sep 13;
   enddo ! end of do ii; 13 Sep 13;
   
   rzmn%ss(ll) = abs(volume) * pi2 * pi2 * quart * third
   
  enddo ! end of do ll ;  8 Jul 15;
  
 !rzmn%ss(0:Lrad) = rzmn%ss(0:Lrad) / rzmn%ss(Lrad) ! normalized volume;  8 Jul 15;
  
  do ll = 1, Lrad
   if( rzmn%ss(ll).lt.rzmn%ss(ll-1) ) then
     write(0,'("bc00aa : input error ; volume not increasing ;")') ; ibc00aa = 3 ; goto 9999
   endif
  enddo

  VV = rzmn%ss(Lrad) ; oV = one / VV
  
  if( allocated(rzmn%Xbc) ) deallocate(rzmn%Xbc)
  if( allocated(rzmn%Xbs) ) deallocate(rzmn%Xbs)
  if( allocated(rzmn%Ybc) ) deallocate(rzmn%Ybc)
  if( allocated(rzmn%Ybs) ) deallocate(rzmn%Ybs)
  
  SALLOCATE(rzmn%Xbc,(0:Lrad,-1:2,1:mn),zero)
  SALLOCATE(rzmn%Ybs,(0:Lrad,-1:2,1:mn),zero)
  SALLOCATE(rzmn%Xbs,(0:Lrad,-1:2,1:mn),zero)
  SALLOCATE(rzmn%Ybc,(0:Lrad,-1:2,1:mn),zero)
  
  SALLOCATE(swork,(0:Lrad,0:3),zero)
  SALLOCATE(rwork,(0:Lrad    ),zero)

  do ii = 1, mn ; mi = rzmn%im(ii) ; mm = min(mi,rzmn%mm)
   
   do ll = 1, Lrad ; lrho = rzmn%ss(ll) * oV ; rm = lrho**mm ! regularization factor is normalized;  8 Jul 15;
    
    rzmn%Xbc(ll,0,ii) = ( rzmn%Rbc(ll,ii) - rzmn%Rbc(0,ii) ) / rm
    rzmn%Xbs(ll,0,ii) = ( rzmn%Rbs(ll,ii) - rzmn%Rbs(0,ii) ) / rm
    rzmn%Ybc(ll,0,ii) = ( rzmn%Zbc(ll,ii) - rzmn%Zbc(0,ii) ) / rm
    rzmn%Ybs(ll,0,ii) = ( rzmn%Zbs(ll,ii) - rzmn%Zbs(0,ii) ) / rm
    
   enddo ! end of do ll; 11 Sep 15;
   
   select case( mi )
   case( 0 )
    rzmn%Xbc(0,0,ii) = zero
    rzmn%Xbs(0,0,ii) = zero
    rzmn%Ybc(0,0,ii) = zero
    rzmn%Ybs(0,0,ii) = zero
   case default
    select case( Lrad )
    case( 1  )
    rzmn%Xbc(0,0,ii) =   rzmn%Xbc(1,0,ii) !  zero-th order extrapolation; 22 Sep 15;
    rzmn%Xbs(0,0,ii) =   rzmn%Xbs(1,0,ii)
    rzmn%Ybc(0,0,ii) =   rzmn%Ybc(1,0,ii)
    rzmn%Ybs(0,0,ii) =   rzmn%Ybs(1,0,ii)
    case( 2: )
    ! first order extrapolation; 22 Sep 15;
    rzmn%Xbc(0,0,ii) = ( rzmn%Xbc(1,0,ii) * rzmn%ss(2) - rzmn%Xbc(2,0,ii) * rzmn%ss(1) ) / ( rzmn%ss(2) - rzmn%ss(1) )
    rzmn%Xbs(0,0,ii) = ( rzmn%Xbs(1,0,ii) * rzmn%ss(2) - rzmn%Xbs(2,0,ii) * rzmn%ss(1) ) / ( rzmn%ss(2) - rzmn%ss(1) )
    rzmn%Ybc(0,0,ii) = ( rzmn%Ybc(1,0,ii) * rzmn%ss(2) - rzmn%Ybc(2,0,ii) * rzmn%ss(1) ) / ( rzmn%ss(2) - rzmn%ss(1) )
    rzmn%Ybs(0,0,ii) = ( rzmn%Ybs(1,0,ii) * rzmn%ss(2) - rzmn%Ybs(2,0,ii) * rzmn%ss(1) ) / ( rzmn%ss(2) - rzmn%ss(1) )
    end select
   end select
   
   ! derivative at origin  ; 11 Sep 15;
   rzmn%Xbc(   0,1,ii) = ( rzmn%Xbc(   1,0,ii) - rzmn%Xbc(     0,0,ii) ) / ( rzmn%ss(   1) - rzmn%ss(     0) )
   rzmn%Xbs(   0,1,ii) = ( rzmn%Xbs(   1,0,ii) - rzmn%Xbs(     0,0,ii) ) / ( rzmn%ss(   1) - rzmn%ss(     0) )
   rzmn%Ybc(   0,1,ii) = ( rzmn%Ybc(   1,0,ii) - rzmn%Ybc(     0,0,ii) ) / ( rzmn%ss(   1) - rzmn%ss(     0) )
   rzmn%Ybs(   0,1,ii) = ( rzmn%Ybs(   1,0,ii) - rzmn%Ybs(     0,0,ii) ) / ( rzmn%ss(   1) - rzmn%ss(     0) )
   
   ! derivative at boundary; 11 Sep 15;
   rzmn%Xbc(Lrad,1,ii) = ( rzmn%Xbc(Lrad,0,ii) - rzmn%Xbc(Lrad-1,0,ii) ) / ( rzmn%ss(Lrad) - rzmn%ss(Lrad-1) ) 
   rzmn%Xbs(Lrad,1,ii) = ( rzmn%Xbs(Lrad,0,ii) - rzmn%Xbs(Lrad-1,0,ii) ) / ( rzmn%ss(Lrad) - rzmn%ss(Lrad-1) )
   rzmn%Ybc(Lrad,1,ii) = ( rzmn%Ybc(Lrad,0,ii) - rzmn%Ybc(Lrad-1,0,ii) ) / ( rzmn%ss(Lrad) - rzmn%ss(Lrad-1) )
   rzmn%Ybs(Lrad,1,ii) = ( rzmn%Ybs(Lrad,0,ii) - rzmn%Ybs(Lrad-1,0,ii) ) / ( rzmn%ss(Lrad) - rzmn%ss(Lrad-1) )
   
   ! spline interpolation; 11 Sep 15;
   call nr00aa( pLrad, rzmn%ss(0:Lrad), rzmn%Xbc(0:Lrad,-1:2,ii), swork(0:Lrad,0:3), rwork(0:Lrad) )
   call nr00aa( pLrad, rzmn%ss(0:Lrad), rzmn%Xbs(0:Lrad,-1:2,ii), swork(0:Lrad,0:3), rwork(0:Lrad) )
   call nr00aa( pLrad, rzmn%ss(0:Lrad), rzmn%Ybc(0:Lrad,-1:2,ii), swork(0:Lrad,0:3), rwork(0:Lrad) )
   call nr00aa( pLrad, rzmn%ss(0:Lrad), rzmn%Ybs(0:Lrad,-1:2,ii), swork(0:Lrad,0:3), rwork(0:Lrad) )
   
  enddo ! end of do ii; 16 Oct 14;
  
  DALLOCATE(rwork)
  DALLOCATE(swork)
  
9998 continue
  
  ibc00aa = 0
  
9999 continue
  
  if( ifail.le. 0 ) write(0,9000) ibc00aa, rzmn%mm, rzmn%ss(0:Lrad)
  
  ifail = ibc00aa
  
9000 format("bc00aa : ifail ="i3" : mm="i3" ; ss="999(f09.04","))
  
  return
  
end subroutine bc00aa

!latex \subroutine{bc00ab}{interpolation of toroidal surfaces; evaluation;}

!latex \bi
!latex \item[1.] This routine evaluates the coordinate transformation defined by \verb+bc00aa+ and provides the metric elements, Jacobian etc.
!latex           It can be called after \verb+bc00aa+.
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : coordinates, bc00ab+   \\ \\
!latex           \verb+type(coordinates)   :: rzmn+   \\ \\ in their source that calls \verb+bc00ab+,
!latex           where \verb+coordinates+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+rzmn+, is arbitrary.
!latex \ei

subroutine bc00ab( rzmn, stz, Lderiv, dRpZ, sqrtg, ifail )
  
  implicit none
  
  type(coordinates)     :: rzmn
  INTEGER               :: Lderiv, ifail
  REAL                  :: stz(1:3), dRpZ(1:3,0:3,0:3), sqrtg
  
  INTEGER               :: ibc00ab, Lrad, pLrad, ii, mi, ni, mm
  REAL                  :: lvol, teta, zeta, VV, oV, lrho, rm, rn, ro, rc, rs, zc, zs
  REAL                  :: drc, drs, dzc, dzs, ddrc, ddrs, ddzc, ddzs
  REAL                  :: arg, carg, sarg, Xbc(-1:2), Xbs(-1:2), Ybc(-1:2), Ybs(-1:2)
  
  dRpZ(1:3,0:3,0:3) = zero ; sqrtg = one ! default intent-out; 22 Sep 15;
  
  Lrad = rzmn%Lrad ; pLrad = Lrad+1 ! internal shorthand; 22 Sep 15;
  
  lvol =         stz(1)
  teta = modulo( stz(2),  pi2  )
  zeta = modulo( stz(3),  pi2  ) 
  
  CHECKINPUT(bc00ab, lvol.lt.zero            , 9999)
  CHECKINPUT(bc00ab, .not.allocated(rzmn%ss) , 9999)
  CHECKINPUT(bc00ab, Lrad.lt.1               , 9999)
  CHECKINPUT(bc00ab, rzmn%mm.lt. 1           , 9999)

  CHECKINPUT(bc00ab, .not.allocated(rzmn%Rbc) , 9999)
  CHECKINPUT(bc00ab, .not.allocated(rzmn%Zbc) , 9999)
  CHECKINPUT(bc00ab, .not.allocated(rzmn%Rbs) , 9999)
  CHECKINPUT(bc00ab, .not.allocated(rzmn%Zbs) , 9999)
  
  VV = rzmn%ss(Lrad) ! radial normalization; outermost volume; 22 Sep 15;

  CHECKINPUT(bc00ab, VV.lt.small             , 9999)
  
  oV = one / VV ; lrho = lvol * oV ! regularization factor is normalized volume; 11 Sep 15;
  
  dRpZ(2,0,0) = zeta ; dRpZ(2,1,0) = zero ; dRpZ(2,2,0) = zero ; dRpZ(2,3,0) =  one ! toroidal angle and derivatives; 22 Sep 15;

! Fourier summation loop; 22 Sep 15;
  
  do ii = 1, rzmn%mn

   mi = rzmn%im(ii) ; ni = rzmn%in(ii) ; arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)
   
   mm = min(mi,rzmn%mm) ! regularization exponent; 22 Sep 15;
   
   call nr00ab( pLrad, rzmn%ss(0:Lrad), rzmn%Xbc(0:Lrad,-1:2,ii), lvol, Xbc(-1:2) ) ! evaluate splines; 22 Sep 15;
   call nr00ab( pLrad, rzmn%ss(0:Lrad), rzmn%Xbs(0:Lrad,-1:2,ii), lvol, Xbs(-1:2) )
   call nr00ab( pLrad, rzmn%ss(0:Lrad), rzmn%Ybc(0:Lrad,-1:2,ii), lvol, Ybc(-1:2) )
   call nr00ab( pLrad, rzmn%ss(0:Lrad), rzmn%Ybs(0:Lrad,-1:2,ii), lvol, Ybs(-1:2) )
   
   select case( mi )
   case( 0   ) ;  rc = rzmn%Rbc(0,ii) +      Xbc(0)
    ;          ;  rs = rzmn%Rbs(0,ii) +      Xbs(0)
    ;          ;  zc = rzmn%Zbc(0,ii) +      Ybc(0)
    ;          ;  zs = rzmn%Zbs(0,ii) +      Ybs(0)
   case( 1:  ) ;  rm = lrho**mm
    ;          ;  rc =                  rm * Xbc(0)
    ;          ;  rs =                  rm * Xbs(0)
    ;          ;  zc =                  rm * Ybc(0)
    ;          ;  zs =                  rm * Ybs(0)
   end select

! can save Fourier harmonics into output; 22 Sep 15;
   
   dRpZ(1,0,0) = dRpZ(1,0,0) + ( + rc * carg + rs * sarg )
!  dRpZ(2,0,0) = zeta
   dRpZ(3,0,0) = dRpZ(3,0,0) + ( + zc * carg + zs * sarg )
   
   if( Lderiv.lt.1 ) cycle      
   
   select case( mi )
   case( 0   ) ; drc =      Xbc(1) ! 22 Sep 15;
    ;          ; drs =      Xbs(1)
    ;          ; dzc =      Ybc(1)
    ;          ; dzs =      Ybs(1)
   case( 1   ) ; drc = rm * Xbc(1) +           oV * Xbc(0) ! 22 Sep 15;
    ;          ; drs = rm * Xbs(1) +           oV * Xbs(0)
    ;          ; dzc = rm * Ybc(1) +           oV * Ybc(0)
    ;          ; dzs = rm * Ybs(1) +           oV * Ybs(0)
   case( 2:  ) ; rn = lrho**(mm-1)
    ;          ; drc = rm * Xbc(1) + mm * rn * oV * Xbc(0) ! 22 Sep 15;
    ;          ; drs = rm * Xbs(1) + mm * rn * oV * Xbs(0)
    ;          ; dzc = rm * Ybc(1) + mm * rn * oV * Ybc(0)
    ;          ; dzs = rm * Ybs(1) + mm * rn * oV * Ybs(0)
   end select                                                                                              

   dRpZ(1,1,0) = dRpZ(1,1,0) + ( + drc * carg + drs * sarg )
   dRpZ(1,2,0) = dRpZ(1,2,0) + ( -  rc * sarg +  rs * carg ) * (+mi)
   dRpZ(1,3,0) = dRpZ(1,3,0) + ( -  rc * sarg +  rs * carg ) * (-ni)
  !dRpZ(2,1,0) = zero
  !dRpZ(2,2,0) = zero
  !dRpZ(2,3,0) =  one
   dRpZ(3,1,0) = dRpZ(3,1,0) + ( + dzc * carg + dzs * sarg )
   dRpZ(3,2,0) = dRpZ(3,2,0) + ( -  zc * sarg +  zs * carg ) * (+mi)
   dRpZ(3,3,0) = dRpZ(3,3,0) + ( -  zc * sarg +  zs * carg ) * (-ni)

   if( Lderiv.lt.2 ) cycle      
   
   select case( mi )
   case( 0   ) ; ddrc =        Xbc(2) ! 22 Sep 15;
    ;          ; ddrs =        Xbs(2)
    ;          ; ddzc =        Ybc(2)
    ;          ; ddzs =        Ybs(2)
   case( 1   ) ; ddrc =   rm * Xbc(2) + 2           * oV * Xbc(1) ! 22 Sep 15;
    ;          ; ddrs =   rm * Xbs(2) + 2           * oV * Xbs(1)
    ;          ; ddzc =   rm * Ybc(2) + 2           * oV * Ybc(1)
    ;          ; ddzs =   rm * Ybs(2) + 2           * oV * Ybs(1)
   case( 2   ) ; ddrc =   rm * Xbc(2) + 2 * mm * rn * oV * Xbc(1) + mm *                 oV * oV * Xbc(0) ! 22 Sep 15;
    ;          ; ddrs =   rm * Xbs(2) + 2 * mm * rn * oV * Xbs(1) + mm *                 oV * oV * Xbs(0)
    ;          ; ddzc =   rm * Ybc(2) + 2 * mm * rn * oV * Ybc(1) + mm *                 oV * oV * Ybc(0)
    ;          ; ddzs =   rm * Ybs(2) + 2 * mm * rn * oV * Ybs(1) + mm *                 oV * oV * Ybs(0)
   case( 3:  ) ; ro = lrho**(mm-2)
    ;          ; ddrc =   rm * Xbc(2) + 2 * mm * rn * oV * Xbc(1) + mm * ( mm-1 ) * ro * oV * oV * Xbc(0) ! 22 Sep 15;
    ;          ; ddrs =   rm * Xbs(2) + 2 * mm * rn * oV * Xbs(1) + mm * ( mm-1 ) * ro * oV * oV * Xbs(0)
    ;          ; ddzc =   rm * Ybc(2) + 2 * mm * rn * oV * Ybc(1) + mm * ( mm-1 ) * ro * oV * oV * Ybc(0)
    ;          ; ddzs =   rm * Ybs(2) + 2 * mm * rn * oV * Ybs(1) + mm * ( mm-1 ) * ro * oV * oV * Ybs(0)
   end select                                                                                              

   dRpZ(1,1,1) = dRpZ(1,1,1) + ( + ddrc * carg + ddrs * sarg )
   dRpZ(1,1,2) = dRpZ(1,1,2) + ( -  drc * sarg +  drs * carg )         * (+mi)
   dRpZ(1,1,3) = dRpZ(1,1,3) + ( -  drc * sarg +  drs * carg )         * (-ni)
   dRpZ(1,2,2) = dRpZ(1,2,2) + ( -   rc * carg -   rs * sarg ) * (+mi) * (+mi)
   dRpZ(1,2,3) = dRpZ(1,2,3) + ( -   rc * carg -   rs * sarg ) * (+mi) * (-ni)
   dRpZ(1,3,3) = dRpZ(1,3,3) + ( -   rc * carg -   rs * sarg ) * (-ni) * (-ni)
  !dRpZ(2,1,1) = zero
  !dRpZ(2,1,2) = zero
  !dRpZ(2,1,3) = zero
  !dRpZ(2,2,2) = zero
  !dRpZ(2,2,3) = zero
  !dRpZ(2,3,3) = zero
   dRpZ(3,1,1) = dRpZ(3,1,1) + ( + ddzc * carg + ddzs * sarg )
   dRpZ(3,1,2) = dRpZ(3,1,2) + ( -  dzc * sarg +  dzs * carg )         * (+mi)
   dRpZ(3,1,3) = dRpZ(3,1,3) + ( -  dzc * sarg +  dzs * carg )         * (-ni)
   dRpZ(3,2,2) = dRpZ(3,2,2) + ( -   zc * carg -   zs * sarg ) * (+mi) * (+mi)
   dRpZ(3,2,3) = dRpZ(3,2,3) + ( -   zc * carg -   zs * sarg ) * (+mi) * (-ni)
   dRpZ(3,3,3) = dRpZ(3,3,3) + ( -   zc * carg -   zs * sarg ) * (-ni) * (-ni)

  enddo ! end of do ii; 16 Oct 14;
  
! Jacobian;
  
  sqrtg = dRpZ(1,0,0) * ( dRpZ(1,2,0) * dRpZ(3,1,0) - dRpZ(1,1,0) * dRpZ(3,2,0) ) ! Jacobian of toroidal coordinates; 02 Jul 13;
  
! complete derivative array; 11 Sep 15;

  dRpZ(1,0,1) = dRpZ(1,1,0) ; dRpZ(2,0,1) = dRpZ(2,1,0) ; dRpZ(3,0,1) = dRpZ(3,1,0)
  dRpZ(1,0,2) = dRpZ(1,2,0) ; dRpZ(2,0,2) = dRpZ(2,2,0) ; dRpZ(3,0,2) = dRpZ(3,2,0)
  dRpZ(1,0,3) = dRpZ(1,3,0) ; dRpZ(2,0,3) = dRpZ(2,3,0) ; dRpZ(3,0,3) = dRpZ(3,3,0)
 !dRpZ(1,1,1) = dRpZ(1,1,1) ; dRpZ(2,1,1) = dRpZ(2,1,1) ; dRpZ(3,1,1) = dRpZ(3,1,1)
  dRpZ(1,2,1) = dRpZ(1,1,2) ; dRpZ(2,2,1) = dRpZ(2,1,2) ; dRpZ(3,2,1) = dRpZ(3,1,2)
  dRpZ(1,3,1) = dRpZ(1,1,3) ; dRpZ(2,3,1) = dRpZ(2,1,3) ; dRpZ(3,3,1) = dRpZ(3,1,3)
 !dRpZ(1,2,2) = dRpZ(1,2,2) ; dRpZ(2,2,2) = dRpZ(2,2,2) ; dRpZ(3,2,2) = dRpZ(3,2,2)
  dRpZ(1,3,2) = dRpZ(1,2,3) ; dRpZ(2,3,2) = dRpZ(2,2,3) ; dRpZ(3,3,2) = dRpZ(3,2,3)
 !dRpZ(1,3,3) = dRpZ(1,3,3) ; dRpZ(2,3,3) = dRpZ(2,3,3) ; dRpZ(3,3,3) = dRpZ(3,3,3)
  
9998 continue
  
  ifail = 0
  
9999 continue
  
  return
  
end subroutine bc00ab

!latex \newpage \section{``toroidal'' subroutines}

!latex In this section are described subroutines that depend on the toroidal coordinates described in \Sec{coordinatetransformation}.

!latex \subroutine{aa00aa}{construct vector potential in toroidal coordinates;}
!latex \ben
!latex \item The vector potential, ${\bf B}=\nabla\times{\bf A}$, may be constructed in toroidal coordinates, $(\r,\t,\z)$, by radial integration.
!latex \item The toroidal coordinates must be provided in a suitable format, and it is suggested that \verb+bc00aa+ be called prior to calling \verb+aa00aa+.
!latex       The magnetic field, ${\bf B}=B^R {\bf e}_R + B^\phi {\bf e}_\phi + B^Z{\bf e}_Z$, is assumed to be given in cylindrical coordinates;
!latex       and the coordinate transformation, ${\bf x}(\r,\t,\z)$, as given in \Eqn{cylindricaltoroidalcoordinates} will be assumed.

!l!tex \item Assuming toroidal symmetry, and introducing `local' Cartesian coordinates
!l!tex       $x\equiv R(\r,\t,\z)-R(0,\t,\z)$ and $y\equiv Z(\r,\t,\z)-Z(0,\t,\z)$, 
!l!tex       \be \begin{array}{ccc} x & = & \alpha \; \r^{1/2} \; \cos\t, \\
!l!tex                              y & = & \beta  \; \r^{1/2} \; \sin\t, \end{array} \label{eq:localCartesian}
!l!tex       \ee
!l!tex       for some constants $\alpha\equiv x_{1,0}$ and $\beta\equiv y_{1,0}$.
!l!tex \item The cylindrical field is assumed to be regular at the axis of the toroidal coordinates,
!l!tex       i.e. that $B^R$, $B^\phi$ and $B^Z$ may be described using a Taylor expansion in locally-defined Cartesian variables.
!l!tex \item The Taylor expansion for a regular function, $f$, is
!l!tex       \be f & = & \sum_{i,j} a_{i,j} x^i j^y.
!l!tex       \ee
!l!tex \item Substituting \Eqn{localCartesian} and repeated use of double-angle formulae, we derive
!l!tex       \be f(\r,\t) & = & \sum_{m} \r^{m/2} f_m(\r) \cos(m\t),
!l!tex       \ee 
!l!tex       where $f_m(\r)$ is a power series in $\r$.
!l!tex \item Assuming $B^R$, $B^\phi$ and $B^Z$ are each regular, the vector transformation, \Eqn{toroidalcylindricalvector}, 
!l!tex       implies that the Fourier harmonics of $\sqrt g B^\r$, $\sqrt g B^\t$ and $\sqrt g B^\z$ obey
!l!tex       \be (\sqrt g B^\r)_m(\r) &=& \left\{ \begin{array}{lcl}  \r^{    1} & f_m(\r), & \mbox{\rm for } m = 0, \\ 
!l!tex                                                                \r^{m/2  } & f_m(\r), & \mbox{\rm for } m > 0, \end{array} \right. \\
!l!tex           (\sqrt g B^\t)_m(\r) &=& \left\{ \begin{array}{lcl}             & g_m(\r), & \mbox{\rm for } m = 0, \\ 
!l!tex                                                                \r^{m/2-1} & g_m(\r), & \mbox{\rm for } m > 0, \end{array} \right. \\
!l!tex           (\sqrt g B^\z)_m(\r) &=& \left\{ \begin{array}{lcl}             & h_m(\r), & \mbox{\rm for } m = 0, \\ 
!l!tex                                                                \r^{m/2\:} & h_m(\r), & \mbox{\rm for } m > 0, \end{array} \right.
!l!tex       \ee
!l!tex       where $f_m(\r)$, $g_m(\r)$ and $h_m(\r)$ are arbitrary polynomials in $\r$.

!latex \item A general representation for the vector potential is ${\bf A} = \bar A_\r \nabla \r + \bar A_\t \nabla \t + \bar A_z \nabla \z$.
!latex       A gauge function, $g(\r,\t,\z)$, may be chosen to simplify this:
!latex       choosing $\partial_\r g = -\bar A_\r$, the $\bar A_\r$ component is cancelled leaving
!latex       ${\bf A} =                            A_\t \nabla \t +      A_\z \nabla \z$.
!latex       The equation ${\bf B}=\nabla\times{\bf A}$ becomes
!latex       \be \begin{array}{lllll} \partial_\t A_\z &-& \partial_\z A_\t &=& \sqrt g B^\r, \\
!latex                                                 &-& \partial_\r A_\z &=& \sqrt g B^\t, \\
!latex                                \partial_\r A_\t & &                  &=& \sqrt g B^\z. \end{array} \label{eq:BiscurlA}
!latex       \ee
!latex \item The Fourier harmonics of the components of the vector potential can be determined by radial integration:
!latex       \be A_{\t,m,n}(\r) & = & A_{\t,m,n}(0) + \int_{0}^{\r} \!\! (\sqrt g B^\z)_{m,n}(\bar \r) \; d\bar\r, \\
!latex           A_{\z,m,n}(\r) & = & A_{\z,m,n}(0) - \int_{0}^{\r} \!\! (\sqrt g B^\t)_{m,n}(\bar \r) \; d\bar\r. 
!latex       \ee 
!latex \item Solving these two equations automatically satifies the second and third equations in \Eqn{BiscurlA}.
!latex       The first equation in \Eqn{BiscurlA} is satisfied if and only if (i) $\nabla\cdot{\bf B}=0$, and (ii) the `integration constants'
!latex       $A_{\t}(0,\t,\z)$ and $A_{\z}(0,\t,\z)$ satisfy $\partial_\t A_\z(0,\t,\z) - \partial_\z A_\t(0,\t,\z) = (\sqrt g B^\r)(0,\t,\z)$.
!l!tex \item The Fourier harmonics of $A_\t$ and $A_\z$ are
!l!tex       \be A_{\t,m}(\r) &=& \left\{ \begin{array}{lcl}  \r^{    1} & a_m(\r), & \mbox{\rm for } m = 0, \\ 
!l!tex                                                        \r^{m/2+1} & a_m(\r), & \mbox{\rm for } m > 0, \end{array} \right. \\
!l!tex           A_{\z,m}(\r) &=& \left\{ \begin{array}{lcl}  \r^{    1} & b_m(\r), & \mbox{\rm for } m = 0, \\ 
!l!tex                                                        \r^{m/2  } & b_m(\r), & \mbox{\rm for } m > 0, \end{array} \right. \\
!l!tex       \ee
!latex       Note that $(\sqrt g B^\r)(0,\t,\z)=0$,
!latex       and this allows the remaining gauge freedom, namely $g(0,\t,\z)$, to be used 
!latex       to set $A_\t(0,\t,\z)=0$ and $A_\z(0,\t,\z)=0$.
!latex \item The `background' toroidal coordinates should usually be provided by \verb+bc00aa+,
!latex       and the domain of integration is $\r\in[0,V]$, where $V\equiv \max(\r)$.
!latex \item The coordinate interpolation is evaluated using \verb+nr00ab+,
!latex       the spline construction to the Fourier harmonics of the vector potential are constructed using \verb+nr00aa+, 
!latex       and the required FFTs are evaluated using \verb+ft02aa+.
!latex \item After calling \verb+aa00aa+, the routine \verb+aa00ab+ may be called
!latex       to evaluate the vector potential at an arbitrary point within the computational domain.
!latex \item     \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : coordinates, vectorpotential, aa00aa+ \\ \\
!latex           \verb+type(coordinates)     :: rzmn+ \\
!latex           \verb+type(vectorpotential) :: Atzmn+ \\ \\ in their source that calls \verb+aa00aa+,
!latex           where \verb+coordinates+ and \verb+vectorpotential+ are derived types (i.e. structures)
!latex           that contains both the required input and output information.
!latex           The variable names, \verb+rzmn+ and \verb+Atzmn+, are arbitrary.
!latex \item     \underline{\bf Required inputs}
!latex \item[  ] \verb+rzmn                : structure ;+
!latex \bi
!latex \item[i.] the background toroidal coordinates, $(\r,\t,\z)$, that will be used;
!latex \item[ii.] it is assumed that before calling \verb+aa00aa+, that the routine \verb+bc00aa+ has been called, and \verb+rzmn+ should be unchanged;
!latex \ei
!latex \item[  ] \verb+Atzmn%Nfp           : integer   ;+
!latex \bi
!latex \item[i.] the field periodicity, which must be the same as the field periodicity of the coordinates;
!latex \ei
!latex \item[  ] \verb+Atzmn%Lrad          : integer   ;+
!latex \bi
!latex \item[i.] the required radial resolution;
!latex \ei
!latex \item[  ] \verb+Atzmn%Mpol          : integer   ;+
!latex \item[  ] \verb+Atzmn%Ntor          : integer   ;+
!latex \bi
!latex \item[i.] the required Fourier resolution;
!latex \ei
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call aa00aa( rzmn, Atzmn, ifail )+
!latex \item     \underline{\bf Outputs}
!latex \item[  ] \verb+Atzmn%gBtc(0:Lrad,-1:2,1:amn)  : real      ;+
!latex \item[  ] \verb+Atzmn%gBts(0:Lrad,-1:2,1:amn)  : real      ;+
!latex \item[  ] \verb+Atzmn%gBzc(0:Lrad,-1:2,1:amn)  : real      ;+
!latex \item[  ] \verb+Atzmn%gBzs(0:Lrad,-1:2,1:amn)  : real      ;+
!latex \bi
!latex \item[i.] the Fourier harmonics of the vector potential;
!latex \ei
!latex \item[  ] \verb+Atzmn%gBstz(1:Ntz,1:3,0:Lrad)  : real      ;+
!latex \bi
!latex \item[i.] the vector potential on a regular grid in real space ;
!latex \ei
!latex \een
  
  subroutine aa00aa( lrzmn, lAtzmn, ifail )
    
    implicit none
    
    type(coordinates)     :: lrzmn
    type(vectorpotential) :: lAtzmn
    
    INTEGER               :: ifail
    
    INTEGER               :: astat, Nfp, jrho, ii, mi, ni, jj, kk, jk, ift02aa
    INTEGER               ::        aLrad, aMpol, aNtor, amn, aNt, aNz, aNtz, paLrad
    INTEGER               ::        cLrad,                                    pcLrad, mm
    
    REAL                  :: rm, rn, rc, rs, zc, zs, drc, drs, dzc, dzs
    REAL                  :: VV, oV, lvol, lrho, teta, zeta, arg, ca, sa

    REAL, allocatable     :: Xbc(:,:), Xbs(:,:), Ybc(:,:), Ybs(:,:), carg(:,:), sarg(:,:)

    INTEGER               :: lfail                    ! input/output from user-supplied bfield; 23 Nov 15;
    REAL                  :: RpZ(1:3), dBRpZ(1:3,0:3) ! input/output from user-supplied bfield; 23 Nov 15;

    REAL                  :: VM(1:3,1:3) ! vector transformation matrix; 23 Nov 15;

    REAL, allocatable     :: Roij(:,:), Zoij(:,:) ! cylindrical coordinates and derivatives; 23 Nov 15;

    REAL                  :: ds, surfaceerror(1:3)

    REAL                  :: stz(1:3), lBstz(1:3)
    
    REAL, allocatable     :: gBstz(:,:) ! reconstructed magnetic field from Fourier vector potential; 23 Nov 15;

    REAL, allocatable     :: even(:), oddd(:), comn(:), simn(:) ! dummy workspace for FFT; 23 Nov 15;

    REAL, allocatable     :: swork(:,:), rwork(:)
    
    iaa00aa = ifail
    
    if( allocated(lAtzmn%ss   ) ) deallocate(lAtzmn%ss   )
    if( allocated(lAtzmn%im   ) ) deallocate(lAtzmn%im   )
    if( allocated(lAtzmn%in   ) ) deallocate(lAtzmn%in   )
    if( allocated(lAtzmn%gBtc ) ) deallocate(lAtzmn%gBtc )
    if( allocated(lAtzmn%gBts ) ) deallocate(lAtzmn%gBts )
    if( allocated(lAtzmn%gBzc ) ) deallocate(lAtzmn%gBzc )
    if( allocated(lAtzmn%gBzs ) ) deallocate(lAtzmn%gBzs )
    if( allocated(lAtzmn%gBstz) ) deallocate(lAtzmn%gBstz)
    
    CHECKINPUT( aa00aa, lrzmn%mn  .le.1           , 9999 )
    CHECKINPUT( aa00aa, lrzmn%Lrad.le.0           , 9999 )
    
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%ss ) , 9999 )
    
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%Rbc) , 9999 )
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%Rbs) , 9999 )
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%Zbc) , 9999 )
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%Zbs) , 9999 )
    
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%Xbc) , 9999 )
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%Xbs) , 9999 )
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%Ybc) , 9999 )
    CHECKINPUT( aa00aa, .not.allocated(lrzmn%Ybs) , 9999 )
    
    CHECKINPUT( aa00aa, lrzmn%mm.lt.1             , 9999 )
    
    CHECKINPUT( aa00aa, lrzmn%ss(lrzmn%Lrad).le.small , 9999 )
    
    cLrad = lrzmn%Lrad ; pcLrad = cLrad + 1 ; VV = lrzmn%ss(cLrad) ; oV = one / VV
    
    CHECKINPUT( aa00aa, lAtzmn%Nfp .le.0 , 9999 )
    CHECKINPUT( aa00aa, lAtzmn%Lrad.le.0 , 9999 )
    CHECKINPUT( aa00aa, lAtzmn%Mpol.lt.1 , 9999 )
    CHECKINPUT( aa00aa, lAtzmn%Ntor.lt.0 , 9999 )
    
    CHECKINPUT( aa00aa, lAtzmn%Nfp .ne. lrzmn%Nfp, 9999 )
    
    aLrad = lAtzmn%Lrad ; paLrad = aLrad + 1 ; aMpol = lAtzmn%Mpol ; aNtor = lAtzmn%Ntor ; Nfp = lAtzmn%Nfp ! internal shorthand; 22 Sep 15;
    
    SALLOCATE( lAtzmn%ss,(0:aLrad),zero )
    
    lAtzmn%ss(0:aLrad) = (/ ( jrho, jrho = 0, aLrad ) /) * VV / aLrad ! radial grid for vector potential; 13 Oct 15;
    
    lAtzmn%Nt = 4 * aMpol ; lAtzmn%Nz = 4 * max( aNtor, 1 )

    aNtz = lAtzmn%Nt * lAtzmn%Nz ; amn = aNtor + 1 + aMpol * ( 2 * aNtor + 1 ) ; lAtzmn%Ntz = aNtz ; lAtzmn%mn = amn
        
    SALLOCATE(lAtzmn%im,(1:amn),0)
    SALLOCATE(lAtzmn%in,(1:amn),0)
    
    call gi00ab( aMpol, aNtor, Nfp, amn, lAtzmn%im(1:amn), lAtzmn%in(1:amn) )
    
    SALLOCATE(lAtzmn%gBtc,(0:aLrad,-1:2,1:amn),zero) ! Fourier harmonics of vector potential; 13 Oct 15;
    SALLOCATE(lAtzmn%gBts,(0:aLrad,-1:2,1:amn),zero)
    SALLOCATE(lAtzmn%gBzc,(0:aLrad,-1:2,1:amn),zero)
    SALLOCATE(lAtzmn%gBzs,(0:aLrad,-1:2,1:amn),zero)
    
    SALLOCATE(Roij,(1:aNtz,0:3),zero)
    SALLOCATE(Zoij,(1:aNtz,0:3),zero)
    
    SALLOCATE(lAtzmn%gBstz,(1:aNtz,1:3,0:aLrad),zero) ! contravariant components of toroidal field on regular grid; 22 Sep 15;
    SALLOCATE(       gBstz,(1:aNtz,1:3        ),zero) ! contravariant components of toroidal field on regular grid; 22 Sep 15;
    
    SALLOCATE(even,(1:amn),zero)
    SALLOCATE(oddd,(1:amn),zero)
    SALLOCATE(comn,(1:amn),zero)
    SALLOCATE(simn,(1:amn),zero)

    SALLOCATE(Xbc,(-1:2,1:lrzmn%mn),zero)
    SALLOCATE(Xbs,(-1:2,1:lrzmn%mn),zero)
    SALLOCATE(Ybc,(-1:2,1:lrzmn%mn),zero)
    SALLOCATE(Ybs,(-1:2,1:lrzmn%mn),zero)    

    SALLOCATE(carg,(1:aNtz,1:lrzmn%mn),zero)
    SALLOCATE(sarg,(1:aNtz,1:lrzmn%mn),zero)
    

    do kk = 0, lAtzmn%Nz-1 ;                            ; zeta = kk * (pi2/Nfp) / lAtzmn%Nz
     do jj = 0, lAtzmn%Nt-1 ; jk = 1 + jj + kk*lAtzmn%Nt ; teta = jj *  pi2      / lAtzmn%Nt
      do ii = 1, lrzmn%mn
       mi = lrzmn%im(ii) ; ni = lrzmn%in(ii)
       arg = mi * teta - ni * zeta ; carg(jk,ii) = cos(arg) ; sarg(jk,ii) = sin(arg)
      enddo
     enddo
    enddo
    
    if( ifail.lt.0 ) then 
      write(0,'("aa00aa : " 10x " : constructing vector potential : Lrad ="i5" ; Mpol ="i4" ; Ntor ="i4" ;")') &
        aLrad, aMpol, aNtor
    endif

    do jrho = 0, aLrad ! calculate magnetic field on regular grid; 19 Nov 14;
     
     lvol = lAtzmn%ss(jrho) ; lrho = lvol * oV ! regularization factor is normalized;  8 Jul 15;
     
     do ii = 1, lrzmn%mn
      call nr00ab( pcLrad, lrzmn%ss(0:cLrad), lrzmn%Xbc(0:cLrad,-1:2,ii), lvol, Xbc(-1:2,ii) ) ! interpolate Fourier harmonics of coordinates;
      call nr00ab( pcLrad, lrzmn%ss(0:cLrad), lrzmn%Xbs(0:cLrad,-1:2,ii), lvol, Xbs(-1:2,ii) )
      call nr00ab( pcLrad, lrzmn%ss(0:cLrad), lrzmn%Ybc(0:cLrad,-1:2,ii), lvol, Ybc(-1:2,ii) )
      call nr00ab( pcLrad, lrzmn%ss(0:cLrad), lrzmn%Ybs(0:cLrad,-1:2,ii), lvol, Ybs(-1:2,ii) )
     enddo
     
     do kk = 0, lAtzmn%Nz-1 ;                            ; zeta = kk * (pi2/Nfp) / lAtzmn%Nz ! zeta is     used to give position to bfield; 23 Nov 15;
      do jj = 0, lAtzmn%Nt-1 ; jk = 1 + jj + kk*lAtzmn%Nt ! teta = jj *  pi2      / lAtzmn%Nt ! teta is not used in this particular loop   ; 23 Nov 15;
       
       Roij(jk,0:3) = zero ! initialize summation; 23 Nov 15;
       Zoij(jk,0:3) = zero
       
       do ii = 1, lrzmn%mn ; mi = lrzmn%im(ii) ; ni = lrzmn%in(ii) ; mm = min(mi,lrzmn%mm)
        
        select case( mi )
        case( 0   ) ;  rc = lrzmn%Rbc(0,ii) +      Xbc(0,ii)
         ;          ;  rs = lrzmn%Rbs(0,ii) +      Xbs(0,ii)
         ;          ;  zc = lrzmn%Zbc(0,ii) +      Ybc(0,ii)
         ;          ;  zs = lrzmn%Zbs(0,ii) +      Ybs(0,ii)
        case( 1:  ) ;  rm = lrho**mm
         ;          ;  rc =                   rm * Xbc(0,ii)
         ;          ;  rs =                   rm * Xbs(0,ii)
         ;          ;  zc =                   rm * Ybc(0,ii)
         ;          ;  zs =                   rm * Ybs(0,ii)
        end select
        
        Roij(jk,0) = Roij(jk,0) + ( +  rc * carg(jk,ii) +  rs * sarg(jk,ii) )
        Zoij(jk,0) = Zoij(jk,0) + ( +  zc * carg(jk,ii) +  zs * sarg(jk,ii) )   
        
        select case( mi )
        case( 0   ) ; drc =      Xbc(1,ii)
         ;          ; drs =      Xbs(1,ii)
         ;          ; dzc =      Ybc(1,ii)
         ;          ; dzs =      Ybs(1,ii)
        case( 1   ) ; drc = rm * Xbc(1,ii) +           oV * Xbc(0,ii)
         ;          ; drs = rm * Xbs(1,ii) +           oV * Xbs(0,ii)
         ;          ; dzc = rm * Ybc(1,ii) +           oV * Ybc(0,ii)
         ;          ; dzs = rm * Ybs(1,ii) +           oV * Ybs(0,ii)
        case( 2:  ) ; rn = lrho**(mm-1)
         ;          ; drc = rm * Xbc(1,ii) + mm * rn * oV * Xbc(0,ii)
         ;          ; drs = rm * Xbs(1,ii) + mm * rn * oV * Xbs(0,ii)
         ;          ; dzc = rm * Ybc(1,ii) + mm * rn * oV * Ybc(0,ii)
         ;          ; dzs = rm * Ybs(1,ii) + mm * rn * oV * Ybs(0,ii)
        end select
        
        Roij(jk,1) = Roij(jk,1) + ( + drc * carg(jk,ii) + drs * sarg(jk,ii) )
        Roij(jk,2) = Roij(jk,2) + ( -  rc * sarg(jk,ii) +  rs * carg(jk,ii) ) * (+mi)
        Roij(jk,3) = Roij(jk,3) + ( -  rc * sarg(jk,ii) +  rs * carg(jk,ii) ) * (-ni)

        Zoij(jk,1) = Zoij(jk,1) + ( + dzc * carg(jk,ii) + dzs * sarg(jk,ii) )
        Zoij(jk,2) = Zoij(jk,2) + ( -  zc * sarg(jk,ii) +  zs * carg(jk,ii) ) * (+mi)
        Zoij(jk,3) = Zoij(jk,3) + ( -  zc * sarg(jk,ii) +  zs * carg(jk,ii) ) * (-ni)
        
       enddo ! end of do ii; 22 Sep 15;

       RpZ(1:3) = (/ Roij(jk,0), zeta, Zoij(jk,0) /) ; itangent = 0

       VM(1,1:3) = (/ - Zoij(jk,2), Roij(jk,3) * Zoij(jk,2) - Roij(jk,2) * Zoij(jk,3), + Roij(jk,2) /) ! vector transformation; 22 Sep 15;
       VM(2,1:3) = (/ + Zoij(jk,1), Roij(jk,1) * Zoij(jk,3) - Roij(jk,3) * Zoij(jk,1), - Roij(jk,1) /)
       VM(3,1:3) = (/   zero      , Roij(jk,2) * Zoij(jk,1) - Roij(jk,1) * Zoij(jk,2),   zero       /)
       
       lfail = -9 ; call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), lfail )

      !if( jrho.eq.0 .and. kk.eq.0 ) then
      !  write(0,'("aa00aa : " 10x " : ("es13.5" ,"es13.5" ,"es13.5" ) : B^(R,p,Z)="es13.5" ,"es13.5" ,"es13.5" ;")') RpZ(1:3), dBRpZ(1:3,0)
      !endif
       
       if( lfail.ne.0 ) then ; iaa00aa = 1 ; goto 9999
       endif
       
       lAtzmn%gBstz(jk,1:3,jrho) = matmul( VM(1:3,1:3), RpZ(1) * dBRpZ(1:3,0) ) ! note inclusion of R factor; 22 Sep 15;

      enddo ! end of do jj; 22 Sep 15;
     enddo ! end of do kk; 22 Sep 15;
     
     call ft02aa( Nfp, lAtzmn%Nt, lAtzmn%Nz, lAtzmn%gBstz(1:aNtz,2,jrho), lAtzmn%gBstz(1:aNtz,3,jrho), & ! not-destructive; 13 Oct 15;
                  amn, lAtzmn%im(1:amn), lAtzmn%in(1:amn), even(1:amn), oddd(1:amn), comn(1:amn), simn(1:amn), ift02aa )
     
     lAtzmn%gBtc(jrho,0,1:amn) = even(1:amn)
     lAtzmn%gBts(jrho,0,1:amn) = oddd(1:amn)
     lAtzmn%gBzc(jrho,0,1:amn) = comn(1:amn)
     lAtzmn%gBzs(jrho,0,1:amn) = simn(1:amn)
          
     if( ifail.le.-1 .and. jrho.gt.0 ) then ! compute "re-construction error" from Fourier harmonics ; 13 Oct 15;
      
      gBstz(1:aNtz,1:3) = zero ; surfaceerror(1:3) = zero
      
      do kk = 0, lAtzmn%Nz-1 ;                            ; zeta = kk * (pi2/Nfp) / lAtzmn%Nz ! loop over real-space grid; 13 Oct 15;
       do jj = 0, lAtzmn%Nt-1 ; jk = 1 + jj + kk*lAtzmn%Nt ; teta = jj *  pi2      / lAtzmn%Nt ! loop over real-space grid; 13 Oct 15;

        do ii = 1, amn ; mi = lAtzmn%im(ii) ; ni = lAtzmn%in(ii) ! sum Fourier harmonics; 13 Oct 15;
         arg = mi * teta - ni * zeta ; ca = cos(arg) ; sa = sin(arg)
         gBstz(jk,2) = gBstz(jk,2) + lAtzmn%gBtc(jrho,0,ii) * ca + lAtzmn%gBts(jrho,0,ii) * sa
         gBstz(jk,3) = gBstz(jk,3) + lAtzmn%gBzc(jrho,0,ii) * ca + lAtzmn%gBzs(jrho,0,ii) * sa
        enddo ! end of do ii ; 13 Oct 15;

       !ds = sqrt( ( Zoij(jk,3) * Roij(jk,2) - Roij(jk,3) * Zoij(jk,2) )**2 + Roij(jk,0)**2 * ( Roij(jk,2)**2 + Zoij(jk,2)**2 ) )
       !surfaceerror(2) = surfaceerror(2) + ds * abs( ( lAtzmn%gBstz(jk,2,jrho) - gBstz(jk,2) ) / lAtzmn%gBstz(jk,2,jrho) )
       !surfaceerror(3) = surfaceerror(3) + ds * abs( ( lAtzmn%gBstz(jk,3,jrho) - gBstz(jk,3) ) / lAtzmn%gBstz(jk,3,jrho) )
        surfaceerror(2) = surfaceerror(2) +      abs( ( lAtzmn%gBstz(jk,2,jrho) - gBstz(jk,2) ) / lAtzmn%gBstz(jk,2,jrho) )
        surfaceerror(3) = surfaceerror(3) +      abs( ( lAtzmn%gBstz(jk,3,jrho) - gBstz(jk,3) ) / lAtzmn%gBstz(jk,3,jrho) )

       enddo ! end of do jj ; 13 Oct 15;
      enddo ! end of do kk ; 13 Oct 15;
      
      surfaceerror(2) = surfaceerror(2) * (pi2/lAtzmn%Nt) * (pi2/Nfp)/lAtzmn%Nz
      surfaceerror(3) = surfaceerror(3) * (pi2/lAtzmn%Nt) * (pi2/Nfp)/lAtzmn%Nz
      
      if( ifail.le.-2 ) then 
        write(0,'("aa00aa : " 10x " : M ="i3" ; N ="i3" ; jrho =",i4," ; <gBt> ="es12.5" ; <gBz> ="es12.5" ;")') &
          aMpol,    aNtor,     jrho,           surfaceerror(2:3)
      endif
      
     endif ! end of if( ifail.le.-1) ; 23 Nov 15;


    enddo ! end of do jrho; 22 Sep 15;

    DALLOCATE(carg)
    DALLOCATE(sarg)

    DALLOCATE(Xbc)
    DALLOCATE(Xbs)
    DALLOCATE(Ybc)
    DALLOCATE(Ybs)

    DALLOCATE(even)
    DALLOCATE(oddd)
    DALLOCATE(comn)
    DALLOCATE(simn)

    DALLOCATE(gBstz) ! just used to calculate the error; 23 Nov 15;
    
! fit splines; 22 Sep 15;
    
    SALLOCATE(swork,(0:aLrad,0:3),zero)
    SALLOCATE(rwork,(0:aLrad    ),zero)
    
    do ii = 1, amn
     call nr00aa( paLrad, lAtzmn%ss(0:aLrad), lAtzmn%gBtc(0:aLrad,-1:2,ii), swork(0:aLrad,0:3), rwork(0:aLrad) ) 
     call nr00aa( paLrad, lAtzmn%ss(0:aLrad), lAtzmn%gBts(0:aLrad,-1:2,ii), swork(0:aLrad,0:3), rwork(0:aLrad) ) 
     call nr00aa( paLrad, lAtzmn%ss(0:aLrad), lAtzmn%gBzc(0:aLrad,-1:2,ii), swork(0:aLrad,0:3), rwork(0:aLrad) ) 
     call nr00aa( paLrad, lAtzmn%ss(0:aLrad), lAtzmn%gBzs(0:aLrad,-1:2,ii), swork(0:aLrad,0:3), rwork(0:aLrad) ) 
    enddo

    DALLOCATE(swork)
    DALLOCATE(rwork)

!    if( ifail.le.-1 ) then
!     do jrho = 0, aLrad ! calculate magnetic field on regular grid; 19 Nov 14;
!      lvol = lAtzmn%ss(jrho)
!      do kk = 0, lAtzmn%Nz-1 ;                            ; zeta = kk * (pi2/Nfp) / lAtzmn%Nz
!       do jj = 0, lAtzmn%Nt-1 ; jk = 1 + jj + kk*lAtzmn%Nt ; teta = jj *  pi2      / lAtzmn%Nt
!        lfail = 1 ; stz(1:3) = (/ lvol, teta, zeta /) ; call aa00ba( lAtzmn, stz(1:3), lBstz(1:3), lfail )
!        reconstructionerror(1) = reconstructionerror(1) + abs( lBstz(1) - lAtzmn%gBstz(jk,1,jrho) )
!       enddo ! end of do   jj; 13 Oct 15;
!      enddo ! end of do   kk; 13 Oct 15;
!     enddo ! end of do jrho; 13 Oct 15;
!    endif

    iaa00aa = 0
    
    DALLOCATE(Roij)
    DALLOCATE(Zoij)

9999 continue
    
    if( ifail.eq. 0 ) write(0,9000) iaa00aa, lAtzmn%Lrad, lAtzmn%Mpol, lAtzmn%Ntor
   !if( ifail.le.-1 ) write(0,9000) iaa00aa, lAtzmn%Lrad, lAtzmn%Mpol, lAtzmn%Ntor, reconstructionerror(1:3) * VV / (aLrad+1)
    
    ifail = iaa00aa
    
9000 format("aa00aa : ifail ="i3" : Lrad ="i6" ; Mpol ="i3" ; Ntor ="i3" ; ":"reconstruction error ="3es13.5" ;")
    
    return
    
  end subroutine aa00aa

!latex \subroutine{aa00ba}{evaluate vector potential (in toroidal coordinates);}
!latex \ben
!latex \item     \underline{\bf Required inputs}
!latex \item[  ] \verb+Atzmn                : structure ;+
!latex \bi
!latex \item[i.] the vector potential returned by \verb+aa00aa+;
!latex \ei
!latex \item[  ] \verb+stz(1:3)            : real      ;+
!latex \bi
!latex \item[i.] the required position;
!latex \ei
!latex \een

  subroutine aa00ba( lAtzmn, stz, gBstz, ifail )
    
    implicit none
    
    type(vectorpotential) :: lAtzmn
    REAL                  ::   stz(1:3)
    REAL                  :: gBstz(1:3)
    
    INTEGER               :: ifail
    
    INTEGER               :: aLrad, paLrad, aMpol, aNtor, Nfp, ii, mi, ni
    REAL                  :: lvol, teta, zeta, gBtc(-1:2), gBts(-1:2), gBzc(-1:2), gBzs(-1:2), arg, carg, sarg
    
    iaa00ba = ifail
    
    CHECKINPUT( aa00ba, lAtzmn%Lrad.le.0 , 9999 )
    CHECKINPUT( aa00ba, lAtzmn%Mpol.lt.1 , 9999 )
    CHECKINPUT( aa00ba, lAtzmn%Ntor.lt.0 , 9999 )
    CHECKINPUT( aa00ba, lAtzmn%Nfp .le.0 , 9999 )

    CHECKINPUT( aa00ba, lAtzmn%mn.le.0, 9999 )

    CHECKINPUT( aa00ba, .not.allocated(lAtzmn%im), 9999 )
    CHECKINPUT( aa00ba, .not.allocated(lAtzmn%in), 9999 )

    CHECKINPUT( aa00ba, .not.allocated(lAtzmn%ss  ), 9999 )

    CHECKINPUT( aa00ba, .not.allocated(lAtzmn%gBtc), 9999 )
    CHECKINPUT( aa00ba, .not.allocated(lAtzmn%gBts), 9999 )
    CHECKINPUT( aa00ba, .not.allocated(lAtzmn%gBzc), 9999 )
    CHECKINPUT( aa00ba, .not.allocated(lAtzmn%gBzs), 9999 )
    
    aLrad = lAtzmn%Lrad ; paLrad = aLrad + 1 ; aMpol = lAtzmn%Mpol ; aNtor = lAtzmn%Ntor ; Nfp = lAtzmn%Nfp ! internal shorthand; 22 Sep 15;
    
    lvol =         stz(1)
    teta = modulo( stz(2),  pi2      )
    zeta = modulo( stz(3),  pi2/Nfp  )
    
    gBstz(1:3) = zero ! initialize summation; 13 Oct 15;
    
    do ii = 1, lAtzmn%mn

     mi = lAtzmn%im(ii) ; ni = lAtzmn%in(ii) ; arg = mi * teta - ni * zeta ; carg = cos(arg) ; sarg = sin(arg)
     
     call nr00ab( paLrad, lAtzmn%ss(0:aLrad), lAtzmn%gBtc(0:aLrad,-1:2,ii), lvol, gBtc(-1:2) )
     call nr00ab( paLrad, lAtzmn%ss(0:aLrad), lAtzmn%gBts(0:aLrad,-1:2,ii), lvol, gBts(-1:2) )
     call nr00ab( paLrad, lAtzmn%ss(0:aLrad), lAtzmn%gBzc(0:aLrad,-1:2,ii), lvol, gBzc(-1:2) )
     call nr00ab( paLrad, lAtzmn%ss(0:aLrad), lAtzmn%gBzs(0:aLrad,-1:2,ii), lvol, gBzs(-1:2) )
     
     gBstz(1) = gBstz(1) - ( - gBtc(-1) * sarg + gBts(-1) * carg ) * (+mi) &
                         - ( - gBzc(-1) * sarg + gBzs(-1) * carg ) * (-ni)
     gBstz(2) = gBstz(2) + (   gBtc( 0) * carg + gBts( 0) * sarg )
     gBstz(3) = gBstz(3) + (   gBzc( 0) * carg + gBzs( 0) * sarg )
     
     if( iaa00ba.lt.0 ) write(0,1000) mi, ni, gBtc(0), gBts(0), gBzc(0), gBzs(0), gBstz(1:3)

1000 format("aa00ba : " 10x " : (m,n)=("i4","i4" ) : gBmn="4es13.5" ; gB="3es13.5" ;")

    enddo ! end of do ii; 13 Oct 15;
    
    iaa00ba = 0
    
9999 continue
    
    if( ifail.le. 0 ) write(0,9000) iaa00ba
    
    ifail = iaa00ba
    
9000 format("aa00ba : ifail ="i3" : ")
    
    return
    
  end subroutine aa00ba
  
!latex \subroutine{ir00aa}{construct irrational flux-surface (incomplete);}
!latex \ben
!latex \item     \underline{\bf Required inputs}
!latex \item[  ] \verb+Atzmn                : structure ;+
!latex \bi
!latex \item[i.] the vector potential returned by \verb+aa00aa+;
!latex \ei
!latex \een

  subroutine ir00aa( lsurface, ifail )
    
    implicit none
    
    type(irrationalsurface ) :: lsurface
    INTEGER                  :: ifail

    INTEGER                  :: astat, Mits, Mpol, Ntor, Nxdof, Npts, Nk, its
    REAL, allocatable        :: wsvd(:), dm(:)
    
    iir00aa = ifail
    
    CHECKINPUT( ir00aa, lsurface%Mits.le.0 , 9999 )
    CHECKINPUT( ir00aa, lsurface%Mpol.le.1 , 9999 )
    CHECKINPUT( ir00aa, lsurface%Ntor.le.0 , 9999 )
    
    Mits = lsurface%Mits ; Mpol = lsurface%Mpol ; Ntor = lsurface%Ntor

    isurface%Nxdof = 2*Mpol + 1 + 2*Mpol + 1 ! #independent degrees of freedom ; 13 Oct 15;

    Nxdof = isurface%Nxdof

    isurface%Npts  = Nxdof                            ! #points along curve to be mapped; 13 Oct 15;

    Npts = isurface%Npts 

    isurface%Nc    = 2 * Npts                                  ! #constraints                    ; 13 Oct 15;
    isurface%Nk    = max(2*Ntor,1)                    ! #discrete toroidal resolution   ; 13 Oct 15;

    Nk = isurface%Nk

    SALLOCATE(isurface%sap,(0:Npts,0:Nk-1),zero) 
    SALLOCATE(isurface%tap,(0:Npts,0:Nk-1),zero)

    SALLOCATE(isurface%xdof,(0:2*Mpol),zero)
    SALLOCATE(isurface%err,(0:2*Npts-1),zero)
    SALLOCATE(isurface%derr,(0:2*Npts-1,0:2*Mpol),zero)


    SALLOCATE(wsvd,(1:Nxdof),zero)
    SALLOCATE(dm,(0:2*Mpol),zero)

    do its = 1, Mits


   !xdof(     0:  Mpol) = irrsurf%sm(0:Mpol) ! pack degrees-of-freedom; 13 Oct 15;
   !xdof(Mpol+1:2*Mpol) = irrsurf%tm(1:Mpol)


    enddo

    DALLOCATE(wsvd)
    DALLOCATE(dm)
    

    DALLOCATE(isurface%xdof)
    DALLOCATE(isurface%err)
    DALLOCATE(isurface%derr)
    
    DALLOCATE(isurface%sap)
    DALLOCATE(isurface%tap)

    iir00aa = 0
    
9999 continue
    
    if( ifail.le. 0 ) write(0,9000) iir00aa
    
    ifail = iir00aa
    
9000 format("ir00aa : ifail ="i3" : ")
    
    return
    
  end subroutine ir00aa
  
!latex \subroutine{qf00aa}{construct quadratic-flux minimizing surface using pseudo fieldline following algorithm;}
!latex \ben
!latex \item     Quadratic-flux minimizing ({\small QFM}) surfaces are constructed using the pseudo-fieldline integration algorithm.
!latex           For a full description of {\small QFM}-surfaces please refer to the recent article by
!latex           Hudson \& Suzuki, \link{dx.doi.org/10.1063/1.4897390}{Phys. Plasmas, 21:102505, 2014} and references therein.
!latex           Only a brief outline will be given here.
!latex \item     The ``pseudo-field'' is defined in toroidal coordinates $(\r,\t,\z)$ as
!latex           \be {\bf B}_\nu \equiv {\bf B} - \nu B^\z {\bf e}_\r,
!latex           \ee
!latex           where $\nu$ is constant along a pseudo fieldline but is otherwise a-priori unknown. 
!latex \item     A $(p,q)$-{\small QFM} surface is a family of $(p,q)$-periodic fieldlines of ${\bf B}_\nu$, 
!latex           along each of which the action-gradient is constant,
!latex           i.e. $\nu=\nu(\alpha)$ where $\alpha$ labels pseudo fieldlines.
!latex \item     The task of finding these periodic pseudo fieldlines is equivalent to finding fixed points of the \mbox{$q$-th} return {\em pseudo}-map,
!latex           $P^q$, which is constructed by integrating along the pseudo-field from an initial point, $(\t_0,\r_0)$, on the \Poincare section $\z=0$,
!latex           around $q$ toroidal periods to arrive at $(\t_q,\r_q)$:
!latex           \be \left(\begin{array}{c} \t_q \\ \r_q\end{array}\right)=P^q\left(\begin{array}{c} \nu \\ \r_0 \end{array}\right), \label{eq:qthpseudomap}
!latex           \ee
!latex           where the dependence on $\t_0$ is suppressed.
!latex \item     Given that $\nu$ is constant along the pseudo-field but that the particular numerical value of $\nu$ is not yet known,
!latex           it is required to find the particular pair $(\nu,\r_0)$ that gives a periodic, integral curve of ${\bf B}_\nu$ at the prescribed angle,
!latex           $\alpha \equiv \t_0$.
!latex \item     The $q$-th return, pseudo tangent-map, $\nabla P^q$, is defined by
!latex           \be \left(\begin{array}{c}
!latex           \delta \t_q \\ \delta \r_q\end{array}\right)=\nabla P^q \cdot \left(\begin{array}{c} \delta \nu \\ \delta \r_0 \end{array}\right),
!latex           \ee
!latex           and can also be determined by pseudo fieldline integration over $\z \in [0,2\pi q]$ by
!latex           \be
!latex           \frac{d}{d \z} \nabla P^q =
!latex           \left(\begin{array}{cc} \partial_{\t} \dot \t & \partial_{\r} \dot \t \\
!latex           \partial_{\t} \dot \r & \partial_{\r} \dot \r \end{array}\right)\cdot \nabla P^q + \left(\begin{array}{cc}0 & 0 \\ 
!latex           1/\sqrt g B^\z & 0 \end{array}\right), \nonumber
!latex           \ee
!latex           with the initial condition
!latex           \be \nabla P^q = \left(\begin{array}{cc} 0 & 0 \\ 0 & 1 \end{array}\right).
!latex           \ee
!latex \item     The pseudo tangent map allows an efficient Newton iterative algorithm for finding fixed points:
!latex           a correction, $(\delta\nu,\delta\r)$, to an initial guess for $(\nu,\r)$ is determined by requiring the pseudo fieldline be periodic,
!latex           \be \left(\begin{array}{c} 
!latex           \t_q \\ \r_q\end{array}\right) + \nabla P^q \cdot \left(\begin{array}{c} \delta \nu \\ 
!latex           \delta \r_0 \end{array}\right) = 
!latex           \left(\begin{array}{c} \t_0+2\pi p \\ \r_0 + \delta \r \end{array}\right). \label{eq:pseudofieldlineintegrationNewtonmethod}
!latex           \ee
!latex \item     For the integrable case, for which there is a true periodic orbit for every value of the poloidal angle, 
!latex           the iterative solution will yield $\nu(\alpha)=0$ for all $\alpha$, and the pseudo field reduces to the true field; 
!latex           and similarly for the non-integrable case: where true periodic fieldlines exist, 
!latex           e.g. at $\alpha_X$ and $\alpha_0$, the solution yields \mbox{$\nu(\alpha_X)=0$} and \mbox{$\nu(\alpha_O)=0$}.
!latex \item     The pseudo fieldline with $\t_0 = \alpha$ serves as an initial guess for the pseudo fieldline with $\t_0 = \alpha + d\alpha$, 
!latex           and we may trace out the entire pseudo surface by varying $\t_0$.
!latex           Note that given that the $(p,q)$ pseudo fieldlines have length equal to $2\pi q$, 
!latex           it is not required to construct the pseudo curves over the range $\alpha \in[0,2\pi]$: 
!latex           periodicity means that it is only required to construct the curves over the range $\alpha \in[0,2\pi/q]$.
!latex \item     At all angle locations where a periodic {\em true} fieldline does not exist, 
!latex           a periodic {\em pseudo} fieldline can still be constructed by suitably choosing $\nu$.
!latex           Intuitively, we may think of $\nu/\sqrt g$ as the amount of radial field that must be subtracted from the true field 
!latex           to `cancel' the resonant effect of the perturbation, 
!latex           and by so doing create a rational `pseudo' surface as a family of rational pseudo fieldlines.
!latex           The function $\nu(\alpha)$ is sinusoidal.
!latex \item     \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : coordinates, qfminsurface, qf00aa+ \\ \\
!latex           \verb+type(coordinates)  :: rzmn+ \\
!latex           \verb+type(qfminsurface) :: qfms+ \\ \\ in their source that calls \verb+qf00aa+,
!latex           where \verb+coordinates+ and \verb+qfminsurface+ are derived types (i.e. structures)
!latex           that contains both the required input and output information.
!latex           The variable names, \verb+rzmn+ and \verb+qfms+, are arbitrary.
!latex \item     \underline{\bf Required inputs}
!latex \item[  ] \verb+rzmn                : structure ;+
!latex \bi
!latex \item[i.] the background toroidal coordinates, $(\r,\t,\z)$, that will be used;
!latex \item[ii.] it is assumed that before calling \verb+qf00aa+, that the routine \verb+bc00aa+ has been called, and \verb+rzmn+ should be unchanged;
!latex \ei
!latex \item[  ] \verb+qfms%Nfp            : integer ;+
!latex \bi
!latex \item[i.] the field periodicity;
!latex \ei
!latex \item[  ] \verb+qfms%pp             : integer ;+
!latex \item[  ] \verb+qfms%qq             : integer ;+
!latex \bi
!latex \item[i.] the rotational transform, $\iotabar \equiv $ \verb+pp/qq+;
!latex \item[ii.] the ``poloidal'' periodicity, \verb+pp+, must be a multiple of the field periodicity, \verb+Nfp+;
!latex \item[iii.] the integers \verb+pp+ and \verb+qq+ must be positive, and \verb+qq+ cannot be zero;
!latex \item[iv.]  note that the sign of the rotational transform depends on the background coordinates;
!latex             and if the background coordinates are such that the rotational transform is negative, 
!latex             then prior to calling \verb+bc00aa+ and \verb+qf00aa+ the ``handedness'' of the background coordinates must be changed;
!latex \ei
!latex \item[  ] \verb+qfms%Np             : integer ;+
!latex \bi
!latex \item[i.] the required poloidal resolution, e.g. \verb+Np=32+;
!latex \ei
!latex \item[  ] \verb+qfms%Ntor           : integer ;+
!latex \item[  ] \verb+qfms%Mpol           : integer ;+
!latex \bi
!latex \item[i.] the required Fourier resolution, e.g. \verb+Ntor=12+, \verb+Mpol=4+;
!latex \ei
!latex \item[  ] \verb+qfms%odetol         : integer ;+
!latex \bi
!latex \item[i.] the o.d.e. integration tolerance, e.g. \verb+odetol+$=10^{-8}$,
!latex \ei
!latex \item[  ] \verb+qfms%pqtol          : integer ;+
!latex \bi
!latex \item[i.] the required tolerance in locating periodic pseudo-fieldlines, e.g. \verb+odetol+$=10^{-6}$,
!latex \ei
!latex \item[  ] \verb+qfms%rr             : real    ;+
!latex \item[  ] \verb+qfms%nu             : real    ;+
!latex \bi
!latex \item[i.] initial guess for $\r$ and $\nu$ at $\t_0=0$;
!latex \item[ii.] only used if \verb+Lrestart+$=0$; (usually a suitable initial guess for $\nu$ is $\nu=0$).
!latex \ei
!latex \item[  ] \verb+qfms%Lrestart       : integer ;+
!latex \bi
!latex \item[i.] if \verb+Lrestart+$=1$, an initial guess for both $\r_i$ and $\nu_i$ at each $\t_i$ for $i=0, $ \verb+Np+ is required on input;
!latex \ei
!latex \item[  ] \verb+qfms%t(0:qNd,0:qNp) : real    ;+
!latex \item[  ] \verb+qfms%r(0:qNd,0:qNp) : real    ;+
!latex \item[  ] \verb+qfms%n(0:qNd,0:qNp) : real    ;+
!latex \bi
!latex \item[i.] only required on input if \verb+Lrestart.eq.1+;
!latex \item[ii.] usually these will be provided by an earlier call to \verb+qf00aa+;
!latex \ei
!latex \item[  ] \verb+ifail               : integer ;+
!latex \bi
!latex \item[i.] if \verb+ifail =  1+, quiet mode;
!latex \item[i.] if \verb+ifail =  0+, screen output is terse;
!latex \item[i.] if \verb+ifail = -1+, 
!latex \ei
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call qf00aa( rzmn, qfms, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+qfms%qNd : integer+
!latex \bi \item[] \verb+qNd = qq * max(Ntor,1)+;
!latex \ei
!latex \item[  ] \verb+qfms%qNp : integer +
!latex \bi \item[] \verb+qNp = qq * Np+ ;
!latex \ei
!latex \item[  ] \verb+qfms%t(0:qNd,0:qNp) : real    ;+
!latex \item[  ] \verb+qfms%r(0:qNd,0:qNp) : real    ;+
!latex \item[  ] \verb+qfms%n(0:qNd,0:qNp) : real    ;+
!latex \item[] \verb+qfms%mn : integer +
!latex \bi \item[i.] \verb-mn = Ntor + 1 + Mpol * ( 2 * Ntor + 1 )-
!latex \ei
!latex \item[  ] \verb+qfms%im(1:mn) : integer+
!latex \item[  ] \verb+qfms%in(1:mn) : integer+
!latex \bi \item[i.] mode identification;
!latex \ei
!latex \item[  ] \verb+qfms%Rbc(1:mn) : +
!latex \item[  ] \verb+qfms%Rbs(1:mn) : +
!latex \item[  ] \verb+qfms%Zbc(1:mn) : +
!latex \item[  ] \verb+qfms%Zbs(1:mn) : +
!latex \bi \item[i.] Fourier harmonics of {\small QFM}-surface in straight pseudo fieldline angle;
!latex \ei
!latex \item[  ] \verb+qfms%ok : integer+
!latex \bi \item[i.] error flag;
!latex \ei
!latex \item[i.] on output: 
!latex \bi \item[] \verb+ifail=0+ : normal execution;
!latex     \item[] \verb+ifail=1+ : input error;
!latex \ei
!latex \item[6.] Comments:
!latex \bi \item[i.] Presently, this routine uses the cylindrical field provided by \verb+bfield+, and the coordinate transformation provided by \verb+bc00aa+.
!latex               Note that the coordinate transformation involves the (slow) Fourier summation of potentially many harmonics in the coordinate transform,
!latex               and consequently \verb+qf00aa+ can be slow.
!latex     \item[ii.] When \verb+aa00aa+ is complete, I will include an option to use the magnetic vector potential in toroidal coordinates,
!latex                as this will eliminate the need for the coordinate transformation, and an exactly divergence-free numerical representation of the 
!latex                magnetic field will be enabled!
!latex \ei
!latex \een
  
  subroutine qf00aa( lrzmn, lqfms, ifail )
    
    implicit none
    
    type(coordinates)     :: lrzmn
    type(qfminsurface)    :: lqfms
    
    INTEGER, parameter    :: Ndof = 2, Ldfjac = Ndof, Lrwork = Ndof * ( 3 * Ndof + 13 ) / 2
    INTEGER               :: astat, ifail, pp, qq, Nfp, nn, lbp, Lrad, mn, Nd, qNd, qNp, ic05xbf, jz, kz, lz, Lderiv
    REAL                  :: teta, xx(1:Ndof), ff(1:Ndof), df(1:Ldfjac,1:Ndof), angle, pqerr
    REAL                  :: rwork(1:Lrwork), tol, stz(1:3), dRpZ(1:3,0:3,0:3), sqrtg

    INTEGER               :: ii, jj, kk, jk, Nt, Nz, Ntz, iflag, iangle, qNpp, ift02aa, ibc00ab, cunit, idiff
    REAL                  :: xd(1:Ndof,-1:1,-1:1), fd(1:Ndof,-1:1,-1:1), fdiff, alpha, itt(-1:2), irr(-1:2), area
    REAL, allocatable     :: Rjk(:), Zjk(:), tt(:,:), rr(:,:), rw(:), sw(:,:), aa(:,:)

    external              :: qf00ab, qf00ac

    iqf00aa = ifail

    Lrad = lrzmn%Lrad ; mn = lrzmn%mn
    
    CHECKINPUT( qf00aa,mn   .le.1                , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%im ) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%in ) , 9999 )
    CHECKINPUT( qf00aa,Lrad .le.0                , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%ss ) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%Rbc) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%Rbs) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%Zbc) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%Zbs) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%Xbc) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%Xbs) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%Ybc) , 9999 )
    CHECKINPUT( qf00aa,.not.allocated(lrzmn%Ybs) , 9999 )
    
    Nfp = lqfms%Nfp ; pp = lqfms%pp ; qq = lqfms%qq

    CHECKINPUT( qf00aa,       Nfp    .le.0     , 9999 )
    CHECKINPUT( qf00aa,       pp     .lt.0     , 9999 )
    CHECKINPUT( qf00aa,       qq     .le.0     , 9999 )
    CHECKINPUT( qf00aa, lqfms%Np     .le.0     , 9999 )
    CHECKINPUT( qf00aa, lqfms%Ntor   .lt.0     , 9999 )
    CHECKINPUT( qf00aa, lqfms%Mpol   .le.0     , 9999 )
    CHECKINPUT( qf00aa, lqfms%odetol .lt.small , 9999 )
    CHECKINPUT( qf00aa, lqfms%pqtol  .lt.small , 9999 )
   
    qfms%Nfp = Nfp ; qfms%pp = pp ; qfms%qq = qq

    qfms%odetol = lqfms%odetol ; qfms%pqtol = lqfms%pqtol
    
    lbp = pp / Nfp
    if( lbp*Nfp .ne. pp ) then ; write(0,'("qf00aa : input error ; (pp/Nfp)*Nfp.ne.pp ;")') ; iqf00aa = 1 ; goto 9999
    endif
    
    do nn = 1, qq-1
     if( mod( nn * lbp, qq ) .eq. 1 ) exit ! search for n s.t. n p = m q + 1, where n & m are integers; 14 Aug 13;
    enddo
    
    lqfms%nn = nn ; lqfms%mm = ( lqfms%nn * lbp - 1 ) / qq

    qfms%Id = 4 ; qfms%offset = thousand
    
    Nd = qfms%Id * max( lqfms%Ntor, 1 ) ; qNd = qq * Nd ; qNp = qq * lqfms%Np

    if( iqf00aa.le.0 ) then 
      write(0,'("qf00aa :         "2x" : ("i4" ,"i4" ) : nn ="i3" ; mm ="i3" ; qNd ="i6" ; Lrestart ="i2" ;")') &
        pp, qq, lqfms%nn, lqfms%mm, qNd, lqfms%Lrestart
    endif

    lqfms%qNd = qNd ; lqfms%qNp = qNp
     qfms%qNd = qNd ;  qfms%qNp = qNp ; qfms%Nd = Nd
    
    if( lqfms%Lrestart.eq.0 ) then
     if( lqfms%rr.lt.small ) then ; write(0,'("qf00aa : input error ; rr.le.small ;")') ; iqf00aa = 1 ; goto 9999
     endif
     if( allocated(lqfms%t) ) deallocate(lqfms%t)
     if( allocated(lqfms%r) ) deallocate(lqfms%r)
     if( allocated(lqfms%n) ) deallocate(lqfms%n)
     SALLOCATE(lqfms%t,(0:qNd,0:qNp), zero )
     SALLOCATE(lqfms%r,(0:qNd,0:qNp), zero )
     SALLOCATE(lqfms%n,(      0:qNp), zero )
    else
     CHECKINPUT(qf00aa, .not.allocated(lqfms%t) , 9999)
     CHECKINPUT(qf00aa, .not.allocated(lqfms%r) , 9999)
     CHECKINPUT(qf00aa, .not.allocated(lqfms%n) , 9999)
    endif
    
    SALLOCATE(qfms%t,(0:qNd,0:qNp),zero)
    SALLOCATE(qfms%r,(0:qNd,0:qNp),zero)
    SALLOCATE(qfms%n,(      0:qNp),zero)
    
    rzmn%Lrad = Lrad ; rzmn%mn = mn ; rzmn%mm = lrzmn%mm
    
    SALLOCATE(rzmn%im,(1:mn),lrzmn%im(1:mn))
    SALLOCATE(rzmn%in,(1:mn),lrzmn%in(1:mn))
        
    SALLOCATE(rzmn%ss,(0:Lrad),lrzmn%ss(0:Lrad))

    SALLOCATE(rzmn%Rbc,(0:Lrad,1:mn),lrzmn%Rbc(0:Lrad,1:mn))
    SALLOCATE(rzmn%Rbs,(0:Lrad,1:mn),lrzmn%Rbs(0:Lrad,1:mn))
    SALLOCATE(rzmn%Zbc,(0:Lrad,1:mn),lrzmn%Zbc(0:Lrad,1:mn))
    SALLOCATE(rzmn%Zbs,(0:Lrad,1:mn),lrzmn%Zbs(0:Lrad,1:mn))

    SALLOCATE(rzmn%Xbc,(0:Lrad,-1:2,1:mn),lrzmn%Xbc(0:Lrad,-1:2,1:mn))
    SALLOCATE(rzmn%Xbs,(0:Lrad,-1:2,1:mn),lrzmn%Xbs(0:Lrad,-1:2,1:mn))
    SALLOCATE(rzmn%Ybc,(0:Lrad,-1:2,1:mn),lrzmn%Ybc(0:Lrad,-1:2,1:mn))
    SALLOCATE(rzmn%Ybs,(0:Lrad,-1:2,1:mn),lrzmn%Ybs(0:Lrad,-1:2,1:mn))

    Lbfieldok = .true.

    angle = pi2 / qq
    
    do iteta = 0, lqfms%Np-1

     if( iteta.eq.1 ) pause
     
     izeta = 0
     
     teta = iteta * angle / lqfms%Np ! for iteta = 0 value of angle is irrelevant; angle will be calculated below; 14 Aug 13;
     
     if( lqfms%Lrestart.eq.1 ) then ! an initial guess is provided; 14 Aug 13;
      ;            ; xx(1:2) = (/ lqfms%r(izeta,iteta), lqfms%n(iteta) + qfms%offset /)    
     else ! obtain initial guess by extrapolation; 14 Aug 13;
      select case( iteta ) ! extrapolation; 31 Jul 13;
      case( 0  )   ; xx(1  ) =                     lqfms%rr 
       ;           ; xx(  2) =                     lqfms%nu                         + qfms%offset
      case( 1  )   ; xx(1  ) = sum( (/     +1 /) * lqfms%r(izeta,iteta-1:iteta-1) )
       ;           ; xx(  2) = sum( (/     +1 /) * lqfms%n(      iteta-1:iteta-1) ) + qfms%offset
      case( 2: )   ; xx(1  ) = sum( (/ -1, +2 /) * lqfms%r(izeta,iteta-2:iteta-1) )
       ;           ; xx(  2) = sum( (/ -1, +2 /) * lqfms%n(      iteta-2:iteta-1) ) + qfms%offset
      end select
     endif ! end of if( lrestart ) ;
     
     niter(0:1) = 0 ; ff(1:2) = zero ; tol = lqfms%pqtol / hundred

     qfms%t(izeta,iteta) = teta
     
     if( iqf00aa.le.-3 ) write(0,1000) pp, qq, iteta, teta, xx(1), xx(2)-qfms%offset
     
     if( iqf00aa.le.-9 ) then ! compute tangent map using finite-differences;  8 Jul 15;
      do idiff = -2, -4, -1
       fdiff = ten**idiff
       do ii = -1, 1
        do jj = -1, 1
         if( abs(ii)+abs(jj).ne.1 ) cycle
         xd(1:Ndof,ii,jj) = xx(1:Ndof) + (/ ii, jj /) * fdiff * half
         iflag = 1 ; call qf00ab( Ndof, xd(1:Ndof,ii,jj), fd(1:Ndof,ii,jj), df(1:Ldfjac,1:Ndof), Ldfjac, iflag )
        enddo
       enddo
       write(0,'("qf00aa :         "2x" : df="4f23.10" ; fdiff="es10.2" ;")') &
         (/ fd(1, 1, 0)-fd(1,-1, 0), fd(2, 1, 0)-fd(2,-1, 0), fd(1, 0, 1)-fd(1, 0,-1), fd(2, 0, 1)-fd(2, 0,-1) /) &
         / fdiff, fdiff
      enddo
      iflag = 2 ; call qf00ab( Ndof, xx(1:Ndof      ), ff(1:Ndof      ), df(1:Ldfjac,1:Ndof), Ldfjac, iflag )
       write(0,'("qf00aa :         "2x" : df="4f23.10" ;       "  10x "  ;")') df(1:Ldfjac,1:Ndof)
      pause
     endif

#ifdef HAVE_NAG
     ic05xbf = 1 ! NAG; 13 Oct 15;
     call C05PBF( qf00ab, Ndof, xx(1:Ndof), ff(1:Ndof), df(1:Ldfjac,1:Ndof), Ldfjac, tol, rwork(1:Lrwork), Lrwork, ic05xbf )
#else 
     call HYBRJ1( qf00ab, Ndof, xx(1:Ndof), ff(1:Ndof), df(1:Ldfjac,1:Ndof), Ldfjac, tol, ic05xbf, rwork(1:Lrwork), Lrwork )
     if (ic05xbf==1) ic05xbf=0 ! because the MINPACK error code is stupid
#endif
     
     if( ic05xbf.eq.-2 ) then ; xx(1:Ndof) = lxx(1:Ndof) ; ff(1:Ndof) = lff(1:Ndof) ! ic05xbf = 0 ! 13 see ps00bb for definition of iflag = -2 ; Dec 13;
     endif

     pqerr = sum( abs(ff(1:Ndof)) ) / Ndof
     
     lqfms%t(0:qNd,iteta) = qfms%t(0:qNd,iteta)
     lqfms%r(0:qNd,iteta) = qfms%r(0:qNd,iteta)

     lqfms%n(      iteta) = xx(2) - qfms%offset

     lqfms%t(0    ,iteta) = teta
     lqfms%r(0    ,iteta) = xx(1)

     if( iqf00aa.le.-1 ) then
      
      select case( ic05xbf ) !  8 Jul 15;                                                               0123456789 
      case( -2 )
        write(0,1000) pp, qq, iteta, teta, xx(1), xx(2)-qfms%offset, ff(1:2), pqerr, ic05xbf, niter(0:1), "accepts ;"
      case( -1 )
        write(0,1000) pp, qq, iteta, teta, xx(1), xx(2)-qfms%offset, ff(1:2), pqerr, ic05xbf, niter(0:1), "B error ;"
      case(  0 )
        write(0,1000) pp, qq, iteta, teta, xx(1), xx(2)-qfms%offset, ff(1:2), pqerr, ic05xbf, niter(0:1), "        ;"
      case(  1 )
        write(0,1000) pp, qq, iteta, teta, xx(1), xx(2)-qfms%offset, ff(1:2), pqerr, ic05xbf, niter(0:1), "input ! ;"
      case(  2 )
        write(0,1000) pp, qq, iteta, teta, xx(1), xx(2)-qfms%offset, ff(1:2), pqerr, ic05xbf, niter(0:1), "restart ;"
      case(  3 )
        write(0,1000) pp, qq, iteta, teta, xx(1), xx(2)-qfms%offset, ff(1:2), pqerr, ic05xbf, niter(0:1), "xtol !  ;"
      case(  4 )
        write(0,1000) pp, qq, iteta, teta, xx(1), xx(2)-qfms%offset, ff(1:2), pqerr, ic05xbf, niter(0:1), "failed  ;"
      case default
        stop "qf00aa : 0123456789 : illegal ifail returned by C05PBF ;"
      end select
      
     endif
     
1000 format("qf00aa : ", 10x ," : (",i4," ,",i4," ) : "i3" : teta =",f10.07," ; (rho,nu)=(", f20.15," ,",es23.15," ) ; ":,&
   "ff=",2es10.02," ; err="es09.2" ; ic05xbf=",i3," ; its=",2i4," ; ",a9)

    !ic05xbf = 0 ! over-rule all error flags;  8 Jul 15;

     if( pqerr.lt.lqfms%pqtol ) ic05xbf = 0 ! shall over-rule error flag returned by NAG; 10 Oct 13;
          
     if( ic05xbf.ne.0 ) exit ! an error has occured; 
     
     if( iteta.eq.0 ) then
      lqfms%rr = xx(1)
      lqfms%nu = xx(2) - qfms%offset
      FATAL(qf00aa,lqfms%nn*Nd.lt.0 .or. lqfms%nn*Nd.gt.qNd, modulo arithmetic error)
      angle = lqfms%t(lqfms%nn*Nd,iteta) - lqfms%mm * pi2
     endif
     
    enddo ! end of do iteta = 0, Np-1 ! 05 Mar 14;
    
    if( ic05xbf.ne.0 ) then ; iqf00aa = 2 ; goto 9998
    endif

    do iteta = lqfms%Np, qNp ! complete pseudo-surface using periodicity identities; 14 Aug 13;
     
     lqfms%n(iteta) = lqfms%n(iteta-lqfms%Np)
     
     do izeta = 0, qNd-1
      
      jz = izeta + lqfms%nn * Nd ; kz = jz / qNd ; lz = jz - kz * qNd
      
      if    ( lz.lt.  0 ) then ; write(0,'("qf00aa : 0123456789 : ("i4" ,"i4" ) : lz.lt.  0 ;")') pp, qq
      elseif( lz.eq.qNd ) then ; write(0,'("qf00aa : 0123456789 : ("i4" ,"i4" ) : lz.eq.qNd ;")') pp, qq
      elseif( lz.gt.qNd ) then ; write(0,'("qf00aa : 0123456789 : ("i4" ,"i4" ) : lz.gt.qNd ;")') pp, qq
      endif
      
      lqfms%t(izeta,iteta) = lqfms%t(lz,iteta-lqfms%Np) + kz * (pi2/Nfp) * pp - lqfms%mm * pi2
      lqfms%r(izeta,iteta) = lqfms%r(lz,iteta-lqfms%Np)
      
     enddo ! end of do izeta; 14 Aug 13;

    enddo ! end of do iteta; 14 Aug 13;

    if( allocated(lqfms%X) ) deallocate(lqfms%X)
    if( allocated(lqfms%Y) ) deallocate(lqfms%Y)

    SALLOCATE(lqfms%X,(0:qNd-1,0:qNp), zero )
    SALLOCATE(lqfms%Y,(0:qNd-1,0:qNp), zero )
     
    do iteta = 0, qNp ! all poloidal angles; 05 Mar 14;
     do izeta = 0, qNd-1
      stz(1:3) = (/ lqfms%r(izeta,iteta), lqfms%t(izeta,iteta), izeta * (pi2/Nfp) / Nd /)
      Lderiv = 0 ; call bc00ab( rzmn, stz(1:3), Lderiv, dRpZ(1:3,0:3,0:3), sqrtg, ibc00ab )
      lqfms%X(izeta,iteta) = dRpZ(1,0,0)
      lqfms%Y(izeta,iteta) = dRpZ(3,0,0)
     enddo ! end of do iteta; 14 Aug 13;
    enddo ! end of izeta; 05 Mar 14;
    
    Nt = 4 * lqfms%Mpol ; Nz = Nd ; Ntz = Nt * Nz ; mn = lqfms%Ntor + 1 + lqfms%Mpol * ( 2 * lqfms%Ntor + 1 )
    
    lqfms%mn = mn

    if( allocated(lqfms%im) ) deallocate(lqfms%im)
    if( allocated(lqfms%in) ) deallocate(lqfms%in)

    SALLOCATE(lqfms%im,(1:mn),0)
    SALLOCATE(lqfms%in,(1:mn),0)
    
    call gi00ab( lqfms%Mpol, lqfms%Ntor, Nfp, mn, lqfms%im(1:mn), lqfms%in(1:mn) )

    if( allocated(lqfms%Rbc)) deallocate(lqfms%Rbc)
    if( allocated(lqfms%Rbs)) deallocate(lqfms%Rbs)
    if( allocated(lqfms%Zbc)) deallocate(lqfms%Zbc)
    if( allocated(lqfms%Zbs)) deallocate(lqfms%Zbs)

    SALLOCATE(lqfms%Rbc,(1:mn),zero)
    SALLOCATE(lqfms%Rbs,(1:mn),zero)
    SALLOCATE(lqfms%Zbc,(1:mn),zero)
    SALLOCATE(lqfms%Zbs,(1:mn),zero)

!   if( allocated(lqfms%a) ) deallocate(lqfms%a)
    
    SALLOCATE(     aa,(0:qNd,0:qNp), zero )
    
    izeta = 0
    
    do iteta = 0, qNp
     area                 = (pi2/Nfp) / Nd * ( sum( lqfms%t(0:qNd-1,iteta) ) + half * pp * (pi2/Nfp) ) !  8 Jul 15;
     aa(izeta,iteta) = area / (pi2/Nfp) / qq - half * pp * (pi2/Nfp)
    enddo

    do iteta =    qNp, 0, -1
     aa(izeta,iteta) = aa(izeta,iteta) - aa(izeta,0) ! shifting; 27 Nov 13;
    enddo
    
    do iteta = 0, qNp
     do izeta = 1, Nd-1
      aa(izeta,iteta) = aa(0,iteta) + izeta * pp * (pi2/Nfp) / Nd / qq ! implement straight field-line angle; 11 Aug 13;
     enddo
    enddo
    
    SALLOCATE(Rjk,(1:Ntz),zero)
    SALLOCATE(Zjk,(1:Ntz),zero)
    
    qNpp = qNp + 1
    
    SALLOCATE(tt,(0:qNp,-1:2),zero)
    SALLOCATE(rr,(0:qNp,-1:2),zero)
    SALLOCATE(sw,(0:qNp, 0:3),zero)
    SALLOCATE(rw,(0:qNp     ),zero)
    
    do kk = 0, Nz - 1

     tt(0:qNp,0) = lqfms%t(kk,0:qNp)
    !tt(0    ,1) = ( tt(  1,0) - tt(    0,0) ) / ( aa(kk,  1) - aa(kk,    0) )
    !tt(  qNp,1) = ( tt(qNp,0) - tt(qNp-1,0) ) / ( aa(kk,qNp) - aa(kk,qNp-1) )  
     tt(0    ,1) = ( tt(  1,0) - tt(qNp-1,0) + pi2 ) / ( aa(kk,  1) - aa(kk,qNp-1) + pi2 )
     tt(  qNp,1) = ( tt(  1,0) - tt(qNp-1,0) + pi2 ) / ( aa(kk,  1) - aa(kk,qNp-1) + pi2 )

     rr(0:qNp,0) = lqfms%r(kk,0:qNp)
    !rr(0    ,1) = ( rr(  1,0) - rr(    0,0) ) / ( aa(kk,  1) - aa(kk,    0) )
    !rr(  qNp,1) = ( rr(qNp,0) - rr(qNp-1,0) ) / ( aa(kk,qNp) - aa(kk,qNp-1) )
     rr(0    ,1) = ( rr(  1,0) - rr(qNp-1,0) + pi2 ) / ( aa(kk,  1) - aa(kk,qNp-1) + pi2 )
     rr(  qNp,1) = ( rr(  1,0) - rr(qNp-1,0) + pi2 ) / ( aa(kk,  1) - aa(kk,qNp-1) + pi2 )
   
     call nr00aa( qNpp, aa(kk,0:qNp), tt(0:qNp,-1:2), sw(0:qNp,0:3), rw(0:qNp) )
     call nr00aa( qNpp, aa(kk,0:qNp), rr(0:qNp,-1:2), sw(0:qNp,0:3), rw(0:qNp) )
     
     do jj = 0, Nt-1 ; jk = 1 + jj + kk*Nt
      
      alpha = jj * pi2 / Nt
      if( alpha.lt.aa(kk,  0) ) then ; iangle = ( aa(kk,  0) - alpha ) / pi2 + 1 ; alpha = alpha + iangle * pi2
      endif
      if( alpha.gt.aa(kk,qNp) ) then ; iangle = ( alpha - aa(kk,qNp) ) / pi2 + 1 ; alpha = alpha - iangle * pi2
      endif
      
      call nr00ab( qNpp, aa(kk,0:qNp), tt(0:qNp,-1:2), alpha, itt(-1:2) )
      call nr00ab( qNpp, aa(kk,0:qNp), rr(0:qNp,-1:2), alpha, irr(-1:2) )

      stz(1:3) = (/ irr(0), itt(0), kk * (pi2/Nfp) / Nd /)
      Lderiv = 0 ; call bc00ab( rzmn, stz(1:3), Lderiv, dRpZ(1:3,0:3,0:3), sqrtg, ibc00ab )

      Rjk(jk) = dRpZ(1,0,0)
      Zjk(jk) = dRpZ(3,0,0)
      
     enddo

    enddo

!    cunit = 20
!    open(cunit, file="RZ", status="unknown", form="unformatted" )
!    write(cunit) Nt, Nz, Ntz   
!    write(cunit) Rjk
!    write(cunit) Zjk
!    close(cunit)

    DALLOCATE(rw)
    DALLOCATE(sw)
    DALLOCATE(rr)
    DALLOCATE(tt)
    
    ift02aa = 1
    call ft02aa( Nfp, Nt, Nz, Rjk(1:Ntz), Zjk(1:Ntz), & ! not destructive; 13 Oct 15;
  mn, lqfms%im(1:mn), lqfms%in(1:mn), lqfms%Rbc(1:mn), lqfms%Rbs(1:mn), lqfms%Zbc(1:mn), lqfms%Zbs(1:mn), ift02aa )

    DALLOCATE(Rjk)
    DALLOCATE(Zjk)

    DALLOCATE(aa)
  
9998 continue
    
    DALLOCATE(qfms%t)
    DALLOCATE(qfms%r)
    DALLOCATE(qfms%n)
    
    DALLOCATE(rzmn%im)
    DALLOCATE(rzmn%in)

    DALLOCATE(rzmn%ss)
    
    DALLOCATE(rzmn%Rbc)
    DALLOCATE(rzmn%Rbs)
    DALLOCATE(rzmn%Zbc)
    DALLOCATE(rzmn%Zbs)

    DALLOCATE(rzmn%Xbc)
    DALLOCATE(rzmn%Xbs)
    DALLOCATE(rzmn%Ybc)
    DALLOCATE(rzmn%Ybs)
    
    iqf00aa = ic05xbf
    
9999 continue
    
    if( ifail.le. 0 ) write(0,9000) iqf00aa, lqfms%pp, lqfms%qq
    
    ifail = iqf00aa ; lqfms%ok = ifail
    
9000 format("qf00aa : ifail ="i3" : ("i4" ,"i4" ) ;")
    
    return
    
  end subroutine qf00aa

#ifndef HAVE_NAG
! dummy routine passed to odepack functions in place of Jacobian
subroutine du00aa( neq, t, y, ml, mu, pd, nrowpd )
  
  implicit none
  
  INTEGER  :: neq, ml, mu, nrowpd
  REAL :: t, y(*), pd(nrowpd,*)
  
  return
  
end subroutine du00aa
#endif
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

end module oculus

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine qf00ab( Ndof, xx, ff, df, Ldfjac, iflag )

  use oculus, only : zero, one, pi2
  use oculus, only : qfms, izeta, iteta, pf00aa, du00aa, lxx, lff, Lbfieldok, iqf00aa, niter, actiongradient, itangent

  implicit none

  INTEGER            :: Ndof, Ldfjac, iflag
  REAL               :: xx(1:Ndof), ff(1:Ndof), df(1:Ldfjac,1:Ndof)
  
  INTEGER, parameter :: Node = 6, Lrwork = 20 * Node
  
  INTEGER            :: id02bjf, Nfp
  REAL               :: teta, rt(1:Node), zetastart, zetaend, pqerr
  CHARACTER          :: relabs
  
  external           :: qf00ac

#ifdef HAVE_NAG
    REAL                   :: rwork(1:Lrwork)
    external               :: D02BJW
#else
    INTEGER, parameter     :: liw=20, lrw=20+16*Node
    INTEGER                :: iwork(liw)
    REAL                   :: rwork(lrw)
    INTEGER                :: iopt,istate,itask,itol,mf,ii
    REAL                   :: atol,rtol,zetac,zetae
#endif
  
  select case( iflag )
  case( 1 )    ; itangent = 0
  case( 2 )    ; itangent = 1
  end select
  
  relabs = 'D' ; Nfp = qfms%Nfp
  
  izeta = 0 ; actiongradient = xx(2) - qfms%offset ; teta = qfms%t(izeta,iteta)
  
  rt(1:Node) = (/ xx(1), teta, one, zero, zero, zero /) ! set starting point for (pseudo) field-line integration; 31 Jul 13;
  
 !qfms%r(izeta,iteta) = rt(1)
  
  zetastart = zero ; zetaend = zetastart + qfms%qq * (pi2/Nfp)
    
#ifdef HAVE_NAG
  id02bjf = 1 ! NAG; 13 Oct 15;
  call D02BJF( zetastart, zetaend, Node, rt(1:Node), pf00aa, qfms%odetol, relabs, qf00ac, D02BJW, rwork(1:Lrwork), id02bjf )
#else
! set the integrator parameters.
  rwork=0.
  iwork=0
! istate=1 :  indicates the first lsode call
  istate=1
! itask=4 :  normal integration with limited over-shoot (set by
!            rwork(1) in the i_lsode loop
  itask=1
! iopt=1 :  optional inputs are used
  iopt=1
! rwork(6) :  set maximum lsode-internal step size
  iwork(6)=400000
! rwork(7) :  set minimum lsode-internal step size
! iwork(6) :  set maximum lsode-internal steps
! iwork(7) :  set maximum lsode-internal error messages printed
! mf=10 :  non-stiff Adams method of integration
  mf=10
! itol=1 :  indicates absolute tolerance is just a scalar
  itol=1
! rtol :  relative tolerance
  rtol=qfms%odetol
! atol :  absolute tolerance
  atol=qfms%odetol
! initializations for loop
  zetac=zetastart
  zetae=zetac
  call qf00ac(zetae,rt)
  do ii=1,qfms%Nd
    call lsode(pf00aa,(/Node/),rt,zetac,zetae, &
               itol,(/rtol/),(/atol/),itask,istate,iopt, &
               rwork,lrw,iwork,liw,du00aa,mf)
    zetac=zetae
    call qf00ac(zetae,rt)
  enddo
  id02bjf=istate
  if (istate>0) id02bjf=0
#endif

  if( .not.Lbfieldok ) id02bjf = -1
  
  select case( iflag )
  case( 1 )    ; niter(0) = niter(0)+1 ; ff(1  ) = rt(1) -   xx(1)
   ;           ;                       ; ff(2  ) = rt(2) - ( teta + qfms%pp * (pi2/Nfp) )
  case( 2 )    ; niter(1) = niter(1)+1 ; df(1,1) = rt(3) -   one ; df(1,2) = rt(4)
   ;           ;                       ; df(2,1) = rt(5)         ; df(2,2) = rt(6)
  end select

  lxx(1:Ndof) = xx(1:Ndof) ; lff(1:Ndof) = ff(1:Ndof)
  
  pqerr = sum( abs(ff(1:Ndof)) ) / Ndof

  if( iqf00aa.ne.-9 ) then
   if( ( iqf00aa.le.-2 .and. iteta.eq.0 ) .or. iqf00aa.le.-3 ) then
    write(0,1000) qfms%pp, qfms%qq, iteta, teta, xx(1), xx(2)-qfms%offset, ff(1:2), pqerr, id02bjf
   endif
  endif

1000 format("qf00ab : ", 10x ," : (",i4," ,",i4," ) : "i3" : teta =",f10.07," ; (rho,nu)=(", f20.15," ,",es23.15," ) ; ":,&
   "ff=",2es10.02," ; err="es09.2" ; id02bjf=",i3," ;")

  if( iflag.eq.1 .and. pqerr.lt.qfms%pqtol ) iflag = -2 ! user accepts     ; terminate; 13 Dec 13;

  if( id02bjf.ne.0                         ) iflag = -1 ! error has occured; terminate;  5 Jun 13;
  
 izeta = qfms%qNd ; qfms%r(izeta,iteta) = rt(1) ; qfms%t(izeta,iteta) = rt(2) ! ensure end-point is recorded;
  
  return
  
end subroutine qf00ab

subroutine qf00ac( zeta, rt ) ! saves information along pseudo-surf into surf structure; 10 Aug 13;

  use oculus, only : izeta, iteta, pi2, qfms, iqf00aa

  implicit none

  INTEGER, parameter :: Node = 6
  REAL               :: zeta, rt(1:Node)
  
!  if( iqf00aa.le.-4 ) then
!   write(0,'("qf00ac : ", 10x ," : (",i4," ,",i4," ) : "i3" : teta =", 10x  ," ; (rho,nu)=(",  18x   " ,",  23x  ," ) ; izeta="i3" ;")') &
! qfms%pp, qfms%qq, iteta,                           izeta
!  endif
  
  if( izeta.le.qfms%qNd ) then ; qfms%r(izeta,iteta) = rt(1)
   ;                           ; qfms%t(izeta,iteta) = rt(2)
  endif
  
  izeta = izeta + 1 ; zeta = izeta * (pi2/qfms%Nfp) / qfms%Nd
  
  return
  
end subroutine qf00ac

!latex \newpage \section{miscellaneous/auxilliary subroutines}

subroutine nr00aa( NN, xx, dy, yn, wk ) ! numerical recipes in fortran: spline;  8 Jul 15;
  
  use oculus, only : zero, half, one, two, three, six
  
  implicit none
  
  INTEGER :: NN
  REAL    :: xx(1:NN), dy(1:NN,-1:2), yn(1:NN,0:3), wk(1:NN)
  
  INTEGER :: ii, jj
  REAL    :: pp, qn, sigma, un, hh, exponent
  
  ii = 1
  
  dy(ii,2) = - half ; hh = xx(ii+1)-xx(ii)
  wk(ii) = ( three / hh ) * ( (dy(ii+1,0)-dy(ii,0)) / hh - dy(ii,1) )
  
  do ii = 2, NN-1 ; hh = xx(ii+1)-xx(ii)
   sigma = (xx(ii)-xx(ii-1)) / (xx(ii+1)-xx(ii-1))
   pp = sigma * dy(ii-1,2) + two
   dy(ii,2) = ( sigma - one ) / pp
   wk(ii) = ( six * ( (dy(ii+1,0)-dy(ii,0)) / hh &
                    - (dy(ii,0)-dy(ii-1,0)) / (xx(ii)-xx(ii-1)) ) / (xx(ii+1)-xx(ii-1)) - sigma * wk(ii-1) ) / pp
  enddo
  
  qn = half
  un = ( three / (xx(NN)-xx(NN-1)) ) * ( dy(NN,1) - (dy(NN,0)-dy(NN-1,0)) / (xx(NN)-xx(NN-1)) )
  
  ;                   ; dy(NN,2) = ( un - qn * wk(NN-1) ) / ( qn * dy(NN-1,2) + one )
  do ii = NN-1, 1, -1 ; dy(ii,2) = dy(ii,2) * dy(ii+1,2) + wk(ii)
  enddo
  
  dy(1,-1) = zero ! initialize integration; 11 Sep 15;

  do ii = 1, NN-1 ; hh = xx(ii+1)-xx(ii) ! convert to polynomial coefficients;  8 Jul 15;

   yn(ii,0)=( 6*dy(ii,0)*xx(ii+1) &
              -(hh-xx(ii+1))*(-6*dy(ii+1,0)+xx(ii+1)*(dy(ii+1,2)*(2*hh-xx(ii+1))+dy(ii,2)*(hh+xx(ii+1)))))/(6*hh)
   yn(ii,1)=(-6*dy(ii,0)+6*dy(ii+1,0) &
             +(dy(ii,2)+2*dy(ii+1,2))*hh**2-6*dy(ii+1,2)*hh*xx(ii+1)+3*(-dy(ii,2)+dy(ii+1,2))*xx(ii+1)**2)/(6*hh)
   yn(ii,2)=(dy(ii+1,2)*hh+dy(ii,2)*xx(ii+1)-dy(ii+1,2)*xx(ii+1))/(2*hh)
   yn(ii,3)=(-dy(ii,2)+dy(ii+1,2))/(6*hh)

   ;                    ; dy(ii+1,-1) = dy(ii  ,-1)
   do jj = 0, 3
    exponent = jj + one ; dy(ii+1,-1) = dy(ii+1,-1) + yn(ii,jj) * ( xx(ii+1)**exponent -xx(ii)**exponent ) / exponent
   enddo

  enddo

 return

end subroutine nr00aa

subroutine nr00ab( NN, xx, dy, px, py ) ! numerical recipes in fortran: splint;  8 Jul 15;
  
  use oculus, only : zero, one, six
  
  implicit none
  
  INTEGER :: NN
  REAL    :: xx(1:NN), dy(1:NN,-1:2), px, py(-1:2)
  INTEGER :: kk, ku, kl
  REAL    :: AA, BB, dA, dB, hh, lx, y0, y1, y2, y3, Al, Au, Bl, Bu
  
  if( px.gt.xx(NN) ) then
   
   lx = xx(NN) ; ku = NN ; kl = NN-1
   
   hh = xx(ku)-xx(kl)
   
   AA = ( xx(ku) - lx     ) / hh ; BB = ( lx     - xx(kl) ) / hh
   dA =          - one      / hh ; dB =   one               / hh
   
   y0 = AA*dy(kl,0) + BB*dy(ku,0) + ( (  AA**3-AA)    * dy(kl,2) + (  BB**3-BB)    * dy(ku,2) ) * hh**2 / six
   y1 = dA*dy(kl,0) + dB*dy(ku,0) + ( (3*AA**2- 1)*dA * dy(kl,2) + (3*BB**2- 1)*dB * dy(ku,2) ) * hh**2 / six
   y2 =                                  AA           * dy(kl,2) +    BB           * dy(ku,2)
   y3 =                                  dA           * dy(kl,2) +    dB           * dy(ku,2)
   
   py(0) = y0 + y1*(px-lx) + y2*(px-lx)**2 / 2 + y3*(px-lx)**3 / 6
   py(1) =      y1         + y2*(px-lx)        + y3*(px-lx)**2 / 2
   py(2) =                   y2                + y3*(px-lx)
  !py(3) =                                     + y3

   py(-1) = zero ! under construction; 13 Oct 15;

   return
   
  endif
  
  kl = 1 ; ku = NN
  
1 if ( ku-kl .gt. 1 )  then
   kk = ( ku + kl ) / 2
   if( xx(kk).gt.px ) then ; ku = kk
   else                    ; kl = kk
   endif
   goto 1
  endif
  
  hh = xx(ku) - xx(kl)
  
! if (h.eq.0.) fatal error; the input xx must be distinct;
  
  AA = ( xx(ku) - px     ) / hh ; BB = ( px     - xx(kl) ) / hh
  dA =          - one      / hh ; dB =   one               / hh
  
  py( 0) = AA*dy(kl,0) + BB*dy(ku,0) + ( (  AA**3-AA)    * dy(kl,2) + (  BB**3-BB)    * dy(ku,2) ) * hh**2 / six
  py( 1) = dA*dy(kl,0) + dB*dy(ku,0) + ( (3*AA**2- 1)*dA * dy(kl,2) + (3*BB**2- 1)*dB * dy(ku,2) ) * hh**2 / six
  py( 2) =                                  AA           * dy(kl,2) +    BB           * dy(ku,2)
 !py( 3) =                                  dA           * dy(kl,2) +    dB           * dy(ku,2)

  Al = ( xx(ku) - xx(kl) ) / hh ; Au = ( xx(ku) - px     ) / hh
  Bl = ( xx(kl) - xx(kl) ) / hh ; Bu = ( px     - xx(kl) ) / hh

  py(-1) = dy(kl,-1) &
         - ( dy(kl,0) - dy(kl,2) * hh**2 / six ) * ( Au**2 - Al**2 ) * hh / 2 - dy(kl,2) * ( Au**4 - Al**4 ) * hh**3 / 24 &
         + ( dy(ku,0) - dy(ku,2) * hh**2 / six ) * ( Bu**2 - Bl**2 ) * hh / 2 + dy(ku,2) * ( Bu**4 - Bl**4 ) * hh**3 / 24
  
  return
  
end subroutine nr00ab

subroutine gi00ab( Mpol, Ntor, Nfp, mn, im, in )
  
  implicit none
  
  INTEGER, intent(in)  :: Mpol, Ntor, Nfp, mn
  INTEGER, intent(out) :: im(mn), in(mn)
  
  INTEGER              :: imn, mm, nn
  
  imn = 0
  
  do mm = 0, 0
   do nn =     0, Ntor
    imn = imn+1 ; im(imn) = mm ; in(imn) = nn*Nfp
   enddo
  enddo
  
  do mm = 1, Mpol
   do nn = -Ntor, Ntor
    imn = imn+1 ; im(imn) = mm ; in(imn) = nn*Nfp
   enddo
  enddo
  
  return
  
end subroutine gi00ab

subroutine ft02aa( Nfp, Nt, Nz, ijreal, ijimag, mn, im, in, even, oddd, cosi, sini, ift02aa )
  
  use oculus, only : zero, one
  
  implicit none
  
  INTEGER             :: Nfp, Nt, Nz, mn, im(1:mn), in(1:mn), ift02aa
  REAL                :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), even(1:mn), oddd(1:mn), cosi(1:mn), sini(1:mn)
  
  INTEGER             :: Ntz, ii, mm, nn, ic06fuf
  REAL                :: jireal(1:Nt*Nz), jiimag(1:Nt*Nz), trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Nt*Nz), sqrtoNtz
  CHARACTER           :: isr

#ifndef HAVE_NAG
  INTEGER             :: astat, Lwka, Lwkb
  REAL, allocatable   :: wka(:), wkb(:)
  COMPLEX             :: jicomp(1:Nt,1:Nz)
#endif
  
  isr = 'I' ; Ntz = Nt * Nz ; sqrtoNtz = sqrt(one*Ntz)
  
  jireal(1:Ntz) = ijreal(1:Ntz)
  jiimag(1:Ntz) = ijimag(1:Ntz)
  
#ifdef HAVE_NAG
  ic06fuf = 1 ! NAG; 13 Oct 15;
  call C06FUF( Nt, Nz, jireal(1:Ntz), jiimag(1:Ntz), isr, trigm(1:2*Nt), trign(1:2*Nz), trigwk(1:2*Ntz), ic06fuf )
  select case( ic06fuf )
  case( 0 ) ;!write(0,'("ft02aa : " 10x " : ifail="i3" ; success ;                    ")')
  case( 1 ) ; write(0,'("ft02aa : " 10x " : ifail="i3" ; M<1 ;                        ")')
  case( 2 ) ; write(0,'("ft02aa : " 10x " : ifail="i3" ; N<1 ;                        ")')
  case( 3 ) ; write(0,'("ft02aa : " 10x " : ifail="i3" ; init error;                  ")')
  case( 4 ) ; write(0,'("ft02aa : " 10x " : ifail="i3" ; not used ;                   ")')
  case( 5 ) ; write(0,'("ft02aa : " 10x " : ifail="i3" ; inconsistent trigm or trign ;")')
  case( 6 ) ; write(0,'("ft02aa : " 10x " : ifail="i3" ; unexpected ;                 ")')
  end select
#else
  jicomp = CMPLX( RESHAPE( jireal, (/Nt,Nz/) ), RESHAPE( jiimag, (/Nt,Nz/) ) )
  Lwka = 2*(Nt+Nz) + INT(LOG(real(Nt))/LOG(2.)) + INT(LOG(real(Nz))/LOG(2.)) + 8
  SALLOCATE(wka,(1:Lwka),zero)
  call cfft2i( Nt, Nz, wka, Lwka, ic06fuf )
  Lwkb = 2 * Ntz
  SALLOCATE(wkb,(1:Lwkb),zero)
  call cfft2f( 1, Nt, Nz, jicomp(1:Nt,1:Nz), wka, Lwka, wkb, Lwkb, ic06fuf ) ! this is destructive
  DALLOCATE(wka,wkb)
  jireal = real( RESHAPE( jicomp, (/Ntz/) ) )
  jiimag = aimag( RESHAPE( jicomp, (/Ntz/) ) )
#endif

  isr = 'S' ! irrelevant; 02 Jul 13;
  
  jireal(1:Ntz) = jireal(1:Ntz) / sqrtoNtz
  jiimag(1:Ntz) = jiimag(1:Ntz) / sqrtoNtz
  
  even(1:mn) = zero
  oddd(1:mn) = zero
  cosi(1:mn) = zero
  sini(1:mn) = zero
  
  do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / Nfp
   
   if    ( mm.gt.0 .and. nn.gt.0 ) then
    
   !FATALMESS(ft02aa,1+(Nt-mm)+(   nn)*Nt.lt.1 .or. 1+(Nt-mm)+(   nn)*Nt.gt.Ntz, subscript error )
   !FATALMESS(ft02aa,1+(   mm)+(Nz-nn)*Nt.lt.1 .or. 1+(   mm)+(Nz-nn)*Nt.gt.Ntz, subscript error )

    even(ii) =   jireal(1+(Nt-mm)+(   nn)*Nt) + jireal(1+(   mm)+(Nz-nn)*Nt)
    oddd(ii) =   jiimag(1+(Nt-mm)+(   nn)*Nt) - jiimag(1+(   mm)+(Nz-nn)*Nt)
    cosi(ii) =   jiimag(1+(Nt-mm)+(   nn)*Nt) + jiimag(1+(   mm)+(Nz-nn)*Nt)
    sini(ii) = - jireal(1+(Nt-mm)+(   nn)*Nt) + jireal(1+(   mm)+(Nz-nn)*Nt) 
    
   elseif( mm.gt.0 .and. nn.eq.0 ) then
    
   !FATALMESS(ft02aa, 1+(Nt-mm)           .lt.1 .or. 1+(Nt-mm)           .gt.Ntz, subscript error )
   !FATALMESS(ft02aa, 1+(   mm)           .lt.1 .or. 1+(   mm)           .gt.Ntz, subscript error )

    even(ii) =   jireal(1+(Nt-mm)           ) + jireal(1+(   mm)           )
    oddd(ii) =   jiimag(1+(Nt-mm)           ) - jiimag(1+(   mm)           )
    cosi(ii) =   jiimag(1+(Nt-mm)           ) + jiimag(1+(   mm)           )
    sini(ii) = - jireal(1+(Nt-mm)           ) + jireal(1+(   mm)           )
    
   elseif( mm.gt.0 .and. nn.lt.0 ) then
    
   !FATALMESS(ft02aa, 1+(Nt-mm)+(Nz+nn)*Nt.lt.1 .or. 1+(Nt-mm)+(Nz+nn)*Nt.gt.Ntz, subscript error )
   !FATALMESS(ft02aa, 1+(   mm)+(  -nn)*Nt.lt.1 .or. 1+(   mm)+(  -nn)*Nt.gt.Ntz, subscript error )
    
    even(ii) =   jireal(1+(Nt-mm)+(Nz+nn)*Nt) + jireal(1+(   mm)+(  -nn)*Nt)
    oddd(ii) =   jiimag(1+(Nt-mm)+(Nz+nn)*Nt) - jiimag(1+(   mm)+(  -nn)*Nt)
    cosi(ii) =   jiimag(1+(Nt-mm)+(Nz+nn)*Nt) + jiimag(1+(   mm)+(  -nn)*Nt)
    sini(ii) = - jireal(1+(Nt-mm)+(Nz+nn)*Nt) + jireal(1+(   mm)+(  -nn)*Nt) 
    
   elseif( mm.eq.0 .and. nn.gt.0 ) then
    
   !FATALMESS(ft02aa, 1+        (Nz-nn)*Nt.lt.1 .or. 1+        (Nz-nn)*Nt.gt.Ntz, subscript error )
   !FATALMESS(ft02aa, 1+        (   nn)*Nt.lt.1 .or. 1+        (   nn)*Nt.gt.Ntz, subscript error )
    
    even(ii) =   jireal(1+        (Nz-nn)*Nt) + jireal(1+        (   nn)*Nt)
    oddd(ii) = - jiimag(1+        (Nz-nn)*Nt) + jiimag(1+        (   nn)*Nt)
    cosi(ii) =   jiimag(1+        (Nz-nn)*Nt) + jiimag(1+        (   nn)*Nt)
    sini(ii) =   jireal(1+        (Nz-nn)*Nt) - jireal(1+        (   nn)*Nt)
    
   elseif( mm.eq.0 .and. nn.eq.0 ) then
    
    even(ii) =   jireal(1                   )
    oddd(ii) =   zero
    cosi(ii) =   jiimag(1                   )
    sini(ii) =   zero
    
   endif
   
  enddo ! end of do ii = 1, mn; identification of modes in concise format; 02 Jul 13;
  
  return

end subroutine ft02aa

subroutine ft01aa( Nfp, mn , im , in , efmn , ofmn , cfmn , sfmn , Nt , Nz , ijreal , ijimag ) ! isr , trigm , trign , trigwk )

  use oculus, only    : zero, half, one

  implicit none
  INTEGER  , intent(in)    :: Nfp, mn, im(mn), in(mn)
  REAL     , intent(in)    :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn)
  INTEGER  , intent(in)    :: Nt, Nz
  REAL     , intent(out)   :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;

  CHARACTER                :: isr
  REAL                     :: trigm(2*Nt), trign(2*Nz), trigwk(2*Nt*Nz)
  
  INTEGER                  :: Ntz, imn, jj, kk, c06fuffail, c06gcffail, mm, nn

#ifndef HAVE_NAG
  INTEGER                  :: astat, ii, Lwka, Lwkb
  REAL     , allocatable   :: wka(:), wkb(:)
  COMPLEX                  :: ijcomp(1:Nt,1:Nz)
#endif
  
  isr = 'I' ; Ntz = Nt * Nz ; ijreal(1:Ntz) = zero ; ijimag(1:Ntz) = zero

  do imn = 1, mn ; mm = im(imn) ; nn = in(imn) / Nfp
   
   if    ( mm.gt.0 .and. nn.gt.0 ) then
    
    ijreal(1+(Nt-mm)+(   nn)*Nt) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)+(Nz-nn)*Nt) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)+(   nn)*Nt) = (   ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)+(Nz-nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.gt.0 .and. nn.eq.0 ) then
    
    ijreal(1+(Nt-mm)           ) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)           ) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)           ) = ( + ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)           ) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.gt.0 .and. nn.lt.0 ) then
    
    ijreal(1+(Nt-mm)+(Nz+nn)*Nt) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)+(  -nn)*Nt) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)+(Nz+nn)*Nt) = ( + ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)+(  -nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.eq.0 .and. nn.gt.0 ) then
    
    ijreal(1+        (Nz-nn)*Nt) = ( + efmn(imn) + sfmn(imn) ) * half
    ijreal(1+        (   nn)*Nt) = ( + efmn(imn) - sfmn(imn) ) * half
    
    ijimag(1+        (Nz-nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+        (   nn)*Nt) = ( + ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.eq.0 .and. nn.eq.0 ) then
    
    ijreal(1) = efmn(imn)
    ijimag(1) = cfmn(imn)
    
   endif
   
  enddo

  ijreal(1:Ntz) = ijreal(1:Ntz) * sqrt(one*Ntz)
  ijimag(1:Ntz) = ijimag(1:Ntz) * sqrt(one*Ntz)

#ifdef HAVE_NAG
  c06gcffail=0
  call C06GCF( ijimag, Ntz , c06gcffail ) ! NAG; 13 Oct 15;
  c06fuffail=0
  call C06FUF( Nt , Nz , ijreal , ijimag , isr , trigm , trign , trigwk , c06fuffail ) ! NAG; 13 Oct 15;
  c06gcffail=0
  call C06GCF( ijimag , Ntz , c06gcffail ) ! NAG; 13 Oct 15;
#else
  ijcomp = CMPLX( RESHAPE( ijreal, (/Nt,Nz/) ), RESHAPE( -1*ijimag, (/Nt,Nz/) ) )
  Lwka = 2*(Nt+Nz) + INT(LOG(real(Nt))/LOG(2.)) + INT(LOG(real(Nz))/LOG(2.)) + 8
  SALLOCATE(wka,(1:Lwka),zero)
  call cfft2i( Nt, Nz, wka, Lwka, c06gcffail )
  Lwkb = 2 * Ntz
  SALLOCATE(wkb,(1:Lwkb),zero)
  call cfft2f( 1, Nt, Nz, ijcomp(1:Nt,1:Nz), wka, Lwka, wkb, Lwkb, c06gcffail ) ! this is destructive
  DALLOCATE(wka,wkb)
  ijreal = real( RESHAPE( ijcomp, (/Ntz/) ) )
  ijimag = -1 * aimag( RESHAPE( ijcomp, (/Ntz/) ) )
#endif

  return

end subroutine ft01aa
