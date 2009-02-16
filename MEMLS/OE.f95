      Program OE
!     Geophysical parameters with error bar and information content
!     Sensitivity of satellite retrieved brightness temperature at select channels to snow cover on sea ice
!     2D-array(i,j)=(row,column)
! -------------------------------------------------------------------------------
!     Geophysical parameters
!     Variables:
!         Array Regime/Source   Property                        Unit         Constraint
!         y(1)  Air/Reanalysis  2 metre temperature             [K]          Displaces atmospheric profile
!         y(2)                  Total column water vapour       [mm = kg/m2]   -"-
!         y(3)  Snow/Memls      Thickness                       [cm]         Common to numerical layers
!         y(4)                  Temperature                     [K]
!         y(5)                  Density                         [kg/m3]
!         y(6)                  Correlation length              [mm]         Common to numerical layers
!         y(7)  Ice/Memls       Thickness                       [cm]
!         y(8)                  Temperature                     [K]
!         y(9)                  Salinity                        [per mil]    (Model initiation FY = 15, MY top/below = 1/2.5)
!     Prescribed: As in a priori thermodynamic model
!               Snow            Water content                   [%]          Fix: 0 (true until melt in June)
!                               Salinity                        [per mil]    Fix: 0 
!! RASMUS: freeboard OK? Ice    Density                         [kg/m3]      Fix: Only MY above the waterline is different
!                               Water content (volume of brine) [%]          Formula: y(8,9)
!                               Correlation length              [mm]         FY = formula (brine)
!     Constraints
!        Type of ice            FY or MY
!        Ice & snow temperature Air-sea thermodynamic equilibrium (0) or piecewise linear to free interface (1)
!                               Integer number gives the degree of freedoms that replaces y(4,8)
!        Range of validity      Emissivity = [0:1]
!                               2 metre temperature = [-43:-2]C to accomodate Memls
!        Deficiency             Exits after 8-10 iterations because Rttov is not closing all of its files
!     Main:
!     1)  Initialize
!     2)  Jacobian by forward difference
!     3)  Optimal estimation
!     Subroutines:
!     a)  Dsvdcmp               Singular value decomposition
!     b)  Inverse               Inverse of diagonal matrix
!     c)  Rttov_ini             Rttov parameters for instruments selected       (max Nin0)
!     d)  Conform               Convert to Rttov mandatory levels and humidity presentation
!     e)  Conform2              Update the atmospheric profile
!     e)  Memls0                Memls surface emissivity
!     f)  Rttov_fwd             Simulated brightness temperature
!     Library and modules:
!     e)  Rttov8.7
! -------------------------------------------------------------------------------
!     Input files:
!     1   Input.txt             Listing of directory paths
!         Input-Ch.txt          Listing of channels to be included
!         Input-Rttov87.txt     Rttov mandatory levels and limits of validity
!         ifil(4)               Data: atmospheric profile mean temperature [K], humidity [kg/kg]
!                                     a priori covariance and mean
!                               Measurement: covariance
!     Output files:
!     4   Syntes.txt            Interface Memls
!         Input/Memls/Memls-out.txt
!     1,2 fil                   List frequencies and zenith angles for selected instruments
!     2   ofil                  Optimal estimation: Tabularize parameter sensitivity
! -------------------------------------------------------------------------------
!     Version |    Date  | Code     | Author
! -------------------------------------------------------------------------------
!         1   | Jul 2007 | Original | T. M. Schroeder
! -------------------------------------------------------------------------------
      Implicit None
      Character(len=29),dimension(2)      :: fil=(/'Input/Memls/Freq-Memls-in.txt','Input/Memls/Angl-Memls-in.txt'/) ! [GHz, deg]
      Character(len=2)                    :: ice                                 ! Type of ice
      Character(len=100)                  :: ifil(4),ofil,idir,string            ! Directory and file
      Integer,parameter                   :: Npar=9,Nin0=10,Nm=43                ! # Parameters, max.instruments, mandatory levels Rttov
! --- Input
      Integer,parameter                   :: Nch0=35,Nle=40                      ! Rttov max. # channels and actual # profile levels
      Character(len=8),dimension(Nin0)    :: sen=' ',sat=' '                     ! Sensor, platform, polarization
      Character(len=1),dimension(Nin0*Nch0) :: pol='X'
      Integer                             :: scene,Nin,Nle0,i,j,k,l,n,cat        ! Temp.profile, # Instruments and levels in atm.profile
      Integer,dimension(Nin0)             :: num=0                               ! Satellite
      Integer,dimension(Nin0*Nch0)        :: Ch
      Real(kind=8)                        :: Saz,SzeAMSU,fac,eps,x(2)            ! Angle azimuth [0:360] and AMSU zenith [0:40]
! --- Atmosphere
      Real(kind=8),dimension(Nle)         :: P,T,q                               ! Mean atmospheric profile [hPa,K,kg/kg]
      Real(kind=8),dimension(Nm)          :: Pm,Tm0,Tm,mwm0,mwm                  ! Predefined levels, at these levels temperature and water vapour volume mixing ratio [hPa,K,ppmv]
      Real(kind=8),dimension(2,Nm)        :: T0,q0                               ! Limits of validity
! --- Rttov_ini
      Character(len=200),dimension(Nin0)  :: rfil                                ! Rttov coefficient file
      Integer                             :: Nch1,Nch2
      Integer,dimension(Nin0)             :: Nch,Nch3
      Real(kind=8),dimension(Nin0)        :: Szem                                ! Zenith angle [0:90], frequency
      Real(kind=8),dimension(Nin0,Nch0)   :: Sze,frq
      Real(kind=8),dimension(Nin0*Nch0)   :: frq0,frq3,Se,Seinv                  ! Measurement: variance diagonal and inverse
      Real(kind=8),dimension(Nin0*Nch0,Nin0*Nch0) :: Seinv2
! --- Memls and Rttov
      Real(kind=8),dimension(Nin0)        :: Teff                                ! Effective temperature lowest frequency mixed pol
      Real(kind=8),dimension(Nin0,Nch0)   :: emi,Tb,Tb0,Tb1                      ! Emissivity (1=H-pol, 2=V-pol), brightness temperature [K]
! --- Optimal estimation
      Real(kind=8)                        :: H                                   ! Information content
      Real(kind=8),dimension(Npar)        :: y0,y,y1,yJ,stdev,Saw,Sawinv         ! Data (Geophysical parameters): a priori, Jacobian and inverse
      Real(kind=8),dimension(Npar,Npar)   :: Sa,Sau,Sav,Sainv2,J2,Jinv2,IA
      Real(kind=8),dimension(Nin0*Nch0)   :: dT
      Real(kind=8),dimension(Nin0*Nch0,Npar):: M,J1                              ! Jacobian
      eps=Dsqrt(1d-6)
  100 Format (i2,1x,e15.6,2(9(1x,e15.6),5x))
! --- End of header -------------------------------------------------------------

! -------------------------------------------------------------------------------
!     1. Initialize
! -------------------------------------------------------------------------------
!     Load instruments
      Open                (unit=1,file='Input/Input.txt',status='old')
      Read                (1,*)                        scene
      Read                (1,*)                        Nin
      Read                (1,'(10(a8))')               (sen(i),i=1,Nin)
      Read                (1,'(18(a8))')               (sat(i),i=1,Nin)
      Read                (1,*)                        (num(i),i=1,Nin)
      Read                (1,*)                        Saz,SzeAMSU
      Read                (1,*)                        fac
      Do i=1,4
         Read             (1,'(a100)')                 ifil(i)
      End do
      Read                (1,'(a100)')                 idir
      Read                (1,'(a100)')                 ofil
      cat=                Index(idir,' ',back=.false.)-1
      Close               (unit=1)
! --- Fetch channels and their zenith angle
      Call Rttov_ini      (fil,idir,cat,sen,sat,num,rfil,Nin0,Nin,Nch0,Nch,SzeAMSU,Sze,Szem,frq,frq0,Nch1)

! --- AMSR contains V- as well as H-pol for all frequencies
!     Modify channels in printing array accordingly
      Nch2=0
      Nch3=0
      k=1
      Do i=1,Nin
         If (sen(i) == 'amsr') then
            Do j=1,Nch(i)
               frq3(k)=frq(i,j)
               frq3(k+1)=frq(i,j)
               pol(k)='V'
               pol(k+1)='H'
               k=k+2
            End do
            Nch3(i)=2*Nch(i)
         Else
            Do j=1,Nch(i)
               frq3(k)=frq(i,j)
               k=k+1
            End do
            Nch3(i)=Nch(i)
         End if
         Nch2=Nch2+Nch3(i)
      End do
! --- Channels to be sampled
      Ch=0
      Open                (unit=1,file='Input/Rttov/Input-Ch.txt',status='old')
      Read                (1,*)                        (Ch(i),i=1,Nch2)
      Close               (unit=1)

! --- Assemble input
      Nle0=0
      Do i=1,4
         Open             (unit=1,file=ifil(i),status='old')
         Read             (1,*)                        string
         If (i < 3) then
! ---       Data: mean atmospheric profile temperature [K] and humidity [kg/kg]. Error bar left unused
            Do j=1,1000
               If (i == 1) then
                  Read    (1,*,End=10)                 P(j),T(j)
                  Nle0=Nle0+1
               Else
                  Read    (1,*,End=10)                 P(j),q(j)
               End if
            End do
         Else if (i == 3) then
! ---       Data: a priori mean and covariance, which has to be mirrored symmetrically into its upper triangle
!           File name is ending in 'FY.txt' or 'MY.txt' as indication of the type of ice
            cat=          Index(ifil(i),'.',back=.false.)-2
            Read          (ifil(i)(cat:cat+1),*)       ice
            Read          (1,*)                        string,y0
            Read          (1,*)                        string,stdev
            stdev=fac*stdev
            Do j=1,Npar
               Read       (1,*)                        string,(Sa(j,k),k=1,j)
               Do k=1,j
                  Sa(j,k)=Sa(j,k)*stdev(j)*stdev(k)
               End do
            End do
            Do j=1,Npar
               Do k=j+1,Npar
                  Sa(j,k)=Sa(k,j)
               End do
            End do
! ---       SVD-decomposition and inversion
!           Inverse of decomposed matrix Ainv = (u*w*vT)inv = v*winv*uT
            Call Dsvdcmp  (Sa,Npar,Npar,Npar,Npar,Saw,Sav)
            Sau=Sa
            Call Inverse  (Saw,Sawinv,Npar,Npar)
            Sainv2=0d0
            Do j=1,Npar
               Sainv2(j,j)=Sawinv(j)
            End do
            Sainv2=Matmul(Sainv2,Transpose(Sau))
            Sainv2=Matmul(Sav,Sainv2)
         Else
! ---       Measurement: variance (sensitivity and accuracy)
!           Sensors sorted according to Input.txt
            Do j=1,Nch1
               Read       (1,*)                        string,string,x
               Se(j)=x(1)**2d0+x(2)**2d0
            End do
! ---       Inversion
            Seinv2=0d0
            Call Inverse  (Se,Seinv,Nin0*Nch0,Nch1)
            Do j=1,Nch1
               Seinv2(j,j)=Seinv(j)
            End do
         End if
   10    Continue
         Close            (unit=1)
      End do

! --- Atmosphere
!     Rttov levels and range of validity
      Open                (unit=1,file='Input/Rttov/Input-Rttov87.txt',status='old')
      Read                (1,*)                        string
      Do i=1,Nm
         Read             (1,*)                        Pm(i),(T0(j,i),j=1,2),(q0(j,i),j=1,2)
      End do
      Close               (unit=1)
! --- Comply with Rttov mandatory levels. Note the lowermost is not the 2-metre temperature
      Call Conform        (P,q,T,Nle0,Pm,T0,Tm0,q0,mwm0,Nm)

! -------------------------------------------------------------------------------
!     2. Jacobian by forward difference
! -------------------------------------------------------------------------------
      Open                (unit=1,file=ofil,status='new')
      Write               (1,'(a10,50(1x,f6.2))')      'Frq[GHz]',((frq(k,j),j=1,Nch(k)),k=1,Nin)
      Write               (1,*)                        'Iteration#  Information  Parameter     Error bar'
! --- Initial condition set to a priori mean value
      y=y0
      n=0
      H=1d0
      Write               (1,100)                      n,H,y0,stdev

   11 Continue
      n=n+1
      M=0d0
      Do i=1,Npar+1
! ---    Abscissa
         If (i == 1) then
! ---       Reference call
            yJ=y
         Else
! ---       Nudging one parameter at a time with trick to reduce the finite precision error
            If (y(i-1) /= 0d0) then
               yJ(i-1)=y(i-1)*(1d0+eps)
            Else
               yJ(i-1)=eps
            End if
! ---       Reset all the other parameters
            Do j=1,Npar
               If (j /= i-1)                           yJ(j)=y(j)
            End do
         End if

! ---    Physical model
!        Surface emissivity for effective temperature taken at a single frequency (*)
         Call Memls0      (fil,frq0,Nin0,Nin,Nch0,Nch,Nch1,Npar,scene,ice,yJ,emi,Teff)
! ---    Scene of required thermodynamics, not taken into account a priori, adjusts the snow & ice temperature--also initially
         If (i == 1) then
            y(4)=yJ(4)
            y(8)=yJ(8)
            If (n == 1) then
               y0(4)=yJ(4)
               y0(8)=yJ(8)
            End if
         End if
! ---    Constrain the atmospheric profile according to y(1-2)
         Call Conform2    (y0,yJ,Npar,Tm0,Tm,mwm0,mwm,Nm)
!Write (1,'(a5,3(3(a15,f6.2),10x))') 'A','Tair: 1000hPa=',T(Nle0),' 1013.25hPa-J=',Tm(Nm),' 0.1hPa-J=',Tm(1), &
!& ' T2m: y0=',y0(1),' 2m-y=',y(1),' 2m-yJ=',yJ(1), &
!&' Tsnow: y0=',y0(4),' y=',y(4),' yJ=',yJ(4)
! ---    Brightness temperature as gauged by the suite of satellites
         Tb=0d0
!Write (1,*) 'Teff=Tbmax ',(Teff(j),j=1,Nin)
         Do j=1,Nin
!Write (1,'(a15,i2,a8,50(x,f8.3))') 'Emissivity',i,sen(j),(emi(j,k),k=1,2*Nch(j))
            Call Rttov_fwd(j,rfil,yJ(1),Pm(Nm),Nin0,Nin,Nch0,Nch,Saz,Szem,Tm,mwm,Nm,emi,Teff,Tb)
!Write (1,'(a15,i2,a8,50(x,f8.3))') 'Tb ',i,sen(j),(Tb(j,k),k=1,Nch3(j))
         End do

! ---    Compute Jacobian
         If (i == 1) then
            y1=y
            Tb1=Tb
            If (n == 1)                                Tb0=Tb
         Else
            l=0
            Do j=1,Nin
               If (j > 1)                              l=l+Nch(j-1)
! ---          Sampled channels, when nudging of y(4) is not altered to zero by equilibrium (*)
               If (yJ(i-1) /= y1(i-1)) then
                  Do k=1,Nch3(j)
                     If (Ch(l+k) > 0)                  M(l+k,i-1)=(Tb(j,k)-Tb1(j,k))/(yJ(i-1)-y1(i-1))
                  End do
               End if
            End do
         End if
      End do
!     Do i=1,Nch2
!        If (Ch(i) > 0) Write (1,'(a15,f8.3,a1,a15,10(e12.4,1x))')  'Channel[GHz]',frq3(i),pol(i),'dTb/dparameter',(M(i,j),j=1,Npar)
!     End do

! -------------------------------------------------------------------------------
!     3. Optimal estimation
! -------------------------------------------------------------------------------
!     Covariance of estimated parameters is inverse of A = (u*w*vT)inv
!     Has verified result and Ainv*A = E and uT*u = vT*v = E
!     A priori and measurement covariances result in u = v because of symmetry
      J1=0d0
      J1=Matmul(Seinv2,M)
      J2=Matmul(Transpose(M),J1)
      J2=Sainv2+J2
! --- SVD-decomposition and inversion
      Call Dsvdcmp        (J2,Npar,Npar,Npar,Npar,Saw,Sav)
      Sau=J2
      Call Inverse        (Saw,Sawinv,Npar,Npar)
      Jinv2=0d0
      Do j=1,Npar
         Jinv2(j,j)=Sawinv(j)
      End do
      Jinv2=Matmul(Jinv2,Transpose(Sau))
      Jinv2=Matmul(Sav,Jinv2)
! --- Output error bar
      Do i=1,Npar
         stdev(i)=Dsqrt(Jinv2(i,i))
      End do
! --- Information content on parameters considered
      IA=0d0
      IA=Matmul(Jinv2,Sainv2)
      Call Information    (Npar,IA,H)

! --- Perform iteration for the sampled channels [Rodgers, 2004: Eq.(5-8)]
      dT=0d0
      l=0
      Do j=1,Nin
         If (j > 1)                                    l=l+Nch(j-1)
         Do k=1,Nch3(j)
            If (Ch(l+k) > 0)                           dT(l+k)=Tb0(j,k)-Tb(j,k)
         End do
      End do
      dT=Matmul(Seinv2,dT)
      yJ=Matmul(Transpose(M),dT)
      yJ=yJ-Matmul(Sainv2,y-y0)
      y=y+Matmul(Jinv2,yJ)

! --- Output
      Write               (1,100)                      n,H,y,stdev
!! Test for convergence
If (n < 50) Go to 11
      Close               (unit=1)
      Stop
      End