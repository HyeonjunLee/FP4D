!> This module contains functions for evaluating elliptic integrals
module elliptics_mod
  implicit none

  contains

    !> Evaluates elliptic integrals K and E
    subroutine elliptics(k,ellip_K,ellip_E,vpic_ierr)
      implicit none
!$acc routine seq
      real(kind=8) :: pi = 3.14159265358979D0
      real(kind=8) :: k, ellip_K, ellip_E
      integer :: vpic_ierr, k_index
      real(kind=8) :: ell_B, ell_D, BX_star, B0_star, DX_star, D0_star, X, mc, kc
      integer :: iter, Bnp1, Dnp1
      real (kind=8), dimension(12) :: &
         B0to1 = (/ 0.790401413584395132,  0.102006266220019155, 0.039878395558551461, &
               0.021737136375982167,  0.013960979767622058, 0.009892518822669142, &
               0.007484612400663336,  0.005934625664295474, 0.004874249053581664, &
               0.004114606930310886,  0.003550452989196177, 0.003119229959988475 /), &
         D0to1 = (/ 0.800602040206397048,  0.313994477771767757, 0.205913118705551955, &
               0.157744346538923994,  0.130595077319933092, 0.113308474489758567, &
               0.101454199173630195,  0.092918784207297437, 0.086565380148168087, &
               0.081727984665103014,  0.077990665729107038, 0.075080426851268007 /), &
         B1to2 = (/ 0.801024064452844894,  0.110695344529634015, 0.047348746716993718, &
               0.028484367255041423,  0.020277811444003597, 0.015965005853099119, &
               0.013441320273553635,  0.011871565736951440, 0.010868363672485521, &
               0.010231587232710565,  0.009849585546666211, 0.009656606347153765 /), &
         D1to2 = (/ 0.834232667811735098,  0.360495281619098276, 0.262379664114505869, &
               0.223723944518094276,  0.206447811775681053, 0.199809440876486856, &
               0.199667451603795275,  0.204157558868236842, 0.212387467960572375, &
               0.223948914061499360,  0.238708097425597860, 0.256707203545463756 /)
      real (kind=8), dimension(13) :: &
         B2to3 = (/ 0.812597772919920493,  0.121109617945510113, 0.057293376831239877, &
               0.038509451602167328,  0.030783430301775233, 0.027290564934732527, &
               0.025916369289445199,  0.025847203343361799, 0.026740923539348855, &
               0.028464314554825705,  0.030995446237278954, 0.034384369179940976, &
               0.038738002072493936 /), &
         D2to3 = (/ 0.873152581892675550,  0.420622230667770216, 0.344231061559450379, &
               0.331133021818721762,  0.345277285052808412, 0.377945322150393392, &
               0.427378012464553881,  0.494671744307822406, 0.582685115665646201, &
               0.695799207728083165,  0.840018401472533403, 1.023268503573606061, &
               1.255859085136282496 /), &
         B3to4 = (/ 0.825323557983515895, 0.133862116083687790, 0.071011293597988675, &
                    0.054178477417387376, 0.049451744948102993, 0.050222196224107476, &
                    0.054742913171830353, 0.062746257927001699, 0.074669881043476886, &
                    0.091480845177733472, 0.114705092110997824, 0.146571132581439876, &
                    0.190257137333846284 /)
      real (kind=8), dimension(14) :: &
         D3to4 = (/ 0.919027039242097348, 0.501002159288247514, 0.468831270566456863, &
                    0.517714227776400015, 0.620843391317303107, 0.782364393786869723, &
                    1.019114535076102913, 1.359345202748496052, 1.845717302358827942, &
                    2.541071703153920729, 3.537404655208041337, 4.969296002977425930, &
                    7.033822870030031126, 10.02004322503447140 /)
      real (kind=8), dimension(13) :: &
         B4to5 = (/ 0.839479570270612971, 0.149916440306396336, 0.090831935819428835, &
                    0.080347033483341786, 0.085638440500470454, 0.101954725932990372, &
                    0.130574811533616015, 0.176105076358849928, 0.246835164402955447, &
                    0.356424476867718855, 0.527002562230102743, 0.794389634259304750, &
                    1.216762532429718021 /)
      real (kind=8), dimension(16) :: &
         D4to5 = (/ 0.974404366546369673, 0.613246805394160910, 0.671096669502166996, &
                    0.870727620185086140, 1.229542231202690761, 1.826605967544420569, &
                    2.806934530997762740, 4.418789329084028134, 7.083236057478765325, &
                    11.51508812055758294, 18.93151118599927464, 31.41199693820496388, &
                    52.52072945457582854, 88.38485473506529806, 149.5663744939804784, &
                    254.3179084310411743 /)
      real (kind=8), dimension(14) :: &
         B5to6 = (/ 0.855469615156419991, 0.170896072689739584, 0.121335229026948226, &
                    0.128201883574947410, 0.164687281451527560, 0.237418908749381742, &
                    0.369208104716495452, 0.605658733847927717, 1.033705561557812744, &
                    1.818988489363267885, 3.279377651273850938, 6.029888380717536331, &
                    11.26979685557794172, 21.35457785038283450 /)
      real (kind=8), dimension(17) :: &
         D5to6 = (/ 1.043455295115133534, 0.779625721928504850, 1.029742360932067582, &
                    1.622037223411353130, 2.787989531185347620, 5.048381487372069147, &
                    9.463277611943484295, 18.18148994942766790, 35.58098059117916870, &
                    70.63393546191445013, 141.8285800834330593, 287.4487512501321663, &
                    587.1153846499230762, 1207.065435225480616, 2495.588727248664223, &
                    5184.692429394806441, 10817.21333690413275 /)
      real (kind=8), dimension(16) :: &
         B6to7 = (/ 0.873920061848643136, 0.199814057482376946, 0.172769615878015213, &
                    0.228106913284202167, 0.370468141118071220, 0.679271252884820555, &
                    1.348008496681757302, 2.827670976853820704, 6.179468250123914084, &
                    13.93568601034281150, 32.21892928105972203, 76.00696295922610103, &
                    182.3214490877540696, 443.5150764411264816, 1091.854722902838829, &
                    2715.765866403819588 /)
      real (kind=8), dimension(18) :: &
         D6to7 = (/ 1.133678336575733166, 1.048643173729970391, 1.753465041198464516, &
                    3.523182726803385513, 7.749476413813974582, 17.98645005585073306, &
                    43.25591634623261333, 106.6815344540960170, 268.0984865731174340, &
                    683.6241148502898048, 1763.497085219187407, 4592.374753831163809, &
                    12053.44101904888928, 31846.66302074208170, 84621.22135905680802, &
                    225956.4231829078900, 605941.5172817588600, 1631082.599539268321 /)
      real (kind=8), dimension(19) :: &
         B7to8 = (/ 0.895902820924731621, 0.243140003766786662, 0.273081875594105532, &
                    0.486280007533573324, 1.082747437228230918, 2.743445290986452500, &
                    7.555817828670234627, 22.05194082493752427, 67.15640644740229408, &
                    211.2722537881770962, 681.9037843053270682, 2246.956231592536517, &
                    7531.483865999711792, 25608.51260130241579, 88140.74740089604971, &
                    306564.4242098446591, 1076036.077811072194, 3807218.502573632648, &
                    13566382.24422139551 /)
      real (kind=8), dimension(21) :: &
         D7to8 = (/ 1.260612826574911614, 1.548665638082676581, 3.553669411871607615, &
                    9.900444676104398756, 30.32056661745247199, 98.18025865888308915, &
                    329.7710104345570550, 1136.655989742890393, 3993.834335746229798, &
                    14242.72958655527085, 51394.75729168872096, 187246.7029146231521, &
                    687653.0923753899027, 2542385.535653982270, 9453781.219347490272, &
                    35328363.01797091708, 132593262.3833930149, 499544968.1840548215, &
                    1888409347.294438724, 7160267534.478937192, 27223307946.96339622 /)
      real (kind=8), dimension(15) :: &
         B8to85 = (/ 0.915922052601931494, 0.294714252429483394, 0.435776709264636140, &
                     1.067328246493644239, 3.327844118563268085, 11.90406004445092906, &
                     46.47838820224626394, 192.7556002578809477, 835.3356299261900064, &
                     3743.124548343029103, 17219.07731004063094, 80904.60401669850158, &
                     386808.3292751742460, 1876487.670110449342, 9216559.908641567755 /)
      real (kind=8), dimension(18) :: &
         D8to85 = (/ 1.402200569110579095, 2.322205897861749447, 7.462158366466719683, &
                     29.43506890797307903, 128.1590924337895775, 591.0807036911982326, &
                     2830.546229607726377, 13917.76431889392230, 69786.10525163921228, &
                     355234.1420341879635, 1830019.186413931054, 9519610.812032515607, &
                     49920868.75574849454, 263567700.9826023474, 1399645765.120061119, &
                     7469935792.837635005, 40041555958.35610574, 215463066814.4966654 /)
      real (kind=8), dimension(19) :: &
         B85to9 = (/ 0.931906061029524828, 0.348448029538453861, 0.666809178846938248, &
                     2.210769135708128663, 9.491765048913406881, 47.09304791027740853, &
                     255.9200460211233087, 1480.029532675805408, 8954.040904734313578, &
                     56052.48220982686950, 360395.7241626000917, 2367539.415273216078, &
                     15829949.57277684102, 107415809.3278511100, 738058546.0239595692, &
                     5126022002.555101497, 35935340655.02416589, 253988125761.2812212, &
                     1808180007145.359570 /)
      real (kind=8), dimension(21) :: &
         D85to9 = (/ 1.541690112721819084, 3.379176214579645449, 14.94058385670236672, &
                     81.91773929235074881, 497.4900546551479866, 3205.184010234846235, &
                     21457.32237355321926, 147557.0156564174712, 1035045.290185256525, &
                     7371922.334832212125, 53143443.95142401142, 386882347.5795976313, &
                     2839458401.528033778, 20982661229.43898942, 155961775401.7662418, &
                     1165096220419.884791, 8742012983013.913805, 65847254626723.66919, &
                     497679873706243.4393, 3773018634056605.405, 28682631948378196.60 /)
      real (kind=8) , dimension(14) :: &
         elliptic_BX = (/ 0D0, -0.25D0, -3.125D-2, -1.171875D-2, -6.103515625D-3, -3.7384033203125D-3, &
                        -2.5234222412109375D-3, -1.817464828491210938D-3, -1.371212303638458252D-3, &
                        -1.071259612217545509D-3, -8.599834109190851450D-4, -7.055772985040675849D-4, &
                        -5.893174027278291760D-4, -4.995976058381756957D-4 /), &
         elliptic_B0 = (/ 1D0, -0.25D0, 4.6875D-2, 2.34375D-2, 1.352945963541666667D-2, 8.740743001302083333D-3, &
                         6.0962677001953125D-3, 4.489111900329589844D-3, 3.441497071513107845D-3, &
                         2.721402519715151617D-3, 2.205478451662837336D-3, 1.823336289104075164D-3, &
                         1.532467036966436629D-3, 1.305981575356390531D-3 /)
      real (kind=8) , dimension(13) :: &
         elliptic_DX = (/ 0.5D0, -0.125D0, -2.34375D-2, -9.765625D-3, -5.340576171875D-3, -3.36456298828125D-3, &
                        -2.313137054443359375D-3, -1.687645912170410156D-3, -1.285511534661054611D-3, &
                        -1.011745189316570759D-3, -8.169842403731308877D-4, -6.735056031175190583D-4, &
                        -5.647625109475029603D-4 /), &
         elliptic_D0 = (/ -1D0, 0D0, 3.90625D-2, 2.018229166666666667D-2, 1.202901204427083333D-2, &
                         7.941436767578125000D-3, 5.623292922973632813D-3, 4.187006609780447824D-3, &
                         3.237116100665714060D-3, 2.576826204497751499D-3, 2.099504446134290894D-3, &
                         1.743372975543576159D-3, 1.470660484741195621D-3 /)
      integer, dimension(11) :: &
         elliptic_Bnum = (/ 11, 11, 12, 12, 12, 13, 15, 18, 14, 18, 13 /), &
         elliptic_Dnum = (/ 11, 11, 12, 13, 15, 16, 17, 20, 17, 20, 12 /)

      ell_B = 0
      ell_D = 0

      k_index = floor(k*10D0)
      select case(k_index)
        case(0)
            if( k .eq. 0) then
                ell_B = pi*0.25
                ell_D = pi*0.25
            else if((0D0 .lt. k) .and. (k .le. 0.1)) then
                Bnp1 = elliptic_Bnum(1)+1
                Dnp1 = elliptic_Dnum(1)+1
                ell_B = B0to1(Bnp1)
                ell_D = D0to1(Dnp1)
                kc = k-0.05D0
                do iter=1,elliptic_Bnum(1)
                    ell_B = ell_B*kc+B0to1(Bnp1-iter)
                enddo
                do iter=1,elliptic_Dnum(1)
                    ell_D = ell_D*kc+D0to1(Dnp1-iter)
                enddo
            endif
        case(1)
            Bnp1 = elliptic_Bnum(2)+1
            Dnp1 = elliptic_Dnum(2)+1
            ell_B = B1to2(Bnp1)
            ell_D = D1to2(Dnp1)
            kc = k-0.15D0
            do iter=1,elliptic_Bnum(2)
                ell_B = ell_B*kc+B1to2(Bnp1-iter)
            enddo
            do iter=1,elliptic_Dnum(2)
                ell_D = ell_D*kc+D1to2(Dnp1-iter)
            enddo
        case(2)
            Bnp1 = elliptic_Bnum(3)+1
            Dnp1 = elliptic_Dnum(3)+1
            ell_B = B2to3(Bnp1)
            ell_D = D2to3(Dnp1)
            kc = k-0.25D0
            do iter=1,elliptic_Bnum(3)
                ell_B = ell_B*kc+B2to3(Bnp1-iter)
            enddo
            do iter=1,elliptic_Dnum(3)
                ell_D = ell_D*kc+D2to3(Dnp1-iter)
            enddo
        case(3)
          Bnp1 = elliptic_Bnum(4)+1
          Dnp1 = elliptic_Dnum(4)+1
          ell_B = B3to4(Bnp1)
          ell_D = D3to4(Dnp1)
          kc = k-0.35D0
          do iter=1,elliptic_Bnum(4)
              ell_B = ell_B*kc+B3to4(Bnp1-iter)
          enddo
          do iter=1,elliptic_Dnum(4)
              ell_D = ell_D*kc+D3to4(Dnp1-iter)
          enddo
        case(4)
          Bnp1 = elliptic_Bnum(5)+1
          Dnp1 = elliptic_Dnum(5)+1
          ell_B = B4to5(Bnp1)
          ell_D = D4to5(Dnp1)
          kc = k-0.45D0
          do iter=1, elliptic_Bnum(5)
              ell_B = ell_B*kc+B4to5(Bnp1-iter)
          enddo
          do iter=1, elliptic_Dnum(5)
              ell_D = ell_D*kc+D4to5(Dnp1-iter)
          enddo
        case(5)
          Bnp1 = elliptic_Bnum(6)+1
          Dnp1 = elliptic_Dnum(6)+1
          ell_B = B5to6(Bnp1)
          ell_D = D5to6(Dnp1)
          kc = k-0.55D0
          do iter=1, elliptic_Bnum(6)
              ell_B = ell_B*kc+B5to6(Bnp1-iter)
          enddo
          do iter=1, elliptic_Dnum(6)
              ell_D = ell_D*kc+D5to6(Dnp1-iter)
          enddo
        case(6)
          Bnp1 = elliptic_Bnum(7)+1
          Dnp1 = elliptic_Dnum(7)+1
          ell_B = B6to7(Bnp1)
          ell_D = D6to7(Dnp1)
          kc = k-0.65D0
          do iter=1, elliptic_Bnum(7)
              ell_B = ell_B*kc+B6to7(Bnp1-iter)
          enddo
          do iter=1, elliptic_Dnum(7)
              ell_D = ell_D*kc+D6to7(Dnp1-iter)
          enddo
        case(7)
          Bnp1 = elliptic_Bnum(8)+1
          Dnp1 = elliptic_Dnum(8)+1
          ell_B = B7to8(Bnp1)
          ell_D = D7to8(Dnp1)
          kc = k-0.75D0
          do iter=1, elliptic_Bnum(8)
              ell_B = ell_B*kc+B7to8(Bnp1-iter)
          enddo
          do iter=1, elliptic_Dnum(8)
              ell_D = ell_D*kc+D7to8(Dnp1-iter)
          enddo
        case(8)
            if((0.8 .lt. k) .and. (k .le. 0.85)) then
                Bnp1 = elliptic_Bnum(9)+1
                Dnp1 = elliptic_Dnum(9)+1
                ell_B = B8to85(Bnp1)
                ell_D = D8to85(Dnp1)
                kc = k-0.825D0
                do iter=1, elliptic_Bnum(9)
                    ell_B = ell_B*kc+B8to85(Bnp1-iter)
                enddo
                do iter=1, elliptic_Dnum(9)
                    ell_D = ell_D*kc+D8to85(Dnp1-iter)
                enddo
            else if((0.85 .lt. k) .and. (k .le. 0.9)) then
                Bnp1 = elliptic_Bnum(10)+1
                Dnp1 = elliptic_Dnum(10)+1
                ell_B = B85to9(Bnp1)
                ell_D = D85to9(Dnp1)
                kc = k-0.875D0
                do iter=1, elliptic_Bnum(10)
                    ell_B = ell_B*kc+B85to9(Bnp1-iter)
                enddo
                do iter=1, elliptic_Dnum(10)
                    ell_D = ell_D*kc+D85to9(Dnp1-iter)
                enddo
            endif
        case(9)
          Bnp1 = elliptic_Bnum(11)+1
          Dnp1 = elliptic_Dnum(11)+1
          BX_star = elliptic_BX(Bnp1)
          B0_star = elliptic_B0(Bnp1)
          DX_star = elliptic_DX(Dnp1)
          D0_star = elliptic_D0(Dnp1)
          mc = 1D0-k

          do iter=1,elliptic_Bnum(11)
              BX_star = BX_star*mc + elliptic_BX(Bnp1-iter)
              B0_star = B0_star*mc + elliptic_B0(Bnp1-iter)
          enddo
          do iter=1,elliptic_Dnum(11)
              DX_star = DX_star*mc + elliptic_DX(Dnp1-iter)
              D0_star = D0_star*mc + elliptic_D0(Dnp1-iter)
          enddo
          X = -log(mc/16)
          ell_B = (B0_star+BX_star*X)/k
          ell_D = (D0_star+DX_star*X)/k
        case(10)
          ell_B = 0D0
          ell_D = 0D0
          vpic_ierr = 2
        case default
#ifndef _OPENACC
          print *, 'INPUT for elliptic integral is out of range :',k
#endif
          vpic_ierr = 2
      end select

      ellip_K = ell_B+ell_D
      ellip_E = ell_B+(1D0-k)*ell_D
    end subroutine elliptics

    !> ellip_agm & ellip_agm_v
    !! To calculate complete elliptic integrals
    !! Using arithmetic-geometric mean (agm) and modified agm (magm)
    !! ellip_agm by E.S. Yoon - Mar. 09, 2015 : original code
    !! ellip_agm_v by Nathan Wichmann - June, 2015 : vectorized ellip_agm
    !! Note that ellip_agm_v is NOT for OpenACC
    subroutine ellip_agm_v(k, ellip_K, ellip_E, vpic_ierr,n)
      implicit none
      real(kind=8) :: pi = 3.14159265358979D0
      real(kind=8), intent(in),dimension(n)  :: k
      real(kind=8), intent(out),dimension(n) :: ellip_K, ellip_E
      !integer, intent(out) :: vpic_ierr
      integer :: vpic_ierr
      real (kind=8),dimension(n) :: x_n, y_n
      real (kind=8),dimension(n) :: Mx_n, My_n, Mz_n
      real (kind=8) :: x_np1
      real (kind=8) :: Mx_np1, My_np1, rd
      real (kind=8) :: e_tol
      integer :: order
      integer :: n_order
      integer :: n,i,iibeg,iiend
      integer,parameter::vl=8

      n_order = 10
      e_tol = 1e-8

      do iibeg=1,n,vl
        iiend = min(iibeg+vl-1,n)

        do i=iibeg,iiend
            Mx_n(i) = 1D0
            My_n(i) = 1D0 - k(i)    !beta^2
            Mz_n(i) = 0D0
            x_n(i) = 1D0
            y_n(i) = sqrt(My_n(i))  !beta
        enddo

        do order=1, n_order
            !!agm
            do i=iibeg,iiend
                x_np1 = (x_n(i)+y_n(i))*0.5D0
                y_n(i) = sqrt(x_n(i)*y_n(i))
                x_n(i) = x_np1 !update results
                !!magm
                Mx_np1 = (Mx_n(i)+My_n(i))*0.5D0
                rd = sqrt( (Mx_n(i)-Mz_n(i))*(My_n(i)-Mz_n(i)) )
                My_np1 = Mz_n(i)+rd
                Mz_n(i) = Mz_n(i)-rd
                Mx_n(i) = Mx_np1
                My_n(i) = My_np1
            enddo

            if(mod(order,2)==0.and.order>1)then
                if( all((dabs(x_n(iibeg:iiend)-y_n(iibeg:iiend)) .lt. e_tol) .and. &
                        (dabs(Mx_n(iibeg:iiend)-My_n(iibeg:iiend)) .lt. e_tol)) ) exit
            endif
        enddo

        do i=iibeg,iiend
            ellip_K(i) = 0.5D0*pi/x_n(i)
            ellip_E(i) = Mx_n(i)*ellip_K(i)
        enddo
      enddo

      return

    end subroutine

    !> ellip_agm & ellip_agm_v
    !! To calculate complete elliptic integrals
    !! Using arithmetic-geometric mean (agm) and modified agm (magm)
    !! ellip_agm by E.S. Yoon - Mar. 09, 2015 : original code
    !! ellip_agm_v by Nathan Wichmann - June, 2015 : vectorized ellip_agm
    !! Note that ellip_agm_v is NOT for OpenACC
    subroutine ellip_agm(k, ellip_K, ellip_E, vpic_ierr)
      implicit none
!$acc routine seq
#ifdef _OPENACC
      real(kind=8), value  :: k
#else
      real(kind=8), intent(in)  :: k
#endif
      real(kind=8) :: pi = 3.14159265358979D0
      real(kind=8), intent(out) :: ellip_K, ellip_E
      integer, intent(out) :: vpic_ierr
      real (kind=8) :: x_n, y_n
      real (kind=8) :: Mx_n, My_n, Mz_n
      real (kind=8) :: x_np1
      real (kind=8) :: Mx_np1, My_np1, rd
      real (kind=8) :: e_tol
      integer :: order
      integer :: n_order

      n_order = 10
      e_tol = 1e-8

      Mx_n = 1D0
      My_n = 1D0 - k    !beta^2
      Mz_n = 0D0
      x_n = 1D0
      y_n = sqrt(My_n)  !beta

      do order=1, n_order
        !!agm
        x_np1 = (x_n+y_n)*0.5D0
        y_n = sqrt(x_n*y_n)
        x_n = x_np1 !update results
        !!magm
        Mx_np1 = (Mx_n+My_n)*0.5D0
        rd = sqrt( (Mx_n-Mz_n)*(My_n-Mz_n) )
        My_np1 = Mz_n+rd
        Mz_n = Mz_n-rd
        Mx_n = Mx_np1
        My_n = My_np1

!#ifndef _OPENACC
        if( (dabs(x_n-y_n) .lt. e_tol) .and. (dabs(Mx_n-My_n) .lt. e_tol) ) exit
!#endif
      enddo

      ellip_K = 0.5D0*pi/x_n
      ellip_E = Mx_n*ellip_K
      vpic_ierr = order

    end subroutine ellip_agm

end module elliptics_mod
