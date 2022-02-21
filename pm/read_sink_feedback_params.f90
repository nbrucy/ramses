subroutine read_stellar_params()
  use cooling_module, only: mH
  use amr_commons, only: dp, myid
  use pm_commons, only: iseed
  use amr_parameters, only:stellar

  use sink_feedback_parameters
  implicit none

  !------------------------------------------------------------------------
  ! Read stellar object related parameters and perform some 'sanity checks'
  !------------------------------------------------------------------------
  namelist/stellar_params/ nstellarmax, stellar_msink_th, sn_direct, &
                         & imf_index, imf_low, imf_high, &
                         & lt_t0, lt_m0, lt_a, lt_b, &
                         & stf_K, stf_m0, stf_a, stf_b, stf_c, &
                         & hii_w, hii_alpha, hii_c, hii_t, hii_T2 , &
                         & sn_feedback_sink,make_stellar_glob,iseed, &
                         & sn_feedback_cr,fcr, &
                         & mstellarini, &
                         & Tsat, Vsat, sn_r_sat, sn_p, sn_e, sn_mass, &
                         & FB_nsource, FB_on, FB_start, FB_end, FB_sourcetype, &
                         & FB_pos_x, FB_pos_y, FB_pos_z, &
                         & FB_mejecta, FB_energy, FB_thermal, &
                         & FB_radius, FB_r_refine, Vdisp, &
                         & ssm_table_directory, use_ssm

  real(dp):: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp):: msun, Myr, km_s

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Initialise mstellarini (should be zero if not set in the namelist)
  mstellarini = 0d0
  
  ! Read namelist file 
  rewind(1)
  read(1, nml=stellar_params, end=111)
  rewind(1)

  if(nstellarmax <= 0) stellar = .false.

  if(.not. stellar) return

  if(imf_index >= -1.0d0) then
      if(myid == 1) write(*, *) 'imf_alpha should be lower than -1'
      call clean_stop
  end if

  if(imf_low <= 0.0d0 .or. imf_low >= imf_high) then
      if(myid == 1) write(*, *) '0 < imf_low < imf_high has to be respected'
      call clean_stop
  end if

  if(stellar_msink_th <= 0.0d0) then
      if(myid == 1) write(*, *) 'stellar_msink_th should be positive'
      call clean_stop
  end if

  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  ! Convert parameters to code units
  msun = 2d33 / scale_d / scale_l**3
  Myr = 1d6 * 365.25d0 * 86400d0 / scale_t
  km_s = 1d5 / scale_v

  imf_low = imf_low * msun
  imf_high = imf_high * msun
  lt_t0 = lt_t0 * Myr
  lt_m0 = lt_m0 * msun
  stellar_msink_th = stellar_msink_th * msun
  mstellarini = mstellarini * msun
  
  !Careful : convert the parameter for ionising flux in code units
  stf_K = stf_K * scale_t ! K is in s**(-1)
  stf_m0 = stf_m0 * msun 

  hii_alpha = hii_alpha / (scale_l**3 / scale_t) ! alpha is in cm**3 / s
  hii_c = hii_c * km_s

  !Careful: normalised age of the time during which the star is emitting HII ionising flux
  hii_t = hii_t * Myr 
  hii_T2 = hii_T2 / scale_T2
  mH_code = mH / (scale_d * scale_l**3) ! make this useful...

  ! PREVIOUS FEED PARAMS

  !normalise the supernova quantities
  sn_p_ref = sn_p / (scale_d * scale_v * scale_l**3)
  sn_e_ref = sn_e / (scale_d * scale_v**2 * scale_l**3)
  sn_mass_ref = sn_mass / (scale_d * scale_l**3) !1 solar mass ejected

  !normalise Vsat which is assumed to be in KM/S
  Vsat = Vsat * 1.e5 / scale_v

  !normalise Vdisp which is assumed to be in KM/S
  Vdisp = Vdisp * 1.e5 / scale_v
      
  ! Set up the single stellar sources
!added by PH to make it use-able without RT
!#ifdef RT
!      if (use_ssm) call ssm_setup(ssm_table_directory)
!#endif

111 return

end subroutine read_stellar_params