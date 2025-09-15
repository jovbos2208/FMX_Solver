module nrlmsis_cwrap
  use, intrinsic :: iso_c_binding
  use msis_calc,    only: msiscalc
  use msis_constants, only: rp
contains

  subroutine msis_eval_c(day, utsec, z_km, lat_deg, lon_deg, f107a, f107, ap_daily, ap_now, &
                         Tn, dn_tot, dn_n2, dn_o2, dn_o, dn_he, dn_h, dn_ar, dn_n, dn_ao, dn_no) bind(C, name="msis_eval_c")
    implicit none
    real(c_double), value :: day, utsec, z_km, lat_deg, lon_deg, f107a, f107, ap_daily, ap_now
    real(c_double) :: Tn, dn_tot, dn_n2, dn_o2, dn_o, dn_he, dn_h, dn_ar, dn_n, dn_ao, dn_no
    real(kind=rp) :: day_r, utsec_r, z_km_r, lat_r, lon_r, f107a_r, f107_r, ap_r(7)
    real(kind=rp) :: Tn_r, dn_r(10)

    day_r   = real(day,   kind=rp)
    utsec_r = real(utsec, kind=rp)
    z_km_r  = real(z_km,  kind=rp)
    lat_r   = real(lat_deg, kind=rp)
    lon_r   = real(lon_deg, kind=rp)
    f107a_r = real(f107a, kind=rp)
    f107_r  = real(f107,  kind=rp)
    ap_r    = 0.0_rp
    ap_r(1) = real(ap_daily, kind=rp)
    ap_r(2) = real(ap_now,   kind=rp)

    call msiscalc(day_r, utsec_r, z_km_r, lat_r, lon_r, f107a_r, f107_r, ap_r, Tn_r, dn_r)

    Tn     = real(Tn_r,  c_double)
    dn_tot = real(dn_r(1), c_double)
    dn_n2  = real(dn_r(2), c_double)
    dn_o2  = real(dn_r(3), c_double)
    dn_o   = real(dn_r(4), c_double)
    dn_he  = real(dn_r(5), c_double)
    dn_h   = real(dn_r(6), c_double)
    dn_ar  = real(dn_r(7), c_double)
    dn_n   = real(dn_r(8), c_double)
    dn_ao  = real(dn_r(9), c_double)
    dn_no  = real(dn_r(10), c_double)
  end subroutine msis_eval_c

end module nrlmsis_cwrap
