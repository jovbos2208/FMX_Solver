module hwm14_cwrap
  use, intrinsic :: iso_c_binding
  use hwm_interface, only: hwm_14
contains
  subroutine hwm14_eval_c(doy, utsec, alt_km, lat_deg, lon_deg, ap3hr, w_meridional, w_zonal) bind(C, name="hwm14_eval_c")
    implicit none
    integer(c_int), value :: doy
    real(c_double), value :: utsec, alt_km, lat_deg, lon_deg, ap3hr
    real(c_double) :: w_meridional, w_zonal
    real(8) :: wm, wz
    call hwm_14(doy, utsec, alt_km, lat_deg, lon_deg, ap3hr, wm, wz)
    w_meridional = wm
    w_zonal = wz
  end subroutine hwm14_eval_c
end module hwm14_cwrap
