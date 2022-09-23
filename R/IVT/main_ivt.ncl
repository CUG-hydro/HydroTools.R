; notes:
; 1. q&level in the unit of hPa
; 2. the right most dimension should be `lat`, `lon`
; 
undef("hPa2Pa")
function hPa2Pa(x)
begin
  
end

undef("cal_IVT")
function cal_IVT(u, v, q, ps, ptop, outfile, debug)
begin
  ; CONSTANTS
  g    = 9.80665          ; m/s2
  
  ;---Change [kg/kg] to [g/kg]; not necessary: but common units for q
  q    = q*1000
  q@units = "g/kg"
  
  plev = q&level                      ; hPa  
  plev = plev*100                     ; [100000,...,30000] Pa [kg/(m-s2)]
  plev@units = "Pa"

  ;---Layer thickness: ; Pa=>[kg/(m-s2)], (time,level,lat,lon) 
  ;---Mass weighting:  (dp/g) => [Pa/(m/s2)] => (Pa-s2)/m => [kg/(m-s2)][s2/m] =>  (kg/m2)
  ;---Reference: http://www.cgd.ucar.edu/cas/catalog/newbudgets/
  dp   = dpres_plevel_Wrap(plev, ps, ptop, 0) ; Pa; layar thickness 
  dpg  = dp/g
  dpg@long_name = "Layer Mass Weighting"
  dpg@units     = "kg/m2"                     ; dp/g, Pa/(m s-2), reduce to kg m-2

  ;---Moisture flux components at each pressure level
  uq   = u*q                                  ; (:,:,:,:)
  uq@long_name = "Zonal Moisture Flux [uq]"
  uq@units = "["+u@units+"]["+q@units+"]"     ; [m/s][g/kg]     
  copy_VarCoords(u,uq)                        ; (time,level,lat,lon)

  vq   = v*q                                  ; (:,:,:,:)
  vq@long_name = "Meridional Moisture Flux [vq]"
  vq@units = "["+v@units+"]["+q@units+"]" 
  copy_VarCoords(v,vq)                        ; (time,level,lat,lon)


  ;---Integrated mass weighted moisture flux components
  uq_dpg = uq*dpg                ; mass weighted 'uq'; [m/s][g/kg][kg/m2]=>[m/s][g/kg]
  iuq    = dim_sum_n(uq_dpg, 1)
  iuq@long_name = "Integrated Zonal UQ [uq*dpg]" 
  iuq@LONG_NAME = "Sum: Mass Weighted Integrated Zonal Moisture Flux [uq*dpg]" 
  iuq@units     = "[m/s][g/kg]"
  copy_VarCoords(u(:,0,:,:), iuq); (time,lat,lon)
  delete(uq_dpg)

  vq_dpg = vq*dpg                ; mass weighted 'vq'; [m/s][g/kg][kg/m2]=>[m/s][g/kg] 
  ivq    = dim_sum_n(vq_dpg, 1)
  ivq@long_name = "Integrated Meridional VQ [vq*dpg]" 
  ivq@LONG_NAME = "Sum: Mass Weighted Integrated Meridional Moisture Flux [vq*dpg]" 
  ivq@units     = "[m/s][g/kg]"
  copy_VarCoords(v(:,0,:,:), ivq); (time,lat,lon)
  delete(vq_dpg)

  ;---Divergence of moisture flux: uv2dvF => global 'fixed' rectilinear grid

  duvq  = uv2dvF_Wrap(uq, vq)    ; (time,level,lat,lon)
  duvq@long_name = "Divergence of Moisture Flux"
  duvq@units     = "g/(kg-s)"    ; (1/m)*[(m/s)(g/kg)] => [g/(kg-s)]

  ;---Mass weighted integration [sum] of the divergence of moisture flux

  duvq_dpg = duvq*dpg            ;  [g/(kg-s)][kg/m2] => [g/(m2-s)]
  iduvq    = dim_sum_n(duvq_dpg, 1)
  iduvq@long_name = "Integrated Mass Wgt MFC" 
  iduvq@LONG_NAME = "Integrated Mass Weighted Moisture Flux Convergence" 
  iduvq@units     = "g/(m2-s)"
  copy_VarCoords(u(:,0,:,:), iduvq)      ; (time,lat,lon)
  delete(duvq_dpg)

  VIMFC =  iduvq           ; keep meta data                         
  VIMFC = -VIMFC           ; Note the preceding -1 [negative precedes integration] 
  VIMFC@long_name = "VIMFC"

  
  ;*************************************************
  ; Calculate divergence: Use Wrap to include meta data
  ; Calculate divergent wind components; used for graphics 
  ;*************************************************
  div = uv2dvF_Wrap(u,v)                ; u,v ==> divergence; (:,:,:,:)

  ud  = new ( dimsizes(u), typeof(u), "No_FillValue")
  vd  = new ( dimsizes(v), typeof(v), "No_FillValue")
  dv2uvf(div,ud,vd)                     ; divergence ==> divergent components  

  copy_VarCoords(u, ud ) 
  copy_VarCoords(u, vd ) 
  ud@long_name  = "Zonal Divergent Wind"
  ud@units      = u@units
  vd@long_name  = "Meridional Divergent Wind"
  vd@units      = v@units

  ; ## save data to nefcdf
  system("rm -rf IVT.nc")
  fout = addfile("IVT.nc","c")
  fout->uq = uq
  fout->vq = vq
  fout->iuq = iuq
  fout->ivq = ivq
  fout->duvq = duvq
  fout->VIMFC = VIMFC
  fout->ud = ud
  fout->vd = vd
  delete(fout)
end
