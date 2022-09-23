;*************************************************
; https://www.ncl.ucar.edu/Applications/Scripts/mfc_div_1.ncl
; 
; mfc_div_1.ncl
;
; Concepts illustrated:
;   - Read daily mean wind components, humidity and sfc. pressure 
;     from different files
;   - Reorder the input (N==>S) grid order to (S==>N) via NCL syntax  ::-1
;   - Calculate mass weightined layer thickness [units="kg/m2"]
;   - Calculate moisture flux [uq, vq]             
;   - Calculate moisture flux divergence using spherical harmonics              
;   - Integrate the moisture flux divergence using mass weighting 
;   - Plot a number of quantities
;*************************************************
;---Calculate the Horizontal Moisture Flux Convergence [MFC]
;*************************************************
;---High frequency source data: hourly/3hr/6hr/12hr/daily .... NOT monthly values
;---References:
;---http://www.cgd.ucar.edu/cas/catalog/newbudgets/
;---http://tornado.sfsu.edu/geosciences/classes/e260/AtmosphericRivers/Moisture%20Flux.pdf
;---https://www.spc.noaa.gov/publications/banacos/mfc-sls.pdf
;===================================================================
;   Data Source: ESRL Physical Sciences Division
;        https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
;   NCEP Reanalysis data provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, 
;   from their Web site at https://www.esrl.noaa.gov/psd/
;===================================================================
undef("print_var")
procedure print_var(x, name)
begin
  print("[" + name + "]" + "-----")
  printVarSummary(x)                          ; (time,level,lat,lon); m/s
  printMinMax(x,0)
end

undef("make_var")
function make_var(x, units)
begin
  x@units = units
  return(x)
end

"""
Arguments:
- `ptop`:
- `date`: 
- files:
- `indir`:
+ `f_q`: 
+ `f_u`: 
+ `f_v`: 
+ `f_ps`: 

# Input Data:
- `shum`: https://psl.noaa.gov/repository/entry/show?entryid=eb41bda7-441e-490f-b7db-b7113741c90e


"""
date  = 20080715        ; NH summer
debug = True ; If false, all time will be calculated
ptop = make_var(300, "hPar") ; 'shum' upper level                  

;---ESRL: CDC data
indir = "./data/"
f_q = indir + "shum.2008.nc"   ; daily data for current year [366 days]
f_u = indir + "uwnd.2008.nc"
f_v = indir + "vwnd.2008.nc"
f_ps= indir + "pres.sfc.2008.nc"

; ## 01. READ DATA -------------------------------------------------------------
fu   = addfile(f_u ,"r")
fv   = addfile(f_v ,"r")
fq   = addfile(f_q ,"r")
fps  = addfile(f_ps,"r")

;---Time
ymd  = cd_calendar(fu->time, -2)    ; ymd[*]: human readable
nt   = ind(ymd .eq. date)             ; date for plotting and testing

if (.not.debug) then                   ; all times                   
  u    = fu->uwnd(:,{1000:ptop},:,:)  ; m/s, (time,level,lat,lon)
  v    = fv->vwnd(:,{1000:ptop},:,:)
  q    = fq->shum                     ; [kg/kg], 1000-300 levels only
  ps   = fps->pres                    ; Pa=>[kg/(m-s2)], (time,lat,lon)
else                                  ; one time step; keep time dimension [ nt:nt: ]                     
  u    = fu->uwnd(nt:nt,{1000:ptop},:,:); m/s, (time,level,lat,lon)
  v    = fv->vwnd(nt:nt,{1000:ptop},:,:)
  q    = fq->shum(nt:nt,:,:,:)        ; [kg/kg], 1000-300 levels only
  ps   = fps->pres(nt:nt,:,:)         ; Pa=>[kg/(m-s2)], (time,lat,lon)
  nt   = 0                            ; only one time step
end if

;---Divergence function [used later] requires S->N grid order
u    = u(:,:,::-1,:)  
v    = v(:,:,::-1,:)
q    = q(:,:,::-1,:)
ps   = ps(:, ::-1,:)

;---Vertical levels
ptop = ptop*100
ptop@units = "Pa"


PRINT_RAW = True
if (PRINT_RAW) then
  print_var(q)  ; (time,level,lat,lon); g/kg
  print_var(u)  ; (time,level,lat,lon); m/s
  print_var(v)  ; (time,level,lat,lon); m/s
  print_var(ps) ; (time,lat,lon); Pa => kg/(m-s2) 
  print_var(uq) ; (time,level,lat,lon); (m/s)(g/kg)
  print_var(vq) ; (time,level,lat,lon); (m/s)(g/kg)
  print_var(dp) ; (time,level,lat,lon); Pa => kg/(m-s2)
  
  print("-----")
  ; examine layer thickness at selected locations
  print(dp(nt,:,{40},{180}))    ; mid-Pacific
  print(dp(nt,:,{40},{255}))    ; Boulder, CO 
  print("-----")
end if

;---Another way to compute Integrated divergence of moisture flux [iduvq_1]

;;IDUVQ = wgt_vertical_n(duvq, dp, 2, 1) 
;;iduvq_0 = IDUVQ[0]
;;iduvq_0 = iduvq_0/g                    ; complete mass weighting
;;iduvq_0@long_name = "Average Mass Weighted MFC" 
;;iduvq_0@LONG_NAME = "Average Mass Weighted Moisture Flux Convergence" 
;;iduvq_0@units     = "g/(m2-s)"

;;iduvq_1 = IDUVQ[1]                     ; same as iuvq_sum
;;iduvq_1 = iduvq_1/g
;;iduvq_1@long_name = "Integrated MFC" 
;;iduvq_1@LONG_NAME = "Integrated Moisture Flux Convergence" 
;;iduvq_1@units     = "g/(m2-s)"

PRINT_RESULT = True
if (PRINT_RESULT) then
  print_var(iuq) ; (time,lat,lon)
  print_var(ivq) ; (time,lat,lon)
  print_var(duvq); (time,lev,lat,lon)
  print_var(iduvq); (time,lat,lon)
end if


if (PRINT_RESULT) then
  print_var(ud) ; (time,level,lat,lon)
  print_var(vd) ; (time,level,lat,lon)
end if
; uq, vq, iuq, ivq, duvq, VIMFC, ud, vd

;*************************************************
; plot results
;*************************************************    

scl5  = 1e5                                  ; arbitrary: used for nicer plot values
sclab5= "(10~S~-5~N~)"                       ; used later   
SCLAB5= "(10~S~5~N~)"           

scl6  = 1e6  
sclab6= "(10~S~-6~N~)"         
SCLAB6= "(10~S~6~N~)"         

plot := new(2,graphic)

wks   = gsn_open_wks("png","mfc_div")        ; send graphics to PNG file
resd                 = True
resd@cnFillOn        = True                  ; color
resd@cnLinesOn       = False                 ; turn off contour lines

resd@cnLevelSelectionMode = "ManualLevels"   ; set manual contour levels
resd@cnMinLevelValF  = -15.                  ; set min contour level
resd@cnMaxLevelValF  =  15.                  ; set max contour level
resd@cnLevelSpacingF =   1.                  ; set contour spacing
;resd@cnFillPalette   = "cmocean_balance"     ; NCL 6.5.0
resd@cnFillPalette   = "ViBlGrWhYeOrRe"

resd@mpFillOn        = False                 ; turn off map fill
resd@vcRefMagnitudeF = 3.                    ; make vectors larger
resd@vcRefLengthF    = 0.025                 ; reference vector length
resd@vcGlyphStyle    = "CurlyVector"         ; turn on curly vectors
resd@vcMinDistanceF  = 0.010                 ; thin the vectors
resd@vcRefAnnoOrthogonalPosF = -1.0          ; move ref vector up
resd@gsnLeftString   = "Divergent Wind"
resd@gsnScalarContour= True                  ; vectors over contours

LEVP = 700
DIV  = div(nt,{LEVP},:,:)                     ; keep meta data
DIV  = DIV*scl6                               ; nicer numbers                 

resd@tiMainString    = "Divergence and Divergent Winds" 
resd@gsnCenterString = LEVP+"hPa: "+date     
resd@gsnRightString  = sclab6+" "+div@units
dplt = gsn_csm_vector_scalar_map(wks,ud(nt,{LEVP},:,:),vd(nt,{LEVP},:,:),DIV,resd)

;--- Moisture Transport [uq, vq] at a specified pressure level

res                   = True             ; plot mods desired

res@gsnDraw           = False            ; don't draw yet
res@gsnFrame          = False            ; don't advance frame yet

res@cnFillOn          = True             ; turn on color
res@cnLinesOn         = False            ; turn off contour lines
res@cnLineLabelsOn    = False            ; turn off contour lines
res@cnFillPalette     = "ViBlGrWhYeOrRe" ; set White-in-Middle color map
res@lbLabelBarOn      = False            ; turn off individual cb's
res@mpFillOn          = False            ; turn off map fill
; Use a common scale
res@cnLevelSelectionMode = "ManualLevels"; manual set levels so lb consistent
res@cnMaxLevelValF       =  140.0        ; max level
res@cnMinLevelValF       = -res@cnMaxLevelValF     ; min level
res@cnLevelSpacingF      =   10.0        ; contour interval

LEVP    = 700
res@gsnCenterString      = LEVP+"hPa"
plot(0) = gsn_csm_contour_map(wks,uq(nt,{LEVP},:,:),res)
plot(1) = gsn_csm_contour_map(wks,vq(nt,{LEVP},:,:),res)

resP                     = True                ; modify the panel plot
resP@gsnPanelMainString  = date+": Unweighted Moisture Flux Components"
resP@gsnPanelLabelBar    = True                ; add common colorbar
gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot

;--- Integrated Moisture Transport [iuq, ivq]

delete(res@gsnCenterString)              ; not used for this plot
res@cnMaxLevelValF       = 10.0          ; min level
res@cnMinLevelValF       = -res@cnMaxLevelValF     ; min level
res@cnLevelSpacingF      =  0.5          ; contour interval

IUQ     = iuq(nt,:,:)                    ; local array: keep meta data
IUQ     = IUQ/scl5                       ; scale for plot
res@gsnRightString  = SCLAB5+" "+iuq@units
plot(0) = gsn_csm_contour_map(wks,IUQ,res)

IVQ     = ivq(nt,:,:)                    ; local array: keep meta data
IVQ     = IVQ/scl5
res@gsnRightString  = SCLAB5+" "+ivq@units
plot(1) = gsn_csm_contour_map(wks,IVQ,res)

resP@gsnPanelMainString  = date+": Mass Wgt. Component Moisture Flux"
gsn_panel(wks,plot,(/2,1/),resP)               ; now draw as one plot

delete( [/IUQ, IVQ/] )                   ; no longer needed

;---Divergence of Moisture Flux
res@cnMaxLevelValF       = 100.0          ; min level
res@cnMinLevelValF       = -res@cnMaxLevelValF     ; min level
res@cnLevelSpacingF      =  5.0          ; contour interval

LEVP    = 700
DUVQ    = duvq(nt,{LEVP},:,:)                    ; keep meta data
DUVQ    = DUVQ*scl6                              ; scale for plot
res@gsnCenterString = LEVP+"hPa"
res@gsnRightString  = sclab6+" "+duvq@units
plot(0) = gsn_csm_contour_map(wks,DUVQ,res)

LEVP    = 500
DUVQ    = duvq(nt,{LEVP},:,:)                    ; keep meta data
DUVQ    = DUVQ*scl6
res@gsnCenterString = LEVP+"hPa"
res@gsnRightString  = sclab6+" "+duvq@units
plot(1) = gsn_csm_contour_map(wks,DUVQ,res)

resP@gsnPanelMainString  = date+": Divergence of Moisture Flux"
gsn_panel(wks,plot,(/2,1/),resP)                ; now draw as one plot

delete(DUVQ)                                        ; no longer needed
delete([/res@gsnCenterString, res@gsnRightString/]) ; not used in next plot 

;---Integrated Divergence of Moisture Flux Convergence [no scaling]

res@gsnDraw              = True
res@gsnFrame             = True
res@lbLabelBarOn         = True        

;res@cnFillPalette        = "cmp_flux"
res@cnMaxLevelValF       =  0.50                ; min level
res@cnMinLevelValF       = -res@cnMaxLevelValF  ; min level
res@cnLevelSpacingF      =  0.050               ; contour interval
res@tiMainString         = date+": VIMFC"

plt = gsn_csm_contour_map(wks,VIMFC(nt,:,:) ,res)
