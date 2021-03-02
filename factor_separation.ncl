;*************************************************
; NCL Graphics: top.ncl
;*************************************************
; Concepts illustrated:
;   - Using 'short2flt' to unpack 'short' data
;   - Drawing line contours over a satellite map
;   - Changing the view of a satellite map
;************************************************
; This file is loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;************************************************
begin
;************************************************
; read in netCDF file
;************************************************

  climname=new(8,string)
  climname=(/"E280","E400","Eo280","Eo400","Ei280","Ei400","Eoi280","Eoi400"/)

  all_files="../code_deepak/dataCode_for_paper/"+climname+".nc"

  f = addfiles(all_files,"r") 
  ListSetType (f, "join")

  temp = f[:]->TREFHT
  lons = temp&lon
  lats = temp&lat

;  temp=temp-273.15

  gwgt  = f[:]->gw

  weight=gwgt(0,:)

E280=new(dimsizes(temp(0,:,:)),"double")
E400=new(dimsizes(temp(0,:,:)),"double")
Eo280=new(dimsizes(temp(0,:,:)),"double")
Eo400=new(dimsizes(temp(0,:,:)),"double")
Ei280=new(dimsizes(temp(0,:,:)),"double")
Ei400=new(dimsizes(temp(0,:,:)),"double")
Eoi280=new(dimsizes(temp(0,:,:)),"double")
Eoi400=new(dimsizes(temp(0,:,:)),"double")

E280=flt2dble(temp(0,:,:))
E400=flt2dble(temp(1,:,:))
Eo280=flt2dble(temp(2,:,:))
Eo400=flt2dble(temp(3,:,:))
Ei280=flt2dble(temp(4,:,:))
Ei400=flt2dble(temp(5,:,:))
Eoi280=flt2dble(temp(6,:,:))
Eoi400=flt2dble(temp(7,:,:))




nvars=4

dTall=E280
dTco2=E280
dToro=E280
dTice=E280


ntype=6


   wlonind=232
   wlatind=100


uAve_gwgt=new((/ntype,nvars/),"double")

   do t=0, ntype-1


dTall=Eoi400-E280

   if (t .eq. 0) then
typename="Lunt et al (2012)"
typenameshort="lunt"
dTco2=1.0/4.0*((E400-E280)+(Eo400-Eo280)+(Eoi400-Eoi280)+(Ei400-Ei280))
dToro=1.0/4.0*((Eo280-E280)+(Eo400-E400)+(Eoi400-Ei400)+(Eoi280-Ei280))
dTice=1.0/4.0*((Ei280-E280)+(Ei400-E400)+(Eoi400-Eo400)+(Eoi280-Eo280))
   end if

   if (t .eq. 1) then
typename="Linear Sum" 
typenameshort="lins"
dTco2=1.0/6.0*(2*(E400-E280)+(Eo400-Eo280)+2*(Eoi400-Eoi280)+(Ei400-Ei280))
dToro=1.0/6.0*(2*(Eo280-E280)+(Eo400-E400)+2*(Eoi400-Ei400)+(Eoi280-Ei280))
dTice=1.0/6.0*(2*(Ei280-E280)+(Ei400-E400)+2*(Eoi400-Eo400)+(Eoi280-Eo280))
   end if

   if (t .eq. 2) then
typename="Scaled Total"
typenameshort="scat"
dTco2p=1.0/4.0*((E400-E280)+(Eo400-Eo280)+(Eoi400-Eoi280)+(Ei400-Ei280))
dTorop=1.0/4.0*((Eo280-E280)+(Eo400-E400)+(Eoi400-Ei400)+(Eoi280-Ei280))
dTicep=1.0/4.0*((Ei280-E280)+(Ei400-E400)+(Eoi400-Eo400)+(Eoi280-Eo280))
dTallp=dTco2p+dTorop+dTicep
dTco2=dTco2p*(Eoi400-E280)/dTallp
dToro=dTorop*(Eoi400-E280)/dTallp
dTice=dTicep*(Eoi400-E280)/dTallp

   print(" ")
   print("testing:")
   print(lons(wlonind))
   print(lats(wlatind))
   print("Eoi400: "+Eoi400(wlatind,wlonind))
   print("E280: "+E280(wlatind,wlonind))

dTco2pp=1.0/4.0*((E400(wlatind,wlonind)-E280(wlatind,wlonind))+(Eo400(wlatind,wlonind)-Eo280(wlatind,wlonind))+(Eoi400(wlatind,wlonind)-Eoi280(wlatind,wlonind))+(Ei400(wlatind,wlonind)-Ei280(wlatind,wlonind)))
dToropp=1.0/4.0*((Eo280(wlatind,wlonind)-E280(wlatind,wlonind))+(Eo400(wlatind,wlonind)-E400(wlatind,wlonind))+(Eoi400(wlatind,wlonind)-Ei400(wlatind,wlonind))+(Eoi280(wlatind,wlonind)-Ei280(wlatind,wlonind)))
dTicepp=1.0/4.0*((Ei280(wlatind,wlonind)-E280(wlatind,wlonind))+(Ei400(wlatind,wlonind)-E400(wlatind,wlonind))+(Eoi400(wlatind,wlonind)-Eo400(wlatind,wlonind))+(Eoi280(wlatind,wlonind)-Eo280(wlatind,wlonind)))

   dTallpp=dTco2pp+dToropp+dTicepp
   dTallto=Eoi400(wlatind,wlonind)-E280(wlatind,wlonind)

   print("dTco2pp: "+dTco2pp)
   print("dToropp: "+dToropp)
   print("dTicepp: "+dTicepp)
   print("dTallpp: "+dTallpp)
   print("dTallto: "+dTallto)

   my_scaling=dTallto/dTallpp
   print("scaling: "+my_scaling)
   
   print("co2: "+dTco2pp*my_scaling)
   print("oro: "+dToropp*my_scaling)
   print("ice: "+dTicepp*my_scaling)



   end if

   if (t .eq. 3) then
typename="Scaled Total (2)"
typenameshort="sca2"
dTco2p=1.0/4.0*((E400-E280)+(Eo400-Eo280)+(Eoi400-Eoi280)+(Ei400-Ei280))
dTorop=1.0/4.0*((Eo280-E280)+(Eo400-E400)+(Eoi400-Ei400)+(Eoi280-Ei280))
dTicep=1.0/4.0*((Ei280-E280)+(Ei400-E400)+(Eoi400-Eo400)+(Eoi280-Eo280))
dTallp=dTco2p+dTorop+dTicep
R=(Eoi400-E280)-(dTallp)
dTco2=dTco2p + R*dTco2p/dTallp
dToro=dTorop + R*dTorop/dTallp
dTice=dTicep + R*dTicep/dTallp
   end if

   if (t .eq. 4) then
typename="Scaled Total (3)"
typenameshort="sca3"
dTco2p=1.0/4.0*((E400-E280)+(Eo400-Eo280)+(Eoi400-Eoi280)+(Ei400-Ei280))
dTorop=1.0/4.0*((Eo280-E280)+(Eo400-E400)+(Eoi400-Ei400)+(Eoi280-Ei280))
dTicep=1.0/4.0*((Ei280-E280)+(Ei400-E400)+(Eoi400-Eo400)+(Eoi280-Eo280))
dTallp=dTco2p+dTorop+dTicep
Rp=((3.0*Eoi400 + Ei400 + Eoi280  + Eo400)-(3.0*E280 + Ei280 + Eo280 + E400))/4.0
dTco2=dTco2p * (Eoi400-E280)/Rp
dToro=dTorop * (Eoi400-E280)/Rp
dTice=dTicep * (Eoi400-E280)/Rp

print("parts for method 4:")
print("Eoi400: "+Eoi400(wlatind,wlonind))
print("Ei400: "+Ei400(wlatind,wlonind))
print("Eoi280: "+Eoi280(wlatind,wlonind))
print("Eo400: "+Eo400(wlatind,wlonind))

print("Components for method 4: "+(3*Eoi400(wlatind,wlonind) + Ei400(wlatind,wlonind) + Eoi280(wlatind,wlonind)  + Eo400(wlatind,wlonind))+" "+(3*E280(wlatind,wlonind) + Ei280(wlatind,wlonind) + Eo280(wlatind,wlonind) + E400(wlatind,wlonind)))
print("Scaling for method 4: "+Rp(wlatind,wlonind))

   end if

   if (t .eq. 5) then
typename="Scaled Total (abs)"
typenameshort="scaa"
dTco2p=1.0/4.0*((E400-E280)+(Eo400-Eo280)+(Eoi400-Eoi280)+(Ei400-Ei280))
dTorop=1.0/4.0*((Eo280-E280)+(Eo400-E400)+(Eoi400-Ei400)+(Eoi280-Ei280))
dTicep=1.0/4.0*((Ei280-E280)+(Ei400-E400)+(Eoi400-Eo400)+(Eoi280-Eo280))
dTallp=dTco2p+dTorop+dTicep
dTalla=abs(dTco2p)+abs(dTorop)+abs(dTicep)
R=(Eoi400-E280)-(dTallp)
dTco2=dTco2p + R*abs(dTco2p)/dTalla
dToro=dTorop + R*abs(dTorop)/dTalla
dTice=dTicep + R*abs(dTicep)/dTalla
   end if

   print(" ")



    minlat = -90
    maxlat = 90
    minlon = -180
    maxlon = 180

;     my_cmap="WhViBlGrYeOrRe"
;     my_cmap="cmp_b2r"
;     my_cmap="BlueWhiteOrangeRed"
;     my_cmap="posneg_1"
     my_cmap="BlueYellowRed"

     my_max=20
     my_min=-20.0
     my_del=1.0

   do v=0, nvars-1

   if (v .eq. 0) then
     name_var="dTall"    
     my_var=dTall
   end if

   if (v .eq. 1) then
     name_var="dTco2"    
     my_var=dTco2
   end if

   if (v .eq. 2) then
     name_var="dTice"    
     my_var=dTice
   end if

   if (v .eq. 3) then
     name_var="dToro"    
     my_var=dToro
   end if

     my_label=name_var


;************************************************
; plotting parameters
;************************************************
  wks = gsn_open_wks("pdf","plot_"+name_var+"_"+typenameshort)             ; send graphics to PNG file

   if (v .ne. 0 .and. t .ne. 0) then
   delete(cmap)
   delete(res)
   end if
 
 res                            = True       ; plot mods desired

 
  cmap                = read_colormap_file(my_cmap)

;  res@cnFillPalette  = cmap


;---Set some resources for contouring and mapping

  res@gsnMaximize        = True     ; maximize plot in frame
  res@cnFillOn           = True     ; turn on contour fill
  res@cnFillPalette      = cmap(2:,:) ; set color map
  res@cnLinesOn          = False    ; turn off contour lines
  res@cnLineLabelsOn     = False    ; turn off line labels
  res@cnInfoLabelOn      = False    ; turn off info label
  res@lbBoxLinesOn       = False    ; turn off labelbar box lines
;  res@cnFillMode         = "MeshFill"   ; for faster draw

;---Pick "nice" contour levels
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = my_min
  res@cnMaxLevelValF       = my_max
  res@cnLevelSpacingF      = my_del

  res@gsnAddCyclic         = True           ; add longitude cyclic point


  res@cnFillOn = True
  res@cnMissingValFillColor = "LightBlue"

;---Zoom in on map
  res@mpMinLatF            = minlat
  res@mpMaxLatF            = maxlat
  res@mpMinLonF            = minlon
  res@mpMaxLonF            = maxlon
  res@mpCenterLonF         = (res@mpMinLonF + res@mpMaxLonF) / 2.

  res@pmTickMarkDisplayMode       = "Always"

  res@mpProjection               = "CylindricalEquidistant" ; choose map projection

  res@mpCenterLonF               = 0.       ; choose center lon
  res@mpCenterLatF               = 0.         ; choose center lat
;  res@mpSatelliteDistF           = 3.0         ; choose satellite view

  res@tiMainString    = "" ; add title

  res@mpGridAndLimbOn        = True               ; turn on grid lines
  res@mpGridLineDashPattern  = 5                  ; lat/lon lines dashed

  res@gsnLeftStringFontHeightF = 0.01
  res@gsnRightStringFontHeightF = 0.01
  res@gsnLeftString          = typename
  res@gsnRightString          = my_label

  res@gsnFrame     =  False                   ; do not advance the frame
  res@gsnDraw      =  False                   ; do not draw the plot

  res@mpDataBaseVersion    = "LowRes"    ; better map outlines
  res@gsnPaperOrientation = "portrait"
  res@mpOutlineOn = True
  res@lbBoxEndCapStyle = "TriangleBothEnds" 
  res@mpGeophysicalLineThicknessF = 0.5


  cres=True
  cres@cnInfoLabelOn     = False
  cres@gsnLeftString          = ""
  cres@gsnRightString          = ""
  cres@cnLineLabelsOn       = False            ; do not use line labels
  cres@gsnDraw              = False             ; do not draw the plot
  cres@gsnFrame             = False             ; do not advance the frame
  cres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
  cres@cnLevels             = ispan(3,6,3)   ; set the contour levels
  cres@gsnContourLineThicknessesScale = 3
;  cres@cnSmoothingOn = True


  map = gsn_csm_contour_map(wks,my_var(:,:),res)
  map2 = gsn_csm_contour(wks,my_var(:,:),cres) 




;  overlay(map,map2)

  draw(map)
  frame(wks)

  uAve_gwgt(t,v) = wgt_areaave(my_var(:,:), weight, 1.0, 1)
  print(typename+" "+t+": "+name_var+" "+v+" = "+uAve_gwgt(t,v))
  
  print(typename+" max: "+max(my_var(:,:)))
  print(typename+" min: "+min(my_var(:,:)))

  dims = dimsizes(my_var)
  x1d = ndtooned(my_var)      ; convert 2D array to 1D for use in maxind
  inds = ind_resolve(maxind (x1d), dims)    ; convert 1D array back to 2D 
  ilat = inds(0,0)        ; select the latitude index where the X array is at its' maximum  
  ilon = inds(0,1)        ; select the longitude index where the X array is at its' maximum
  print("Maximum value located at index "+ilon+", "+ilat) 
  lat_max = lats(ilat)   ; insert the latitude index into the lat coordinate variable
  lon_max = lons(ilon)   ; insert the longitude index into the lon coordinate variable
  print("Maximum value located at "+lon_max+", "+lat_max) 

  print(" ")

  end do

  print("total: "+sum(uAve_gwgt(t,1:3)))
  print(" ")

  end do

   do t=1, ntype-1
   do v=0, nvars-1
   perc=100.0*(uAve_gwgt(t,v)-uAve_gwgt(0,v))/uAve_gwgt(0,v)
   print(t+": "+v+" = "+perc) 

   end do
   end do

end
