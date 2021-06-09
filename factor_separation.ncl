;*************************************************
; run using: 
; ncl factor_separation.ncl
;*************************************************

begin

do_testing_bound=0
do_testing_maxval=1
do_avg=1
do_plots=0
do_perc=0


; Read in netcdf files

  climname=new(8,string)
  climname=(/"E280","E400","Eo280","Eo400","Ei280","Ei400","Eoi280","Eoi400"/)

  all_files="../code_deepak/dataCode_for_paper/"+climname+".nc"

  f = addfiles(all_files,"r") 
  ListSetType (f, "join")

  temp = f[:]->TREFHT
  lons = temp&lon
  lats = temp&lat

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

ntype=4

nx=288
ny=192

dTco2pmax=new((/ny,nx/),"double")
dTco2pmin=new((/ny,nx/),"double")

uAve_gwgt=new((/ntype,nvars/),"double")


dTall=Eoi400-E280

dTco2p=1.0/4.0*((E400-E280)+(Eo400-Eo280)+(Eoi400-Eoi280)+(Ei400-Ei280))
dTorop=1.0/4.0*((Eo280-E280)+(Eo400-E400)+(Eoi400-Ei400)+(Eoi280-Ei280))
dTicep=1.0/4.0*((Ei280-E280)+(Ei400-E400)+(Eoi400-Eo400)+(Eoi280-Eo280))
dTallp=dTco2p+dTorop+dTicep

if (do_testing_bound .eq. 1) then
do j=0,ny-1
do i=0,nx-1
dTco2pmax(j,i)=max((/(E400(j,i)-E280(j,i)),(Eo400(j,i)-Eo280(j,i)),(Eoi400(j,i)-Eoi280(j,i)),(Ei400(j,i)-Ei280(j,i))/))
dTco2pmin(j,i)=min((/(E400(j,i)-E280(j,i)),(Eo400(j,i)-Eo280(j,i)),(Eoi400(j,i)-Eoi280(j,i)),(Ei400(j,i)-Ei280(j,i))/))
end do
end do
end if


   do t=0, ntype-1



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
typename="Scaled Residual"
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

   if (t .eq. 3) then
typename="Scaled Residual (relative)"
typenameshort="scat"
dTco2=dTco2p*(Eoi400-E280)/dTallp
dToro=dTorop*(Eoi400-E280)/dTallp
dTice=dTicep*(Eoi400-E280)/dTallp
   end if



if (do_testing_bound .eq. 1) then

print("Start Testing")

do j=0,ny-1
do i=0,nx-1

if ( dTco2(j,i) .lt. dTco2pmin(j,i) ) then
print("UNBOUNDED -ve: "+t+" "+i+" "+j+" "+dTco2(j,i)+" "+dTco2pmax(j,i)+" "+dTco2pmin(j,i))
end if

if ( dTco2(j,i) .gt. dTco2pmax(j,i) ) then
print("UNBOUNDED +ve: "+t+" "+i+" "+j+" "+dTco2(j,i)+" "+dTco2pmax(j,i)+" "+dTco2pmin(j,i))
end if

end do
end do

print("End Testing")

end if

 
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

  if (do_avg .eq. 1) then
  uAve_gwgt(t,v) = wgt_areaave(my_var(:,:), weight, 1.0, 1)
  print(typename+" "+t+": "+name_var+" "+v+" = "+uAve_gwgt(t,v))
  end if

  if (do_testing_maxval .eq. 1) then
  print(typename+" max: "+max(my_var(:,:)))
  print(typename+" min: "+min(my_var(:,:)))
  dims = dimsizes(my_var)
  x1d = ndtooned(my_var)      ; convert 2D array to 1D for use in maxind
  inds = ind_resolve(maxind (x1d), dims)    ; convert 1D array back to 2D 
  ilat = inds(0,0)        ; select the latitude index where the X array is at its maximum  
  ilon = inds(0,1)        ; select the longitude index where the X array is at its maximum
  print("Maximum value located at index "+ilon+", "+ilat) 
  lat_max = lats(ilat)   ; insert the latitude index into the lat coordinate variable
  lon_max = lons(ilon)   ; insert the longitude index into the lon coordinate variable
  print("Maximum value located at "+lon_max+", "+lat_max) 
  end if




  if (do_plots .eq. 1) then

   minlat = -90
    maxlat = 90
    minlon = -180
    maxlon = 180

     my_cmap="BlueYellowRed"

     my_max=20
     my_min=-20.0
     my_del=1.0

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


;---Set some resources for contouring and mapping
  res@gsnMaximize        = True     ; maximize plot in frame
  res@cnFillOn           = True     ; turn on contour fill
  res@cnFillPalette      = cmap(2:,:) ; set color map
  res@cnLinesOn          = False    ; turn off contour lines
  res@cnLineLabelsOn     = False    ; turn off line labels
  res@cnInfoLabelOn      = False    ; turn off info label
  res@lbBoxLinesOn       = False    ; turn off labelbar box lines

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

  map = gsn_csm_contour_map(wks,my_var(:,:),res)
  map2 = gsn_csm_contour(wks,my_var(:,:),cres) 

  draw(map)
  frame(wks)

  end if

  end do

  if (do_avg .eq. 1) then
  print("total: "+sum(uAve_gwgt(t,1:3)))
  print(" ")
  end if


  end do

   if (do_perc .eq. 1) then
   print("for comment in paper about percentages:")
   do t=1, ntype-1
   do v=0, nvars-1
   perc=100.0*(uAve_gwgt(t,v)-uAve_gwgt(0,v))/uAve_gwgt(0,v)
   print(t+": "+v+" = "+perc) 
   end do
   end do
   end if

end
