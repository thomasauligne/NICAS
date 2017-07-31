load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;model = "arp"
;xp = "6B60"
model = "aro"
xp = "5840"
;suffix = "_deformed"
;suffix = "_resol-10"
;suffix = "_heterogeneous"
suffix = "_homogeneous_sqrt"
levplot = 1
ilplot = 1

; Load data :
; -----------

data = addfile("../../data/" + model + "/" + xp + "/" + model + "_" + xp + "_dirac_" + sprinti("%0.3i",levplot) + ".nc","r")
longitude = data->longitude
latitude = data->latitude
var = data->hr(ilplot-1,:,:)
var = where(var.lt.1.0e-12,1.0e-12,var)
if (model.eq."aro") then
   lon = longitude
   lat = latitude
   hr = var
end if
if (model.eq."arp") then
   dims = dimsizes(longitude)
   nlat = dims(0)
   nlon = dims(1)
   nc0 = num(.not.ismissing(longitude))
   lon = new(nc0,double)
   lat = new(nc0,double)
   hr = new(nc0,double)
   ic0 = 0
   do ilon=0,nlon-1
     do ilat=0,nlat-1
       if (.not.ismissing(longitude(ilat,ilon))) then
         lon(ic0) = longitude(ilat,ilon)
         lat(ic0) = latitude(ilat,ilon)
         hr(ic0) = var(ilat,ilon)
         ic0 = ic0+1
       end if
     end do
   end do
end if

; Resources :
; -----------

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnMaximize = True
res@gsnAddCyclic = False

res@cnFillOn = True
res@cnFillMode = "AreaFill"
res@trGridType = "TriangularMesh"
res@cnMonoFillPattern = True
res@cnMonoFillColor = False
res@lbLabelBarOn = True
res@lbOrientation = "vertical"
res@cnInfoLabelOn = False
res@cnLineLabelsOn = False
res@cnLinesOn = False
res@cnNoDataLabelOn = False
res@cnMissingValFillColor = 0
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF = 0.05
res@cnMaxLevelValF = 1.0
res@cnMinLevelValF = 0.0

res@mpOutlineOn = True
res@mpDataBaseVersion = "MediumRes"
res@mpGridLatSpacingF = 20.0
res@mpDataSetName = "Earth..4"
res@mpOutlineBoundarySets = "Geophysical"
res@mpLandFillColor = -1
res@mpGridAndLimbDrawOrder = "PreDraw"
res@mpPerimOn = False
res@mpGreatCircleLinesOn = True
if (model.eq."aro") then
   res@mpProjection = "LambertConformal"
   res@mpGridAndLimbOn = False
   res@mpLimitMode = "LatLon"
   res@mpMinLonF = min(lon)
   res@mpMaxLonF = max(lon)-3.0
   res@mpMinLatF = min(lat)
   res@mpMaxLatF = max(lat)
end if
if (model.eq."arp") then
   res@mpProjection = "WinkelTripel"
   res@mpGridAndLimbOn = True
   res@mpGridLineDashPattern = 2
   res@mpGridLineColor = -1
   res@mpGridSpacingF = 45.0
   res@mpCenterLonF = 0.0
   res@mpCenterLatF = 0.0
end if

wks_epsi = "epsi"

; Display :
; ---------

res@sfXArray = lon
res@sfYArray = lat

system("mkdir -p ../fig")

output = "../fig/" + model + "_" + xp + "_dirac" + "_" + sprinti("%0.3i",levplot) + suffix
wks = gsn_open_wks(wks_epsi,output)
gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

plot = gsn_csm_contour_map(wks,hr,res)

draw(plot)
frame(wks)
delete(wks)

system("epstopdf " + output + ".epsi;rm -f " + output + ".epsi")

end
