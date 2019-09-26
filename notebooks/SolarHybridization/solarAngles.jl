function solarAngles(h,tzone,lat,long)
  #=============================================================================
   This function calculates the various solar angles and returns the incidence
   angle of direct rays on a horizontal surface (with single-axis tracking)
   # Written by: M. D. Stuber, March 28, 2018, rev. May 15, 2019 (Julia 1.0)
   # This code is used in the paper: Stuber (2018), DOI: 10.3390/pr6070076
   INPUT
    h: hour of year (1:8760)
    tzone: time-zone
    long: longitude
    lat: latitude
   OUTPUT
    ia: solar incidence angle
  =============================================================================#
  d = floor(h/24+1)# hour to day number conversion
  b=360/365*(d-81)
  et = 9.87*sin(2*b*pi/180)-7.53*cos(b*pi/180)-1.5*sin(b*pi/180) # equation of time
  TC = 4*(15*tzone-long)+et # solar time calculation
  ha = 15*(h+TC/60-12) # hour angle
  da = 23.45*sin(b*pi/180) # declination angle
  za = (acos(sin(lat*pi/180)*sin(da*pi/180)+cos(lat*pi/180)*cos(da*pi/180).*cos(ha*pi/180)))*180/pi # zenith angle
  se = (asin(cos(za*pi/180)))*180/pi # solar elevation angle
  # azimuth angle calculation
  if (cos(ha*pi/180)>tan(da*pi/180)/tan(lat*pi/180))
    aa = (asin(sin(ha*pi/180)*cos(da*pi/180)/cos(se*pi/180)))*180/pi
  elseif (ha<=-1e-16)
    aa = (-pi+abs(asin(sin(ha*pi/180)*cos(da*pi/180)/cos(se*pi/180))))*180/pi
  else
    aa = (pi-asin(sin(ha*pi/180)*cos(da*pi/180)/cos(se*pi/180)))
  end
  if se>0 sun = 1 else sun = 0 end # is the sun up?
  # calculate the mirror tilt angle for 1-axis tracking
  if sun==0 mt = -90 else mt = atan(tan(za*pi/180)*cos((90-aa)*pi/180))*180/pi end
  ia = acos(sqrt(cos(za*pi/180)^2+cos(da*pi/180)^2*sin(ha*pi/180)^2))*180/pi # incidence angle
  return ia # return the incidence angle
end;
