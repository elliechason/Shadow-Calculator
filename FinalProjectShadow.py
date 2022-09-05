def shadow(y, m, day, dd, sigma, phi, height, e):
   y, m, day, sigma, phi, height, e = float(y), float(m), float(day), float(sigma),\
                                       float(phi), float(height), float(e)
   if m == 1 or m == 2:
       y -= 1
       m += 12
   d = dd + day
   b = 2-int(y/100)+int(int(y/100)/4)
   jd = int(365.25*(y+4716))+int(30.6001*(m+1))+d+b-1524.5
   t = y -2000
   delta_t = 63.86+(0.3345*t)-(0.060374*(t**2))+(0.0017275*(t**3))+\
             (0.000651814*(t**4))+(0.00002373599*(t**5))

   jde = jd + (delta_t/86400)
   jc = (jd-2451545)/36525
   jce = (jde-2451545)/36525
   jme = jce/10

   from EPT import EPT, B_EPT, R_EPT, nutation
   from math import pi, atan, sin, cos, tan, atan2, asin, radians
   #l is Earth heliocentric longitude
   l_all = EPT(jme)
   l0, l1, l2, l3, l4, l5 = l_all[0], l_all[1], l_all[2], l_all[3], l_all[4], l_all[5]
   l_rad = (l0+(l1*jme)+(l2*(jme**2))+(l3*(jme**3))+(l4*(jme**4))+(l5*(jme**5)))/(10**8)
   l_unlimit = l_rad * 180 / pi
   if l_unlimit > 0:
       l = l_unlimit%360
   elif l_unlimit < 0:
       l = 360-(abs(l_unlimit)%360)

   #b is Earth heliocentric latitude
   b_all = B_EPT(jme)
   b0, b1 = b_all[0], b_all[1]
   b_rad = (b0+(b1*jme))/(10**8)
   b = b_rad * 180 / pi

   #r is Earth radius vector
   r_all = R_EPT(jme)
   r0, r1, r2, r3, r4 = r_all[0], r_all[1], r_all[2], r_all[3], r_all[4]
   r = (r0+(r1*jme)+(r2*(jme**2))+(r3*(jme**3))+(r4*(jme**4)))/(10**8)

   #theta is geocentric longitude
   theta_unlimit = l + 180
   if theta_unlimit > 0:
       theta = theta_unlimit%360
   elif theta_unlimit < 0:
       theta = 360-(abs(theta_unlimit)%360)
   #beta is geocentric latitude
   beta = b*(-1)

   #x0 is mean elongation of moon from sun
   x0 = 297.85036+445267.111480*jce-0.0019142*(jce**2)+((jce**3)/189474)
   #x1 is mean anomaly of sun
   x1 = 357.52772+35999.050340*jce-.0001603*(jce**2)-((jce**3)/300000)
   #x2 is mean anomaly of moon
   x2 = 134.96298+477198.867398*jce+.0086972*(jce**2)+((jce**3)/56250)
   #x3 is moon's argument of latitude
   x3 = 93.27191+483202.017538*jce-.0036825*(jce**2)+((jce**3)/327270)
   #x4 is longitude of ascending node of mean moon orbit
   x4 = 125.04452-1934.136261*jce+.0020708*(jce**2)+((jce**3)/450000)
   #Nutation in longitude and obliquity
   nu_vals = nutation(jce, x0, x1, x2, x3, x4)
   delta_psi = nu_vals[0]
   delta_epsilon = nu_vals[1]
   u = jme/10
   epsilon0 = 84381.448-(4680.93*u)-(1.55*(u**2))+(1999.25*(u**3))-\
              (51.38*(u**4))-(249.67*(u**5))-(39.05*(u**6))+(7.12*(u**7))+\
              (27.87*(u**8))+(5.79*(u**9))+(2.45*(u**10))
   #true obliquity of ecliptic in degrees
   epsilon = (epsilon0/3600)+delta_epsilon
   #delta_tao is aberration correction in degrees
   delta_tao = 20.4898/(3600*r)
   #lamb is apparent sun longitude in degrees
   lamb = theta + delta_psi + delta_tao

   v0_unlimit = 280.46061837+360.98564736629*(jd-2451545)+0.000387933*(jc**2)-((jc**3)/38710000)
   if v0_unlimit > 0:
       v0 = v0_unlimit%360
   elif v0_unlimit < 0:
       v0 = 360-(abs(v0_unlimit)%360)
   #v is apparent sidereal time at Greenwich in degrees
   v = v0+(delta_psi*cos(radians(epsilon)))

   numer = sin(radians(lamb))*cos(radians(epsilon))-tan(radians(beta))*sin(radians(epsilon))
   denom = cos(radians(lamb))
   alpha_rad = atan2(numer, denom)
   alpha_unlimit = alpha_rad*180/pi
   #alpha is geocentric sun right ascension
   if alpha_unlimit > 0:
       alpha = alpha_unlimit%360
   elif alpha_unlimit < 0:
       alpha = 360-(abs(alpha_unlimit)%360)

   #delta is geocentric sun declination in degrees
   delta_rad = asin((sin(radians(beta))*cos(radians(epsilon)))+(cos(radians(beta))*sin(radians(epsilon))*sin(radians(lamb))))

   delta = delta_rad*180/pi
   #h is observer local hour angle in degrees
   h_unlimit = v+sigma-alpha
   if h_unlimit > 0:
       h = h_unlimit%360
   elif h_unlimit < 0:
       h = 360-(abs(h_unlimit)%360)
   #xi is equatorial parallax in degrees
   xi = 8.794/(3600*r)
   u = atan(0.99664719*tan(radians(phi)))
   x = cos(u)+((e/6378140)*cos(radians(phi)))
   y = (0.99664719*sin(u))+((e/6378140)*sin(radians(phi)))
   numer2 = -1*x*sin(radians(xi))*sin(radians(h))
   denom2 = cos(radians(delta))-(x*sin(radians(xi))*cos(radians(h)))
   delta_alpha_rad = atan2(numer2, denom2)
   delta_alpha = delta_alpha_rad*180/pi
   #alpha_p is topocentric sun right ascension in degrees
   alpha_p = alpha+delta_alpha

   numer3 = (sin(radians(delta))-(y*sin(radians(xi))))*cos(radians(delta_alpha))
   denom3 = cos(radians(delta))-(x*(sin(radians(xi)))*cos(radians(h)))
   #delta_p is topocentric sun declination in degrees
   delta_p_rad = atan2(numer3, denom3)
   delta_p = delta_p_rad*180/pi

   #h_p is topocentric local hour angle in degrees
   h_p = h - delta_alpha

   elev_ang_rad = asin((sin(radians(phi))*sin(radians(delta_p)))+\
             (cos(radians(phi))*cos(radians(delta_p))*cos(radians(h_p))))
   elev_ang = elev_ang_rad*180/pi
   topo_zenith = 90-elev_ang
   #Gamma is topocentric astronomers azimuth angle in degrees
   numer4 = sin(radians(h_p))
   denom4 = (cos(radians(h_p))*sin(radians(phi)))-(tan(radians(delta_p))*cos(radians(phi)))
   Gamma_rad = atan2(numer4, denom4)
   Gamma_unlimit = Gamma_rad*180/pi
   if Gamma_unlimit > 0:
       Gamma = Gamma_unlimit%360
   elif Gamma_unlimit < 0:
       Gamma = 360-(abs(Gamma_unlimit)%360)
   #Phi is topocentric azimuth angle for navigation and radiation
   Phi_unlimit = Gamma + 180
   if Phi_unlimit > 0:
       Phi = Phi_unlimit%360
   elif Phi_unlimit < 0:
       Phi = 360-(abs(Phi_unlimit)%360)

   shad_len = height/tan(radians(elev_ang))
   if elev_ang <0:
       shad_len = 0
   #north, south, east, west = 0, 90, 180, 270
   if Phi == 0:
       shad_angle = 'south'
   elif Phi <= 45:
       shad_angle = str(Phi)+' degrees west of south'
   elif Phi > 45 and Phi<=135:
       if Phi == 90:
           shad_angle = 'west'
       elif Phi < 90:
           shad_angle = str(0-Phi)+' degrees south of west'
       else:
          shad_angle = str(Phi-90)+' degrees north of west'
   elif Phi>135 and Phi<=225:
       if Phi == 180:
           shad_angle = 'north'
       elif Phi < 180:
           shad_angle = str(180-Phi)+' degrees west of north'
       else:
          shad_angle = str(Phi-180)+' degrees east of north'   
   elif Phi > 225 and Phi<=315:
       if Phi == 270:
           shad_angle = 'east'
       elif Phi < 270:
           shad_angle = str(270-Phi)+' degrees north of east'
       else:
          shad_angle = str(Phi-270)+' degrees south of east'
   else:
       shad_angle = str(360-Phi) + ' degrees east of south'

   return [shad_len, elev_ang, shad_angle, Phi]


sig = float(input("Please enter your geographical longitude (positive for east, negative for west)\n"))
phii = float(input("Please enter your geographical latitude (positive for north, negtive for south)\n"))
obj_height = input("Please enter your height in meters. If you would like the default value of\n1.7 meters, hit enter.\n")
if obj_height == '':
    obj_height = 1.7
else:
    obj_height = float(obj_height)
elev = float(input("Please enter your elevation in meters\n"))

choice = input("Enter 'Now' to calculate the sun's position for the current time and date,\
or enter 'Custom' to calculate the position for any time and date.\n")
from time import gmtime, strftime
if choice.capitalize() == 'Now':
   month = int(strftime("%m", gmtime()))
   day_month = int(strftime("%d", gmtime()))
   year = int(strftime("%Y", gmtime()))
   time = strftime("%H:%M:%S", gmtime())

elif choice.capitalize() == 'Custom':
   year = int(input("What is the year? \n"))
   month = input("What is the month? \n")
   months = {"January": 1, "February": 2, "March": 3, "April": 4, "May": 5,"June": 6\
             ,"July": 7, "August": 8, "September": 9, "October": 10, "November":11\
             , "December": 12}
   try:
       typee = int(month)
       month = typee
   except ValueError:
       month = months[month.capitalize()]
   day_month = int(input("What is the day of the month? \n"))
   time = input("What is the universal time? Seperate hour, minute and seconds with colons (e.g.\
   One thirty and thirty seconds PM is written as 1:30:30PM \n")    

time = time.split(":")
print(time)
hour = int(time[0])
minute = int(time[1])
sec_info = list(time[2])
sec = int(''.join(sec_info[:2]))

if choice.capitalize() == 'Custom':
   if sec_info[2]==' ':
       det = sec_info[3]
   else:
       det = sec_info[2]
   if det.upper() == 'A':
       if hour == 12:
           hour = 0
   elif det.upper() == 'P' and hour != 12:
       hour += 12
SECS_IN_DAY = 24*3600
totsec = (hour*3600) + (minute*60) + sec
deciday = totsec/SECS_IN_DAY

ans_info = shadow(year, month, day_month, deciday, sig, phii, obj_height, elev)
output = "Elevation: {} degrees\nAzimuth: {} degrees Eastward from North\nDirection of Shadow: {}\nShadow Length: {} meters".format(ans_info[1], ans_info[3], ans_info[2], ans_info[0])
if float(ans_info[1])<0:
   output += "\nNighttime"
print(output)
from math import cos, sin, radians
import numpy as np
#import ffmpeg
from matplotlib import pyplot as plt
from IPython.display import Image
from celluloid import Camera
fig = plt.figure()
camera = Camera(fig)

timespan = np.linspace(0, 1, 100)
acc=0
sunriseset = []
e_ang = []
for i in timespan:

   arg = i
   shadoww = shadow(year, month, day_month, arg, sig, phii, obj_height, elev)
   shad_len = float(shadoww[0])
   elev_angle = float(shadoww[1])
   shad_ang_polar = float(shadoww[3])
   totsecs = arg*SECS_IN_DAY
   totmins = int(totsecs//60)
   secs = int(totsecs%60)
   hrs = int(totmins//60)
   mins = int(totmins%60)
   if abs(abs(elev_angle)-.8333) <= 1:
       sunriseset == sunriseset.append("{}:{}:{}".format(hrs, mins, secs))
       e_ang == e_ang.append(elev_angle)
   if elev_angle >= 0:
       acc += 1
       xcoord = shad_len*cos(radians(90-shad_ang_polar))
       ycoord = shad_len*sin(radians(90-shad_ang_polar))
       if acc== 1:
           newtext = "Shadow at Time: {}:{}:{}".format(hrs, mins, secs)
           plt.plot([0, xcoord], [0, ycoord], 'k', label="Shadow")
           plt.text(xcoord+10, ycoord+10, newtext)
       else:
           newtext = "Shadow at Time: {}:{}:{}".format(hrs, mins, secs)
           plt.plot([0, xcoord], [0, ycoord], 'k')
           plt.text(xcoord+10, ycoord+10, newtext)
       camera.snap()
   else:
       newtext = "Nightime: {}:{}:{}".format(hrs, mins, secs)
       plt.text(0, 40, newtext)
       camera.snap()
title = "Shadow length and direction for 24 hrs\nStarting at 12:00am on {}/{}/{}".format(month, day_month, year)
plt.title(title)
title = "Shadow length and direction for 24 hrs\nStarting at 00:00:00 on {}/{}/{}".\
       format(month, day_month, year)
plt.title(title)
plt.xlabel("Bird's Eye View - Positive is East, Negative is West")
plt.ylabel("Bird's Eye View\nPositive is North, Negative is South")
plt.legend(loc="upper right")
animation = camera.animate()
animation.save('Shadow.mp4', writer = 'ffmpeg')
