# Shadow-Calculator

How To Use This Code

This program takes the users geographical longitude and latitude, height in meters
with a default value of 1.7, elevation in meters,  current or
customized time in Universal Time (UT) or Greenwich Mean Time (GMT)*, and current or customized
date. More on each input is described below:

    Geographical latitude and longitude: Note that this asks for geographic
    coordinates instead of geodectic (both are commonly used). North and East are
    positive, South and West are negative

    Height: This is entered in meters, and users can hit enter to use 1.7 meters
    (average height of human)

    Elevation: This should be taken above sea level

    Time: If the user enters "Now", the current UT time will be entered as the time.
    If the user would prefer to enter a different time, they may do so by typing "Custom",
    which will lead to prompts asking for the year, month (enter either month name or number,
    with Janaury = 1, February = 2, etc), day of the month (should be a numerical value),
    and time (it uses a 12 hr clock, so AM or PM must be specified)

The program will return the topographic elevation angle (measured from horizon to sun's
geometric center), the topographic azimuth angle for astronomers**(measured eastward from north), the shadow angle
along with cardinal direction, and shadow length.

*The difference between UT and GMT is miniscule, never more than .9 seconds
**Note that the azimuth for astronomers and the azimuth for navigators is the same, they just have
different reference points (astronomers is Eastward from North, navigators is Westward from South)

How My Project Works

I mainly followed Ibrahim Reda and Afshin Andreas's work in their guide "Solar Position Algorithm for
Solar Radiation Applications". Because calculation were used upon large datatables, I created text files
containing the Earth Periodic Terms for L, B, and R, along with the Periodic Terms in Nutation in Logitude
and Obliquity. These textfiles are accompanied by functions that perform the necessary operations on them,
EPT for L, B_EPT for B, R_EPT for R, and nutation for nutation. First, the function shadow calculates
the Julian Date, then it goes on to calculate different astronomical terms, all of which are labelled
in the program. I deviated from "Solar Position Algorithm for Solar Radiation Applications" when it came
to calculating delta_t, which is the difference between Terrestrial Time (TT) and UT or GMT. This value is only measured
by observation, so I used NASA's polynomial approximation of
 t = y -2000
 delta_t = 63.86+(0.3345*t)-(0.060374*(t**2))+(0.0017275*(t**3))+(0.000651814*(t**4))+(0.00002373599*(t**5))
where y is the year. I also used basic trigonometry to calculate shadow length as height/tan(angle of elevation).
Angle of elevation is negative at night, so length will be 0 when this occurs instead of negative. Shadow angle is simply
opposite to the azimuth angle.

    

    
