import matplotlib.pyplot as plt
import numpy as np
import datetime
import pandas as pd

# 기타 참조?
# https://pyorbital.readthedocs.io/en/feature-moon-phase/index.html#

# Astronomical Coordinate System
# https://docs.astropy.org/en/stable/coordinates/index.html#convenience-methods
# Solar System Ephemerides
# https://docs.astropy.org/en/stable/coordinates/solarsystem.html
# Skycoord Concept
# https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html
# https://9574m.tistory.com/30

import astropy

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body_barycentric, get_body
from astropy.time import Time
from astropy import units as u
import math

from astropy.coordinates import solar_system_ephemeris


def get_altaz(lat, lon, time):
    solar_system_ephemeris.set('de432s')
    loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
    # t = Time("2020-01-01 00:00")
    # t = Time.strptime(time, '%Y-%m-%d %H:%M:%S')
    t = Time(time)
    print(t)
    # calc utc offset
    utc_offset = 9*u.hour
    t = t - utc_offset

    sun = get_body('sun', t, loc)
    result = str(sun.transform_to(astropy.coordinates.AltAz(obstime=t,location=loc)))
    # result = sun.transform_to(astropy.coordinates.AltAz(obstime=t,location=loc))
    # alt = result.alt
    # az = result.az
    # print(result)
    result = result.split('(')[-1].replace(')>', '').split(', ')
    alt = float(result[1])
    az = float(result[0])
    return alt, (az+180)%360
    # print(calc_theta(alt, az, 45, 190))



# 2023.01.12 용권이 일사각 구하기 코드
def get_tau(alt, type, zenith):

    A = alt
    type_dict = {
        'tropical':[0.95, 0.98, 1.02],
        'midlatitude_summer': [0.97, 0.99, 1.02],
        'subarctic_summer': [0.99, 0.99, 1.01],
        'midlatitude_winter': [1.03, 1.01, 1.00]
    }
    r_0 = type_dict[type][0]
    r_1 = type_dict[type][1]
    r_k = type_dict[type][2]
    theta_z = zenith


    a_0 = (0.4237 - 0.00821(6-A)**2)*r_0
    a_1 = (0.5055 + 0.00595(6.5-A)**2)*r_1
    k = (0.2711 + 0.01858(2.5-A)**2)*r_k
    tau_b = a_0 + a_1 * math.exp(-k/math.cos(math.radians(theta_z)))
    return tau_b

def get_solar_time(lon, B, time):
    # Standard meridian 계산
    Lst = (int((lon - 7.5) / 15)+1)*15
    print(lon, Lst)
    # Lst = 135
    print(B)
    E = 229.2*(0.000075 + 0.001868 * math.cos(math.radians(B)) - 0.032077 * math.sin(math.radians(B)) - 0.014615 * math.cos(math.radians(2*B)) - 0.04089 * math.sin(math.radians(2*B)))
    print(E)
    delta_min = 4 * (Lst - lon) + E # L_loc = lon
    print(delta_min)
    solar_time = time + datetime.timedelta(minutes=delta_min)
    return solar_time

def get_hour_angle(solar_time):
    return -180 + 15 * (solar_time.hour + solar_time.minute/60 + solar_time.second/3600)

def get_hourly_beam_radiation(n, lat, lon, altitude, climate, time, face_az, slope, reflectance, solar_time = False):
    coeff_dict = {
    }
    print(n)
    B = (n-1)*360/365
    G_sc = 1367 # W/m^2
    # G_on = G_sc * (1.100110 + 0.034221*math.cos(math.radians(B)) + 0.001280*math.sin(math.radians(B)) + 0.000719*math.cos(math.radians(2*B)) + 0.000077*math.sin(math.radians(2*B)))
    G_on = G_sc * (1 + 0.033 * math.cos(math.radians(360*n/365)))

    # solar time 계산
    if solar_time == True:
        solar_time = time
    else:
        solar_time = get_solar_time(lon, B, time)
    print(solar_time)

    # solar time으로 hour angle 계산
    hour_angle = get_hour_angle(solar_time)

    print(hour_angle)

    # declination 계산 (delta)
    declination = 23.45 * math.sin(math.radians(360 * (284 + n) / 365))
    # declination = (180/math.pi) * (0.006918 - 0.399912 * math.cos(math.radians(B)) + 0.070257 * math.sin(math.radians(B)) - 0.006758 * math.cos(math.radians(2*B)) + 0.000907 * math.sin(math.radians(2*B)) - 0.002697 * math.cos(math.radians(3*B)) + 0.00148 * math.sin(math.radians(3*B)))
    print('declination')
    print(declination)

    # incidence angle 계산 (theta) - 지금은 안 씀
    # incidence_angle = math.degrees(math.acos(
    #     math.sin(math.radians(declination)) * math.sin(math.radians(lat)) * math.cos(math.radians(slope))
    #     - math.sin(math.radians(declination)) * math.cos(math.radians(lat)) * math.sin(math.radians(slope)) * math.cos(math.radians(azimuth_surface))
    #     + math.cos(math.radians(declination)) * math.cos(math.radians(lat)) * math.cos(math.radians(slope)) * math.cos(math.radians(hour_angle))
    #     + math.cos(math.radians(declination)) * math.sin(math.radians(lat)) * math.sin(math.radians(slope)) * math.cos(math.radians(azimuth_surface)) * math.cos(math.radians(hour_angle))
    #     + math.cos(math.radians(declination)) * math.sin(math.radians(slope)) * math.sin(math.radians(azimuth_surface)) * math.sin(math.radians(hour_angle))
    # ))

    # solar zenith angle 계산 (slope가 없다고 가정했을 때의 incidence angle)
    zenith_angle = math.degrees(math.acos(
      math.cos(math.radians(lat)) * math.cos(math.radians(declination)) * math.cos(math.radians(hour_angle)) + math.sin(math.radians(lat)) * math.sin(math.radians(declination))
    ))
    # horizontal surface
    cos_zen = math.cos(math.radians(lat)) * math.cos(math.radians(declination)) * math.cos(math.radians(hour_angle)) + math.sin(math.radians(lat)) * math.sin(math.radians(declination))
    if cos_zen < 0:
        coeff_dict['n'] = n
        coeff_dict['B'] = B
        coeff_dict['G_on'] = G_on
        coeff_dict['solar_time'] = solar_time
        coeff_dict['hour_angle'] = hour_angle
        coeff_dict['declination'] = declination
        coeff_dict['cos_zen'] = cos_zen
        coeff_dict['I_slope_beam'] = 0
        coeff_dict['I_slope_diffuse'] = 0
        coeff_dict['I_slope_reflect'] = 0
        coeff_dict['I_slope_total'] = 0
        coeff_dict['R_b_ave'] = 0
        coeff_dict['start_hour_angle'] = 0
        coeff_dict['end_hour_angle'] = 0

        return 0, 0, coeff_dict
    # vertical surface
    # eq 1.6.4 (지금은 안 만듬)

    print('zenith angle')
    print(cos_zen)

    # G_on에 tau를 곱해서 직달일사량을 구함
    if altitude > 2500:
        print('altitude exceeds range (2.5km)')
        raise IOError
    a_0_star = 0.4237 - 0.00821 * (6 - altitude/1000)**2
    a_1_star = 0.5055 + 0.00595 * (6.5 - altitude/1000) ** 2
    k_star = 0.2711 + 0.01858 * (2.5 - altitude/1000) ** 2
    print(a_0_star, a_1_star, k_star)
    a_coeff_dict = {
        'tropical':(0.95, 0.98, 1.02),
        'midlat_summer': (0.97, 0.99, 1.02),
        'subarc_summer': (0.99, 0.99, 1.01),
        'midlat_winter': (1.03, 1.01, 1.00)
    }
    a_0 = a_0_star * a_coeff_dict[climate][0]
    a_1 = a_1_star * a_coeff_dict[climate][1]
    k = k_star * a_coeff_dict[climate][2]
    print(a_0, a_1, k)
    tau_b = a_0 + a_1 * math.exp(-k/cos_zen)
    print('tau_b')
    print(tau_b)
    print('G_on')
    print(G_on)


    G_cnb = G_on * tau_b
    # G_cb = G_on * math.cos(math.radians(incidence_angle)) * tau_b * math.cos(math.radians(zenith_angle))
    # 직달일사 계산
    G_cb = G_cnb * cos_zen

    print('Gcb')
    print(G_cb)
    # 확산일사 tau_d 계산
    tau_d = 0.271-0.294*tau_b
    G_cd = G_on * cos_zen * tau_d

    # sunrise angle, sunset angle 계산 : 1.6.5
    w = math.degrees(math.acos(
        -(math.sin(math.radians(lat))*math.sin(math.radians(declination))/math.cos(math.radians(lat))/math.cos(math.radians(declination)))
    ))
    sunrise_angle = -w
    sunset_angle = w

    # 경사면의 경우 계산, 해당 없으면 0 출력
    if slope == 0:
        # 해당없는 부분 출력
        I_slope_beam = 0
        I_slope_diffuse = 0
        I_slope_reflect = 0
        I_slope_total = 0
        R_b_ave = 0
        start_hour_angle = 0
        end_hour_angle = 0
    else:
        # 평균 직달일사비 경사계수 R_b_ave 구하기 (1시간 평균값)
        start_hour_angle = get_hour_angle(solar_time - datetime.timedelta(minutes=30))
        end_hour_angle = get_hour_angle(solar_time + datetime.timedelta(minutes=30))
        # 일사/일출 시간대의 hour angle 계산
        if sunrise_angle > start_hour_angle and sunrise_angle < end_hour_angle:
            start_hour_angle = sunrise_angle
        elif sunset_angle > start_hour_angle and sunset_angle < end_hour_angle:
            end_hour_angle = sunset_angle
        hourangle_diff = end_hour_angle - start_hour_angle
        hourangle_sin_diff = math.sin(math.radians(end_hour_angle)) - math.sin(math.radians(start_hour_angle))
        hourangle_cos_diff = math.cos(math.radians(end_hour_angle)) - math.cos(math.radians(start_hour_angle))

        print(start_hour_angle, end_hour_angle)

        a = (
                (math.sin(math.radians(declination))*math.sin(math.radians(lat))*math.cos(math.radians(slope))
                - math.sin(math.radians(declination))*math.cos(math.radians(lat))*math.sin(math.radians(slope))*math.cos(math.radians(face_az)))
                *(1/180)*math.pi*hourangle_diff
                + (math.cos(math.radians(declination))*math.cos(math.radians(lat))*math.cos(math.radians(slope))
                + math.cos(math.radians(declination))*math.sin(math.radians(lat))*math.sin(math.radians(slope))*math.cos(math.radians(face_az)))
                *hourangle_sin_diff
                - math.cos(math.radians(declination))*math.sin(math.radians(slope))*math.sin(math.radians(face_az))
                *hourangle_cos_diff
             )
        b = (
                math.cos(math.radians(lat))*math.cos(math.radians(declination))*hourangle_sin_diff + math.sin(math.radians(lat))*math.sin(math.radians(declination))*(1/180)*math.pi*hourangle_diff
            )
        R_b_ave = a/b
        print('R_b')
        print(R_b_ave)
        print(a, b)
        print(G_cb, G_cd)

        I_slope_beam = G_cb * R_b_ave
        I_slope_diffuse = G_cd * ((1 + math.cos(math.radians(slope)))/2)
        I_slope_reflect = (G_cb + G_cd) * reflectance * ((1 - math.cos(math.radians(slope)))/2)

        # 직달일사가 음수인 경우 (해가 판때기 뒤편에 있는 경우) 0으로 바꿈
        if I_slope_beam < 0:
            I_slope_beam = 0

        I_slope_total = I_slope_beam + I_slope_diffuse + I_slope_reflect


    coeff_dict['I_slope_beam'] = I_slope_beam
    coeff_dict['I_slope_diffuse'] = I_slope_diffuse
    coeff_dict['I_slope_reflect'] = I_slope_reflect
    coeff_dict['I_slope_total'] = I_slope_total
    coeff_dict['R_b_ave'] = R_b_ave
    coeff_dict['start_hour_angle'] = start_hour_angle
    coeff_dict['end_hour_angle'] = end_hour_angle



    coeff_dict['n'] = n
    coeff_dict['B'] = B
    coeff_dict['G_on'] = G_on
    coeff_dict['solar_time'] = solar_time
    coeff_dict['hour_angle'] = hour_angle
    coeff_dict['declination'] = declination
    coeff_dict['cos_zen'] = cos_zen

    print(coeff_dict)

    return G_cb, G_cd, coeff_dict

def get_vertical_insolation(lat, lon, day, hour, zenith, cos_zen):
    import math
    longitude = lon
    latitude = lat
    doy = day

    # # Constants
    # latitude = float(input("Enter latitude: "))
    # doy = int(input("Enter day of year: "))
    # hour = int(input("Enter hour of day: "))

    # Perez Model parameters
    a = [0.409, 0.5016, 0.6309]
    b = [-0.0608, -0.0388, -0.0040]
    c = [0.0007, 0.0021, 0.0027]
    d = [-0.0685, -0.0406, -0.0193]

    # Solar time calculation
    longitude_correction = 4 * (longitude - 15 * (hour - 12))
    solar_time = hour + 0.25 * (longitude_correction - latitude)

    # Zenith angle calculation
    # zenith = math.acos(math.sin(latitude) * math.sin(math.pi / 180 * solar_time) + math.cos(latitude) * math.cos(
    #     math.pi / 180 * solar_time) * math.cos(math.pi / 180 * solar_time))

    if math.cos(zenith) < 0:
        return 0
    print('cos_zennnn')
    print(math.degrees(zenith))
    print(math.cos(zenith))

    # Perez model calculation
    diffuse_fraction = a[0] + b[0] * math.exp(-c[0] / math.cos(zenith)) + d[0] * math.exp(
        -1.41 / math.pow(math.cos(zenith), 1.4))
    direct_fraction = a[1] + b[1] * math.exp(-c[1] / math.cos(zenith)) + d[1] * math.exp(
        -1.41 / math.pow(math.cos(zenith), 1.4))
    sky_diffuse_fraction = a[2] + b[2] * math.exp(-c[2] / math.cos(zenith)) + d[2] * math.exp(
        -1.41 / math.pow(math.cos(zenith), 1.4))

    print(diffuse_fraction, direct_fraction)

    # Extraterrestrial radiation
    Ra = 1367 * (1 + 0.033 * math.cos(2 * math.pi * doy / 365))
    # I_on = 1367 * (1 + 0.033 * math.cos(math.radians(360 * n / 365)))

    # Clear sky solar irradiance on a vertical surface
    I = Ra * diffuse_fraction + Ra * direct_fraction * math.sin(zenith)

    # print("Hourly Insolation on a vertical surface (Perez Model): {:.2f} W/m^2".format(I))
    return I

def calc_earth_insolation():
    df_dict = {
        'time': [],
        'beam': [],
        'diff': [],
        'total': [],
        'n': [],
        'B': [],
        'G_on': [],
        'solar_time': [],
        'hour_angle': [],
        'declination': [],
        'cos_zen': [],
        'altitude': [],
        'azimuth': [],
        'I_slope_beam': [],
        'I_slope_diffuse': [],
        'I_slope_reflect': [],
        'I_slope_total': [],
        'R_b_ave': [],
        'start_hour_angle': [],
        'end_hour_angle': []
    }

    # 조건값 입력
    lat = 37.5665
    lon = 126.9780
    alt = 0
    slope = 90
    reflectance = 0.6
    face_az = 0
    climate = 'midlat'
    solar_time = False
    basetime = datetime.datetime(2021, 1, 1, 0, 0, 0)

    # climate는 아래 값 중 하나로 입력
    # 'tropical', 'midlat', 'subarc_summer'

    time = basetime
    for i in range(8760):
    # basetime = datetime.datetime(2021, 3, 4, 6, 30, 0)
    # time = basetime
    # for i in range(7):
        n = ((time-basetime).days) + 1
        print(time)
        print(n)
        alt, az = get_altaz(lat, lon, str(time))
        df_dict['altitude'].append(alt)
        df_dict['azimuth'].append(az)
        # 여름/겨울 처리
        if climate == 'midlat':
            if time.month >= 4 and time.month <= 10:
                climate = 'midlat_summer'
            else:
                climate = 'midlat_winter'
        # get_hourly_beam_radiation(n, lat, lon, altitude, climate, time, face_az, slope, reflectance)
        beam, diff, coeff_dict = get_hourly_beam_radiation(n, lat, lon, alt, climate, time, face_az, slope, reflectance, solar_time=solar_time)

        total = beam + diff
        df_dict['time'].append(time)
        df_dict['beam'].append(beam)
        df_dict['diff'].append(diff)
        df_dict['total'].append(total)
        for key in coeff_dict.keys():
            df_dict[key].append(coeff_dict[key])


        day = n
        hour = time.hour
        # print(df_dict['cos_zen'])
        zenith = math.acos(df_dict['cos_zen'][-1])

        # df_dict['insolation_vertical'].append(get_vertical_insolation(lat, lon, day, hour, zenith, 0))

        time = time + datetime.timedelta(hours=1)

    for key in df_dict.keys():
        print(key)
        print(len(df_dict[key]))

    res = pd.DataFrame(df_dict)
    res.to_excel('./result_V12.xlsx')
    res.to_csv('./result_V12.csv', encoding='ANSI')
    raise IOError



    time = datetime.datetime(2021, 8, 22, 11, 30, 0)
    # (n, slope, lat, lon, altitude, climate, azimuth_surface, time)
    print(get_hourly_beam_radiation(234, 43.071934640533165, 89.40120883118382, 270, 'midlat_summer', time))
    raise IOError

    def calc_theta(alt, az, tilt, rotate):
        res = math.degrees(math.acos(math.cos(math.radians(az-rotate))/math.cos(math.radians(alt-tilt))))
        return res

    # 43.071934640533165, -89.40120883118382





    raise IOError










#

def norm_deg(deg):
    deg = deg % 360
    if deg < 0:
        deg = deg + 360
    return deg

def print_dict(dict_var, name = None):
    # 스피드업 비활성화
    i = 0
    # print('=================================================')
    # if name != None:
    #     print('values for ' + name)
    # for key in dict_var.keys():
    #     print(str(key) + '\t:\t' + str(dict_var[key]))
    # print('=================================================')

def make_periodic_terms_nutation():
    # periodic_terms_copied = [{"arg":[0,0,0,0,1],"sin":[-171996,-174.2],"cos":[92025,8.9]},{"arg":[-2,0,0,2,2],"sin":[-13187,-1.6],"cos":[5736,-3.1]},{"arg":[0,0,0,2,2],"sin":[-2274,-0.2],"cos":[977,-0.5]},{"arg":[0,0,0,0,2],"sin":[2062,0.2],"cos":[-895,0.5]},{"arg":[0,1,0,0,0],"sin":[1426,-3.4],"cos":[54,-0.1]},{"arg":[0,0,1,0,0],"sin":[712,0.1],"cos":[-7,0]},{"arg":[-2,1,0,2,2],"sin":[-517,1.2],"cos":[224,-0.6]},{"arg":[0,0,0,2,1],"sin":[-386,-0.4],"cos":[200,0]},{"arg":[0,0,1,2,2],"sin":[-301,0],"cos":[129,-0.1]},{"arg":[-2,-1,0,2,2],"sin":[217,-0.5],"cos":[-95,0.3]},{"arg":[-2,0,1,0,0],"sin":[-158,0],"cos":[0,0]},{"arg":[-2,0,0,2,1],"sin":[129,0.1],"cos":[-70,0]},{"arg":[0,0,-1,2,2],"sin":[123,0],"cos":[-53,0]},{"arg":[2,0,0,0,0],"sin":[63,0],"cos":[0,0]},{"arg":[0,0,1,0,1],"sin":[63,0.1],"cos":[-33,0]},{"arg":[2,0,-1,2,2],"sin":[-59,0],"cos":[26,0]},{"arg":[0,0,-1,0,1],"sin":[-58,-0.1],"cos":[32,0]},{"arg":[0,0,1,2,1],"sin":[-51,0],"cos":[27,0]},{"arg":[-2,0,2,0,0],"sin":[48,0],"cos":[0,0]},{"arg":[0,0,-2,2,1],"sin":[46,0],"cos":[-24,0]},{"arg":[2,0,0,2,2],"sin":[-38,0],"cos":[16,0]},{"arg":[0,0,2,2,2],"sin":[-31,0],"cos":[13,0]},{"arg":[0,0,2,0,0],"sin":[29,0],"cos":[0,0]},{"arg":[-2,0,1,2,2],"sin":[29,0],"cos":[-12,0]},{"arg":[0,0,0,2,0],"sin":[26,0],"cos":[0,0]},{"arg":[-2,0,0,2,0],"sin":[-22,0],"cos":[0,0]},{"arg":[0,0,-1,2,1],"sin":[21,0],"cos":[-10,0]},{"arg":[0,2,0,0,0],"sin":[17,-0.1],"cos":[0,0]},{"arg":[2,0,-1,0,1],"sin":[16,0],"cos":[-8,0]},{"arg":[-2,2,0,2,2],"sin":[-16,0.1],"cos":[7,0]},{"arg":[0,1,0,0,1],"sin":[-15,0],"cos":[9,0]},{"arg":[-2,0,1,0,1],"sin":[-13,0],"cos":[7,0]},{"arg":[0,-1,0,0,1],"sin":[-12,0],"cos":[6,0]},{"arg":[0,0,2,-2,0],"sin":[11,0],"cos":[0,0]},{"arg":[2,0,-1,2,1],"sin":[-10,0],"cos":[5,0]},{"arg":[2,0,1,2,2],"sin":[-8,0],"cos":[3,0]},{"arg":[0,1,0,2,2],"sin":[7,0],"cos":[-3,0]},{"arg":[-2,1,1,0,0],"sin":[-7,0],"cos":[0,0]},{"arg":[0,-1,0,2,2],"sin":[-7,0],"cos":[3,0]},{"arg":[2,0,0,2,1],"sin":[-7,0],"cos":[3,0]},{"arg":[2,0,1,0,0],"sin":[6,0],"cos":[0,0]},{"arg":[-2,0,2,2,2],"sin":[6,0],"cos":[-3,0]},{"arg":[-2,0,1,2,1],"sin":[6,0],"cos":[-3,0]},{"arg":[2,0,-2,0,1],"sin":[-6,0],"cos":[3,0]},{"arg":[2,0,0,0,1],"sin":[-6,0],"cos":[3,0]},{"arg":[0,-1,1,0,0],"sin":[5,0],"cos":[0,0]},{"arg":[-2,-1,0,2,1],"sin":[-5,0],"cos":[3,0]},{"arg":[-2,0,0,0,1],"sin":[-5,0],"cos":[3,0]},{"arg":[0,0,2,2,1],"sin":[-5,0],"cos":[3,0]},{"arg":[-2,0,2,0,1],"sin":[4,0],"cos":[0,0]},{"arg":[-2,1,0,2,1],"sin":[4,0],"cos":[0,0]},{"arg":[0,0,1,-2,0],"sin":[4,0],"cos":[0,0]},{"arg":[-1,0,1,0,0],"sin":[-4,0],"cos":[0,0]},{"arg":[-2,1,0,0,0],"sin":[-4,0],"cos":[0,0]},{"arg":[1,0,0,0,0],"sin":[-4,0],"cos":[0,0]},{"arg":[0,0,1,2,0],"sin":[3,0],"cos":[0,0]},{"arg":[0,0,-2,2,2],"sin":[-3,0],"cos":[0,0]},{"arg":[-1,-1,1,0,0],"sin":[-3,0],"cos":[0,0]},{"arg":[0,1,1,0,0],"sin":[-3,0],"cos":[0,0]},{"arg":[0,-1,1,2,2],"sin":[-3,0],"cos":[0,0]},{"arg":[2,-1,-1,2,2],"sin":[-3,0],"cos":[0,0]},{"arg":[0,0,3,2,2],"sin":[-3,0],"cos":[0,0]},{"arg":[2,-1,0,2,2],"sin":[-3,0],"cos":[0,0]}]
    periodic_terms_copied = [{"arg":[0,0,0,0,1],"sin":[-171996,-174.2],"cos":[92025,8.9]},{"arg":[-2,0,0,2,2],"sin":[-13187,-1.6],"cos":[5736,-3.1]},{"arg":[0,0,0,2,2],"sin":[-2274,-0.2],"cos":[977,-0.5]},{"arg":[0,0,0,0,2],"sin":[2062,0.2],"cos":[-895,0.5]},{"arg":[0,1,0,0,0],"sin":[1426,-3.4],"cos":[54,-0.1]},{"arg":[0,0,1,0,0],"sin":[712,0.1],"cos":[-7,0]},{"arg":[-2,1,0,2,2],"sin":[-517,1.2],"cos":[224,-0.6]},{"arg":[0,0,0,2,1],"sin":[-386,-0.4],"cos":[200,0]},{"arg":[0,0,1,2,2],"sin":[-301,0],"cos":[129,-0.1]},{"arg":[-2,-1,0,2,2],"sin":[217,-0.5],"cos":[-95,0.3]},{"arg":[-2,0,1,0,0],"sin":[-158,0],"cos":[0,0]},{"arg":[-2,0,0,2,1],"sin":[129,0.1],"cos":[-70,0]},{"arg":[0,0,-1,2,2],"sin":[123,0],"cos":[-53,0]},{"arg":[2,0,0,0,0],"sin":[63,0],"cos":[0,0]},{"arg":[0,0,1,0,1],"sin":[63,0.1],"cos":[-33,0]},{"arg":[2,0,-1,2,2],"sin":[-59,0],"cos":[26,0]},{"arg":[0,0,-1,0,1],"sin":[-58,-0.1],"cos":[32,0]},{"arg":[0,0,1,2,1],"sin":[-51,0],"cos":[27,0]},{"arg":[-2,0,2,0,0],"sin":[48,0],"cos":[0,0]},{"arg":[0,0,-2,2,1],"sin":[46,0],"cos":[-24,0]},{"arg":[2,0,0,2,2],"sin":[-38,0],"cos":[16,0]},{"arg":[0,0,2,2,2],"sin":[-31,0],"cos":[13,0]},{"arg":[0,0,2,0,0],"sin":[29,0],"cos":[0,0]},{"arg":[-2,0,1,2,2],"sin":[29,0],"cos":[-12,0]},{"arg":[0,0,0,2,0],"sin":[26,0],"cos":[0,0]},{"arg":[-2,0,0,2,0],"sin":[-22,0],"cos":[0,0]},{"arg":[0,0,-1,2,1],"sin":[21,0],"cos":[-10,0]},{"arg":[0,2,0,0,0],"sin":[17,-0.1],"cos":[0,0]},{"arg":[2,0,-1,0,1],"sin":[16,0],"cos":[-8,0]},{"arg":[-2,2,0,2,2],"sin":[-16,0.1],"cos":[7,0]},{"arg":[0,1,0,0,1],"sin":[-15,0],"cos":[9,0]},{"arg":[-2,0,1,0,1],"sin":[-13,0],"cos":[7,0]},{"arg":[0,-1,0,0,1],"sin":[-12,0],"cos":[6,0]},{"arg":[0,0,2,-2,0],"sin":[11,0],"cos":[0,0]},{"arg":[2,0,-1,2,1],"sin":[-10,0],"cos":[5,0]},{"arg":[2,0,1,2,2],"sin":[-8,0],"cos":[3,0]},{"arg":[0,1,0,2,2],"sin":[7,0],"cos":[-3,0]},{"arg":[-2,1,1,0,0],"sin":[-7,0],"cos":[0,0]},{"arg":[0,-1,0,2,2],"sin":[-7,0],"cos":[3,0]},{"arg":[2,0,0,2,1],"sin":[-7,0],"cos":[3,0]},{"arg":[2,0,1,0,0],"sin":[6,0],"cos":[0,0]},{"arg":[-2,0,2,2,2],"sin":[6,0],"cos":[-3,0]},{"arg":[-2,0,1,2,1],"sin":[6,0],"cos":[-3,0]},{"arg":[2,0,-2,0,1],"sin":[-6,0],"cos":[3,0]},{"arg":[2,0,0,0,1],"sin":[-6,0],"cos":[3,0]},{"arg":[0,-1,1,0,0],"sin":[5,0],"cos":[0,0]},{"arg":[-2,-1,0,2,1],"sin":[-5,0],"cos":[3,0]},{"arg":[-2,0,0,0,1],"sin":[-5,0],"cos":[3,0]},{"arg":[0,0,2,2,1],"sin":[-5,0],"cos":[3,0]},{"arg":[-2,0,2,0,1],"sin":[4,0],"cos":[0,0]},{"arg":[-2,1,0,2,1],"sin":[4,0],"cos":[0,0]},{"arg":[0,0,1,-2,0],"sin":[4,0],"cos":[0,0]},{"arg":[-1,0,1,0,0],"sin":[-4,0],"cos":[0,0]},{"arg":[-2,1,0,0,0],"sin":[-4,0],"cos":[0,0]},{"arg":[1,0,0,0,0],"sin":[-4,0],"cos":[0,0]},{"arg":[0,0,1,2,0],"sin":[3,0],"cos":[0,0]},{"arg":[0,0,-2,2,2],"sin":[-3,0],"cos":[0,0]},{"arg":[-1,-1,1,0,0],"sin":[-3,0],"cos":[0,0]},{"arg":[0,1,1,0,0],"sin":[-3,0],"cos":[0,0]},{"arg":[0,-1,1,2,2],"sin":[-3,0],"cos":[0,0]},{"arg":[2,-1,-1,2,2],"sin":[-3,0],"cos":[0,0]},{"arg":[0,0,3,2,2],"sin":[-3,0],"cos":[0,0]},{"arg":[2,-1,0,2,2],"sin":[-3,0],"cos":[0,0]}]
    temp_arr = []
    for item in periodic_terms_copied:
        temp_arr.append([item['arg'], item['sin'], item['cos']])
    return temp_arr

    # periodic terms for longitude and distance of the moon
    # periodic_terms_lon_dist = [
    #     [
    # [0, 0, 1, 0], 6288774, -20905355],
    # [[2, 0, -1, 0], 1274027, -3699111],
    # [[2, 0, 0, 0], 658314, -2955968],
    # [[0, 0, 2, 0], 213618, -569925],
    # [[0, 1, 0, 0], -185116, 48888],
    # [[0, 0, 0, 2], -114332, -3149],
    # [[2, 0, -2, 0], 58793, 246158],
    # [[2, -1, -1, 0], 57066, -152138],
    # [[2, 0, 1, 0], 53322, -170733],
    # [[2, -1, 0, 0], 45758, -204586],
    # [[0, 1, -1, 0], -40923, -129620],
    # [[1, 0, 0, 0], -34720, 108743],
    # [[0, 1, 1, 0], -30383, 104755],
    # [[2, 0, 0, -2], 15327, 10321],
    # [[0, 0, 1, 2], -12528, 0],
    # [[0, 0, 1, -2], 10980, 79661],
    # [[4, 0, -1, 0], 10675, -34782],
    # [[0, 0, 3, 0], 10034, -23210],
    # [[4, 0, -2, 0], 8548, -21636],
    # [[2, 1, -1, 0], -7888, 24208],
    # [[2, 1, 0, 0], -6766, 30824],
    # [[1, 0, -1, 0], -5163, -8379],
    # [[1, 1, 0, 0], 4987, -16675],
    # [[2, -1, 1, 0], 4036, -12831],
    # [[2, 0, 2, 0], 3994, -10445],
    # [[4, 0, 0, 0], 3861, -11650],
    # [[2, 0, -3, 0], 3665, 14403],
    # [[0, 1, -2, 0], -2689, -7003],
    # [[2, 0, -1, 2], -2602, 0],
    # [[2, -1, -2, 0], 2390, 10056],
    # [[1, 0, 1, 0], -2348, 6322],
    # [[2, -2, 0, 0], 2236, -9884],
    #
    # [[0, 1, 2, 0], -2120, 5751],
    # [[0, 2, 0, 0], -2069, 0],
    # [[2, -2, -1, 0], 2048, -4950],
    # [[2, 0, 1, -2], -1773, 4130],
    # [[2, 0, 0, 2], -1595, 0],
    # [[4, -1, -1, 0], 1215, -3958],
    # [[0, 0, 2, 2], -1110, 0],
    # [[3, 0, -1, 0], -892, 3258],
    # [[2, 1, 1, 0], -810, 2616],
    # [[4, -1, -2, 0], 759, -1897],
    # [[0, 2, -1, 0], -713, -2117],
    # [[2, 2, -1, 0], -700, 2354],
    # [[2, 1, -2, 0], 691, 0],
    # [[2, -1, 0, -2], 596, 0],
    # [[4, 0, 1, 0], 549, -1423],
    # [[0, 0, 4, 0], 537, -1117],
    # [[4, -1, 0, 0], 520, -1571],
    # [[1, 0, -2, 0], -487, -1739],
    # [[2, 1, 0, -2], -399, 0],
    # [[0, 0, 2, -2], -381, -4421],
    # [[1, 1, 1, 0], 351, 0],
    # [[3, 0, -2, 0], -340, 0],
    # [[4, 0, -3, 0], 330, 0],
    # [[2, -1, 2, 0], 327, 0],
    # [[0, 2, 1, 0], -323, 1165],
    # [[1, 1, -1, 0], 299, 0],
    # [[2, 0, 3, 0], 294, 0],
    # [[2, 0, -1, -2], 0, 8752]
    # ]

def make_periodic_terms():
    periodic_terms_copied = [{"arg":[0,0,1,0],"sin":6288774,"cos":-20905355},{"arg":[2,0,-1,0],"sin":1274027,"cos":-3699111},{"arg":[2,0,0,0],"sin":658314,"cos":-2955968},{"arg":[0,0,2,0],"sin":213618,"cos":-569925},{"arg":[0,1,0,0],"sin":-185116,"cos":48888},{"arg":[0,0,0,2],"sin":-114332,"cos":-3149},{"arg":[2,0,-2,0],"sin":58793,"cos":246158},{"arg":[2,-1,-1,0],"sin":57066,"cos":-152138},{"arg":[2,0,1,0],"sin":53322,"cos":-170733},{"arg":[2,-1,0,0],"sin":45758,"cos":-204586},{"arg":[0,1,-1,0],"sin":-40923,"cos":-129620},{"arg":[1,0,0,0],"sin":-34720,"cos":108743},{"arg":[0,1,1,0],"sin":-30383,"cos":104755},{"arg":[2,0,0,-2],"sin":15327,"cos":10321},{"arg":[0,0,1,2],"sin":-12528,"cos":0},{"arg":[0,0,1,-2],"sin":10980,"cos":79661},{"arg":[4,0,-1,0],"sin":10675,"cos":-34782},{"arg":[0,0,3,0],"sin":10034,"cos":-23210},{"arg":[4,0,-2,0],"sin":8548,"cos":-21636},{"arg":[2,1,-1,0],"sin":-7888,"cos":24208},{"arg":[2,1,0,0],"sin":-6766,"cos":30824},{"arg":[1,0,-1,0],"sin":-5163,"cos":-8379},{"arg":[1,1,0,0],"sin":4987,"cos":-16675},{"arg":[2,-1,1,0],"sin":4036,"cos":-12831},{"arg":[2,0,2,0],"sin":3994,"cos":-10445},{"arg":[4,0,0,0],"sin":3861,"cos":-11650},{"arg":[2,0,-3,0],"sin":3665,"cos":14403},{"arg":[0,1,-2,0],"sin":-2689,"cos":-7003},{"arg":[2,0,-1,2],"sin":-2602,"cos":0},{"arg":[2,-1,-2,0],"sin":2390,"cos":10056},{"arg":[1,0,1,0],"sin":-2348,"cos":6322},{"arg":[2,-2,0,0],"sin":2236,"cos":-9884},{"arg":[0,1,2,0],"sin":-2120,"cos":5751},{"arg":[0,2,0,0],"sin":-2069,"cos":0},{"arg":[2,-2,-1,0],"sin":2048,"cos":-4950},{"arg":[2,0,1,-2],"sin":-1773,"cos":4130},{"arg":[2,0,0,2],"sin":-1595,"cos":0},{"arg":[4,-1,-1,0],"sin":1215,"cos":-3958},{"arg":[0,0,2,2],"sin":-1110,"cos":0},{"arg":[3,0,-1,0],"sin":-892,"cos":3258},{"arg":[2,1,1,0],"sin":-810,"cos":2616},{"arg":[4,-1,-2,0],"sin":759,"cos":-1897},{"arg":[0,2,-1,0],"sin":-713,"cos":-2117},{"arg":[2,2,-1,0],"sin":-700,"cos":2354},{"arg":[2,1,-2,0],"sin":691,"cos":0},{"arg":[2,-1,0,-2],"sin":596,"cos":0},{"arg":[4,0,1,0],"sin":549,"cos":-1423},{"arg":[0,0,4,0],"sin":537,"cos":-1117},{"arg":[4,-1,0,0],"sin":520,"cos":-1571},{"arg":[1,0,-2,0],"sin":-487,"cos":-1739},{"arg":[2,1,0,-2],"sin":-399,"cos":0},{"arg":[0,0,2,-2],"sin":-381,"cos":-4421},{"arg":[1,1,1,0],"sin":351,"cos":0},{"arg":[3,0,-2,0],"sin":-340,"cos":0},{"arg":[4,0,-3,0],"sin":330,"cos":0},{"arg":[2,-1,2,0],"sin":327,"cos":0},{"arg":[0,2,1,0],"sin":-323,"cos":1165},{"arg":[1,1,-1,0],"sin":299,"cos":0},{"arg":[2,0,3,0],"sin":294,"cos":0},{"arg":[2,0,-1,-2],"sin":0,"cos":8752}]
    temp_arr = []
    for item in periodic_terms_copied:
        temp_arr.append([item['arg'], item['sin'], item['cos']])

    periodic_terms_copied_2 = [{"arg":[0,0,0,1],"sin":5128122},{"arg":[0,0,1,1],"sin":280602},{"arg":[0,0,1,-1],"sin":277693},{"arg":[2,0,0,-1],"sin":173237},{"arg":[2,0,-1,1],"sin":55413},{"arg":[2,0,-1,-1],"sin":46271},{"arg":[2,0,0,1],"sin":32573},{"arg":[0,0,2,1],"sin":17198},{"arg":[2,0,1,-1],"sin":9266},{"arg":[0,0,2,-1],"sin":8822},{"arg":[2,-1,0,-1],"sin":8216},{"arg":[2,0,-2,-1],"sin":4324},{"arg":[2,0,1,1],"sin":4200},{"arg":[2,1,0,-1],"sin":-3359},{"arg":[2,-1,-1,1],"sin":2463},{"arg":[2,-1,0,1],"sin":2211},{"arg":[2,-1,-1,-1],"sin":2065},{"arg":[0,1,-1,-1],"sin":-1870},{"arg":[4,0,-1,-1],"sin":1828},{"arg":[0,1,0,1],"sin":-1794},{"arg":[0,0,0,3],"sin":-1749},{"arg":[0,1,-1,1],"sin":-1565},{"arg":[1,0,0,1],"sin":-1491},{"arg":[0,1,1,1],"sin":-1475},{"arg":[0,1,1,-1],"sin":-1410},{"arg":[0,1,0,-1],"sin":-1344},{"arg":[1,0,0,-1],"sin":-1335},{"arg":[0,0,3,1],"sin":1107},{"arg":[4,0,0,-1],"sin":1021},{"arg":[4,0,-1,1],"sin":833},{"arg":[0,0,1,-3],"sin":777},{"arg":[4,0,-2,1],"sin":671},{"arg":[2,0,0,-3],"sin":607},{"arg":[2,0,2,-1],"sin":596},{"arg":[2,-1,1,-1],"sin":491},{"arg":[2,0,-2,1],"sin":-451},{"arg":[0,0,3,-1],"sin":439},{"arg":[2,0,2,1],"sin":422},{"arg":[2,0,-3,-1],"sin":421},{"arg":[2,1,-1,1],"sin":-366},{"arg":[2,1,0,1],"sin":-351},{"arg":[4,0,0,1],"sin":331},{"arg":[2,-1,1,1],"sin":315},{"arg":[2,-2,0,-1],"sin":302},{"arg":[0,0,1,3],"sin":-283},{"arg":[2,1,1,-1],"sin":-229},{"arg":[1,1,0,-1],"sin":223},{"arg":[1,1,0,1],"sin":223},{"arg":[0,1,-2,-1],"sin":-220},{"arg":[2,1,-1,-1],"sin":-220},{"arg":[1,0,1,1],"sin":-185},{"arg":[2,-1,-2,-1],"sin":181},{"arg":[0,1,2,1],"sin":-177},{"arg":[4,0,-2,-1],"sin":176},{"arg":[4,-1,-1,-1],"sin":166},{"arg":[1,0,1,-1],"sin":-164},{"arg":[4,0,1,-1],"sin":132},{"arg":[1,0,-1,-1],"sin":-119},{"arg":[4,-1,0,-1],"sin":115},{"arg":[2,-2,0,1],"sin":107}]
    temp_arr_2 = []
    for item in periodic_terms_copied_2:
        temp_arr_2.append([item['arg'], item['sin']])

    periodic_terms_lon_dist = temp_arr
    periodic_terms_lat = temp_arr_2
    return periodic_terms_lon_dist, periodic_terms_lat

def process_decimal(number):
    if 'e' in number:
        temp_res = number.split('e')
        return float(temp_res[0]) * 10 ** (int(temp_res[1]))
    else:
        return float(number)

def get_axial_element(barycentric_obj):
    string_data = str(barycentric_obj)
    res = string_data.split(',')

    res[0] = process_decimal(res[0][1:])
    res[1] = process_decimal(res[1][1:])
    res[2] = process_decimal(res[2][1:-4])
    print(res)
    return res





# # solar_system_ephemeris.set('jpl')
# solar_system_ephemeris.set('de432s')
#
# print(solar_system_ephemeris.bodies)
#
# # moon_pos = get_body_barycentric('moon', astropy.time.Time.now())
# moon_pos = get_body_barycentric('moon', astropy.time.Time.now())
# sun_pos = get_body_barycentric('sun', astropy.time.Time.now())
# earth_pos = get_body_barycentric('earth', astropy.time.Time.now())
#
def haversine(lat_1, lon_1, lat_2, lon_2):
    # haversine formula로 두 점 사이의 중심각을 구함
    # 원래는 두 점 사이의 거리를 구하는 식인데, 잘라내서 응용함
    # https://velog.io/@revimal/Mathematics-Haversine-Formula
    def haversine(rad):
        versine = 1 - math.cos(rad)
        return versine/2

    def archav(val):
        return 2*math.asin(math.sqrt(val))

    val = haversine(math.radians(lat_2 - lat_1)) + math.cos(math.radians(lat_1)) * math.cos(math.radians(lat_2)) * haversine(math.radians(lon_2 - lon_1))
    theta = math.degrees(archav(val))
    return theta

def get_azimuth_latlon(lat_1, lon_1, lat_2, lon_2):
    # https://spiralmoon.tistory.com/entry/Algorithm-%EC%A7%80%EA%B5%AC%EC%97%90%EC%84%9C-%EB%91%90-%EC%A0%90-%EC%82%AC%EC%9D%B4%EC%9D%98-%EB%B0%A9%EC%9C%84%EA%B0%81-%EA%B5%AC%ED%95%98%EA%B8%B0
    # y = math.sin(math.radians(lon_2 - lon_1)) * math.cos(math.radians(lat_2))
    # x = math.cos(math.radians(lat_1)) * math.sin(math.radians(lat_2)) - math.sin(math.radians(lat_1)) * math.cos(math.radians(lat_2)) * math.cos(math.radians(lon_2 - lon_1))
    # theta = math.atan2(y, x)
    # bearing = math.degrees(theta) % 360
    # return bearing
    # https://m.blog.naver.com/PostView.naver?isHttpsRedirect=true&blogId=ezpy_&logNo=220129817873
    from math import pow, degrees, radians, atan2
    from scipy import cos, sin, arctan, sqrt, arctan2
    Lat1, Lat2 = radians(lat_1), radians(lat_2)
    Long1, Long2 = radians(lon_1), radians(lon_2)
    y = sin(Long2 - Long1) * cos(Lat2)
    x = cos(Lat1) * sin(Lat2) - sin(Lat1) * cos(Lat2) * cos(Long2 - Long1)
    return degrees(atan2(y, x))

def calc_jd(year, month, day):
    # Chapter 7 참조
    # import calendar
    # day_of_month = calendar.monthrange(year, month)[1]

    if month <= 2:
        year = year - 1
        month = month + 12
    a = int(year/100)

    # Gregorian Calendar의 개정 시행 이전 날자일 경우 (Day in julian calendar), b를 0으로 둔다
    # https://namu.wiki/w/%EA%B7%B8%EB%A0%88%EA%B3%A0%EB%A6%AC%EB%A0%A5
    if year < 1582:
        b = 0
    elif year == 1582 and month < 10:
        b = 0
    elif year == 1582 and month == 10 and day < 5:
        b = 0
    else:
        b = 2 - a + int(a/4)

    jd = int(365.25*(year + 4716)) + int(30.6001*(month + 1)) + day + b - 1524.53
    return round(jd, 1)

import decimal

def calc_t(jde):
    # Chapter 47에 사용될 century단위의 JDE를 계산
    return round((jde-2451545)/36525, 12)

def rev_calc_jde(t):
    # t로 jde 를 역산
    return t*36525 + 2451545

def calc_rst(D, M, Mp, F, E, mean_ascending_node, t):
    # D : Mean elongation of the moon
    # M : Sun's mean anomaly
    # M` : Moon's mean anomaly (Mp)
    # F : Moon's argument of latitude

    K_1 = 119.75 + 131.849 * t
    K_2 = 72.56 + 20.186 * t

    rho = (-0.02752 * math.cos(math.radians(Mp))
           -0.02245 * math.sin(math.radians(F))
           + 0.00684 * math.cos(math.radians(Mp-2*F))
           - 0.00293 * math.cos(math.radians(2*F))
           - 0.00085 * math.cos(math.radians(2*F-2*D))
           - 0.00054 * math.cos(math.radians(Mp-2*D))
           - 0.00020 * math.sin(math.radians(Mp+F))
           - 0.00020 * math.cos(math.radians(Mp+2*F))
           - 0.00020 * math.cos(math.radians(Mp-F))
           + 0.00014 * math.cos(math.radians(Mp+2*F-2*D))
           )

    sigma = (-0.02816 * math.sin(math.radians(Mp))
           + 0.02244 * math.cos(math.radians(F))
           - 0.00682 * math.sin(math.radians(Mp-2*F))
           - 0.00279 * math.sin(math.radians(2*F))
           - 0.00083 * math.sin(math.radians(2*F-2*D))
           + 0.00069 * math.sin(math.radians(Mp-2*D))
           + 0.00040 * math.cos(math.radians(Mp+F))
           - 0.00025 * math.sin(math.radians(2*Mp))
           - 0.00023 * math.sin(math.radians(Mp+2*F))
           + 0.00020 * math.cos(math.radians(Mp-F))
           + 0.00019 * math.sin(math.radians(Mp-F))
           + 0.00013 * math.sin(math.radians(Mp + 2 * F - 2 * D))
           - 0.00010 * math.cos(math.radians(Mp - 3 * F))
           )

    tau = (0.02520 * E * math.sin(math.radians(M))
             + 0.00473 * math.sin(math.radians(2*Mp - 2*F))
             - 0.00467 * math.sin(math.radians(Mp))
             + 0.00396 * math.sin(math.radians(K_1))
             + 0.00276 * math.sin(math.radians(2*Mp - 2*D))
             + 0.00196 * math.sin(math.radians(mean_ascending_node))
             - 0.00183 * math.cos(math.radians(Mp - F))
             + 0.00115 * math.sin(math.radians(Mp - 2*D))
             - 0.00096 * math.sin(math.radians(Mp - D))
             + 0.00046 * math.sin(math.radians(2*F - 2*D))
             - 0.00039 * math.sin(math.radians(Mp - F))
             - 0.00032 * math.sin(math.radians(Mp - M - D))
             + 0.00027 * math.sin(math.radians(2*Mp - M - 2*D))
             + 0.00023 * math.sin(math.radians(K_2))
             - 0.00014 * math.sin(math.radians(2*D))
             + 0.00014 * math.cos(math.radians(2*Mp - 2*F))
             - 0.00012 * math.sin(math.radians(Mp - 2*F))
             - 0.00012 * math.sin(math.radians(2*Mp))
             + 0.00011 * math.sin(math.radians(2*Mp - 2*M - 2*D))
             )

    return rho, sigma, tau

def calc_nutation(t):
    periodic_terms_nutation = make_periodic_terms_nutation()

    # D : Mean elongation of the moon from the sun
    # M : Sun's mean anomaly
    # M` : Moon's mean anomaly
    # F : Moon's argument of latitude

    lunar_mean_elongation_nutation = 297.85036 + 445267.111480*t - 0.0019142*t**2 + t**3 / 189474    # moon_from_sun
    solar_mean_anomaly_nutation = 357.52772 + 35999.050340 * t - 0.0001603 * t**2 - t**3/300000 # on_earth
    lunar_mean_anomaly_nutation = 134.96298 + 477198.867398*t + 0.0086972 * t**2 + t**3 / 56250
    lunar_argument_of_latitude_nutation = 93.27191 + 483202.017538 * t - 0.0036825 * t**2 + t**3 / 327270

    # longitude of the ascending node of the Moon's mean orbit on the ecliptic (measured from the mean equinox of the date)
    lunar_mean_orbit_longitude = 125.04452 - 1934.136261 * t + 0.0020708 * t**2 + t**3 / 450000 # sigma

    lunar_mean_elongation_nutation = lunar_mean_elongation_nutation % 360
    solar_mean_anomaly_nutation = solar_mean_anomaly_nutation % 360
    lunar_mean_anomaly_nutation = lunar_mean_anomaly_nutation % 360
    lunar_argument_of_latitude_nutation = lunar_argument_of_latitude_nutation % 360
    lunar_mean_orbit_longitude = lunar_mean_orbit_longitude % 360

    nutation_longitude_arr = []
    nutation_obliquity_arr = []
    for item in periodic_terms_nutation:
        arg = item[0][0] * lunar_mean_elongation_nutation + item[0][1] * solar_mean_anomaly_nutation + item[0][2] * lunar_mean_anomaly_nutation + item[0][3] *lunar_argument_of_latitude_nutation + item[0][4] * lunar_mean_orbit_longitude
        nutation_longitude_arr.append((item[1][0] + item[1][1]) * math.sin(math.radians(arg)))
        nutation_obliquity_arr.append((item[2][0] + item[2][1]) * math.cos(math.radians(arg)))

    # print(sum(nutation_longitude_arr)/10000/3600)
    nutation_longitude = sum(nutation_longitude_arr)/10000/3600  # unit : deg
    nutation_obliquity = sum(nutation_obliquity_arr)/10000/3600  # unit : deg

    # L = 280.4665 + 36000.7698 * t
    # L_prime = 218.3165 + 481267.8813 * t
    #
    # nutation_longitude = -17.20/3600 * math.sin(math.radians(lunar_mean_orbit_longitude) + 1.32/3600 * math.sin(math.radians(2*L)) - 0.23/3600 * math.sin(math.radians(2*L_prime)) + 0.21/3600 * math.sin(math.radians(2*lunar_mean_orbit_longitude)) )

    # mean obliquity of the ecliptic
    # mean_obliquity = 23 + 26/60 + 21.448/3600 - 46.8150/3600 * t - 0.00059/3600 * t**2 + 0.001813/3600 * t**3 # old formula
    u = t/100
    mean_obliquity = 23 + 26/60 + 21.448/3600 - 4680.93/3600 * u - 1.55 * u**2 + 1999.25 * u**3 - 51.38 * u**4 - 249.67 * u**5 - 39.05 * u**6 + 7.12 * u**7 + 27.87 * u**8 + 5.79 * u**9 + 2.45 * u**10

    true_obliquity = mean_obliquity + nutation_obliquity

    attr_dict = {
        't': t,
        'lunar_mean_elongation_nutation': lunar_mean_elongation_nutation,
        'solar_mean_anomaly_nutation': solar_mean_anomaly_nutation,
        'lunar_mean_anomaly_nutation': lunar_mean_anomaly_nutation,
        'lunar_argument_of_latitude_nutation': lunar_argument_of_latitude_nutation,
        'lunar_mean_orbit_longitude': lunar_mean_orbit_longitude,
        'nutation_longitude': nutation_longitude,
        'nutation_obliquity': nutation_obliquity,
        'mean_obliquity': mean_obliquity,
        'true_obliquity': true_obliquity
    }
    print_dict(attr_dict, 'nutation')
    return nutation_longitude, nutation_obliquity, mean_obliquity, true_obliquity

def calc_lunar_attr(t):
    periodic_terms_lon_dist, periodic_terms_lat = make_periodic_terms()

    # L` : Mean longitude of the moon
    # D : Mean elongation of the moon
    # M : Sun's mean anomaly
    # M` : Moon's mean anomaly
    # F : Moon's argument of latitude

    lunar_mean_longitude = 218.3164477 + 481267.88123421*t - 0.0015786 * t**2 + t**3 / 538841 - t**4 / 65194000
    lunar_mean_elongation = 297.8501921 + 445267.1114034*t - 0.0018819*t**2 + t**3 / 545868 -  t**4 / 113065000
    solar_mean_anomaly = 357.5291092 + 35999.0502909*t - 0.0001536*t**2 + t**3 / 24490000
    lunar_mean_anomaly = 134.9633964 + 477198.8675055*t + 0.0087414*t**2 + t**3 / 69699 - t**4 / 14712000
    lunar_argument_of_latitude = 93.2720950 + 483202.0175233*t - 0.0036539*t**2 - t**3 / 3526000 + t**4 / 863310000

    lunar_mean_longitude = round(lunar_mean_longitude % 360, 6)
    lunar_mean_elongation = round(lunar_mean_elongation % 360, 6)
    solar_mean_anomaly = round(solar_mean_anomaly % 360, 6)
    lunar_mean_anomaly = round(lunar_mean_anomaly % 360, 6)
    lunar_argument_of_latitude = round(lunar_argument_of_latitude % 360, 6)


    # a = 6288774*math.sin(lunar_mean_anomaly) + 3958*math.sin(math.radians(119.75 + 131.849*t)) + 1962*math.sin(lunar_mean_longitude-lunar_argument_of_latitude) + 318*math.sin(math.radians(53.09 + 479264.290*t))
    # b = -20905355*math.cos(lunar_mean_anomaly)

    a_1 = round((119.75 + 131.849*t) % 360, 2)
    a_2 = round((53.09 + 479264.290*t) % 360, 2)
    a_3 = round((313.45 + 481266.484*t) % 360, 2)
    e = round(1 - 0.002516*t - 0.0000074*t**2, 6)


    # t = -0.077221081451
    #
    # lunar_mean_longitude = 134.290182
    # lunar_mean_elongation = 113.842304
    # solar_mean_anomaly = 97.643514
    # lunar_mean_anomaly = 5.150833
    # lunar_argument_of_latitude = 219.889721
    #
    # a_1 = 109.57
    # a_2 = 123.78
    # a_3 = 229.53
    # e = 1.000194



    longitude_arr = []
    distance_arr = []
    latitude_arr = []

    for item in periodic_terms_lon_dist:
        if item[0][1] == 0:
            longitude_arr.append(item[1] * math.sin(math.radians(lunar_mean_elongation*item[0][0] + solar_mean_anomaly*item[0][1] + lunar_mean_anomaly*item[0][2] + lunar_argument_of_latitude*item[0][3])))
            distance_arr.append(item[2] * math.cos(math.radians(lunar_mean_elongation * item[0][0] + solar_mean_anomaly * item[0][1] + lunar_mean_anomaly * item[0][2] + lunar_argument_of_latitude * item[0][3])))
        else:
            longitude_arr.append(item[1] * e**abs(item[0][1]) * math.sin(math.radians(lunar_mean_elongation * item[0][0] + solar_mean_anomaly * item[0][1] + lunar_mean_anomaly * item[0][2] + lunar_argument_of_latitude * item[0][3])))
            distance_arr.append(item[2] * e**abs(item[0][1]) * math.cos(math.radians(lunar_mean_elongation * item[0][0] + solar_mean_anomaly * item[0][1]+ lunar_mean_anomaly * item[0][2] + lunar_argument_of_latitude * item[0][3])))

    for item in periodic_terms_lat:
        if item[0][1] == 0:
            latitude_arr.append(item[1] * math.sin(math.radians(lunar_mean_elongation * item[0][0] + solar_mean_anomaly * item[0][1] + lunar_mean_anomaly * item[0][2] + lunar_argument_of_latitude * item[0][3])))
        else:
            latitude_arr.append(item[1] * e ** abs(item[0][1]) * math.sin(math.radians(lunar_mean_elongation * item[0][0] + solar_mean_anomaly * item[0][1] + lunar_mean_anomaly * item[0][2] + lunar_argument_of_latitude * item[0][3])))

    longitude_arr.append(3958*math.sin(math.radians(a_1)))
    longitude_arr.append(1962*math.sin(math.radians(lunar_mean_longitude - lunar_argument_of_latitude)))
    longitude_arr.append(318*math.sin(math.radians(a_2)))

    latitude_arr.append(-2235*math.sin(math.radians(lunar_mean_longitude)))
    latitude_arr.append(382*math.sin(math.radians(a_3)))
    latitude_arr.append(175*math.sin(math.radians(a_1 - lunar_argument_of_latitude)))
    latitude_arr.append(175*math.sin(math.radians(a_1 + lunar_argument_of_latitude)))
    latitude_arr.append(127*math.sin(math.radians(lunar_mean_longitude - lunar_mean_anomaly)))
    latitude_arr.append(-115*math.sin(math.radians(lunar_mean_longitude + lunar_mean_anomaly)))

    # lambda : geocentric longitude
    # beta : geocentric latitude
    # delta : center distance between earth and moon
    # pi : equatorial horizontal parallax

    attr_lambda = lunar_mean_longitude + sum(longitude_arr) / 1000000   # unit : deg
    attr_beta = sum(latitude_arr) / 1000000 # unit : deg
    attr_delta = 385000.56 + round(sum(distance_arr), 0) / 1000   # unit : kilometers
    attr_pi = math.degrees(math.asin(6378.14/attr_delta))

    attr_dict = {
        't':t,
        'lunar_mean_longitude':lunar_mean_longitude,
        'lunar_mean_elongation':lunar_mean_elongation,
        'solar_mean_anomaly':solar_mean_anomaly,
        'lunar_mean_anomaly':lunar_mean_anomaly,
        'lunar_argument_of_latitude':lunar_argument_of_latitude,
        'a_1':a_1,
        'a_2':a_2,
        'a_3':a_3,
        'E':e,
        'sum_longitude':round(sum(longitude_arr), 0),
        'sum_latitude':round(sum(latitude_arr), 0),
        'sum_distance':round(sum(distance_arr), 0),
        'lambda':round(attr_lambda, 6),
        'beta':round(attr_beta, 6),
        'delta':round(attr_delta, 6),
        'pi':round(attr_pi, 6)
            }

    # get nutation
    nutation_longitude, nutation_obliquity, mean_obliquity, true_obliquity = calc_nutation(t)

    # apparent_lambda : lambda + nutation in longitude (Chapter #22)
    # apparent_lambda = attr_lambda + nutation_longitude
    apparent_lambda = attr_lambda + 0.004610

    # apparent right ascension and declination
    apparent_ascension = math.degrees(math.atan((math.sin(math.radians(apparent_lambda)) * math.cos(math.radians(true_obliquity)) - math.tan(math.radians(attr_beta)) * math.sin(math.radians(true_obliquity)))/math.cos(math.radians(apparent_lambda))))
    apparent_declination = math.degrees(math.asin(math.sin(math.radians(attr_beta))*math.cos(math.radians(true_obliquity)) + math.cos(math.radians(attr_beta)) * math.sin(math.radians(true_obliquity)) * math.sin(math.radians(apparent_lambda))))
    # 탄젠트의 주기가 짧아서, 양수로 강제전환 시켜 줌
    if apparent_ascension < 0:
        apparent_ascension = 180 + apparent_ascension

    # 스피드업 비활성화
    # print(apparent_ascension, apparent_declination)

    # Lunar node and lunar perigee
    lunar_mean_ascending_node = 125.0445479 - 1934.1362891*t + 0.0020754 * t**2 + t**3 / 467441 - t**4 / 60616000
    lunar_mean_perigee = 83.3532465 + 4069.0137287 * t - 0.0103200 * t**2 - t**3 / 80053 + t**4 / 18999000

    lunar_true_ascending_node = lunar_mean_ascending_node - 1.4979*math.sin(2*math.radians(lunar_mean_elongation-lunar_argument_of_latitude)) - 0.1500 * math.sin(math.radians(solar_mean_anomaly)) - 0.1226 * math.sin(2*math.radians(lunar_mean_elongation)) - 0.1176 * math.sin(2*math.radians(lunar_argument_of_latitude)) - 0.0801*math.sin(2*math.radians(lunar_mean_anomaly-lunar_argument_of_latitude))

    attr_dict['apparent_lambda'] = apparent_lambda
    attr_dict['apparent_ascension'] = apparent_ascension
    attr_dict['apparent_declination'] = apparent_declination
    attr_dict['lunar_mean_ascending_node'] = lunar_mean_ascending_node
    attr_dict['lunar_mean_perigee'] = lunar_mean_perigee
    attr_dict['lunar_true_ascending_node'] = lunar_true_ascending_node
    print_dict(attr_dict, 'moon_position')

    return attr_dict


def calc_solar_attr(t):

    solar_mean_longitude = 280.46646 + 36000.76983 * t + 0.0003032 * t**2
    solar_mean_anomaly = 357.52911 + 35999.05029 * t - 0.0001537 * t**2
    earth_orbit_eccentricity = 0.016708634 - 0.000042037 * t - 0.0000001267 * t**2

    solar_mean_longitude = solar_mean_longitude % 360
    solar_mean_anomaly = solar_mean_anomaly % 360
    e = earth_orbit_eccentricity


    solar_center = (1.914602 - 0.004817 * t - 0.000014 * t**2) * math.sin(math.radians(solar_mean_anomaly)) + (0.019993 - 0.000101 * t) * math.sin(2 * math.radians(solar_mean_anomaly)) + 0.000289 * math.sin(3*math.radians(solar_mean_anomaly))

    solar_true_longitude = solar_mean_longitude + solar_center
    solar_true_anomaly = solar_mean_anomaly + solar_center

    solar_radius_vector = (1.000001018*(1-e**2)) / (1 + e * math.cos(math.radians(solar_true_anomaly)))

    # can be upgraded for high accuracy with the methodology in Chapter.26 (지금은 안 함!)
    solar_mean_ascending_node = 125.04 - 1934.136 * t
    solar_apparent_longitude = solar_true_longitude - 0.00569 - 0.00478 * math.sin(math.radians(solar_mean_ascending_node))

    # get nutation
    nutation_longitude, nutation_obliquity, mean_obliquity, true_obliquity = calc_nutation(t)

    solar_right_ascension = math.degrees(math.atan((math.cos(math.radians(true_obliquity)) * math.sin(math.radians(solar_true_longitude)))/math.cos(math.radians(solar_true_longitude))))
    solar_right_declination = math.degrees(math.asin(math.sin(math.radians(true_obliquity)) * math.sin(math.radians(solar_true_longitude))))
    solar_apparent_ascension = math.degrees(math.atan((math.cos(math.radians(true_obliquity + 0.00256*math.cos(math.radians(solar_mean_ascending_node)))) * math.sin(math.radians(solar_apparent_longitude)))/math.cos(math.radians(solar_apparent_longitude))))
    solar_apparent_declination = math.degrees(math.asin(math.sin(math.radians(true_obliquity + 0.00256*math.cos(math.radians(solar_mean_ascending_node)))) * math.sin(math.radians(solar_apparent_longitude))))

    # right_ascension과 apparent_ascension은 아크탄젠트여서 보정해 주어야 함. (지금은 귀찮으니까 빼자...)

    attr_dict = {
        't':t,
        'solar_mean_longitude':solar_mean_longitude,
        'solar_mean_anomaly':solar_mean_anomaly,
        'earth_orbit_eccentricity(e)':e,
        'solar_center':solar_center,
        'solar_true_longitude':solar_true_longitude,
        'solar_true_anomaly':solar_true_anomaly,
        'solar_radius_vector':solar_radius_vector,
        'solar_mean_ascending_node':solar_mean_ascending_node,
        'solar_apparent_longitude':solar_apparent_longitude,
        'true_obliquity':true_obliquity,
        'solar_right_ascension':solar_right_ascension,
        'solar_right_declination':solar_right_declination,
        'solar_apparent_ascension':solar_apparent_ascension,
        'solar_apparent_declination':solar_apparent_declination
            }
    print_dict(attr_dict, 'sun_position')

    return attr_dict

def calc_selenographic_solar_pos(t, lat, lon):
    lunar_attr = calc_lunar_attr(t)
    solar_attr = calc_solar_attr(t)
    lambda_zero = solar_attr['solar_apparent_longitude']   # 이것이 geocentric인지 확인이 필요함
    lambda_nosuffix = lunar_attr['apparent_lambda']
    beta_nosuffix = lunar_attr['beta']
    solar_radius_vector = solar_attr['solar_radius_vector'] * u.AU.to(u.kilometer)

    heliocentric_lambda = lambda_zero + 180 + lunar_attr['delta']/solar_radius_vector * 57.296 * math.cos(math.radians(beta_nosuffix)) * math.sin(math.radians(lambda_zero - lambda_nosuffix))
    heliocentric_beta = beta_nosuffix * lunar_attr['delta']/solar_radius_vector

    # 스피드업 비활성화
    # print(lambda_zero)
    # print(lunar_attr['delta'])
    # print(solar_radius_vector)

    heliocentric_lambda = round(heliocentric_lambda % 360, 6)
    heliocentric_beta = round(heliocentric_beta % 360, 6)
    # print(heliocentric_lambda, heliocentric_beta)

    # Proces eq.53.1 with substitution (lambda -> lambda_Heliocentric, beta -> beta_Heliocentric)
    # get nutation
    nutation_longitude, nutation_obliquity, mean_obliquity, true_obliquity = calc_nutation(t)
    # get lunar inclination (value adopted by the International Astronomical Union)
    lunar_inclination = 1.54242 # I

    W = heliocentric_lambda - nutation_longitude - lunar_attr['lunar_mean_ascending_node']
    A = math.degrees(math.atan(
        (math.sin(math.radians(W)) * math.cos(math.radians(heliocentric_beta)) * math.cos(math.radians(lunar_inclination)) - math.sin(math.radians(heliocentric_beta)) * math.sin(math.radians(lunar_inclination))) / (math.cos(math.radians(W)) * math.cos(math.radians(heliocentric_beta)))
    ))
    W = norm_deg(W)
    tan_a = (math.sin(math.radians(W)) * math.cos(math.radians(heliocentric_beta)) * math.cos(math.radians(lunar_inclination)) - math.sin(math.radians(heliocentric_beta)) * math.sin(math.radians(lunar_inclination))) / (math.cos(math.radians(W)) * math.cos(math.radians(heliocentric_beta)))

    test_dict = {}
    test_dict['A'] = A
    test_dict['W'] = W
    test_dict['beta'] = heliocentric_beta
    test_dict['I'] = lunar_inclination
    test_dict['A_1'] = math.sin(math.radians(W))
    test_dict['A_2'] = math.cos(math.radians(heliocentric_beta))
    test_dict['A_3'] = math.cos(math.radians(lunar_inclination))
    test_dict['A_4'] = math.sin(math.radians(lunar_inclination))
    test_dict['A_5'] = math.cos(math.radians(W))
    test_dict['A_6'] = math.cos(math.radians(heliocentric_beta))
    # 아크탄젠트 처리 코드
    if math.cos(math.radians(W)) < 0:
        inverse = True
    else:
        inverse = False

    atan_a = A

    A = norm_deg(A)

    # [애매] 일단 A를 3사분면에 위치시킨다
    # if A < 180:
    #     A = A + 180


    l_prime_zero = A - lunar_attr['lunar_argument_of_latitude']
    b_prime_zero = math.degrees(math.asin(
        (-math.sin(math.radians(W)) * math.cos(math.radians(heliocentric_beta)) * math.sin(math.radians(lunar_inclination))) - (math.sin(math.radians(heliocentric_beta)) * math.cos(math.radians(lunar_inclination)))
    ))
    # print(b_prime_zero)

    # calculating physical librations

    # D : Mean elongation of the moon
    # M : Sun's mean anomaly
    # M` : Moon's mean anomaly
    # F : Moon's argument of latitude

    lunar_mean_elongation = lunar_attr['lunar_mean_elongation']
    solar_mean_anomaly = lunar_attr['solar_mean_anomaly']
    lunar_mean_anomaly = lunar_attr['lunar_mean_anomaly']
    lunar_argument_of_latitude= lunar_attr['lunar_argument_of_latitude']
    E = lunar_attr['E']

    rho, sigma, tau = calc_rst(lunar_mean_elongation, solar_mean_anomaly, lunar_mean_anomaly, lunar_argument_of_latitude, E, lunar_attr['lunar_mean_ascending_node'], t)

    # 달 중심으로 구하는 것이 아닌 부분, 테스트용으로 임시 사용한 뒤 주석 처리하였음 (2023.01.10)
    # 해당 부분은 Chapter.53에 해당됨.
    # W = lunar_attr['apparent_lambda'] - nutation_longitude - lunar_attr['lunar_mean_ascending_node']
    # A = math.degrees(math.atan(
    #     (math.sin(math.radians(W)) * math.cos(math.radians(lunar_attr['beta'])) * math.cos(
    #         math.radians(lunar_inclination)) - math.sin(math.radians(lunar_attr['beta'])) * math.sin(
    #         math.radians(lunar_inclination))) / (math.cos(math.radians(W)) * math.cos(math.radians(lunar_attr['beta'])))
    # ))
    # print(A)
    # l_prime = A - lunar_attr['lunar_argument_of_latitude']
    # b_prime = math.degrees(math.asin(
    #     (-math.sin(math.radians(W)) * math.cos(math.radians(lunar_attr['beta'])) * math.sin(
    #         math.radians(lunar_inclination))) - (
    #                 math.sin(math.radians(lunar_attr['beta'])) * math.cos(math.radians(lunar_inclination)))
    # ))
    # # print(b_prime)
    #
    #
    # print(rho, sigma, tau)
    # A = norm_deg(A+180)
    # W = norm_deg(W)
    # l_doubleprime = -tau + (rho*math.cos(math.radians(A)) + sigma * math.sin(math.radians(A))) * math.tan(math.radians(b_prime))
    # b_doubleprime = sigma * math.cos(math.radians(A)) - rho * math.sin(math.radians(A))
    #
    # l = l_prime + l_doubleprime
    # b = b_prime + b_doubleprime
    # print(lunar_attr['lunar_mean_ascending_node'])
    # print(lunar_attr['apparent_lambda'], lunar_attr['beta'])
    # print(W, A)
    # print(l_doubleprime, b_doubleprime)
    # print(l_prime, b_prime)
    # print(l, b)
    # raise IOError

    # Proces eq.53.2 with substitution (b_prime -> b_prime_zero)
    l_doubleprime_zero = -tau + (rho*math.cos(math.radians(A)) + sigma * math.sin(math.radians(A))) * math.tan(math.radians(b_prime_zero))
    b_doubleprime_zero = sigma * math.cos(math.radians(A)) - rho * math.sin(math.radians(A))

    l_zero = l_prime_zero + l_doubleprime_zero
    b_zero = b_prime_zero + b_doubleprime_zero


    selenographic_colongitude_sun = (90-l_zero) % 360
    if selenographic_colongitude_sun < 0:
        selenographic_colongitude_sun = selenographic_colongitude_sun + 360
    # print(selenographic_colongitude_sun)
    # selenographic_colongitude_sun = 360 - selenographic_colongitude_sun

    # l_zero = norm_deg(l_zero)
    # b_zero = norm_deg(b_zero)
    # l_prime_zero = norm_deg(l_prime_zero)
    # b_prime_zero = norm_deg(b_prime_zero)
    # l_doubleprime_zero = norm_deg(l_doubleprime_zero)
    # b_doubleprime_zero = norm_deg(b_doubleprime_zero)

    # print(W, A)
    # print(selenographic_colongitude_sun)
    # print(l_zero, b_zero)
    # print(l_prime_zero, b_prime_zero)
    # print(l_doubleprime_zero, b_doubleprime_zero)

    selenographic_lon = lon
    selenographic_lat = lat

    # print(selenographic_lon, selenographic_lat)

    # test_dict = {}

    # calculate altitude h
    h = math.degrees(math.asin(
        math.sin(math.radians(b_zero))*math.sin(math.radians(selenographic_lat)) + math.cos(math.radians(b_zero)) * math.cos(math.radians(selenographic_lat)) * math.sin(math.radians(selenographic_colongitude_sun + selenographic_lon))
    ))


    # 값반전 테스트
    test_dict['raw_value'] = math.sin(math.radians(b_zero))*math.sin(math.radians(selenographic_lat)) + math.cos(math.radians(b_zero)) * math.cos(math.radians(selenographic_lat)) * math.sin(math.radians(selenographic_colongitude_sun + selenographic_lon))
    test_dict['asin'] = math.asin(
        math.sin(math.radians(b_zero))*math.sin(math.radians(selenographic_lat)) + math.cos(math.radians(b_zero)) * math.cos(math.radians(selenographic_lat)) * math.sin(math.radians(selenographic_colongitude_sun + selenographic_lon))
    )
    test_dict['a'] = math.sin(math.radians(b_zero)) * math.sin(math.radians(selenographic_lat))
    test_dict['b'] = math.cos(math.radians(b_zero)) * math.cos(math.radians(selenographic_lat)) * math.sin(math.radians(selenographic_colongitude_sun + selenographic_lon))
    test_dict['b1'] = math.cos(math.radians(b_zero))
    test_dict['b2'] = math.cos(math.radians(selenographic_lat))
    test_dict['b3'] = math.sin(math.radians(selenographic_colongitude_sun + selenographic_lon))
    test_dict['sel_colon_sun'] = selenographic_colongitude_sun
    test_dict['sel_lon'] = selenographic_lon
    test_dict['lpz'] = l_prime_zero
    test_dict['ldpz'] = l_doubleprime_zero
    test_dict['tan_a'] = tan_a
    test_dict['atan_a'] = atan_a


    # -90~-180, 90~180의 아크사인 값 조절
    # if h > 90:
    #     h = h - 90
    # elif h < -90:
    #     h = h + 90
    # print(h)

    # 태양-달의 거리를 구함
    # 1. haversine formula로 사이각을 구한다
    geocentric_angle_between_sun_moon = haversine(0, solar_attr['solar_true_longitude'], lunar_attr['beta'], lunar_attr['lambda'])
    # 2. 제 2 코사인법칙으로 거리를 구한다
    distance_sun_moon = math.sqrt(solar_radius_vector**2 + lunar_attr['delta']**2 - 2*solar_radius_vector * lunar_attr['delta'] * math.cos(math.radians(geocentric_angle_between_sun_moon)))

    # subsolar point는 l_zero와 b_zero임
    # haversine formula로 사이각을 구한다
    lunar_centerangle = haversine(b_zero, l_zero, selenographic_lat, selenographic_lon)
    # 제 2 코사인법칙으로 최종적으로 관측위치상에서의 거리를 구한다
    moon_radius = 1736  # from moon fact sheet - https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
    distance_sun_moon_onpoint = math.sqrt(distance_sun_moon**2 + moon_radius**2 - 2*distance_sun_moon * moon_radius * math.cos(math.radians(lunar_centerangle)))
    # print(lunar_centerangle)
    # 사인법칙으로 관측위치의 고도를 구한다
    i = 180 - math.degrees(math.asin(distance_sun_moon/distance_sun_moon_onpoint * math.sin(math.radians(lunar_centerangle))))
    # -90~-180, 90~180의 아크사인 값 조절
    # if i > 90:
    #     i = i - 90
    # elif i < -90:
    #     i = i + 90

    # 위/경도로 방위도 구한다 (azimuth)
    # solar_azimuth_onpoint = get_azimuth_latlon(b_zero, l_zero, selenographic_lat, selenographic_lon)

    # subsolar point로 방위를 구하는 다른 식
    # https://en.wikipedia.org/wiki/Solar_azimuth_angle
    # s_x = math.cos(math.radians(l_zero - selenographic_lon))
    # s_y = math.cos(math.radians(selenographic_lat))*math.sin(math.radians(b_zero)) - math.sin(math.radians(selenographic_lat))*math.cos(math.radians(b_zero))*math.cos(math.radians(l_zero - selenographic_lon))
    #
    # s_x = math.cos(math.radians(b_zero)) * math.sin(math.radians(l_zero - selenographic_lon))
    # s_y = math.cos(math.radians(selenographic_lat))*math.sin(math.radians(b_zero)) - math.sin(math.radians(selenographic_lat))*math.cos(math.radians(b_zero))*math.cos(math.radians(l_zero - selenographic_lon))
    #
    # solar_azimuth_onpoint = math.degrees(math.atan2(s_x, s_y))

    def solar_azimuth_nasa(latitude_observation, longitude_observation, latitude_subsolar, longitude_subsolar):
        S_x = math.cos(math.radians(latitude_subsolar)) * math.sin(math.radians(longitude_subsolar - longitude_observation))
        S_y = math.cos(math.radians(latitude_observation)) * math.sin(math.radians(latitude_subsolar)) - math.sin(math.radians(latitude_observation)) * math.cos(math.radians(latitude_subsolar)) * math.cos(math.radians(longitude_subsolar - longitude_observation))
        S_z = math.sin(math.radians(latitude_observation)) * math.sin(math.radians(latitude_subsolar)) - math.cos(math.radians(latitude_observation)) * math.cos(math.radians(latitude_subsolar)) * math.cos(math.radians(longitude_subsolar - longitude_observation))
        az = math.atan2(-S_x, -S_y)
        return (math.degrees(az) + 360) % 360

    def solar_azimuth(latitude_observation, longitude_observation, latitude_subsolar, longitude_subsolar):
        print(latitude_observation, longitude_observation, latitude_subsolar, longitude_subsolar)
        latitude_observation = np.radians(latitude_observation)
        longitude_observation = np.radians(longitude_observation)
        latitude_subsolar = np.radians(latitude_subsolar)
        longitude_subsolar = np.radians(longitude_subsolar)

        azimuth = np.arctan2(np.sin(longitude_subsolar - longitude_observation),
                             np.cos(latitude_observation) * np.cos(latitude_subsolar) * np.cos(
                                 longitude_subsolar - longitude_observation) + np.sin(latitude_observation) * np.sin(
                                 latitude_subsolar))
        azimuth = np.degrees(azimuth)
        azimuth = (azimuth + 360) % 360
        azimuth = (450 - azimuth) % 360
        return azimuth

    solar_azimuth_onpoint = solar_azimuth_nasa(selenographic_lat, selenographic_lon, b_zero, l_zero)

    # 위/경도를 구형좌표계로 놓고 구해보기
    # def spheric_to_cartesian(r, lon, lat):
    #     return (r*math.cos(math.radians(lon))*math.cos(math.radians(lat)), r*math.sin(math.radians(lon))*math.cos(math.radians(lat)), r*math.sin(math.radians(lat)))
    # # 방향만 구하면 되니 일단 1로 놓는다
    # (x_sub, y_sub, z_sub) = spheric_to_cartesian(1, l_zero, b_zero)
    # (x_obs, y_obs, z_obs) = spheric_to_cartesian(1, selenographic_lon, selenographic_lat)
    #
    #

    # 스피드업 비활성화
    # print('final result')
    # print(i, h, solar_azimuth_onpoint, distance_sun_moon_onpoint)

    if inverse == True:
        h = h * -1

    return i, h, solar_azimuth_onpoint, distance_sun_moon_onpoint, test_dict

def calc_insolation_lesi(i, moon_sun_distance, earth_sun_distance):
    # s varies from 1361.8 to 1369.2 [W/m^2]
    # textbook (but earth case) suggests 1377 [W/m^2]
    # s = 1339
    s = 1361.8
    distance = moon_sun_distance/(earth_sun_distance*u.AU.to(u.kilometer))
    # return s/distance**2
    return s * math.cos(math.radians(i)) / distance**2, s/distance**2

def calc_solar_earth_distance(t):
    # M : Mean anomaly of the Sun
    # C : Sun's equation of the center
    # v : True anomaly of the Sun
    M = 357.52911+ 35999.05029*t - 0.0001537*(t**2)
    C = (
            (1.914602 - 0.004817*t - 0.000014*(t**2))*math.sin(math.radians(M))
            + (0.019993 - 0.000101*t) * math.sin(math.radians(2*M))
            + 0.000289 * math.sin(math.radians(3*M))
        )
    v = M + C

    # e : eccetricity of the Earth's orbit
    # R : Sun's radius vector / Distance between the centers of the Sun and the Earth
    e = 0.016708634 - 0.000042037*t - 0.0000001267*(t**2)
    R = 1.000001018*(1-e**2)/(1+e*math.cos(math.radians(v)))
    return R

# jd = astropy.time.Time('1977-04-26T00:00:00').jd1
# print(jd)

# t = calc_t(calc_jd(1992, 10, 13))
# t = calc_t(calc_jd(1992, 4, 12))
# t = calc_t(calc_jd(1987, 4, 10))
# calc_solar_attr(t)
# calc_lunar_attr(t)

def calc_point(time, lat, lon):
    void, h, az, dist, test_dict = calc_selenographic_solar_pos(time, lat, lon)
    earth_sun_distance = calc_solar_earth_distance(t)
    insol = calc_insolation_lesi((90 - h), dist, earth_sun_distance)
    if insol < 0:
        insol = 0
    res_dict = {
        'h': h,
        'incidence': 90 - h,
        'az': az,
        'dist': dist,
        'insolation': insol
    }
    return res_dict
    return res_dict['insolation']

def process_tilt_cond(insolation, slope, face_az, az, incidence):
    # 1. 기존에 쓰던 코드
    # # 먼저 수직경사에 대해 정사영을 1번 수행한다.
    # # 경사는 0~180도로 가정
    # angular_diff = abs((90 - slope) - altitude) # 180 - (slope + 90)
    # # 일차 결과물 도출
    # res = insolation * math.cos(math.radians(angular_diff))
    #
    # # 수평방위차 (Azimuth)에 대해 정사영을 1번 더 수행한다.
    # # 정남향을 0으로 정한다.
    # if az >= face_az - 90 and az <= face_az + 90:
    #     angular_diff = abs(face_az - az)
    #     res = res * math.cos(math.radians(angular_diff))
    # else:
    #     res = 0

    # 2. 지구 모델에 맞추어 대응하도록 짠 코드
    zenith = 90 - incidence
    cos_theta = math.cos(math.radians(zenith)) * math.cos(math.radians(slope)) + math.sin(math.radians(zenith)) * math.sin(math.radians(slope)) * math.cos(math.radians(az-face_az))
    # theta = math.degrees(math.acos(cos_theta))
    R_b = cos_theta/math.cos(math.radians(zenith))
    print('--------------------')
    print(cos_theta)
    print(math.cos(math.radians(zenith)))
    print(R_b)

    res = insolation * R_b

    # # 3. Tilt면의 Incidence 계산 코드
    # def calculate_incidence_angle(phi, azimuth, beta, surface_azimuth):
    #     """
    #     Calculates the incidence angle on a tilted surface from the horizontal incidence and azimuth.
    #
    #     Args:
    #     phi -- solar altitude angle in radians
    #     azimuth -- solar azimuth angle in radians (measured clockwise from North)
    #     beta -- tilt angle of the surface in radians
    #     surface_azimuth -- azimuth angle of the surface in radians (measured clockwise from North), default is 0
    #
    #     Returns:
    #     The incidence angle on the tilted surface in radians.
    #     """
    #     phi = math.radians(phi)
    #     beta = math.radians(beta)
    #     azimuth = math.radians(azimuth)
    #     surface_azimuth = math.radians(surface_azimuth)
    #
    #     alpha = azimuth - surface_azimuth
    #     theta = math.acos(math.cos(phi) * math.cos(alpha) * math.cos(beta) + math.sin(phi) * math.sin(beta))
    #     return math.degrees(theta)
    #
    # tilt_incidence = calculate_incidence_angle(incidence, az, slope, face_az)
    return res, math.degrees(math.acos(cos_theta))

def calc_reflected_insolation(insolation, slope, reflectance):
    view_factor = (1 - math.cos(math.radians(slope))) / 2
    return insolation * reflectance * view_factor

def calc_duration(stt, edt, lat, lon, iter=1000, reflectance = 0.3, slope=0, face_az=0):
    stt = rev_calc_jde(stt)
    edt = rev_calc_jde(edt)
    base = stt
    res_arr = []
    for i in range(iter):
        t = calc_t(base + i * (edt-stt)/iter)
        void, h, az, dist, test_dict = calc_selenographic_solar_pos(t, lat, lon)
        # h가 음수인 경우(zenith angle이 90을 초과할 때)
        # if h >= 90:
        #     h = 0
        #     incidence = 0
        # else:
        #     incidence = 90-h
        earth_sun_distance = calc_solar_earth_distance(t)
        insol, normal_insol = calc_insolation_lesi((90-h), dist, earth_sun_distance)
        # face_az = az
        tilt_insol, tilt_incidence = process_tilt_cond(insol, slope, face_az, az, h)
        if insol < 0:
            insol = 0
            normal_insol = 0
            tilt_insol = 0
        if tilt_insol < 0:
            tilt_insol = 0
        res_arr.append({
            'time': rev_calc_jde(t),
            'i': i,
            'h': h,
            'incidence':90-h,
            'tilted_incidence':tilt_incidence,
            'az': az,
            'dist': dist,
            'insolation': insol,
            'tilted_insolation': tilt_insol,
            'normal_insolation': normal_insol,
            'test_raw_val': test_dict['raw_value'],
            'test_asin': test_dict['asin'],
            'test_a': test_dict['a'],
            'test_b': test_dict['b'],
            'test_b1': test_dict['b1'],
            'test_b2': test_dict['b2'],
            'test_b3': test_dict['b3'],
            'sel_colon_sun': test_dict['sel_colon_sun'],
            'sel_lon': test_dict['sel_lon'],
            'lpz': test_dict['lpz'],
            'ldpz': test_dict['ldpz'],
            'tan_a': test_dict['tan_a'],
            'atan_a': test_dict['atan_a'],
            'A_1': test_dict['A_1'],
            'A_2': test_dict['A_2'],
            'A_3': test_dict['A_3'],
            'A_4': test_dict['A_4'],
            'A_5': test_dict['A_5'],
            'A_6': test_dict['A_6'],
            'W': test_dict['W'],
            'beta': test_dict['beta'],
            'I': test_dict['I'],
            'A': test_dict['A']
        })

    import pandas as pd
    import matplotlib.pyplot as plt


    data = pd.DataFrame(res_arr)
    return data['tilted_insolation'].values.tolist(), data['tilted_incidence'].values.tolist(), data['az'].values.tolist(), data['normal_insolation'].values.tolist()
    #
    # data.to_excel('./lunar_result2.xlsx')
    plt.scatter(data.index, data['h'], label='h')
    # plt.scatter(data.index, data['incidence'], label='incidence')
    # plt.scatter(data.index, data['sel_colon_sun'], label='sel_colon_sun')
    # plt.scatter(data.index, data['sel_lon'], label='sel_lon')
    # plt.scatter(data.index, data['test_b3'], label='b3')
    # plt.scatter(data.index, data['lpz'], label='lpz')
    # plt.scatter(data.index, data['ldpz'], label='ldpz')
    plt.scatter(data.index, data['tan_a'], label='tan_a')
    plt.scatter(data.index, data['atan_a'], label='atan_a')
    # plt.scatter(data.index, data['A_1'], label='A_1')
    # plt.scatter(data.index, data['A_2'], label='A_2')
    # plt.scatter(data.index, data['A_3'], label='A_3')
    # plt.scatter(data.index, data['A_4'], label='A_4')
    plt.scatter(data.index, data['A_5']*100, label='A_5')
    plt.scatter(data.index, data['A_6']*100 , label='A_6')
    plt.scatter(data.index, data['W'], label='W')
    # plt.scatter(data.index, data['beta'], label='beta')
    # plt.scatter(data.index, data['I'], label='I')
    plt.scatter(data.index, data['A'], label='A')
    plt.legend()
    plt.title('selenographic solar altitude angle')
    plt.show()
    plt.clf()
    plt.scatter(data.index, data['insolation'], label='insolation')
    plt.title('Lunar-surface effective solar insolation')
    plt.show()
    plt.clf()


# jd = calc_jd(2022, 1, 15)
# t = calc_t(jd)
# print(jd)
# print(t)
# print(calc_point(t, 33, 75))
#
# raise IOError


#
# t_3 = calc_t(calc_jd(2023, 1, 15))
# calc_selenographic_solar_pos(t_3, 0, 0)
# raise IOError

# t = calc_t(calc_jd(2022, 1, 1))
# t_2 = calc_t(calc_jd(2023, 12, 31))
# res = calc_duration(t, t_2, 0, 0, iter=3000, reflectance=0.3)
# raise IOError

# t = calc_t(calc_jd(2022, 5, 21))
#
# print(t)
# print(calc_point(t, 45, 45))
#
# print(calc_jd(2022, 5, 21))
#
# raise IOError

lat_arr = [-90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
lat_arr = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
# lat_arr = [-60, 60]
day_arr = list(range(1, 366))
# print(day_arr)
# raise IOError
res_arr = []
for i in lat_arr:
    t = calc_t(calc_jd(2023, 1, 1))
    t_2 = calc_t(calc_jd(2023, 12, 31))
    res, alt, az, n_res = calc_duration(t, t_2, i, 0, iter=24*365, reflectance=0.12, slope=90, face_az=90)
    res_df = pd.DataFrame({'res_tilt':res, 'incidence_tilt':alt, 'az':az, 'normal_res':n_res})
    res_df.to_excel('./year_result_wall(east).xlsx')
    raise IOError
    res_arr.append(res)

for i in lat_arr:
    idx = lat_arr.index(i)
    plt.plot(day_arr, res_arr[idx], label=i)

plt.legend()
plt.show()


# 3D 그래프 플로팅 코드
# fig = plt.figure(figsize=(8, 3))
# ax1 = fig.add_subplot(121, projection='3d')
# # ax2 = fig.add_subplot(122, projection='3d')
#
# _xx, _yy = np.meshgrid(day_arr, lat_arr)
# x, y = _xx.ravel(), _yy.ravel()
#
# # t = calc_t(calc_jd(2023, 1, 1))
#
# # top = calc_point(calc_t(calc_jd(2023, 1, 1) + x), y, 0)
# top = np.transpose(np.array(res_arr))
# top = top.ravel()
# print(top)
# bottom = np.zeros_like(top)
# # print(x)
# # print(y)
# print(len(day_arr))
# print(len(lat_arr))
# print(np.shape(top))
# print(np.ndim(top))
#
# width = 1
# depth = 10
#
# ax1.bar3d(x, y, bottom, width, depth, top, shade=True)
# ax1.set_title('Shaded')
# plt.show()
raise IOError

t = calc_t(calc_jd(2023, 1, 1))
t_2 = calc_t(calc_jd(2023, 12, 31))
calc_duration(t, t_2, 45, 45, iter=3650, reflectance = 0.3)




raise IOError
print(t)







raise IOError





















print(astropy.time.Time('1992-04-12T00:00:00').jd1)
jd = astropy.time.Time('1992-04-12T00:00:00').jd1
D = jd-2451545.0
g = 357.529 + 0.98560028*D
q = (280.459 + 0.98564736*D) % 360
print(q)

print(g)
print(math.radians(g))
L = q + 1.915*math.sin(math.radians(g)) + 0.0020 * math.sin(math.radians(2*g))
print(L)
print('=================================================================')
lambda_h = L + 180 + 368438.72865299/process_decimal('1.49971545e+08') * 57.296 * math.cos(math.radians(13.73591361)) * math.sin(math.radians(L - 134.79900662))
print(lambda_h)

print(get_body('sun', astropy.time.Time('1992-04-12T00:00:00')))
print(get_body('moon', astropy.time.Time('1992-04-12T00:00:00')))

raise IOError


print(moon_pos)
moon_pos = get_axial_element(moon_pos)
print(sun_pos)
sun_pos = get_axial_element(sun_pos)
print(earth_pos)
earth_pos = get_axial_element(earth_pos)

# selenographic lon/lat (입력 값)
moon_lon = 0
moon_lat = 0

# sub_solar point lon/lat 계산
# 태양과 달의 좌표를 기준으로 구한다
# 문헌에서 언급된 길이단위는 AU라서, 이걸 km로 환산해야 한다.
moon_radius = process_decimal('1.161363636e-5')*u.AU.to(u.kilometer)
print(moon_radius)
# 달 반지름이 1737.4로 나오니 얼추 맞는 것으로 보임 (계산결과값 1737.3752705400987)

def get_subsolar_point(sun_pos, target_pos, target_radius):
    # Cartesian coord로 구한다.
    vec = [sun_pos[0]-target_pos[0], sun_pos[1]-target_pos[1], sun_pos[2]-target_pos[2]]
    # print(vec)
    dist = math.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    # print(dist)
    res = [target_pos[0] + vec[0] * target_radius / dist,
           target_pos[1] + vec[1] * target_radius / dist,
           target_pos[2] + vec[2] * target_radius / dist,
           ]
    return res
    # res-target pos 로 달의 반지름이 다시 나오는지를 검산
    # print(math.sqrt((target_pos[0]-res[0])**2 + (target_pos[1]-res[1])**2 + (target_pos[2]-res[2])**2))
    # 1737.3752705400987 -> 1737.3752705390932
    # 대략 0.0000000010055 (km)의 오차.... (0.001mm임)

subsolar_point = get_subsolar_point(sun_pos, moon_pos, moon_radius)



