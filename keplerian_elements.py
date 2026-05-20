# Mobin Ravan

import numpy as np
import matplotlib.pyplot as plt

# Earth's gravitational constant (WGS84)
GM = 398600.4418  # km^3/s^2


#----------- Part 1 Proj: Kepler elements 

# data for S1 (INTELSAT-20) - position vector & velocity vector 
r1_x = 16575.5879
r1_y = -38775.5254
r1_z = 5.3421

v1_x = 2.8271
v1_y = 1.2079
v1_z = -0.0008

# data for S2 (MONILIYA3-50) - position vector & velocity vector 
r2_x = -495.4135
r2_y = 8636.4268
r2_z = -3136.0422
 
v2_x = -4.0917
v2_y = 5.5089
v2_z = 4.9534

# function to calculate magnitude of a vector
def norm3(x, y, z):
    return np.sqrt(x*x + y*y + z*z)

# function to calculate dot product of two vectors
def dot3(x1, y1, z1, x2, y2, z2):
    return x1*x2 + y1*y2 + z1*z2

# function to calculate cross product of two vectors (returns a new vector)
def cross3(x1, y1, z1, x2, y2, z2):
    cx = y1*z2 - z1*y2
    cy = z1*x2 - x1*z2
    cz = x1*y2 - y1*x2
    return cx, cy, cz

# main function to calculate Keplerian elements from position and velocity
def kepler_from_rv(rx, ry, rz, vx, vy, vz, GM):
    # 1 - specific angular momentum vector (h)
    hx, hy, hz = cross3(rx, ry, rz, vx, vy, vz)
    h = norm3(hx, hy, hz)
    
    # 2 - ascending node vector (n) = Z cross (h)
    Zx, Zy, Zz = 0, 0, 1
    nx, ny, nz = cross3(Zx, Zy, Zz, hx, hy, hz)
    n = norm3(nx, ny, nz)
    
    # 3 - inclination angle (i)
    i_rad = np.arccos(hz / h)
    
    # 4 - right ascension of ascending node (Omega)
    if n != 0:
        Omega_rad = np.arccos(nx / n)
        if ny < 0:
            Omega_rad = 2*np.pi - Omega_rad
    else:
        Omega_rad = 0
    
    # 5 - semi-major axis (a) from energy equation
    r = norm3(rx, ry, rz)
    v = norm3(vx, vy, vz)
    E_specific = (v*v)/2 - GM/r
    a = -GM / (2 * E_specific)
    
    # 6 - eccentricity (e) from h and a
    e = np.sqrt(1 - (h*h) / (GM * a))
    
    # 7 - eccentricity vector (e_vec) pointing toward perigee
    # formula: e_vec = ( (v^2 - GM/r) * r_vec - (r_vec . v_vec) * v_vec ) / GM
    r_dot_v = dot3(rx, ry, rz, vx, vy, vz)
    ev_x = ((v*v - GM/r) * rx - r_dot_v * vx) / GM
    ev_y = ((v*v - GM/r) * ry - r_dot_v * vy) / GM
    ev_z = ((v*v - GM/r) * rz - r_dot_v * vz) / GM
    
    # 8 - argument of perigee (omega)
    if n != 0:
        cos_omega = (nx*ev_x + ny*ev_y + nz*ev_z) / (n * e)
        # to avoid small numerical errors
        if cos_omega > 1:
            cos_omega = 1
        if cos_omega < -1:
            cos_omega = -1
        omega_rad = np.arccos(cos_omega)
        if ev_z < 0:
            omega_rad = 2*np.pi - omega_rad
    else:
        omega_rad = 0
    
    # 9 - true anomaly (nu)
    cos_nu = (ev_x*rx + ev_y*ry + ev_z*rz) / (e * r)
    if cos_nu > 1:
        cos_nu = 1
    if cos_nu < -1:
        cos_nu = -1
    nu_rad = np.arccos(cos_nu)
    if r_dot_v < 0:
        nu_rad = 2*np.pi - nu_rad
    
    # 10 - eccentric anomaly (E) and mean anomaly (M) (calculated M to show its value for better analysis)
    cos_E = (e + np.cos(nu_rad)) / (1 + e * np.cos(nu_rad))
    if cos_E > 1:
        cos_E = 1
    if cos_E < -1:
        cos_E = -1
    E_rad = np.arccos(cos_E)
    if nu_rad > np.pi:
        E_rad = 2*np.pi - E_rad
    M_rad = E_rad - e * np.sin(E_rad)
    
    # convert radians to degrees
    i_deg = i_rad * 180.0 / np.pi
    Omega_deg = Omega_rad * 180.0 / np.pi
    omega_deg = omega_rad * 180.0 / np.pi
    nu_deg = nu_rad * 180.0 / np.pi
    M_deg = M_rad * 180.0 / np.pi
    
    return a, e, i_deg, Omega_deg, omega_deg, M_deg, nu_deg

# calculation for S1 (INTELSAT-20)
a1, e1, i1, Omega1, omega1, M1, nu1 = kepler_from_rv(r1_x, r1_y, r1_z, v1_x, v1_y, v1_z, GM)
print("-"*22)
print("Satellite INTELSAT-20")
print("="*30)
print("a (km)          :", round(a1, 6))
print("e              :", round(e1, 6))
print("i (deg)         :", round(i1, 6))
print("Omega (deg)      :", round(Omega1, 6))
print("Omega(AP) (deg)    :", round(omega1, 6))
print("M (deg)          :", round(M1, 6))
print("nu (deg)      :", round(nu1, 6))

# calculation for satellite 2 (MONILIYA3-50)
a2, e2, i2, Omega2, omega2, M2, nu2 = kepler_from_rv(r2_x, r2_y, r2_z, v2_x, v2_y, v2_z, GM)

print("\n" + "."*29)
print("Satellite MONILIYA3-50")
print("="*30)
print("a (km)       :", round(a2, 6))
print("e              :", round(e2, 6))
print("i (deg)         :", round(i2, 6))
print("Omega (deg)      :", round(Omega2, 6))
print("Omega(AP) (deg)    :", round(omega2, 6))
print("M (deg)        :", round(M2, 6))
print("nu (deg)        :", round(nu2, 6))

#-----
#--------------- Part 2 Proj: Draw orbit ellipse
def draw_orbit(a, e, nu_deg, title, sat_name):
    nu_rad = nu_deg * np.pi / 180.0
    
    # calculate ellipse points for angles from 0 to 360 degrees
    theta_list = []
    r_list = []
    for t in range(0, 361):
        rad = t * np.pi / 180.0
        r_t = a * (1 - e*e) / (1 + e * np.cos(rad))
        x_t = r_t * np.cos(rad)
        y_t = r_t * np.sin(rad)
        theta_list.append(rad)
        r_list.append(r_t)
    
    # convert to arrays for plotting
    x_orbit = []
    y_orbit = []
    for i in range(len(theta_list)):
        r_tmp = a * (1 - e*e) / (1 + e * np.cos(theta_list[i]))
        x_orbit.append(r_tmp * np.cos(theta_list[i]))
        y_orbit.append(r_tmp * np.sin(theta_list[i]))
    
    # satellite position at given time
    r_sat = a * (1 - e*e) / (1 + e * np.cos(nu_rad))
    x_sat = r_sat * np.cos(nu_rad)
    y_sat = r_sat * np.sin(nu_rad)
    
    # ellipse center relative to focus (Earth)
    center_x = -a * e
    center_y = 0
    
    # plotting
    plt.figure(figsize=(8, 8))
    plt.plot(x_orbit, y_orbit, 'b-', linewidth=2)
    plt.plot(0, 0, 'go', markersize=10)
    plt.plot(x_sat, y_sat, 'ro', markersize=8)
    plt.plot(center_x, center_y, 'ks', markersize=5)
    plt.plot([0, center_x], [0, center_y], 'k--', linewidth=0.8, alpha=0.5)
    
    plt.xlabel('x (km) - Perifocal frame')
    plt.ylabel('y (km) - Perifocal frame')
    plt.title(title)
    plt.legend(['Orbit', 'Earth (focus)', sat_name + ' position', 'Ellipse center', 'Line of apsides'])
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    plt.show()

# draw orbit for S1
draw_orbit(a1, e1, nu1, 'INTELSAT-20 Orbit (GEO)', 'INTELSAT-20')

# draw orbit for S2
draw_orbit(a2, e2, nu2, 'MONILIYA3-50 Orbit', 'MONILIYA3-50')

#----------------------------
#----------------------- Part 3: Hill's elements for INTELSAT-20
print("\n" + "-"*37)
print("Hill's elements for INTELSAT-20 (GEO)")
print("-"*37)

# calculate r (distance)
r_sat = a1 * (1 - e1*e1) / (1 + e1 * np.cos(nu1 * np.pi/180.0))
print("r (km)          :", round(r_sat, 6))

# radial velocity (r_dot)
h_val = np.sqrt(GM * a1 * (1 - e1*e1))
nu_rad = nu1 * np.pi / 180.0
r_dot = (e1 * h_val * np.sin(nu_rad)) / (a1 * (1 - e1*e1))
print("r_dot (km/s)   :", "{:.6e}".format(r_dot))

# argument of latitude u = omega + nu
omega_rad = omega1 * np.pi / 180.0
u_rad = omega_rad + nu_rad
print("u (rad)        :", round(u_rad, 6), " ", "(", round(u_rad*180/np.pi, 6), "deg)")

# Omega (right ascension of ascending node)
Omega_rad = Omega1 * np.pi / 180.0
print("Omega (rad)     :", round(Omega_rad, 6), " ", "(", round(Omega_rad*180/np.pi, 6), "deg)")

# G = specific angular momentum
G_val = h_val
print("G (km²/s)      :", round(G_val, 6))

# H = G * cos(i) (component along Z axis)
i_rad = i1 * np.pi / 180.0
H_val = G_val * np.cos(i_rad)
print("H (km²/s)       :", round(H_val, 6))



#----- part 3 additional: proof of Geostationary Orbit (GEO) for INTELSAT-20
#---
print("\n" + "-"*35)
print("Proof: INTELSAT-20 is in Geostationary Orbit (GEO)")
print("\n" + "-"*35)

# step by step using Kepler's 3rd law
# Kepler's 3rd law: T = 2*pi * sqrt(a^3 / GM)

a_geo = a1
GM_earth = GM

a_3 = a_geo * a_geo * a_geo

a_3_div_GM = a_3 / GM_earth

sqrt_part = np.sqrt(a_3_div_GM)

T_seconds = 2 * np.pi * sqrt_part

# convert seconds to hours
T_hours = T_seconds / 3600

print("\n")
print("According to Kepler's 3rd law: T = 2π √(a³/GM)")
print("Orbital period of INTELSAT-20 =", round(T_seconds, 2), "seconds =", round(T_hours, 4), "hours")

# check if T_hours is approximately 24 or not
if T_hours >= 23.9 and T_hours <= 24.1:
    print("\n")
    print("Since the orbital period is about 24 hours (≈", round(T_hours, 2), "h),")
    print("   INTELSAT-20 is in a Geostationary Orbit (GEO).")
else:
    print("\n")
    print("Orbital period is not 24 hours, therefore it is not GEO.")