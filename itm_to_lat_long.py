import csv
from math import *


############################################################################

# Meridian Arc

############################################################################

def arcmer(a, equad, lat1, lat2):
    b = a * sqrt(1 - equad)
    n = (a - b) / (a + b)
    a0 = 1. + ((n ** 2) / 4.) + ((n ** 4) / 64.)
    a2 = (3. / 2.) * (n - ((n ** 3) / 8.))
    a4 = (15. / 16.) * ((n ** 2) - ((n ** 4) / 4.))
    a6 = (35. / 48.) * (n ** 3)
    s1 = a / (1 + n) * (a0 * lat1 - a2 * sin(2. * lat1) + a4 * sin(4. * lat1) - a6 * sin(6. * lat1))
    s2 = a / (1 + n) * (a0 * lat2 - a2 * sin(2. * lat2) + a4 * sin(4. * lat2) - a6 * sin(6. * lat2))
    return s2 - s1


###############################################################################
#
# Transverse Mercator Inverse Projection
#
###############################################################################

def xy2geo(m, p, a, equad, lat0, lon0):
    lat0 = radians(lat0)
    lon0 = radians(lon0)
    sigma1 = p
    fil = lat0 + sigma1 / (a * (1 - equad))
    deltafi = 1

    while deltafi > 0.0000000001:
        sigma2 = arcmer(a, equad, lat0, fil)
        RO = a * (1 - equad) / ((1 - equad * (sin(fil) ** 2)) ** (3. / 2.))
        deltafi = (sigma1 - sigma2) / RO
        fil = fil + deltafi

    N = a / sqrt(1 - equad * (sin(fil)) ** 2)
    RO = a * (1 - equad) / ((1 - equad * (sin(fil) ** 2)) ** (3. / 2.))
    t = tan(fil)
    psi = N / RO
    lat = fil - (t / RO) * ((m ** 2) / (2. * N)) + (t / RO) * ((m ** 4) / (24. * (N ** 3))) * (
                -4. * (psi ** 2) - 9. * psi * (1. - t ** 2) + 12. * (t ** 2)) - (t / RO) * (
                      m ** 6 / (720. * (N ** 5))) * (
                      8. * (psi ** 4) * (11. - 24. * (t ** 2)) - 12. * (psi ** 3) * (21. - 71. * (t ** 2)) + 15. * (
                          psi ** 2) * (15. - 98. * (t ** 2) + 15. * (t ** 4)) + 180. * psi * (
                                  5. * (t ** 2) - 3. * (t ** 4)) - 360. * (t ** 4)) + (t / RO) * (
                      (m ** 8) / (40320. * (N ** 7))) * (1385. + 3633. * (t ** 2) + 4095. * (t ** 4) + 1575. * (t ** 6))
    lon = (m / (N)) - ((m ** 3) / (6. * (N ** 3))) * (psi + 2. * (t ** 2)) + ((m ** 5) / (120. * (N ** 5))) * (
                -4. * (psi ** 3) * (1. - 6. * (t ** 2)) + (psi ** 2) * (9. - 68. * (t ** 2)) + 72. * psi * (
                    t ** 2) + 24. * (t ** 4)) - ((m ** 7) / (5040. * (N ** 7))) * (
                      61. + 662. * (t ** 2) + 1320. * (t ** 4) + 720. * (t ** 6))
    lon = lon0 + lon / cos(fil)
    lat = degrees(lat)
    lon = degrees(lon)

    return lat, lon


#############################################################################

# Irish Transverse Mercator - Inverse

#############################################################################

def itm2geo(x, y):
    # GRS-80
    a = 6378137.
    equad = 0.00669437999
    # Natural Origin
    lat0 = 53.5
    lon0 = -8.
    k0 = 0.999820
    p = (y - 750000.) / k0
    m = (x - 600000.) / k0
    lat, lon = xy2geo(m, p, a, equad, lat0, lon0)
    return lat, lon


#############################################################################

# Process CSV by iterating through a list of ITM's and returning latitude and longitude

#############################################################################

def process_csv(input_file, output_file):
    try:
        with open(input_file, newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            with open(output_file, 'w', newline='', encoding='utf-8') as csvfile_out:
                fieldnames = ['x', 'y', 'latitude', 'longitude']
                writer = csv.DictWriter(csvfile_out, fieldnames=fieldnames)
                writer.writeheader()

                for row in reader:
                    x_str = row['x'].strip()
                    y_str = row['y'].strip()

                    if x_str and y_str:  # Check if both x and y are not empty
                        try:
                            x = float(x_str)
                            y = float(y_str)
                            lat, lon = itm2geo(x, y)
                            writer.writerow({'x': x, 'y': y, 'latitude': lat, 'longitude': lon})
                        except ValueError as e:
                            print(f"Skipping row with invalid data {row}: {e}")
                    else:
                        print(f"Skipping row with missing data: {row}")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except UnicodeDecodeError as e:
        print(f"Encoding error: {e}")


if __name__ == '__main__':
    input_file = '/Desktop/coordinates.csv'
    output_file = '/Desktop/converted_coordinates.csv'
    process_csv(input_file, output_file)

#############################################################################


