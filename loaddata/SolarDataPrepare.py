__author__ = 'freak'
from os import listdir
import csv
import math
import codecs

power_data_path = '/usr/claudio_solar_data/Arizona/DA/DPV/'
radiation_sites_file = '/usr/claudio_solar_data/Arizona/SolarSitesArizona.csv'
power_rad_file = '/usr/claudio_solar_data/Arizona/powr_rad_site_combi.csv'

power_data_files = listdir(power_data_path)
pr = csv.writer(codecs.open(power_rad_file, 'wb', 'utf-16'), delimiter=',')



for powfile in power_data_files:
    site_details = powfile.split('_')
    powlat = float(site_details[1])
    powlong = float(site_details[2])
    dist_site = {}

    rad_sites = csv.reader(codecs.open(radiation_sites_file, 'rU', 'utf-16'), delimiter=',')
    for rad in rad_sites:
        id = int(rad[0])
        radlat = float(rad[5])
        radlong = float(rad[6])
        dist_site[id] = math.sqrt(math.pow((powlat - radlat), 2) + math.pow((powlong - radlong), 2))

    #print(dist_site)
    site = min(dist_site, key=dist_site.get)
    #print(site)
    pr.writerow([powfile, powlat, powlong, site])


if __name__ == "__main__":
    pass