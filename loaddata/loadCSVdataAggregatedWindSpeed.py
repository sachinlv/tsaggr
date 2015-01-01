__author__ = 'freak'

import MySQLdb
from os import listdir
import csv

db = MySQLdb.connect(host="localhost",
                     user="sachin",
                      passwd="password",
                      db="eastwind")
cur = db.cursor()

file_dir = "/usr/wind_data/"
file_list = listdir(file_dir)


#for file in file_list[301:600]:
for file in file_list[310:400]:
    if(file.endswith('.CSV')):
        file_path = file_dir + file
        print('loading file:: ' + file_path)

        table_name = file.split('.')[0]
        create_table = "CREATE TABLE IF NOT EXISTS aggr_onshore_" + \
                       table_name + \
                       "( mesdt INT NULL, avgwind FLOAT NULL, power NVARCHAR(8000));"
        cur.execute(create_table)
        db.commit()
        tbl_data_count = "SELECT COUNT(*) FROM aggr_onshore_" + table_name
        cur.execute(tbl_data_count)
        count = cur.fetchone()
        if count[0] == 0:
            with open(file_path, 'rb') as f:
                i = 3
                while(i): #skip first 3 lines of csv file
                    next(f)
                    i -= 1

                csvread = csv.reader(f)

                count = 144
                avg_wind_speed = 0.0
                power = ''
                date = 0
                for row in csvread:
                    if(count):
                        date = int(row[0])
                        avg_wind_speed = avg_wind_speed + float(row[2])
                        power = power + row[3] + ','
                        count -= 1
                    elif not count:
                        avg_wind_speed=avg_wind_speed/144
                        power = power[:-1]
                        load_data = "INSERT INTO aggr_onshore_" + table_name +"(mesdt, avgwind, power)"+\
                            " VALUES (%d, %f, '%s') " % (date, avg_wind_speed, power)
                        cur.execute(load_data)
                        db.commit()
                        count = 144
                        power = ''
                        avg_wind_speed = 0
                        date = 0

        else:
            print("File already loaded")
