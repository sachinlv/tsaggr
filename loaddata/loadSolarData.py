__author__ = 'freak'
import MySQLdb
from os import listdir

db = MySQLdb.connect(host="localhost",
                     user="sachin",
                      passwd="password",
                      db="solar")
cur = db.cursor()

file_dir = "/usr/solar_data/specific_sites3/"
#file_list = listdir(file_dir)
file_list = ['specific_site30.csv',
             'specific_site31.csv',
             'specific_site32.csv',
             'specific_site33.csv',
             'specific_site34.csv',
             'specific_site35.csv',
             'specific_site36.csv',
             'specific_site37.csv',
             'specific_site38.csv',
             'specific_site39.csv']

file_list.sort()
for file in file_list[0:10]:
    if(file.endswith('.csv')):
        file_path = file_dir + file
        print('loading file:: ' + file_path)
        table_name = file.split('.')[0]
        create_table = "CREATE TABLE IF NOT EXISTS solar_" + \
                       table_name + \
                       "( mesdt INT NULL, mestim INT NULL, rad FLOAT NULL, pow FLOAT NULL);"
        cur.execute(create_table)
        db.commit()
        tbl_data_count = "SELECT COUNT(*) FROM solar_" + table_name
        cur.execute(tbl_data_count)
        count = cur.fetchone()
        if count[0] == 0:
            load_data = "LOAD DATA INFILE '" + file_path + "' INTO TABLE solar_" + table_name + \
                        " FIELDS TERMINATED BY ',' ENCLOSED BY '" + \
                        '"' + \
                        "' LINES TERMINATED BY '\\n' IGNORE 1 LINES;"

            print(load_data)
            cur.execute(load_data)
            db.commit()
        else:
            print("File already loaded")

