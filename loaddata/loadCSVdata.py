__author__ = 'freak'
import MySQLdb
from os import listdir

db = MySQLdb.connect(host="localhost",
                     user="sachin",
                      passwd="password",
                      db="eastwind")
cur = db.cursor()

file_dir = "/usr/wind_data/"
#file_list = listdir(file_dir)
file_list = ['SITE_00538.CSV',
             'SITE_00366.CSV',
             'SITE_00623.CSV',
             'SITE_00418.CSV',
             'SITE_00627.CSV',
             'SITE_00532.CSV',
             'SITE_00499.CSV',
             'SITE_00571.CSV',
             'SITE_03247.CSV',
             'SITE_00622.CSV']

file_list.sort()
for file in file_list[0:9]:
    if(file.endswith('.CSV')):
        file_path = file_dir + file
        print('loading file:: ' + file_path)
        table_name = file.split('.')[0]
        create_table = "CREATE TABLE IF NOT EXISTS onshore_" + \
                       table_name + \
                       "( mesdt INT NULL, mestim INT NULL, spd FLOAT NULL, pow FLOAT NULL);"
        cur.execute(create_table)
        db.commit()
        tbl_data_count = "SELECT COUNT(*) FROM onshore_" + table_name
        cur.execute(tbl_data_count)
        count = cur.fetchone()
        if count[0] == 0:
            load_data = "LOAD DATA INFILE '" + file_path + "' INTO TABLE onshore_" + table_name + \
                        " FIELDS TERMINATED BY ',' ENCLOSED BY '" + \
                        '"' + \
                        "' LINES TERMINATED BY '\\n' IGNORE 3 LINES;"

            print(load_data)
            cur.execute(load_data)
            db.commit()
        else:
            print("File already loaded")
