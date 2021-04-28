import os
import sqlite3
import re
from sqlite3 import Error
from utils import *
import constants

def sanitizeSQL(s):
    return re.sub(r'[^A-Za-z0-9 ]+', '', s)

def split(s):
    return [x for x in s]

def newConnection(database):
    con = None
    try:
        con = sqlite3.connect(database)

    except Error as e:
        print(e)
    return con

def newTable(connection, repeat_style):
    cols = ["id INT NOT NULL PRIMARY KEY", "length INTEGER"]
    for col in constants.DATACOLS:
        cols.append(f"{col} REAL")
    mon_cols_list = []
    rs = list(set(split(repeat_style)))
    rs.sort(reverse=True)
    for r in rs:
        element = " ".join([sanitizeSQL(r), 'BLOB'])
        cols.insert(1, element)
    sql_cols = ", ".join(cols)
    sql = f"CREATE TABLE IF NOT EXISTS results ({sql_cols});"
    try:
        cur = connection.cursor()
        cur.execute(sql)
    except Error as e:
        print(e)

def insertData(connection, data, repeat_style):
    data = listelementsToString(data)
    rs = list(set(split(repeat_style)))
    rs.sort(reverse=True)
    cols = constants.DATACOLS[:]
    cols.insert(0, "length")
    for r in rs:
        cols.insert(0, r)
    cols.insert(0,"id")
    sql = f"INSERT INTO results ({', '.join(cols)}) VALUES ({', '.join(data)});"
    cur = connection.cursor()
    cur.execute(sql)
    connection.commit()

def updateData(connection, data, id):
    setter = []
    for idx, col in enumerate(constants.DATACOLS):
        setter.append(f"{col} = {str(data[idx])}")
    sql = f"UPDATE results SET {' ,'.join(setter)} WHERE id = {str(id)};"
    cur = connection.cursor()
    cur.execute(sql)
    connection.commit()
