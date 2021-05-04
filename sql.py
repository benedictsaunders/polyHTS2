import os
import sqlite3
import re
from sqlite3 import Error
from utils import *
import constants

def sanitizeSQL(s):
    return re.sub(r'[^A-Za-z0-9_.]', '', s)

def split(s):
    return [x for x in s]

def newConnection(database):
    con = None
    try:
        con = sqlite3.connect(database)

    except Error as e:
        print(e)
    return con

def newTable(connection, tableName, repeat_style):
    cols = ["id INT NOT NULL PRIMARY KEY", "length INTEGER"]
    for col in constants.DATACOLS:
        cols.append(f"{col} REAL")
    for idx, col in enumerate(cols):
        if "smiles REAL" in col:
            cols[idx] = "smiles BLOB"
    mon_cols_list = []
    rs = list(set(split(repeat_style)))
    rs.sort(reverse=True)
    for r in rs:
        element = " ".join([sanitizeSQL(r), 'BLOB'])
        cols.insert(1, element)
    sql_cols = ", ".join(cols)
    sql = f"CREATE TABLE IF NOT EXISTS {sanitizeSQL(tableName)} ({sql_cols});"
    try:
        cur = connection.cursor()
        cur.execute(sql)
    except Error as e:
        print(e)

def insertData(connection, data, tableName, repeat_style):
    data = listelementsToString(data)
    rs = list(set(split(repeat_style)))
    rs.sort(reverse=True)
    cols = constants.DATACOLS[:]
    cols.insert(0, "length")
    for r in rs:
        cols.insert(0, r)
    cols.insert(0,"id")
    sql = f"INSERT INTO {sanitizeSQL(tableName)} ({', '.join(cols)}) VALUES ({', '.join(data)});"
    cur = connection.cursor()
    cur.execute(sql)
    connection.commit()

def updateData(connection, tableName, data, id):
    print(constants.DATACOLS)
    print(data)
    setter = []
    for idx, col in enumerate(constants.DATACOLS):
        setter.append(f"{col} = {str(data[idx])}")
    sql = f"UPDATE {sanitizeSQL(tableName)} SET {' ,'.join(setter)} WHERE id = {str(id)};"
    cur = connection.cursor()
    cur.execute(sql)
    connection.commit()
