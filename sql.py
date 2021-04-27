import os
import sqlite3
from sqlite3 import Error

def newConnection(database):
    con = None
    try:
        con = sqlite3.connect(database)
        print(sqlite3.version)
    except Error as e:
        print(e)
    return con

def newTable(connection, makeTable):
    try:
        cur = connection.cursor()
        cur.execute(makeTable)
    except Error as e:
        print(e)

def insertData(connection, entry):
    sql = """ INSERT INTO ? (name,field2,field3)
                VALUES(?,?,?) """
    cur = connection.cursor()
    cur.execute(sql, entry)
    connection.commit()
    return cur.lastrowid()

test_table = """ CREATE TABLE IF NOT EXISTS test (
                    id integer PRIMARY KEY,
                    name text NOT NULL,
                    field2 text,
                    field3 text
                ); """
entry = ('test', 'name_field', 'abc', '123')

conn = newConnection('test.db')
newTable(conn, test_table)
lri = insertData(conn, entry)
