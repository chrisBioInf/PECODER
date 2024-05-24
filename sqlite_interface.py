#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:45:16 2024

@author: christopher
"""



db_columns = {
    'table_ecod': 'ecod', 
    'pdb' :'VARCHAR(10)',
    'chain' : 'VARCHAR(2)',
    'number' : 'INT',
    'pdb_range' : 'VARCHAR(50)',
    't_id': 'VARCHAR(10)',
    'ecod_domain_id' : 'VARCHAR(10)',
    'ecod_a_name': 'VARCHAR(75)',
    'ecod_x_name': 'VARCHAR(75)',
    'ecod_h_name': 'VARCHAR(75)',
    'ecod_t_name': 'VARCHAR(75)',
    'ecod_f_name': 'VARCHAR(75)',
    'ecod_ligand': 'VARCHAR(75)',
    }


interactions_columns = {
    'table_ecod': 'ecod', 
    'pdb' :'VARCHAR(10)',
    'chain' : 'VARCHAR(2)',
    'number' : 'INT',
    'residue': 'VARCHAR(3)',
    'distance': 'FLOAT',
    'cofactor': 'VARCHAR(10)',
    'pdb_range' : 'VARCHAR(50)',
    'ecod_domain_id' : 'VARCHAR(10)',
    'ecod_a_name': 'VARCHAR(75)',
    'ecod_x_name': 'VARCHAR(75)',
    'ecod_h_name': 'VARCHAR(75)',
    'ecod_t_name': 'VARCHAR(75)',
    'ecod_f_name': 'VARCHAR(75)',
    'ecod_ligand': 'VARCHAR(75)',
    }


def create_table_if_not_existing(connection):
    cursor = connection.cursor()
    cursor.execute("DROP TABLE IF EXISTS ecod")
    cursor.execute("""CREATE TABLE {table_ecod} (
        pdb {pdb},
        chain {chain},
        number {number},
        pdb_range {pdb_range},
        t_id {t_id},
        ecod_domain_id {ecod_domain_id},
        ecod_a_name {ecod_a_name}, 
        ecod_x_name {ecod_x_name}, 
        ecod_h_name {ecod_h_name}, 
        ecod_t_name {ecod_t_name},
        ecod_f_name {ecod_f_name},
        ecod_ligand {ecod_ligand}
        )""".format(**db_columns))
    
    cursor.execute("DROP TABLE IF EXISTS interactions")
    cursor.execute("""CREATE TABLE interactions (
        pdb {pdb},
        chain {chain},
        number {number},
        residue {residue},
        distance {distance},
        cofactor {cofactor},
        pdb_range {pdb_range},
        ecod_domain_id {ecod_domain_id},
        ecod_a_name {ecod_a_name}, 
        ecod_x_name {ecod_x_name}, 
        ecod_h_name {ecod_h_name}, 
        ecod_t_name {ecod_t_name},
        ecod_f_name {ecod_f_name},
        ecod_ligand {ecod_ligand}
        )""".format(**interactions_columns))
    
    connection.commit()
    return cursor.lastrowid
    

def add_rows_to_table(connection, table, df):
    cursor = connection.cursor()
    val_string = ",".join(['?' for x in df.columns])
    sql = "INSERT INTO {table} VALUES ({value_string})".format(table=table,
                                                              value_string=val_string)
    for row in list(df.itertuples(index=False, name=None)):
        cursor.execute(sql, row)
    
    connection.commit()
    return cursor
    

