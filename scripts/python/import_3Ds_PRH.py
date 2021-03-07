#!/usr/bin/env python

"""
Annotate and import solved 3D structures into MongoDB
"""

import sys, os, math, datetime
from pyrna.task import Task
from pyrna import parsers
from pyrna.computations import Rnaview
from bson.objectid import ObjectId
from pymongo import MongoClient
from pyrna.db import PDB as db_PDB
from pyrna.db import RNA3DHub as db_RNA3D
import urllib


import requests, json

def query_by_dict(dict_v):
    json_str = json.dumps(dict_v, indent=0).replace("\n", "")
    url = "https://search.rcsb.org/rcsbsearch/v1/query?json={:s}".format(json_str)
    response = requests.get(url)
    return response

class Client():
    def __init__(self):
        pass


    def query_all_RNAs(self):
        # see: https://search.rcsb.org/search-attributes.html
        # rcsb_entry_info.polymer_entity_count_RNA
        #   Should be > 0 for RNA
        # rcsb_entry_info.polymer_entity_count_DNA
        #   Should be zero for no DNA
        # rcsb_entry_info.polymer_entity_count_protein
        #   Should b e zero for no protein
        # rcsb_entry_info.polymer_entity_count_nucleic_acid_hybrid
        #   should be zero for no RNA
        json_dict = {
              "query": {
                  "type": "group",
                  "logical_operator": "and",
                  "nodes": [
                      {
                        "type": "terminal",
                        "service": "text",
                          # should have > 0 RNA
                          "parameters": {
                              "attribute": "rcsb_entry_info.polymer_entity_count_RNA",
                              "operator": "greater",
                              "value": 0
                          }
                    },
                      {
                          # should have = 0 DNA
                          "type": "terminal",
                          "service": "text",
                          "parameters": {
                              "attribute": "rcsb_entry_info.polymer_entity_count_DNA",
                              "operator": "equals",
                              "value": 0
                          }},
                      {
                          # should have = 0 protein
                          "type": "terminal",
                          "service": "text",
                          "parameters": {
                              "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                              "operator": "equals",
                              "value": 0
                          }},

                      ]
              },
            # get everything
            "request_options": {
                "return_all_hits": True
            },
              "return_type": "entry"
            }
        response = query_by_dict(json_dict)
        assert response.ok
        response_json = json.loads(response.text)
        ids = \
            sorted(set([i["identifier"] for i in response_json['result_set']]))
        return ids


def single_3ds_cluster_or_None(db,pdb_old_version,pdb_id,rnaview):
    try:
        to_parse = pdb_old_version.get_entry(pdb_id)
    except urllib.error.HTTPError as e:  # se the first pdb_id in the list of ids making a cluster
        print(e)
        print("No annotation for %s" % pdb_id)
        to_parse = None
    if to_parse is not None:
        to_annotate_list =parsers.parse_pdb(to_parse)
        for ts in to_annotate_list :
            try:
                ss = None
                if annotate:
                    ss, ts = rnaview.annotate(ts, canonical_only=canonical_only)
                save(db, ss, ts, pdb_id, limit)
            except Exception as e:
                print(e)
                print("No annotation for %s" % pdb_id)
                save(db, None, ts, pdb_id, limit)
    return to_parse

def import_3Ds(db_host = 'localhost', db_port = 27017, rna3dhub = False, canonical_only = True, annotate = False, limit = 5000):
    client = MongoClient(db_host, db_port)
    db_name = ""

    if rna3dhub:
        db_name = "RNA3DHub"
    else:
        rna3dHub = None
        db_name = "PDB"

    db = client[db_name]
    rnaview = Rnaview()

    if not rna3dhub:
        pdb = Client()
        pdb_ids = pdb.query_all_RNAs()
        print("%i 3Ds to process"%len(pdb_ids))
        pdb_old_version = db_PDB()
        n = len(pdb_ids)
        for i,pdb_id in enumerate(pdb_ids):
            print("Running versus {:s} ({:d}/{:d})".format(pdb_id,i+1,n))
            db_element = db['tertiaryStructures'].find_one({'source':"db:pdb:%s"%pdb_id})
            if db_element:
                continue
            parsed = parsers.parse_pdb(pdb_old_version.get_entry(pdb_id))
            for j,ts in enumerate(parsed):
                try:
                    ss = None
                    if annotate:
                        ss, ts = rnaview.annotate(ts, canonical_only = canonical_only)
                    save(db, ss, ts, pdb_id, limit)
                    print("\tRecover {:s} / {:d}".format(pdb_id,j))
                except Exception as e:
                    print(e)
                    print ("\tNo annotation for %s"%pdb_id)
                    save(db, None, ts, pdb_id, limit)
    else:
        rna3dHub = db_RNA3D()
        clusters = rna3dHub.get_clusters()
        print ("%i 3Ds to process"%len(clusters))
        pdb_old_version = db_PDB()
        for cluster in clusters['pdb-ids']:
            if len(cluster) == 0:
                continue
            all_ids = sorted(set([cluster_to_split.split('|')[0]
                                  for cluster_to_split in cluster]))
            for pdb_id in all_ids:
                only_None_if_failed = \
                    single_3ds_cluster_or_None(db, pdb_old_version, pdb_id, rnaview)
                if only_None_if_failed is not None:
                    # found an example for this cluster
                    break

def save(db, secondary_structure, tertiary_structure, pdbId, limit):
    if db['junctions'].count() >= limit:
        print ("Limit of %i junctions reached"%limit)
        sys.exit()

    tertiary_structure.source="db:pdb:%s"%pdbId

    if secondary_structure:

        computation = {
            'inputs': [tertiary_structure._id+"@tertiaryStructures"],
            'outputs': [secondary_structure._id+"@secondaryStructures"],
            'tool': "tool:rnaview:N.A.",
            'date': str(datetime.datetime.now())
        }

        if secondary_structure.rna == tertiary_structure.rna:
            ncRNA = {
                '_id': secondary_structure.rna._id,
                'source': secondary_structure.rna.source,
                'name': secondary_structure.rna.name,
                'sequence': secondary_structure.rna.sequence,
            }
            if not db['ncRNAs'].find_one({'_id':ncRNA['_id']}):
                db['ncRNAs'].insert_one(ncRNA)
        else:
            ncRNA = {
                '_id': secondary_structure.rna._id,
                'source': secondary_structure.rna.source,
                'name': secondary_structure.rna.name,
                'sequence': secondary_structure.rna.sequence,
            }
            if not db['ncRNAs'].find_one({'_id':ncRNA['_id']}):
                db['ncRNAs'].insert_one(ncRNA)
            ncRNA = {
                '_id': tertiary_structure.rna._id,
                'source': tertiary_structure.rna.source,
                'name': tertiary_structure.rna.name,
                'sequence': tertiary_structure.rna.sequence,
            }
            if not db['ncRNAs'].find_one({'_id':ncRNA['_id']}):
                db['ncRNAs'].insert_one(ncRNA)

        secondary_structure.find_junctions()

        ss_descr = {
            '_id': secondary_structure._id,
            'source': secondary_structure.source,
            'name': secondary_structure.name,
            'rna': secondary_structure.rna._id+"@ncRNAs"
        }

        helices_descr = []
        for helix in secondary_structure.helices:
            helix_desc = {
                'name': helix['name'],
                'location': helix['location']
            }
            if "interactions" in helix:
                interactions_descr = []
                for interaction in helix['interactions']:
                    interactions_descr.append({
                        'orientation': interaction['orientation'],
                        'edge1': interaction['edge1'],
                        'edge2': interaction['edge2'],
                        'location': interaction['location']
                    })
                helix_desc['interactions'] = interactions_descr

            helices_descr.append(helix_desc)

        ss_descr['helices'] = helices_descr

        single_strands_descr = []
        for single_strand in secondary_structure.single_strands:
            single_strands_descr.append({
                'name': single_strand['name'],
                'location': single_strand['location']
            })

        ss_descr['singleStrands'] = single_strands_descr

        tertiary_interactions_descr = []
        for tertiary_interaction in secondary_structure.tertiary_interactions:
            tertiary_interactions_descr.append({
                'orientation': tertiary_interaction['orientation'],
                'edge1': tertiary_interaction['edge1'],
                'edge2': tertiary_interaction['edge2'],
                'location': tertiary_interaction['location']
            })

        ss_descr['tertiaryInteractions'] = tertiary_interactions_descr

        db['secondaryStructures'].insert_one(ss_descr)

    ncRNA = {
        '_id': tertiary_structure.rna._id,
        'source': tertiary_structure.rna.source,
        'name': tertiary_structure.rna.name,
        'sequence': tertiary_structure.rna.sequence,
    }
    if not db['ncRNAs'].find_one({'_id':ncRNA['_id']}):
        db['ncRNAs'].insert_one(ncRNA)

    ts_descr = {
        '_id': tertiary_structure._id,
        'source': tertiary_structure.source,
        'name': tertiary_structure.name,
        'rna': tertiary_structure.rna._id+"@ncRNAs",
        'numbering-system': tertiary_structure.numbering_system
    }

    residues_descr = {}
    keys=[]
    for k in tertiary_structure.residues:
        keys.append(k)

    keys.sort() #the absolute position are sorted

    for key in keys:
        atoms = tertiary_structure.residues[key]['atoms']

        atoms_descr = []

        for atom in atoms:
            atoms_descr.append({
                'name': atom['name'],
                'coords': atom['coords']
            })
        residues_descr[str(key)] = {
            'atoms': atoms_descr
        }

    ts_descr['residues'] = residues_descr

    if not db['tertiaryStructures'].find_one({'_id':ts_descr['_id']}):
        db['tertiaryStructures'].insert_one(ts_descr)

        if secondary_structure:

            for junction in secondary_structure.junctions:
                junction_descr = {
                    '_id': str(ObjectId()),
                    'molecule': secondary_structure.rna._id+"@ncRNAs",
                    'tertiary-structure': {
                        'id':tertiary_structure._id+'@tertiaryStructures',
                        'source': tertiary_structure.source
                    },
                    'description': junction['description'],
                    'location': junction['location']
                }
                computation['outputs'].append(junction_descr['_id']+"@junctions")

                db['junctions'].insert_one(junction_descr)

            db['computations'].insert_one(computation)

if __name__ == '__main__':
    db_host = 'localhost'
    db_port = 27017
    rna3dhub = False
    canonical_only = False
    annotate = False
    limit = 5000

    if "-h" in sys.argv:
        print ("Usage: ./import_3Ds.py [-p x] [-mh x] [-mp x] [-l x] [-rna3dhub] [-canonical_only] [-annotate]")
        print ('- mh: the mongodb host (default: localhost)\n')
        print ('- mp: the mongodb port (default: 27017)\n')
        print ('- l: limit of junctions to be stored (default: 5000)\n')
        print ('- rna3dhub: use the 3D structures from the non-redundant set\n')
        print ('- canonical_only: a secondary structure is made with canonical base-pairs only')
        print ('- annotate: annotate each 3D structure imported')
        sys.exit(-1)

    if "-mh" in sys.argv:
        db_host = sys.argv[sys.argv.index("-mh")+1]
    if "-mp" in sys.argv:
        db_port = int(sys.argv[sys.argv.index("-mp")+1])
    if "-l" in sys.argv:
        limit = int(sys.argv[sys.argv.index("-l")+1])
    rna3dhub =  "-rna3dhub" in sys.argv
    canonical_only =  "-canonical_only" in sys.argv
    annotate = "-annotate" in sys.argv

    import_3Ds(db_host = db_host, db_port = db_port, rna3dhub = rna3dhub, canonical_only = canonical_only, annotate = annotate, limit = limit)
