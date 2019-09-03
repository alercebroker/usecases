import requests
import pandas as pd
from pandas.io.json import json_normalize
from IPython.display import HTML
from astropy.table import Table, Column
#from xml.etree import ElementTree
#import numpy as np
#from io import BytesIO
#from PIL import Image
#import base64

class AlerceParseError(Exception):
    pass

class AlerceOidError(Exception):
    pass

class AlerceOidError(Exception):
    pass

class AlerceAPI(object):
    'ALeRCE API python wrapper class'
    
    # init
    def __init__(self, **kwargs):
                
        self.ztf_url = "http://ztf.alerce.online"
        if "ztf_url" in kwargs.keys():
            self.ztf_url = kwargs["ztf_url"]
            
        self.catsHTM_url = "http://catshtm.alerce.online"
        if "catsHTM_url" in kwargs.keys():
            self.catsHTM_url = kwargs["catsHTM_url"]
            
        self.oid = ""


    def query(self, params):
        """Query the ALeRCE API to get matching objects into a pandas dataframe.

        The current fields to query the db are the following:

        {
        total: number (if not set the total is counted and the query is slower),
        records_per_pages: number (default 20),
        page: number (default 1),
        sortBy: string columnName (default nobs),
        query_parameters:{
           filters:{
              //ZTF Object id
              oid: "ZTFXXXXXX",
              //Number of detections
              nobs: { 
                  min: int
                  max: int
              },
              //Late Classifier (Random Forest)
              classrf: ["CEPH","DSCT","EB","LPV","RRL","SNe","Other"] or int,
              pclassrf: float [0-1],
              //Early Classifier (Stamp Classifier)
              classearly: ["AGN","SN","VS","asteroid","bogus"] or int,
              pclassearly: float [0-1],
              },
              //Coordinate based search (RA,DEC) and Search Radius.
              coordinates:{
                ra: float degrees,
                dec: float degrees,
                sr: float degrese
                },
              dates:{
              //First detection (Discovery date)
              firstmjd: {
                 min: float mjd,
                 max: float mjd
                }
              }
           }
        }
        The response contains the following columns in a pandas dataframe:
        
        {
        "total": int,
        "num_pages": int,
        "page": int,
        "result": {
           <ObjectId>: <ObjectStats>
           }
        }

"""
        
        # show api results
        r = requests.post(url = "%s/query" % self.ztf_url, json = params) 
        df = pd.DataFrame(r.json())
        query_results = json_normalize(df.result)
        query_results.set_index('oid', inplace=True)
        return query_results

    def get_sql(self, params):
        'get sql from given json query parameters'
        
        r = requests.post(url = "%s/get_sql" % self.ztf_url, json = params)
        return r.content

    def get_detections(self, oid):
        'get detections given oid as pandas dataframe'
        
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_detections" % self.ztf_url, json = params) 
        df = pd.DataFrame(r.json())
        detections = json_normalize(df.result.detections)
        detections.sort_values(by=['mjd'], inplace=True)
        detections.set_index('candid', inplace=True)
        return detections

    
    def get_non_detections(self, oid):
        '# get non detections given oid as pandas  dataframe'
        
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_non_detections" % self.ztf_url, json = params)
        try:
            df = pd.DataFrame(r.json())
        except AlerceParseError:
            print("ERROR (get_non_detections): could not convert API output to dataframe")
        non_detections = json_normalize(df.result.non_detections)
        if non_detections.shape[0] > 0:
            non_detections.sort_values(by=['mjd'], inplace=True)
            non_detections.set_index('mjd', inplace=True)
        return non_detections

    def get_stats(self, oid, verbose=False):
        'get stats given oid as pandas dataframe'
        
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_stats" % self.ztf_url, json = params)

        try:
            df = pd.DataFrame(r.json())
        except AlerceParseError:
            print("ERROR (get_non_detections): could not convert API output to dataframe")
        
        stats = json_normalize(df.result.stats)
        stats.set_index('oid', inplace=True)
        
        # save some features
        if verbose:
            print("Setting ra, dec, firstMJD")
        self.oid = oid
        self.meanra = stats.meanra
        self.meandec = stats.meandec
        self.firstMJD = stats.firstmjd
        
        return stats

    def get_probabilities(self, oid, doearly=True, dolate=True):
        'get probabilities given oid as pandas dataframe (late or early)'
        
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_probabilities" % self.ztf_url, json = params)

        if doearly:
            try:
                early = json_normalize(r.json()["result"]["probabilities"]["early_classifier"])
                early.set_index("oid", inplace=True)
            except AlerceParseError:
                print("ERROR (get_probabilities): could not convert early probabilities API output to dataframe")
            
        if dolate:
            try:
                late = json_normalize(r.json()["result"]["probabilities"]["random_forest"])
                late.set_index("oid", inplace=True)
                return early, late
            except AlerceParseError:
                print("ERROR (get_probabilities): could not convert early probabilities API output to dataframe")

        result = {}
        if doearly:
            result["early_probabilities"] = early
        if dolate:
            result["late_probabilities"] = early

        return result

    def get_features(self, oid):
        'get features given oid as pandas dataframe'
    
        #oid
        params = {
            "oid": oid
        }
        
        # show api results
        r = requests.post(url = "%s/get_features" % self.ztf_url, json = params) 
        features = json_normalize(r.json())
        features.set_index('oid', inplace=True)
        return features
    
    def catsHTM_conesearch(self, oid, catalog_name, radius):
        'get catsHTM crossmatch given oid, catalog_name (all is also allowed) and search radius in arcsec'
        
        # get ra, dec
        if oid != self.oid:
            self.get_stats(oid)
        
        params = {
            "catalog": "%s" % catalog_name,
            "ra": "%f" % self.meanra,
            "dec": "%f" % self.meandec,
            "radius": "%f" % radius
        }
        if catalog_name != "all":
            result = requests.get(url = "%s/conesearch" % self.catsHTM_url, params = params)
        else:
            result = requests.get(url = "%s/conesearch_all" % self.catsHTM_url, params = params)
        
        votables = {}

        try:
            if "catalogs" not in result.json().keys():
                return
        except:
            return
        
        for idx, r in enumerate(result.json()["catalogs"]):

            key = list(r.keys())[0]
            if r[key] == {}:
                continue
            t = Table()
            for field in r[key].keys():
                data = r[key][field]['values'] #list(map(lambda x: x["value"], r[key][field]))
                t.add_column(Column(data, name=field))
                t[field].unit = r[key][field]['units']
            t["cat_name"] = Column(["catsHTM_%s" % key], name="cat_name")
            votables[key] = t
        
        return votables
    

    def catsHTM_redshift(self, oid, radius):
        'get redshift given oid using catsHTM crossmatch'

        # get ra, dec
        if oid != self.oid:
            self.get_stats(oid)
                
        # search data in catsHTM
        xmatches = self.catsHTM(oid, "all", radius)
        
        # check whether redshift exists
        for xmatch in xmatches:
            for catname in xmatch.keys():
                for col in xmatch[catname].keys():
                    if col == 'z':
                        print("Found redshift in catalog %s: %f" % (catname, xmatch[catname][col]["value"]))
                        return float(xmatch[catname][col]["value"])

        # check whether photometric redshift exists
        for xmatch in xmatches:
            for catname in xmatch.keys():
                for col in xmatch[catname].keys():
                    if col in ['zphot',  'z_phot']:
                        print("Found photometric redshift in catalog %s: %f" % (catname, xmatch[catname][col]["value"]))
                        return float(xmatch[catname][col]["value"])
                    
        return

    def plot_stamp(self, oid, candid=None):
        'plot stamp in a notebook given oid. It uses IPython HTML.'
        
        # if candid is None, get minimum candid
        if candid is None:
            candid = min(self.get_detections(oid).index)
            
        science = "http://avro.alerce.online/get_stamp?oid=%s&candid=%s&type=science&format=png" % (oid, candid)
        images="""
        <div>ZTF oid: %s, candid: %s</div>
        <div>&emsp;&emsp;&emsp;&emsp;&emsp;
        Science
        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 
        Template
        &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; 
        Difference
        <div class="container">
        <div style="float:left;width:20%%"><img src="%s"></div>
        <div style="float:left;width:20%%"><img src="%s"></div>
        <div style="float:left;width:20%%"><img src="%s"></div>
        </div>
        """ % (oid, candid, science, science.replace("science", "template"), science.replace("science", "difference"))
        display(HTML(images))
