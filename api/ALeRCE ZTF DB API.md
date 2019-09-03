# ALeRCE ZTF DB API
This API gives access to ALeRCE Annotated ZTF API, with corrected magnitudes, objects statistics and the object data.

The access point to the API is [ztf.alerce.online], querying the API is done with **http POST requests**  with parameters in **JSON** format.

## Query the DB
### ztf.alerce.online/query
This will query the DB objects and get the stats for the matching objects.

The current fields to query the db are the following:

```
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
        sr: float degrees
    },
    dates:{
        //First detection (Discovery date)
        firstmjd: {
         min: float mjd,
         max: float mjd
        }
    }
}
```

The response contains the following:
```
{
    "total": int,
    "num_pages": int,
    "page": int,
    "result": {
        <ObjectId>: <ObjectStats>
    }
    
}
```
Where ObjectId is the ZTF id and ObjectStats is a JSON document with the stats.

#### Example
Get the recent objects
```
curl -X POST ztf.alerce.online/query -d @- << EOF
{
   "query_parameters":{
       "dates":{
           "firstmjd":{
               "min": 58682
           }
       }
   } 
}
EOF
```

Get the last 100 objects classified as SNe
```
curl -X POST ztf.alerce.online/query -d @- << EOF
{
    "records_per_pages":100,
    "query_parameters": {
      "filters": {"classearly": 2},
      "dates":{
        "firstmjd":{
            "min":58680
        }
      }
    },
    "sortBy": "pclassearly",
    "total":100
}
EOF
```



### ztf.alerce.online/get_sql
Get the SQL Query done to the ZTF DB using the same parameters as **/query**

### Example
```
curl -X POST ztf.alerce.online/get_sql -d @- << EOF
{
   "query_parameters":{
       "dates":{
           "firstmjd":{
               "min": 58682
           }
       }
   } 
}
EOF
```


## Query an object
To get an specific ZTF Object information.
### ztf.alerce.online/get_detection
Get all the detections for an object.

```
curl -X POST ztf.alerce.online/get_detections -d @- << EOF
{
   "oid":"ZTF18abbvavt"
}
EOF
```



### ztf.alerce.online/get_non_detection
Get non detections of an object.

```
curl -X POST ztf.alerce.online/get_non_detections -d @- << EOF
{
   "oid":"ZTF18abbvavt"
}
EOF
```

### ztf.alerce.online/get_stats
Get stats for an object.
```
curl -X POST ztf.alerce.online/get_stats -d @- << EOF
{
   "oid":"ZTF18abbvavt"
}
EOF
```

### ztf.alerce.online/get_probabilities
Get probabilities of the models for an object.
```
curl -X POST ztf.alerce.online/get_probabilities -d @- << EOF
{
   "oid":"ZTF18abbvavt"
}
EOF
```

### ztf.alerce.online/get_features
Get features computed for an object (list of features pending).
```
curl -X POST ztf.alerce.online/get_features -d @- << EOF
{
   "oid":"ZTF18abbvavt"
}
EOF
```
