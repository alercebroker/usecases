# ALeRCE AVRO/Stamps Access API

API to get an specified avro/stamp from alert files of an ZTF Object. This API uses **GET requests** parameters.

## avro.alerce.online/get_avro
Get an specified avro
```
    curl http://avro.alerce.online/get_stamp?oid=ZTF19abguqsi&candid=928474644915015004
```

Parameters:
- oid: ZTF Object ID.
- candid: ZTF Alert ID.

## avro.alerce.online/get_stamp
Get an specific stamp
```
 curl http://avro.alerce.online/get_stamp?oid=ZTF19abguqsi&candid=928474644915015004&format=png&type=science
```
Parameters:
- oid: ZTF Object ID.
- candid: ZTF Alert ID.
- format: ["png","fits"]
- type: ["science","template","difference"]
