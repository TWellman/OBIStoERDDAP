OBIS-USA Source Data erddap
===========================
This repository contains files required to place the OBIS-USA source data CSVs stored in ScienceBase behind an erddap
instance.  It will include any required custom ERDDAP configuration, as well as the harvester that reads source
data from ScienceBase and stages it for use by the ERDDAP server.

Running ERDDAP
--------------
The production server can be found at https://erddap-obis.digitalcrust.org/erddap.

For development, we are using the ERDDAP Docker container from Axiom, found at https://github.com/axiom-data-science/docker-erddap.

To run the container:

````
docker run \
    -d \
    -p 80:8080 \
    -p 443:8443 \
    -v /path/to/your/erddap/content:/opt/tomcat/content/erddap \
    -v /path/to/your/erddap/bigParentDirectory:/erddapData \
    --name erddap \
    axiom/docker-erddap
````

To get a command shell on the container once it is running:

````
docker exec -it erddap bash
````

ERDDAP will be available at http://localhost/erddap/index.html

Some configuration borrowed from the ERDDAP running at http://54.237.216.54:8080/erddap/index.html

Running obis.py
---------------
The obis.py script downloads CSVs and ZIP files marked as "Final Processed Source" from the OBIS collection in
ScienceBase, creates netCDF v3 files from them, and then creates a datasets.xml to configure ERDDAP to serve
them.  The obis.py script has help which can be viewed with the --help option.

Provisional Software Disclaimer
---------------
Under USGS Software Release Policy, the software codes here are considered preliminary, not released officially, and posted to this repo for informal sharing among colleagues.

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.