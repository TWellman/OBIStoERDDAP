docker run -i -t -p 80:8080 -p 443:8443 -v /Users/jllong/src/obisSourceERDDAP/content/erddap:/opt/tomcat/content/erddap -v /Users/jllong/etc/erddap/erddap_data:/erddapData --name erddap axiom/docker-erddap

