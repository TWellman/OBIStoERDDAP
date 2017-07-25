#!/usr/local/bin/python3
"""
  This script downloads final source CSV files from the OBIS-USA collection, converts them to
  netCDF and stages them for consumption by ERDDAP.  It also generates the datasets.xml ERDDAP
  configuration file to serve those netCDF files as data sets.

  CSVs can be served from ERDDAP.  However, from the ERDDAP EDDTableFromAsciiFiles documentation:
      ASCII files are easy to work with, but they are not an efficient way to store/retreive data.
      For greater efficiency, save the files as NetCDF v3 .nc files (with one dimension, "row",
      shared by all variables) instead.
"""
from pysb import SbSession
from datetime import datetime, timedelta, timezone
from xml.dom import minidom
import xml.etree.ElementTree  as ET
import numpy
import netCDF4
import messytables
import csv
import os
import sys
import glob
import traceback
import decimal
import re
import getopt

opt_fields = [
    'collectionid=', 'sourcedir=', 'erddapdir=', 'serverdir=', 'only_csv', 'only_netcdf',
    'only_datasets_xml', 'virtual_datasets', 'window=', 'rowsperfile=', 'overwrite',
    'sample=', 'verbose', 'help']
def usage():
    print("""
    %s [options]

Options
-------------------
    --collectionid=<id> The ScienceBase Item ID of the collection to process.
    --sourcedir=<path> Directory in which to place source CSV fetch_csvs (defaults to ./source_data).
    --erddapdir=<path> Directory in which to place netCDF files (defaults to ./erddap_data/nc).
    --server_dir=<path> Directory in which netCDF files will be placed on the ERDDAP server.
      Defaults to /etddapData/nc.
    --only_csv Download CSVs from ScienceBase (no netCDF or datasets.xml creation).
    --only_netcdf Create netCDF files from existing CSVs.
    --only_datasets_xml Create datasets.xml from existing netCDF files.
    --virtual_datasets Create a virtual dataset for all netCDF files created from a single CSV.
    --window=<n> Use the first n rows to guess CSV column datatypes (defaults to 500).
    --rowsperfile=<n> Split the CSV data into netCDF files containing n rows each.
    --overwrite Overwrite existing files.  Default is to skip files that exist locally.
    --sample=<n> Only process the first n rows of the CSV (for testing purposes).
    --verbose Additional logging.
""" % (sys.argv[0]))

#
# ScienceBase Item ID of the OBIS-USA source data collection
#
collection_id = '579b64c6e4b0589fa1c98118'

#
# Directory in which to place downloaded source CSVs from ScienceBase
#
source_data_dir = './source_data'

#
# Directory in which to place converted netCDF files for use in ERDDAP
#
erddap_data_dir = './erddap_data/nc'

#
# Whether to fetch metadata from ScienceBase.  If this is True and fetch_csvs
# is false, it will only fetch metadata.
#
fetch_metadata = True

#
# Whether to fetch source files from ScienceBase.
#
fetch_csvs = True

#
# Whether to create netCDF files from the sources
#
create_netcdf_files = True

#
# Whether to create a datasets.xml from the netCDF files found in the erddap_data_dir
#
create_datasets_xml = True

#
# Whether to create a single virtual dataset from multiple netCDFs that were created from a single csv
#
create_virtual_datasets = False

#
# Whether to turn on verbose logging
#
verbose = False

#
# Size in number of rows of the sample.  Only takes effect if sample is greater than zero.
#
sample_size = 0

#
# Number of rows to use to guess column type
#
window = 500

#
# Max number of rows per netCDF file.  If greater than zero, multiple files will be created if necessary.
# Zero Value signals to create a single netCDF file for the dataset.
#
rows_per_file = 0

#
# Whether to overwrite an existing CSV or netCDF.  If set to False, processing of the input file will be skipped.
#
overwrite = False

#
# Path to the netCDF files on the destination server, for the datasets.xml file.
#
server_nc_directory = '/erddapData/nc'

#
# retrieve_source_files
#    sb:  ScienceBase pysb session
#    collection_id:  ScienceBase Item ID of the OBIS-USA source data collection
#    source_data_dir:  Directory in which to place downloaded source data
#    download:  Flag whether to download the file.  Download if True, otherwise only retrieve metadata
# Finds all CSVs marked as "Final Processed Source" in the given collection, retrieves metadata
# for them, and, optionally, downloads the file.
#
# Returns: Dict containing ScienceBase metadata keyed by source file name.
#
def retrieve_source_files(sb, collection_id, source_data_dir, download):
    global overwrite

    source_files = {}
    original_sources_id = collection_id
    results = sb.find_items({
        'ancestors': original_sources_id,
        'fields': 'title, body, files',
        'max': 1000})
    while results and 'items' in results:
        for item in results['items']:
            if 'files' in item:
                for item_file in item['files']:
                    if 'title' in item_file and 'final processed source' in item_file['title'].lower():
                        print('%s(%s): %s' % (item['title'], item['id'], item_file['name']))
                        if download:
                            if not overwrite and os.path.isfile(os.path.join(source_data_dir, item_file['name'])):
                                print('Skipping %s' % (item_file['name']))
                            else:
                                sb.download_file(item_file['url'], item_file['name'], source_data_dir)
                        source_files[item_file['name']] = item
        results = sb.next(results)
    return source_files

#
# csv_to_nc
#     file_name: CSV (or ZIP containing a CSV) file to parse
#     source_data_dir: Directory containing CSV source files
#     metadata: ScienceBase item metadata for the file
#     erddap_data_dir: Directory in which to write the output netCDF file
#     rows_per_file
#
# Convert the given CSV to a netCDF with a "row" dimension.
#
def csv_to_nc(file_name, source_data_dir, metadata, erddap_data_dir, rows_per_file):
    global overwrite

    base_name, extension = os.path.splitext(file_name)
    full_file_name = os.path.join(source_data_dir, file_name)
    if 'csv' in extension.lower() or 'zip' in extension.lower():
        fpath = os.path.join(erddap_data_dir, base_name + '*.nc')
        file_exists = len(glob.glob(os.path.join(erddap_data_dir, base_name + '*.nc'))) > 0

        if not overwrite and file_exists:
            print('Not writing netCDF for %s, already exists' % (full_file_name))
            return
        (row_set, types, fh) = process_csv(full_file_name, extension)
        try:
            (variables, variable_data, num_rows) = normalize_data_for_netcdf(row_set, types)
            offset = 0
            if rows_per_file > 0:
                page_size = rows_per_file
            else:
                page_size = num_rows
            while offset < num_rows:
                create_netcdf(erddap_data_dir, base_name, metadata, variables, variable_data, offset, page_size, num_rows)
                offset += page_size
        except:
            raise
        finally:
            if fh:
                fh.close()
    else:
        print('Unknown file %s type (%s), skipping' % (full_file_name, extension))

#
# process_csv
#     full_file_name: Path and filename of the csv
#     extension: File extension.  Can be "zip"" or "csv."
#
# Process the CSV, guessing column data types.
#
# Returns: The proceseed row_set, a dictionary of column types keyed by column name, and the
# open file handle (use this to close the file once the row_set is processed).
#
def process_csv(full_file_name, extension):
    global window
    row_set = []
    types = []
    fh = None
    print('Processing %s' % (full_file_name))
    try:
        fh = open(full_file_name, 'rb')
        # Parse the CSV into a table set.  CSVs only have one table, so grab the first.
        if 'csv' in extension.lower():
            table_set = messytables.CSVTableSet(fh, window=window)
        elif 'zip' in extension.lower():
            table_set =  messytables.zip.ZIPTableSet(fh, window=window)

        row_set = table_set.tables[0]

        # Grab the header information, and set the iterator past it.
        offset, headers = messytables.headers_guess(row_set.sample)
        row_set.register_processor(messytables.headers_processor(headers))
        row_set.register_processor(messytables.offset_processor(offset + 1))

        # Guess column types and tell the row set to apply these types
        types = messytables.type_guess(row_set.sample, strict=True)
        row_set.register_processor(messytables.types_processor(types, strict=False))        
    except:        
        traceback.print_exc(file=sys.stdout)
        print('An error occurred, skipping %s' % (full_file_name))
        if fh:
            fh.close()
        row_set = types = []
        fh = None
    return (row_set, types, fh)
    

#
# normalize_data_for_netcdf
#     row_set: Messytables Rowset to process_csv
#     types: Dictionary of column types keyed by column NameError
#
# Processes the Messytables Rowset, normalizing the values in each cell for the column data type.
#
# Returns: List of variables containing interesting data, a dictionary of variable data keyed bytes
# column type, and the number of rows in the CSV.
#
def normalize_data_for_netcdf(row_set, types):
    global sample_size
    # Iterate the data and store for pre-processing
    variables = []
    variable_data = {}
    num_rows = 0

    for row in row_set:
        num_rows += 1
        for column_index, cell in enumerate(row):
            if cell.column not in variable_data:
                # Special case:  ERDDAP requires that altitude and depth both be numeric
                if cell.column.lower() in ['altitude', 'depth']:
                    types[column_index] = messytables.types.DecimalType()
                variables.append((cell.column, types[column_index]))
                variable_data[cell.column] = []
            variable_data[cell.column].append(nc_value(cell, types[column_index]))
        if sample_size > 0 and num_rows == sample_size:
            break

    interesting_variables = []
    for variable, variable_type in variables:
        if has_interesting_values(variable_type, variable_data[variable]):
            if isinstance(variable_type, messytables.types.StringType):
                strlen = len(max(variable_data[variable], key=len))
                new_data = []
                for val in variable_data[variable]:
                    # Strip out non-ascii characters
                    ve = val.encode('ascii', errors='ignore').decode('ascii')
                    # Convert string to fixed length character array
                    new_data.append(netCDF4.stringtoarr(ve, strlen))
                variable_data[variable] = new_data
            interesting_variables.append((variable, variable_type))

    return (interesting_variables, variable_data, num_rows)

#
# create_netcdf
#     erddap_data_dir: Directory in which to write the netCDF files
#     base_name: File base name, used as the dataset name
#     metadata: ScienceBase metadata for the dataset
#     variables: List of tuples containing column names and types
#     variable_data: Dictionary of processed variable data, keyed by column name
#     offset: Index of starting row
#     page_size: Number of rows to write to this netCDF file
#     num_rows: Total number of rows in the data
#
# Creates a netCDF file from the processed CSV data.  Accepts an offset and page_size to allow multiple
# files to be written for each dataset.
#
def create_netcdf(erddap_data_dir, base_name, metadata, variables, variable_data, offset, page_size, num_rows):
    global verbose

    if num_rows > page_size:
        base_name = "%s_part%d" % (base_name, offset / page_size + 1)
    nc_fname = os.path.join(erddap_data_dir, base_name + ".nc")
    end_row = min(offset + page_size, num_rows)
    print("Creating %s for rows %d to %d" % (nc_fname, offset + 1, end_row))

    try:
        rootgrp = netCDF4.Dataset(nc_fname, "w", format="NETCDF3_CLASSIC")
        rootgrp.set_fill_on()

        # Set metadata attributes
        if metadata:
            attrs = {}
            if 'title' in metadata:
                attrs['title'] = metadata['title']
            if 'body' in metadata:
                attrs['summary'] = metadata['body']
            if 'link' in metadata:
                attrs['infourl'] = metadata['link']['url']
            rootgrp.setncatts(attrs)

        # Create netCDF dimensions for number of rows, and the size of each string field
        row = rootgrp.createDimension('row', end_row - offset)
        i = 0
        dims = {}
        uninteresting_vars = []

        for variable, variable_type in variables:
            if isinstance(variable_type, messytables.types.StringType):
                # Create a dimension for the character array length
                i += 1
                dim_name = "STRING%d" % (i)
                dims[variable] = dim_name
                # Use max length of all data so sizes stay constant across multiple netCDF files
                rootgrp.createDimension(dim_name, max([len(max(variable_data[variable], key=len)), 1]))

        # Create netCDF variables for each column.  We have to do this in a second pass, because all dimensions must
        # exist before we create variables.
        for variable, variable_type in variables:
            variable_data_slice = variable_data[variable][offset:end_row]
            if verbose:
                print('variable(%s)(%s)' % (variable, variable_type))
            if isinstance(variable_type, messytables.types.StringType):
                nc_var = rootgrp.createVariable(variable,nc_type(variable_type),('row',dims[variable],))
                nc_var[:] = numpy.asarray(variable_data_slice)
            else:
                nc_var = rootgrp.createVariable(variable,nc_type(variable_type),('row',))
                nc_var[:] = variable_data_slice
        rootgrp.close()
    except:
        print("Error creating netCDF %s, removing" % (nc_fname))
        traceback.print_exc(file=sys.stdout)
        os.remove(nc_fname)
#
# has_interesting_values
#     variable_type: Variable type
#     variable_data: Variable data
#
# Return whether the variable data contains data besides fill vaules
#
def has_interesting_values(variable_type, variable_data):
    ret = False
    if not isinstance(variable_type, messytables.types.BoolType):
        for var_value in variable_data:
            if var_value is not netCDF4.default_fillvals[nc_type(variable_type)]:
                ret = True
                break
    return ret

#
# nc_value
#     cell:  Messytables cell
#     column_type:  Type of the column the cell is in
#
# Return the given value in a format suitable for netCDF4
# TODO: Add nc attribute to netCDF for date mask.
#
def nc_value(cell, column_type):
    ret = netCDF4.default_fillvals[nc_type(column_type)]
    try:
        if isinstance(column_type, messytables.types.DateType):
            ret = cell.value.date().isoformat()
        elif isinstance(column_type, messytables.types.StringType):
            if cell.value not in ['NA']:
                ret = str(cell.value).strip()
            else:
                ret = ''
        elif isinstance(column_type, messytables.types.IntegerType):
            ret = int(float(cell.value))
        elif isinstance(column_type, messytables.types.DecimalType):
            ret = float(cell.value)
        elif isinstance(column_type, messytables.types.BoolType):
            ret = str(cell.value)
        elif cell.value:
            ret = str(cell.value)
    except:
        if verbose and cell.value is not None:
            print('filling %s(%s) column with %s for %s' % (cell.column, str(column_type), str(ret), str(cell.value)))
    finally:
        return ret

#
# nc_type
#     orig_type:  Type from messytables to convert to netCDF type
#
# Returns the netCDF type corresponding to the given messytables type
#
def nc_type(orig_type):
    ret = 'S1'
    if isinstance(orig_type, messytables.types.IntegerType):
        ret = 'i4'
    elif isinstance(orig_type, messytables.types.DecimalType):
        ret = 'f8'
    elif isinstance(orig_type, messytables.types.DateType):
        ret = 'S1'
    elif isinstance(orig_type, messytables.types.StringType):
        ret = 'S1'
    elif isinstance(orig_type, messytables.types.BoolType):
        ret = 'S1'
    return ret

#
# print_object
#     o: Object to print
# Prints the input object as comma-separated values (useful for debugging)
#
def print_object(o):
    print(', '.join("%s: %s" % item for item in vars(o).items()))

#
# write_datasets_xml
#     erddap_data_dir: Data containing the netCDF files
#     config_dir: Directory to which to write the datasets.xml file
#     create_virtual_datasets: Whether to combine multiple netCDF files for the
#         same dataset into one virtual ERDDAP dataset.
#
# Write the datasets.xml file for the given netCDF files.
#
def write_datasets_xml(erddap_data_dir, config_dir, create_virtual_datasets):
    datasets_xml_root = create_datasets_xml_root()
    processed = []
    r = re.compile(r'^(.*)_part\d+$')
    for file_name in glob.glob(os.path.join(erddap_data_dir, '*.nc')):
        base_name, extension = os.path.splitext(os.path.basename(file_name))
        m = re.match(r, base_name)
        if create_virtual_datasets and m:
            dataset_name = m.group(1)
        else:
            dataset_name = base_name
        if dataset_name not in processed:
            print('Writing dataset %s to datasets.xml' % (dataset_name))
            add_dataset(datasets_xml_root, dataset_name, file_name, server_nc_directory, create_virtual_datasets)
            processed.append(dataset_name)

    xmlstr = minidom.parseString(ET.tostring(datasets_xml_root)).toprettyxml(indent="   ", encoding="UTF-8")
    with open(os.path.join(config_dir, "datasets.xml") , "wb") as f:
        f.write(xmlstr)

#
# create_datasets_xml_root
#
# Create the DOM and add top level elements for the datasets.xml file
#
def create_datasets_xml_root():
    root = ET.Element('erddapDatasets')

    comment = ET.Comment('Generated from OBIS sources in ScienceBase')
    root.append(comment)

    convertToPublicSourceUrl = ET.SubElement(root, 'convertToPublicSourceUrl')

    requestBlacklist = ET.SubElement(root, 'requestBlacklist')
    requestBlacklist.text = '...'

    subscriptionEmailBlacklist = ET.SubElement(root, 'subscriptionEmailBlacklist')
    subscriptionEmailBlacklist.text = '...'

    user = ET.Comment('<user username="..." password="..." roles="..." />')
    root.append(user)
    return root

#
# add_dataset
#     root: Root DOM for the datasets.xml
#     dataset_name: Name of the dataset
#     file_name: Name of the netCDF file to add to the datasets.xml
#     server_nc_directory: netCDF file path on the ERDDAP server
#     create_virtual_datasets: Whether to combine multiple netCDF files for the
#         same dataset into one virtual ERDDAP dataset.
#
# Add a dataset to the datasets.xml
#
def add_dataset(root, dataset_name, file_name, server_nc_directory, create_virtual_datasets):
    #
    # Root elements
    #
    dataset = ET.SubElement(root, 'dataset', attrib={'type': 'EDDTableFromMultidimNcFiles', 'datasetID': dataset_name, 'active': 'true'})
    reloadEveryNMinutes = ET.SubElement(dataset, 'reloadEveryNMinutes')
    reloadEveryNMinutes.text = '10080'
    updateEveryNMillis = ET.SubElement(dataset, 'updateEveryNMillis')
    updateEveryNMillis.text = '10000'
    fileDir = ET.SubElement(dataset, 'fileDir')
    fileDir.text = server_nc_directory
    fileNameRegex = ET.SubElement(dataset, 'fileNameRegex')
    if create_virtual_datasets:
        fileNameRegex.text = dataset_name + ".*\.nc"
    else:
        fileNameRegex.text = dataset_name + "\.nc"
    recursive = ET.SubElement(dataset, 'recursive')
    recursive.text = 'false'
    pathRegex = ET.SubElement(dataset, 'pathRegex')
    pathRegex.text= '.*'
    metadataFrom = ET.SubElement(dataset, 'metadataFrom')
    metadataFrom.text = 'last'
    preExtractRegex = ET.SubElement(dataset, 'preExtractRegex')
    postExtractRegex = ET.SubElement(dataset, 'postExtractRegex')
    extractRegex = ET.SubElement(dataset, 'extractRegex')
    columnNameForExtract = ET.SubElement(dataset, 'columnNameForExtract')
    removeMVRows = ET.SubElement(dataset, 'removeMVRows')
    removeMVRows.text = 'true'
    sortFilesBySourceNames = ET.SubElement(dataset, 'sortFilesBySourceNames')
    fileTableInMemory = ET.SubElement(dataset, 'fileTableInMemory')
    fileTableInMemory.text = 'false'
    accessibleViaFiles = ET.SubElement(dataset, 'accessibleViaFiles')
    accessibleViaFiles.text = 'false'

    #
    # Get global attributes from netCDF
    #
    nc_file = netCDF4.Dataset(file_name, 'r')
    nc_attrs = nc_file.ncattrs()
    #
    # Top level attributes
    #
    addAttributes = ET.SubElement(dataset, 'addAttributes')
    cdm_data_type = ET.SubElement(addAttributes, 'att', attrib={'name': 'cdm_data_type'})
    cdm_data_type.text = 'Other'
    Conventions = ET.SubElement(addAttributes, 'att', attrib={'name': 'Conventions'})
    Conventions.text = 'COARDS, CF-1.6, ACDD-1.3'
    infoUrl = ET.SubElement(addAttributes, 'att', attrib={'name': 'infoUrl'})
    if 'infourl' in nc_attrs:
        infoUrl.text = nc_file.getncattr('infourl')
    else:
        infoUrl.text = "Local source"
    institution = ET.SubElement(addAttributes, 'att', attrib={'name': 'institution'})
    institution.text = 'US Geological Survey'
    keywords = ET.SubElement(addAttributes, 'att', attrib={'name': 'keywords'})
    # Set the keywords text once the variable names are known
    dataset_license = ET.SubElement(addAttributes, 'att', attrib={'name': 'license'})
    dataset_license.text = 'This USGS product is considered to be in the U.S. public domain'
    sourceUrl = ET.SubElement(addAttributes, 'att', attrib={'name': 'sourceUrl'})
    if 'infourl' in nc_attrs:
        sourceUrl.text = nc_file.getncattr('infourl')
    else:
        sourceUrl.text = 'Local source'
    standard_name_vocabulary = ET.SubElement(addAttributes, 'att', attrib={'name': 'standard_name_vocabulary'})
    standard_name_vocabulary.text = 'CF Standard Name Table v29'
    #subsetVariables = ET.SubElement(addAttributes, 'att', attrib={'name': 'subsetVariables'})
    # Set the keywords text once the variable names are known
    summary = ET.SubElement(addAttributes, 'att', attrib={'name': 'summary'})
    if 'summary' in nc_attrs:
        summary.text = nc_file.getncattr('summary')
    else:
        summary.text = 'No summary info available'
    title = ET.SubElement(addAttributes, 'att', attrib={'name': 'title'})
    if 'title' in nc_attrs:
        title.text = nc_file.getncattr('title')
    else:
        title.text = 'Data from a local source'

    #
    # Variables
    #
    var_names = []
    altitude_and_depth = 'altitude' in nc_file.variables and 'depth' in nc_file.variables

    for nc_var in nc_file.variables:
        var_name = str(nc_var)
        dataVariable = ET.SubElement(dataset, 'dataVariable')
        sourceName = ET.SubElement(dataVariable, 'sourceName')
        sourceName.text = var_name
        destinationName = ET.SubElement(dataVariable, 'destinationName')
        # Special case: ERDDAP is strict about variables named "time""
        # Special case: ERDDAP does not allow a dataset to contain both "altitude" and "depth"
        if (var_name == 'time') or (altitude_and_depth and var_name in ['altitude', 'depth']):
            var_name = 'source_' + var_name
        destinationName.text = var_name
        dataType = ET.SubElement(dataVariable, 'dataType')
        dataType.text = erddap_datatype(nc_file.variables[nc_var].datatype)
        add_field_attrs(dataVariable, nc_var)
        var_names.append(var_name)

    # Now we can set keywords and subsetVariables
    keywords.text = ','.join(var_names)
    #subsetVariables.text = ",".join(var_names)

#
# add_field_attrs
#     dataVariable: dataVariable element in the datasets.xml DOM to which to add attributes
#     nc_var: netCDF variable for which to add attributes
#
# Add the appropriate attributes for the given variable.
#
def add_field_attrs(dataVariable, nc_var):
    addAttributes = ET.SubElement(dataVariable, 'addAttributes')

    var_name = str(nc_var)
    ioos_category = get_ioos_category(var_name)
    units = None
    standard_name = None
    if ioos_category == 'Location':
        if var_name.lower() in ['latitude', 'longitude']:
            coordinate_axis_type = var_name[:3].title()
            if coordinate_axis_type =='Lat':
                axis = 'Y'
                units = 'degrees_north'
            else:
                axis = 'X'
                units = 'degrees_east'
        elif var_name.lower() in ['altitude', 'depth']:
            units = 'm'
        standard_name = var_name.lower()

    ioos_category_att = ET.SubElement(addAttributes, 'att', attrib={'name': 'ioos_category'})
    ioos_category_att.text = ioos_category
    long_name_att = ET.SubElement(addAttributes, 'att', attrib={'name': 'long_name'})
    long_name_att.text = var_name.replace('_', ' ').title()
    if standard_name:
        standard_name_att = ET.SubElement(addAttributes, 'att', attrib={'name': 'standard_name'})
        standard_name_att.text = standard_name
    if units:
        units_att = ET.SubElement(addAttributes, 'att', attrib={'name': 'units'})
        units_att.text = units

#
# get_ioos_category
#     var_name: Name of the variable
#
# Determine ioos_category for the given variable.
#
# The current valid values in ERDDAP are Bathymetry, Biology, Bottom Character, Colored Dissolved Organic Matter,
# Contaminants, Currents, Dissolved Nutrients, Dissolved O2, Ecology, Fish Abundance, Fish Species, Heat Flux,
# Hydrology, Ice Distribution, Identifier, Location, Meteorology, Ocean Color, Optical Properties, Other, Pathogens,
# pCO2, Phytoplankton Species, Pressure, Productivity, Quality, Salinity, Sea Level, Statistics, Stream Flow,
# Surface Waves, Taxonomy, Temperature, Time, Total Suspended Matter, Unknown, Wind, Zooplankton Species,
# and Zooplankton Abundance.
#
def get_ioos_category(var_name):
    ret = 'Unknown'
    if var_name.lower() in ['altitude', 'depth', 'latitude', 'longitude']:
        ret = 'Location'
    return ret

#
# erddap_datatype
#     nc_dtype: dtype of netCDF variable
# Returns the ERDDAP mapping of the given netCDF variable
#
ERDDAP_CHAR = 'char'
ERDDAP_BYTE = 'byte'
ERDDAP_UBYTE = 'byte'
ERDDAP_SHORT = 'short'
ERDDAP_USHORT = 'short'
ERDDAP_INT = 'int'
ERDDAP_UINT = 'int'
ERDDAP_INT64 = 'long'
ERDDAP_UINT64 = 'long'
ERDDAP_FLOAT = 'float'
ERDDAP_DOUBLE = 'double'
ERDDAP_STRING = 'String'
def erddap_datatype(nc_dtype):
    ret = None
    if nc_dtype == 'S1':
        ret = ERDDAP_STRING
    elif nc_dtype == 'c':
        ret = ERDDAP_CHAR
    elif nc_dtype == 'i1' or nc_dtype == 'b' or nc_dtype == 'B':
        ret = ERDDAP_BYTE
    elif nc_dtype == 'u1':
        ret = ERDDAP_UBYTE
    elif nc_dtype == 'i2' or nc_dtype == 'h':
        ret = ERDDAP_SHORT
    elif nc_dtype == 'u2':
        ret = ERDDAP_USHORT
    elif nc_dtype == 'i4' or nc_dtype == 'i' or nc_dtype == 'l':
        ret = ERDDAP_INT
    elif nc_dtype == 'u4':
        ret = ERDDAP_UINT
    elif nc_dtype == 'i8':
        ret = ERDDAP_INT64
    elif nc_dtype == 'u8':
        ret = ERDDAP_UINT64
    elif nc_dtype == 'f4' or nc_dtype == 'f':
        ret = ERDDAP_FLOAT
    elif nc_dtype == 'f8' or nc_dtype == 'd':
        ret = ERDDAP_DOUBLE
    return ret

#################################################
# main
#################################################
#
# Process commandline arguments
#
try:
    opts, args = getopt.getopt(sys.argv[1:], '', opt_fields)
except getopt.GetoptError as err:
    print(str(err))
    usage()
    exit(2)

for o, a in opts:
    if o in '--collectionid':
        collection_id = a
    elif o in '--sourcedir':
        source_data_dir = a
    elif o in '--erddapdir':
        erddap_data_dir = a
    elif o in '--serverdir':
        server_nc_directory = a
    elif o in '--only_csvs':
        fetch_metadata = False
        fetch_csvs = True
        create_netcdf_files = False
        create_datasets_xml = False
    elif o in '--only_netcdf':
        fetch_metadata = True
        fetch_csvs = False
        create_netcdf_files = True
        create_datasets_xml = False
    elif o in '--only_datasets_xml':
        fetch_metadata = False
        fetch_csvs = False
        create_netcdf_files = False
        create_datasets_xml = True
    elif o in '--virtual_datasets':
        create_virtual_datasets = True
    elif o in '--window':
        window = int(a)
    elif o in '--rowsperfile':
        rows_per_file = int(a)
    elif o in '--overwrite':
        istest = True
    elif o in '--sample':
        sample_size = int(a)
    elif o in '--verbose':
        verbose = True
    elif o in '--help':
        usage()
        exit(0)
    else:
        assert False, 'unhandled option'

#
# Retrieve the source data CSVs from ScienceBase
#
if fetch_metadata or fetch_csvs:
    sb = SbSession() # No login, to ensure only public data is processed
    sources = retrieve_source_files(sb, collection_id, source_data_dir, fetch_csvs)
else:
    sources = {}

#
# Create netCDF files from the source CSVs
#
if create_netcdf_files:
    for file_name in os.listdir(source_data_dir):
        if file_name in sources:
            metadata = sources[file_name]
        else:
            metadata = None
        csv_to_nc(file_name, source_data_dir, metadata, erddap_data_dir, rows_per_file)

#
# Create the datasets.xml for the netCDF files
#
if create_datasets_xml:
    write_datasets_xml(erddap_data_dir, "./", create_virtual_datasets)
