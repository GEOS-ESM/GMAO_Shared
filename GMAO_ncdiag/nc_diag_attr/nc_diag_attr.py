# nc_diag_attr

from netCDF4 import Dataset, getlibversion
import netCDF4
import argparse
import sys
import traceback
import numpy

try:
    import ujson as json
except:
    import json

# Version information
__version__ = "0.9b"

VERSION_STR = 'nc_diag_attr v' + __version__ + "\n\n" + \
    "Using the following library/runtime versions:\n" + \
    ("    netcdf4-python v%s\n" % netCDF4.__version__) + \
    ("    NetCDF v%s\n" % getlibversion()) + \
    ("    HDF5 v%s\n" % netCDF4.__hdf5libversion__) + \
    ("    Python v%s\n" % sys.version.split("\n")[0].strip())

# CLI Arguments
global args

def parse_cli_args():
    global args
    parser = argparse.ArgumentParser( #prog='ipush',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Tool to add/modify global and variable attributes for NetCDF files",
    version = VERSION_STR)
    
    disable_group = parser.add_mutually_exclusive_group()
    
    parser.add_argument("-V", "--verbose", 
        dest="verbose", action="store_true", default=False,
        help = "enable verbose output")
    parser.add_argument("-p", "--pretty", 
        dest="pretty_output", action="store_true", default=False,
        help = "enable colorful, pretty output - don't enable if logging")
    
    disable_group.add_argument("-ng", "--no-global",
        dest="global_attributes", action="store_false", default=True,
        help = "disable global attribute adding/modifying")
    disable_group.add_argument("-nv", "--no-var", 
        dest="var_attributes", action="store_false", default=True,
        help = "disable variable attribute adding/modifying")
    
    parser.add_argument("-rc", metavar = "RESOURCE_FILE", dest="resource_file",
        help = "input JSON resource file name with attributes to write", required = True)
    parser.add_argument("nc4_files", help = "NetCDF4 files to apply attributes to", nargs="+")
    
    args = parser.parse_args()

def error_msg(msg):
    global args
    if args.pretty_output:
        print("\033[31m **   ERROR: %s\033[0m" % msg)
    else:
        print(" **   ERROR: %s" % msg)

def warning_msg(msg):
    global args
    if args.verbose:
        if args.pretty_output:
            print("\033[33m ** WARNING: %s\033[0m" % msg)
        else:
            print(" ** WARNING: %s" % msg)

def info_msg(msg):
    global args
    if args.verbose:
        if args.pretty_output:
            print("\033[34m **    INFO: %s\033[0m" % msg)
        else:
            print(" **    INFO: %s" % msg)

global current_line
current_line = ""

# ANSI line updater - if enabled!
def line_msg(msg):
    global args, current_line
    if args.pretty_output:
        # Move cursor to beginning:
        sys.stdout.write("\r")
        # Erase the current line
        sys.stdout.write(len(current_line) * " ")
        # Backspace back to the beginning (we could use \r here...)
        sys.stdout.write(len(current_line) * "\b")
        # Print new message
        sys.stdout.write(msg)
        # Go back to beginning
        sys.stdout.write(len(msg) * "\b")
        # Flush output - if not flushed, output may not show up
        sys.stdout.flush()
        # Set new current line
        current_line = msg
    else:
        print(msg)

def line_msg_done():
    global args, current_line
    if args.verbose and args.pretty_output:
        # Move down from current line and erase current line buffer
        sys.stdout.write("\n")
        sys.stdout.flush()
        current_line = ""

global entry_num, entry_total, entry_str
def init_counter(total_ele, entry):
    global entry_num, entry_total, entry_str
    if args.verbose:
        entry_num = 0
        entry_total = total_ele
        entry_str = entry

def progress_counter(filename):
    global entry_num, entry_total, entry_str
    if args.verbose:
        entry_num += 1
        line_msg("%s %i/%i: %s" % (entry_str, entry_num, entry_total, filename))

def main():
    # Parse arguments
    parse_cli_args()
    
    # Sanity checks
    
    # Check to make sure that the JSON resource file exists!
    try:
        resource_file_fh = open(args.resource_file, "r")
    except IOError:
        error_msg("Resource file '%s' is not accessible or does not exist!" % args.resource_file)
        exit(1)
    
    # Check to make sure that the JSON resource file is valid!
    try:
        resource_data = json.loads(resource_file_fh.read())
    except KeyboardInterrupt:
        info_msg("CTRL-C detected, exiting.")
        exit(0)
    except:
        error_msg("Resource file '%s' is not a valid JSON file!" % args.resource_file)
        print(traceback.format_exc())
        exit(1)
    
    # Close file - we got the data already!
    resource_file_fh.close()
    
    # Print verbose version information
    if args.verbose:
        info_msg("Using following versions:")
        info_msg("    netcdf4-python v%s" % netCDF4.__version__)
        info_msg("    NetCDF v%s" % getlibversion())
        info_msg("    HDF5 v%s" % netCDF4.__hdf5libversion__)
        info_msg("    Python v%s\n" % sys.version.split("\n")[0].strip())
        info_msg("Reading and validating NetCDF4 files...")
    
    # Check to make sure the NetCDF4 files are legitimate!
    nc4_files_root = []
    
    init_counter(len(args.nc4_files), "Reading/verifying file")
    for nc4_file in args.nc4_files:
        try:
            open(nc4_file, "r").close()
        except KeyboardInterrupt:
            info_msg("CTRL-C detected, exiting.")
            exit(0)
        except IOError:
            error_msg("The NetCDF4 file '%s' does not exist!" % nc4_file)
            exit(1)
        
        progress_counter(nc4_file)
        
        try:
            rootgrp = Dataset(nc4_file, "a", format="NETCDF4")
            nc4_files_root.append({ "file" : nc4_file, "group" : rootgrp })
        except KeyboardInterrupt:
            info_msg("CTRL-C detected, exiting.")
            exit(0)
        except:
            error_msg("'%s' is not a valid NetCDF4 file!" % nc4_file)
            exit(1)
    
    line_msg_done()
    
    # Global attributes
    if args.global_attributes:
        # Check if we have a global attributes entry in the resource file
        if not "global_attributes" in resource_data:
            warning_msg("Resource file '%s' does not have any global attributes, skipping." % args.resource_file)
        else:
            # Initialize our counter
            init_counter(len(nc4_files_root), "Applying global attributes to file")
            
            for nc4_entry in nc4_files_root:
                # Update progress counter
                progress_counter(nc4_entry["file"])
                
                for global_attr_key in resource_data["global_attributes"]:
                    global_attr_val = resource_data["global_attributes"][global_attr_key]
                    
                    # We need to convert unicode to ASCII
                    if type(global_attr_val) == unicode:
                        global_attr_val = str(global_attr_val)
                    
                    # BUG fix - NetCDF really, really, REALLY does not like
                    # 64-bit integers. We forcefully convert the value to a 
                    # 32-bit signed integer, with some help from numpy!
                    if type(global_attr_val) == int:
                        global_attr_val = numpy.int32(global_attr_val)
                    
                    setattr(nc4_entry["group"], global_attr_key, global_attr_val)
            line_msg_done()
    
    # Variable attributes
    if args.var_attributes:
        # Check if we have a variable attributes entry in the resource file
        if not "variable_attributes" in resource_data:
            warning_msg("Resource file '%s' does not have any variable attributes, skipping." % args.resource_file)
        else:
            # Initialize our counter
            init_counter(len(nc4_files_root), "Applying variable attributes to file")
            
            for nc4_entry in nc4_files_root:
                # Update progress counter
                progress_counter(nc4_entry["file"])
                
                # Iterate through all of our var_attr variables
                for var in resource_data["variable_attributes"]:
                    if var in nc4_entry["group"].variables.keys():
                        for var_attr_key in resource_data["variable_attributes"][var]:
                            var_attr_val = resource_data["variable_attributes"][var][var_attr_key]
                            var_attr_key = str(var_attr_key)
                            
                            # We need to convert unicode to ASCII
                            if type(var_attr_val) == unicode:
                                var_attr_val = list(str(var_attr_val))
                            
                            # BUG fix - NetCDF really, really, REALLY does not like
                            # 64-bit integers. We forcefully convert the value to a 
                            # 32-bit signed integer, with some help from numpy!
                            if type(var_attr_val) == int:
                                var_attr_val = numpy.int32(var_attr_val)
                            
                            setattr(nc4_entry["group"].variables[var], var_attr_key, var_attr_val)
                    else:
                        warning_msg("Can't find variable %s in file %s!" % (var, nc4_entry["file"]))
            line_msg_done()
    
    # Close everything
    init_counter(len(nc4_files_root), "Saving changes to file")
    for nc4_entry in nc4_files_root:
        progress_counter(nc4_entry["file"])
        nc4_entry["group"].close()
    
    line_msg_done()
    
    info_msg("Attribute appending complete!")

if __name__ == "__main__":
    main()
