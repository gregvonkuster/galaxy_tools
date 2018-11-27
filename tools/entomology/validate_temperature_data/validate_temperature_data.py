#!/usr/bin/env python
import argparse
import datetime
import decimal
import re
import shutil
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--data_type', dest='data_type', default=None, help='Temperature data type, normals or actuals')
parser.add_argument('--input_actuals', dest='input_actuals', default=None, help='Daily actuals temperature data')
parser.add_argument('--input_normals', dest='input_normals', default=None, help='30 year normals temperature data')
parser.add_argument('--output', dest='output', help='Output dataset'),
args = parser.parse_args()

ACTUALS_HEADER = "LATITUDE,LONGITUDE,DATE,DOY,TMIN,TMAX"
NORMALS_HEADER = "stationid,latitude,longitude,elev_m,name,st,mmdd,doy,tmin,tmax"

def add_error_msg(accumulated_msgs, msg):
    return "%s\n%s" % (accumulated_msgs, msg)


def empty_value(line_no, label, accumulated_msgs):
    return add_error_msg(accumulated_msgs, "The required %s value is missing on line %d." % (label, line_no))


def stop_error(msg):
    sys.exit(msg)


def validate_date_string(line_no, date_string, accumulated_msgs):
    try:
        datetime.datetime.strptime(date_string, '%Y-%m-%d')
        return accumulated_msgs
    except ValueError:
        return add_error_msg(accumulated_msgs, "Line %d contains an incorrect date format (%s must be YYYY-MM-DD)." % (line_no, date_string))


def validate_decimal(line_no, decimal_string, accumulated_msgs, label):
    try:
        decimal.Decimal(decimal_string)
        return accumulated_msgs
    except Exception:
        return add_error_msg(accumulated_msgs, "Line %d contains an incorrect %s decimal value (%s)." % (line_no, label, decimal_string))


def validate_integer(line_no, integer_string, accumulated_msgs, label):
    if integer_string.isdigit():
        return accumulated_msgs
    return add_error_msg(accumulated_msgs, "Line %d contains an incorrect %s integer value (%s)." % (line_no, label, integer_string))


def validate_mmdd(line_no, mmdd, accumulated_msgs):
    try:
        datetime.datetime.strptime(mmdd, '%m-%d')
        return accumulated_msgs
    except ValueError:
        # Handle Feb 29.
        items = mmdd.split("-")
        try:
            month = int(items[0])
            day = int(items[1])
            if month == 2 and day == 29:
                return accumulated_msgs
        except Exception:
            # Error message accumulated below.
            pass
        return add_error_msg(accumulated_msgs, "Line %d contains an incorrect date format (%s must be mm-dd)." % (line_no, mmdd))


accumulated_msgs = ""
last_doy = 0
# Parse the input file, skipping the header, and validating
# that each data line consists of 31 comma-separated items.
if args.data_type == "normals":
    input_file = args.input_normals
    num_normals_rows = 0
else:
    input_file = args.input_actuals
with open(input_file, "r") as ih:
    for i, line in enumerate(ih):
        line = line.rstrip("\r\n")
        items = line.split(",")
        if args.data_type == "normals":
            num_normals_rows += 1
            if i == 0:
                if line != NORMALS_HEADER:
                    accumulated_msgs = add_error_msg(accumulated_msgs, "The header is invalid, must be %s" % NORMALS_HEADER)
                continue
            if i > 367:
                accumulated_msgs = add_error_msg(accumulated_msgs, "The input file contains more than 367 lines (must be 1 header line and 366 data lines).")
                stop_error(accumulated_msgs)
            if len(items) != 10:
                accumulated_msgs = add_error_msg(accumulated_msgs, "Line %d contains %s columns, (must be 10)." % (i, len(items)))
                stop_error(accumulated_msgs)
            stationid = items[0].strip()
            if len(stationid) == 0:
                accumulated_msgs = empty_value(i, "stationid", accumulated_msgs)
            latitude = items[1].strip()
            accumulated_msgs = validate_decimal(i, latitude, accumulated_msgs, "latitude")
            longitude = items[2].strip()
            accumulated_msgs = validate_decimal(i, longitude, accumulated_msgs, "longitude")
            elev_m = items[3].strip()
            accumulated_msgs = validate_decimal(i, elev_m, accumulated_msgs, "elev_m")
            name = items[4].strip()
            if len(name) == 0:
                accumulated_msgs = empty_value(i, "name", accumulated_msgs)
            st = items[5].strip()
            if len(st) == 0:
                accumulated_msgs = empty_value(i, "st", accumulated_msgs)
            mmdd = items[6].strip()
            accumulated_msgs = validate_mmdd(i, mmdd, accumulated_msgs)
            doy = items[7].strip()
            accumulated_msgs = validate_integer(i, doy, accumulated_msgs, "doy")
            # Make sure the DOY values are consecutive.
            try:
                if int(doy) != (last_doy + 1):
                    accumulated_msgs = add_error_msg(accumulated_msgs, "Line %d contains a DOY (%s) that is not conexcutive (previous DOY is %d)." % (i, doy, last_doy))
                    stop_error(accumulated_msgs)
                else:
                    last_doy += 1
            except Exception:
                # The error for an invalid integer was captured above.
                pass
            tmin = items[8].strip()
            accumulated_msgs = validate_decimal(i, tmin, accumulated_msgs, "tmin")
            tmax = items[9].strip()
            accumulated_msgs = validate_decimal(i, tmax, accumulated_msgs, "tmax")
        else:
            if i == 0:
                if line != ACTUALS_HEADER:
                    accumulated_msgs = add_error_msg(accumulated_msgs, "The header is invalid, must be %s" % ACTUALS_HEADER)
                continue
            if i > 367:
                accumulated_msgs = add_error_msg(accumulated_msgs, "The input file contains more than 367 lines (must be 1 header line and no more than 366 data lines).")
                stop_error(accumulated_msgs)
            if len(items) != 6:
                accumulated_msgs = add_error_msg(accumulated_msgs, "Line %d contains %s columns, (must be 6)." % (i, len(items)))
                stop_error(accumulated_msgs)
            latitude = items[0].strip()
            accumulated_msgs = validate_decimal(i, latitude, accumulated_msgs, "LATITUDE")
            longitude = items[1].strip()
            accumulated_msgs = validate_decimal(i, longitude, accumulated_msgs, "LONGITUDE")
            date_string = items[2].strip()
            accumulated_msgs = validate_date_string(i, date_string, accumulated_msgs)
            doy = items[3].strip()
            accumulated_msgs = validate_integer(i, doy, accumulated_msgs, "doy")
            # Make sure the DOY values are consecutive.
            if i == 1:
                last_doy = int(doy)
            else:
                try:
                    if int(doy) != (last_doy + 1):
                        accumulated_msgs = add_error_msg(accumulated_msgs, "Line %d contains a DOY (%s) that is not conexcutive (previous DOY is %d)." % (i, doy, last_doy))
                        stop_error(accumulated_msgs)
                    else:
                        last_doy += 1
                except Exception:
                    # The error for an invalid integer was captured above.
                    pass
            tmin = items[4].strip()
            accumulated_msgs = validate_decimal(i, tmin, accumulated_msgs, "tmin")
            tmax = items[5].strip()
            accumulated_msgs = validate_decimal(i, tmax, accumulated_msgs, "tmax")
    if args.data_type == "normals" and num_normals_rows != 367:
        accumulated_msgs = add_error_msg(accumulated_msgs, "The input file contains %d rows, (must be 367)." % num_normals_rows)

if len(accumulated_msgs) > 0:
    stop_error(accumulated_msgs)

shutil.copyfile(input_file, args.output)
