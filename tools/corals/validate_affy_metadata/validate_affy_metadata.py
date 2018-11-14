#!/usr/bin/env python
"""
Validate the metadata file associated with Affymetrix 96 well plate data.
"""
import argparse
import datetime
import decimal
import re
import shutil
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--input', dest='input', help='Metadata file for Affymetrix 96 well plate data')
parser.add_argument('--output', dest='output', help='Output dataset'),
args = parser.parse_args()

EMAIL_MAX_LEN = 255
VALID_EMAIL_RE = re.compile("[^@]+@[^@]+\.[^@]+")


def add_error_msg(accumulated_msgs, msg):
    return "%s\n%s" % (accumulated_msgs, msg)


def empty_value(line_no, label, accumulated_msgs):
    return add_error_msg(accumulated_msgs, "The required %s value is missing on line %d." % (label, line_no))


def stop_error(msg):
    sys.exit(msg)


def validate_date_string(line_no, date_string, accumulated_msgs):
    try:
        datetime.datetime.strptime(date_string, '%y/%m/%d')
        return accumulated_msgs
    except ValueError:
        return add_error_msg(accumulated_msgs, "Line %d contains an incorrect date format (%s must be YY/MM/DD)." % (line_no, date_string))


def validate_decimal(line_no, decimal_string, accumulated_msgs):
    try:
        decimal.Decimal(decimal_string)
        return accumulated_msgs
    except Exception:
        return add_error_msg(accumulated_msgs, "Line %d contains an incorrect decimal value (%s)." % (line_no, decimal_string))


def validate_email(line_no, email, accumulated_msgs):
    if not (VALID_EMAIL_RE.match(email)):
        return add_error_msg(accumulated_msgs, "Line %d contains an invalid email address (%s).  " % (line_no, email))
    elif len(email) > EMAIL_MAX_LEN:
        return add_error_msg(accumulated_msgs, "Line %d contains an email address (%) that is longer than the maximum length, %d characters." % (line_no, email))
    return accumulated_msgs


accumulated_msgs = ""
# Parse the input file, skipping the header, and validating
# that each data line consists of 31 comma-separated items.
with open(args.input, "r") as ih:
    for i, line in enumerate(ih):
        if i == 0:
            # Skip the header.
            continue
        line = line.rstrip("\r\n")
        if i > 97:
            accumulated_msgs = add_error_msg(accumulated_msgs, "The input file contains more than 96 data lines.")
            stop_error(accumulated_msgs)
        items = line.split(",")
        if len(items) != 31:
            accumulated_msgs = add_error_msg(accumulated_msgs, "Line %d contains %s columns, (must be 31)." % (i, len(items)))
            stop_error(accumulated_msgs)
        # Required.
        sample_id = items[0]
        if len(sample_id) == 0:
            accumulated_msgs = empty_value(i, "sample_id", accumulated_msgs)
        # Required and validated.
        date_entered_db = items[1]
        accumulated_msgs = validate_date_string(i, date_entered_db, accumulated_msgs)
        # Required.
        user_specimen_id = items[2]
        if len(user_specimen_id) == 0:
            accumulated_msgs = empty_value(i, "user_specimen_id", accumulated_msgs)
        # Optional.
        duplicate_sample = items[3]
        # Optional.
        matching_samples = items[4]
        # Optional.
        field_call = items[5]
        # Optional.
        bcoral_genet_id = items[6]
        # Optional.
        bsym_genet_id = items[7]
        # Required.
        reef = items[8]
        if len(reef) == 0:
            accumulated_msgs = empty_value(i, "reef", accumulated_msgs)
        # Required.
        region = items[9]
        if len(region) == 0:
            accumulated_msgs = empty_value(i, "region", accumulated_msgs)
        # Required and validated.
        latitude = items[10]
        accumulated_msgs = validate_decimal(i, latitude, accumulated_msgs)
        # Required and validated.
        longitude = items[11]
        accumulated_msgs = validate_decimal(i, longitude, accumulated_msgs)
        # Optional.
        geographic_origin = items[12]
        # Optional.
        sample_location = items[13]
        # Optional.
        latitude_outplant = items[14]
        # Optional.
        longitude_outplant = items[15]
        # Optional.
        depth = items[16]
        # Optional.
        dist_shore = items[17]
        # Optional.
        disease_resist = items[18]
        # Optional.
        bleach_resist = items[19]
        # Optional.
        mortality = items[20]
        # Optional.
        tle = items[21]
        # Optional.
        spawning = items[22]
        # Required.
        collector = items[23]
        if len(collector) == 0:
            accumulated_msgs = empty_value(i, "collector", accumulated_msgs)
        # Required.
        org = items[24]
        if len(org) == 0:
            accumulated_msgs = empty_value(i, "org", accumulated_msgs)
        # Required and validated.
        collection_date = items[25]
        accumulated_msgs = validate_date_string(i, date_entered_db, accumulated_msgs)
        # Required and validated.
        contact_email = items[26]
        accumulated_msgs = validate_email(i, contact_email, accumulated_msgs)
        # Required.
        seq_facility = items[27]
        if len(seq_facility) == 0:
            accumulated_msgs = empty_value(i, "seq_facility", accumulated_msgs)
        # Optional.
        array_version = items[28]
        # Optional.
        data_sharing = items[29]
        # Optional.
        data_hold = items[30]

if len(accumulated_msgs) > 0:
    stop_error(accumulated_msgs)

shutil.copyfile(args.input, args.output)
