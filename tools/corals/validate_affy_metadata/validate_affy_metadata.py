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
VALID_EMAIL_RE = re.compile(r'[^@]+@[^@]+\.[^@]+')


def add_error_msg(accumulated_msgs, msg):
    return "%s\n%s" % (accumulated_msgs, msg)


def empty_value(line_no, label, accumulated_msgs):
    return add_error_msg(accumulated_msgs, "The required %s value is missing on line %d." % (label, line_no))


def stop_error(msg):
    sys.exit(msg)


def string_as_boolean_string(string):
    if str(string).lower() in ['true', 'yes', 'on', '1']:
        return 'True'
    else:
        return 'False'


def validate_date_string(line_no, date_string, column, accumulated_msgs):
    if len(date_string) == 0:
        return accumulated_msgs
    try:
        datetime.datetime.strptime(date_string, '%Y-%m-%d')
        return accumulated_msgs
    except ValueError:
        return add_error_msg(accumulated_msgs, "Line %d contains an incorrect date format (%s must be YYYY-MM-DD) for column %s." % (line_no, date_string, column))


def validate_decimal(line_no, decimal_string, column, accumulated_msgs):
    try:
        decimal.Decimal(decimal_string)
        return accumulated_msgs
    except Exception:
        return add_error_msg(accumulated_msgs, "Line %d contains an incorrect decimal value (%s) for column %s." % (line_no, decimal_string, column))


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
        # Keep 1-based line value for error messages.
        line_no = i + 1
        line = line.rstrip("\r\n")
        items = line.split("\t")
        if len(items) != 32:
            accumulated_msgs = add_error_msg(accumulated_msgs, "Line %d contains %s columns, (must be 32)." % (line_no, len(items)))
            stop_error(accumulated_msgs)
        # Required and validated.
        # Required.
        user_specimen_id = items[0]
        if len(user_specimen_id) == 0:
            accumulated_msgs = empty_value(line_no, "user_specimen_id", accumulated_msgs)
        # Optional.
        field_call = items[1]
        # Optional.
        bcoral_genet_id = items[2]
        # Optional.
        bsym_genet_id = items[3]
        # Required.
        reef = items[4]
        if len(reef) == 0:
            accumulated_msgs = empty_value(line_no, "reef", accumulated_msgs)
        # Required.
        region = items[5]
        if len(region) == 0:
            accumulated_msgs = empty_value(line_no, "region", accumulated_msgs)
        # Required and validated.
        latitude = items[6]
        accumulated_msgs = validate_decimal(line_no, latitude, "latitude", accumulated_msgs)
        # Required and validated.
        longitude = items[7]
        accumulated_msgs = validate_decimal(line_no, longitude, "longitude", accumulated_msgs)
        # Optional.
        geographic_origin = items[8]
        # Optional.
        colony_location = items[9]
        depth = items[10]
        # If depth has a value, then it must be a decimal.
        if len(depth) > 0:
            accumulated_msgs = validate_decimal(line_no, depth, "depth", accumulated_msgs)
        # Optional.
        disease_resist = items[11]
        # Optional.
        bleach_resist = items[12]
        # Optional.
        mortality = items[13]
        # Optional.
        tle = items[14]
        # Optional.
        spawning = string_as_boolean_string(items[15])
        # Required.
        collector_last_name = items[16]
        if len(collector_last_name) == 0:
            accumulated_msgs = empty_value(line_no, "collector_last_name", accumulated_msgs)
        # Required.
        collector_first_name = items[17]
        if len(collector_first_name) == 0:
            accumulated_msgs = empty_value(line_no, "collector_first_name", accumulated_msgs)
        # Required.
        org = items[18]
        if len(org) == 0:
            accumulated_msgs = empty_value(line_no, "org", accumulated_msgs)
        # Required and validated.
        collection_date = items[19]
        accumulated_msgs = validate_date_string(line_no, collection_date, "collection_date", accumulated_msgs)
        # Required and validated.
        contact_email = items[20]
        accumulated_msgs = validate_email(line_no, contact_email, accumulated_msgs)
        # Required.
        seq_facility = items[21]
        if len(seq_facility) == 0:
            accumulated_msgs = empty_value(line_no, "seq_facility", accumulated_msgs)
        # Optional.
        array_version = items[22]
        # Optional.
        public = string_as_boolean_string(items[23])
        # Optional.
        public_after_date = items[24]
        accumulated_msga = validate_date_string(line_no, public_after_date, "public_after_date", accumulated_msgs)
        # Required and validated.
        sperm_motility = items[25]
        accumulated_msgs = validate_decimal(line_no, sperm_motility, "sperm_motility", accumulated_msgs)
        # Required and validated.
        healing_time = items[26]
        accumulated_msgs = validate_decimal(line_no, healing_time, "healing_time", accumulated_msgs)
        # Optional.
        dna_extraction_method = items[27]
        # Optional.
        dna_concentration = items[28]
        # If dna_concentration has a value, then it must be decimal.
        if len(dna_concentration) > 0:
            accumulated_msgs = validate_decimal(line_no, dna_concentration, "dna_concentration", accumulated_msgs)
        # Optional.
        registry_id = items[29]
        # Optional.
        result_folder_name = items[30]
        # Optional.
        plate_barcode = items[31]


if len(accumulated_msgs) > 0:
    stop_error(accumulated_msgs)

shutil.copyfile(args.input, args.output)
