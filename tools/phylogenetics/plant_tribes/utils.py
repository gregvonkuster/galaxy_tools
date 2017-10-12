import os
import shutil
import subprocess
import sys

FSTDERR = 'stderr.txt'
FSTDOUT = 'stdout.txt'


def check_execution_errors(rc, fstderr, fstdout):
    if rc != 0:
        fh = open(fstdout, 'rb')
        out_msg = fh.read()
        fh.close()
        fh = open(fstderr, 'rb')
        err_msg = fh.read()
        fh.close()
        msg = '%s\n%s\n' % (str(out_msg), str(err_msg))
        stop_err(msg)


def get_response_buffers():
    fstderr = os.path.join(os.getcwd(), FSTDERR)
    fherr = open(fstderr, 'wb')
    fstdout = os.path.join(os.getcwd(), FSTDOUT)
    fhout = open(fstdout, 'wb')
    return fstderr, fherr, fstdout, fhout


def move_directory_files(source_dir, destination_dir, copy=False, remove_source_dir=False):
    source_directory = os.path.abspath(source_dir)
    destination_directory = os.path.abspath(destination_dir)
    if not os.path.isdir(destination_directory):
        os.makedirs(destination_directory)
    for dir_entry in os.listdir(source_directory):
        source_entry = os.path.join(source_directory, dir_entry)
        if copy:
            shutil.copy(source_entry, destination_directory)
        else:
            shutil.move(source_entry, destination_directory)
    if remove_source_dir:
        os.rmdir(source_directory)


def run_command(cmd):
    fstderr, fherr, fstdout, fhout = get_response_buffers()
    proc = subprocess.Popen(args=cmd, stderr=fherr, stdout=fhout, shell=True)
    rc = proc.wait()
    # Check results.
    fherr.close()
    fhout.close()
    check_execution_errors(rc, fstderr, fstdout)


def stop_err(msg):
    sys.exit(msg)
