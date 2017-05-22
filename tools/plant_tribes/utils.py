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


def move_directory_files(source_dir, destination_dir, copy=False):
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


def write_html_output(output, title, dir):
    with open(output, 'w') as fh:
        dir_items = sorted(os.listdir(dir))
        # Directories can only contain either files or directories, but not both.
        item_path = os.path.join(dir, dir_items[0])
        if os.path.isdir(item_path):
            header = 'Directories'
        else:
            header = 'Datasets'
        fh.write('<html><head><h3>%s: %d items</h3></head>\n' % (title, len(dir_items)))
        fh.write('<body><p/><table cellpadding="2">\n')
        fh.write('<tr><th>%s</th></tr>\n' % header)
        for index, fname in enumerate(dir_items):
            if index % 2 == 0:
                bgcolor = '#D8D8D8'
            else:
                bgcolor = '#FFFFFF'
            link = '<a href="%s" type="text/plain">%s</a>\n' % (fname, fname)
            fh.write('<tr bgcolor="%s"><td>%s</td></tr>\n' % (bgcolor, link))
        fh.write('</table></body></html>\n')
