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
        fh.write('<html><head><h3>%s: %d files</h3></head>\n' % (title, len(os.listdir(dir))))
        fh.write('<body><p/><table cellpadding="2">\n')
        fh.write('<tr><th>Size</th><th>Name</th></tr>\n')
        for index, fname in enumerate(sorted(os.listdir(dir))):
            if index % 2 == 0:
                bgcolor = '#D8D8D8'
            else:
                bgcolor = '#FFFFFF'
            try:
                size = str(os.path.getsize(os.path.join(dir, fname)))
            except:
                size = 'unknown'
            link = '<a href="%s" type="text/plain">%s</a>\n' % (fname, fname)
            fh.write('<tr bgcolor="%s"><td>%s</td><td>%s</td></tr>\n' % (bgcolor, size, link))
        fh.write('</table></body></html>\n')
